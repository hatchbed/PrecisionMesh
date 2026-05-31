// CGAL must be included before any OCCT header (OCCT defines a conflicting Handle macro).
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <precision_mesh/step_tessellation.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <list>
#include <unordered_map>
#include <unordered_set>

#include <spdlog/spdlog.h>

#include <BRep_Tool.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>

namespace PMP = CGAL::Polygon_mesh_processing;


// ── Inverted-face repair (A2-detect-fix) ─────────────────────────────────────
// Some valid sub-faces with extreme UV-parameter anisotropy are mis-triangulated by
// BRepMesh (folded / crossed triangles).  We detect triangles whose winding opposes
// the face's outward normal and rebuild the face's connectivity from its EXISTING
// boundary + interior nodes via a constrained Delaunay triangulation on the best-fit
// plane.  Vertex positions are never changed — only connectivity — so all geometry
// invariants (vertices on exact edges/faces) hold trivially.  CGAL's exact predicates
// are robust to the parameter anisotropy that defeats BRepMesh.

typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K> CRVb;
typedef CGAL::Triangulation_face_base_with_info_2<int, K>           CRFbInfo;
typedef CGAL::Constrained_triangulation_face_base_2<K, CRFbInfo>    CRFb;
typedef CGAL::Triangulation_data_structure_2<CRVb, CRFb>           CRTDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CRTDS, CGAL::Exact_predicates_tag> CRCDT;

// Standard CGAL nesting-level domain marking: a face is inside iff its level is odd.
static void cr_mark_domains(CRCDT& cdt) {
    for (auto f = cdt.all_faces_begin(); f != cdt.all_faces_end(); ++f)
        f->info() = -1;
    std::list<CRCDT::Face_handle> queue;
    CRCDT::Face_handle inf = cdt.infinite_face();
    inf->info() = 0;
    queue.push_back(inf);
    while (!queue.empty()) {
        CRCDT::Face_handle fh = queue.front();
        queue.pop_front();
        for (int i = 0; i < 3; i++) {
            CRCDT::Face_handle nb = fh->neighbor(i);
            if (nb->info() != -1) continue;
            nb->info() = fh->info() + (cdt.is_constrained(CRCDT::Edge(fh, i)) ? 1 : 0);
            queue.push_back(nb);
        }
    }
}

static std::array<double,3> cr_normalize(std::array<double,3> v) {
    double len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (len > 1e-300) { v[0]/=len; v[1]/=len; v[2]/=len; }
    return v;
}

// Outward (solid) normal of the face: surface normal at the UV centre, flipped for
// TopAbs_REVERSED faces — matching the convention BRepMesh's triangles already use.
static std::array<double,3> cr_face_outward_normal(const TopoDS_Face& face) {
    Standard_Real u1,u2,v1,v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);
    std::array<double,3> n = {0,0,1};
    Handle(Geom_Surface) s = BRep_Tool::Surface(face);
    if (!s.IsNull()) {
        gp_Pnt p; gp_Vec du, dv;
        s->D1(0.5*(u1+u2), 0.5*(v1+v2), p, du, dv);
        gp_Vec c = du.Crossed(dv);
        n = {c.X(), c.Y(), c.Z()};
    }
    n = cr_normalize(n);
    if (face.Orientation() == TopAbs_REVERSED) { n[0]=-n[0]; n[1]=-n[1]; n[2]=-n[2]; }
    return n;
}

static std::array<double,3> cr_tri_normal(const std::vector<Point>& vb,
                                          size_t a, size_t b, size_t c) {
    double e1[3] = {(double)vb[b][0]-(double)vb[a][0],
                    (double)vb[b][1]-(double)vb[a][1],
                    (double)vb[b][2]-(double)vb[a][2]};
    double e2[3] = {(double)vb[c][0]-(double)vb[a][0],
                    (double)vb[c][1]-(double)vb[a][1],
                    (double)vb[c][2]-(double)vb[a][2]};
    return {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
}

// Boundary constraint segments as (vertex_buffer index) pairs.  Reconstructed from each
// edge's 3D CURVE + vertex projection, NOT from PolygonOnTriangulation node indices:
// for shared/sewn edges BRepMesh can index the polygon against the wrong triangulation
// (crossed seams) or drop wire corners entirely (T-junctions) — both corruptions we are
// repairing.  Edge endpoint vertices (from the wire) are always included, creating a new
// mesh vertex in `vb` for any corner BRepMesh dropped.  Orientation is irrelevant
// (a-b == b-a as a constraint).
static std::vector<std::pair<size_t,size_t>> cr_boundary_segments(
        const TopoDS_Face& face, const std::unordered_map<int,size_t>& vmap,
        std::vector<Point>& vb, int seg_idx) {
    std::vector<size_t> verts;
    verts.reserve(vmap.size());
    for (const auto& kv : vmap) verts.push_back(kv.second);

    // Tolerance for "vertex lies on this edge", relative to the face's 3D extent.
    double lo[3] = {1e300,1e300,1e300}, hi[3] = {-1e300,-1e300,-1e300};
    for (size_t vi : verts)
        for (int j = 0; j < 3; j++) {
            double cc = (double)vb[vi][j];
            lo[j] = std::min(lo[j], cc); hi[j] = std::max(hi[j], cc);
        }
    double diag = std::sqrt((hi[0]-lo[0])*(hi[0]-lo[0]) +
                            (hi[1]-lo[1])*(hi[1]-lo[1]) +
                            (hi[2]-lo[2])*(hi[2]-lo[2]));
    double tol = std::max(1e-9, 1e-6 * diag);
    double tol2 = tol * tol;

    // Match a 3D point to an existing candidate vertex, or append a new one (recovering
    // wire corners BRepMesh dropped).  `pool` grows so a corner shared by two edges is
    // created once and reused.
    std::vector<size_t> pool = verts;
    auto dist2 = [&](size_t vi, const gp_Pnt& p) {
        double dx=(double)vb[vi][0]-p.X(), dy=(double)vb[vi][1]-p.Y(), dz=(double)vb[vi][2]-p.Z();
        return dx*dx + dy*dy + dz*dz;
    };
    auto match_or_create = [&](const gp_Pnt& p) -> size_t {
        size_t best = 0; double bd = 1e300; bool any = false;
        for (size_t vi : pool) { double d = dist2(vi, p); if (d < bd) { bd = d; best = vi; any = true; } }
        if (any && bd <= tol2) return best;
        vb.push_back({p.X(), p.Y(), p.Z()});
        size_t ni = vb.size() - 1;
        pool.push_back(ni);
        return ni;
    };

    std::vector<std::pair<size_t,size_t>> segs;
    int edge_idx = 0;
    for (TopExp_Explorer ee(face, TopAbs_EDGE); ee.More(); ee.Next(), edge_idx++) {
        TopoDS_Edge edge = TopoDS::Edge(ee.Current());
        Standard_Real f, l;
        Handle(Geom_Curve) crv = BRep_Tool::Curve(edge, f, l);
        if (crv.IsNull()) {
            spdlog::debug("    seg {} bseg: edge {} has no 3D curve", seg_idx, edge_idx);
            continue;
        }
        // Vertices lying on this edge, ordered by curve parameter.
        std::vector<std::pair<double,size_t>> on_edge;
        for (size_t vi : verts) {
            gp_Pnt p((double)vb[vi][0], (double)vb[vi][1], (double)vb[vi][2]);
            GeomAPI_ProjectPointOnCurve proj(p, crv, f, l);
            if (proj.NbPoints() > 0 && proj.LowerDistance() <= tol)
                on_edge.push_back({proj.LowerDistanceParameter(), vi});
        }
        // Always include the edge's two endpoint vertices; create them if BRepMesh
        // dropped the corner (otherwise the quad loses a node → T-junction / hole).
        TopoDS_Vertex evs[2] = { TopExp::FirstVertex(edge), TopExp::LastVertex(edge) };
        for (const TopoDS_Vertex& ev : evs) {
            if (ev.IsNull()) continue;
            gp_Pnt p = BRep_Tool::Pnt(ev);
            bool present = false;
            for (const auto& pr : on_edge) if (dist2(pr.second, p) <= tol2) { present = true; break; }
            if (present) continue;
            size_t vi = match_or_create(p);
            GeomAPI_ProjectPointOnCurve proj(p, crv, f, l);
            double param = (proj.NbPoints() > 0) ? proj.LowerDistanceParameter() : f;
            on_edge.push_back({param, vi});
        }
        std::sort(on_edge.begin(), on_edge.end());
        for (size_t i = 0; i + 1 < on_edge.size(); i++)
            if (on_edge[i].second != on_edge[i+1].second)
                segs.push_back({on_edge[i].second, on_edge[i+1].second});
    }
    return segs;
}

// Outcome of an inverted-face repair attempt (for global stats).
enum class CrBail { Ok, FewSegments, NotPlanar, SelfIntersect };

// Rebuild a face's triangulation from its existing nodes. Returns {} (caller keeps the
// original) if the boundary can't be recovered or the face isn't near-planar enough.
static std::vector<std::array<size_t,3>> cr_retriangulate_face(
        const TopoDS_Face& face, const Handle(Poly_Triangulation)& tri,
        const std::unordered_map<int,size_t>& vmap,
        std::vector<Point>& vb, const std::array<double,3>& n_out, int seg_idx,
        CrBail& reason) {
    reason = CrBail::Ok;
    auto segs = cr_boundary_segments(face, vmap, vb, seg_idx);  // may append corner verts to vb
    if (segs.size() < 3) {
        reason = CrBail::FewSegments;
        return {};
    }

    // Vertex set = interior nodes (vmap) ∪ boundary nodes from segments (which may
    // include corner vertices just created in vb to replace ones BRepMesh dropped).
    std::vector<size_t> verts;
    {
        std::unordered_set<size_t> s;
        for (const auto& kv : vmap) s.insert(kv.second);
        for (const auto& sg : segs) { s.insert(sg.first); s.insert(sg.second); }
        verts.assign(s.begin(), s.end());
    }
    if (verts.size() < 3) { reason = CrBail::FewSegments; return {}; }

    double c[3] = {0,0,0};
    for (size_t vi : verts) { c[0]+=(double)vb[vi][0]; c[1]+=(double)vb[vi][1]; c[2]+=(double)vb[vi][2]; }
    double inv = 1.0 / verts.size();
    c[0]*=inv; c[1]*=inv; c[2]*=inv;

    // Best-fit-plane basis from the axis least aligned with the face outward normal.
    std::array<double,3> nrm = cr_normalize(n_out);
    if (nrm[0]==0 && nrm[1]==0 && nrm[2]==0) { reason = CrBail::NotPlanar; return {}; }
    int mi = 0;
    if (std::abs(nrm[1]) < std::abs(nrm[mi])) mi = 1;
    if (std::abs(nrm[2]) < std::abs(nrm[mi])) mi = 2;
    std::array<double,3> ax = {0,0,0}; ax[mi] = 1.0;
    double da = ax[0]*nrm[0]+ax[1]*nrm[1]+ax[2]*nrm[2];
    std::array<double,3> e1 = cr_normalize({ax[0]-da*nrm[0], ax[1]-da*nrm[1], ax[2]-da*nrm[2]});
    std::array<double,3> e2 = { nrm[1]*e1[2]-nrm[2]*e1[1],
                                nrm[2]*e1[0]-nrm[0]*e1[2],
                                nrm[0]*e1[1]-nrm[1]*e1[0] };

    auto plane2d = [&](size_t vi) -> std::pair<double,double> {
        double dv[3] = {(double)vb[vi][0]-c[0],(double)vb[vi][1]-c[1],(double)vb[vi][2]-c[2]};
        return { dv[0]*e1[0]+dv[1]*e1[1]+dv[2]*e1[2], dv[0]*e2[0]+dv[1]*e2[1]+dv[2]*e2[2] };
    };

    double max_dev = 0, ext = 1e-30;
    for (size_t vi : verts) {
        double dv[3] = {(double)vb[vi][0]-c[0],(double)vb[vi][1]-c[1],(double)vb[vi][2]-c[2]};
        max_dev = std::max(max_dev, std::abs(dv[0]*nrm[0]+dv[1]*nrm[1]+dv[2]*nrm[2]));
        auto uv = plane2d(vi);
        ext = std::max(ext, std::max(std::abs(uv.first), std::abs(uv.second)));
    }
    bool planar = (max_dev <= 0.1 * ext);

    // Near-planar faces triangulate on the best-fit PLANE — robust even when the
    // surface's UV parameterization is folded/non-injective over the sub-face (as on
    // some subdivided faces, where the cap and far edges run in opposite U directions
    // so the seams become crossing diagonals in UV).  Only genuinely curved faces fall
    // back to the (normalized) UV domain.
    bool use_uv = false;
    std::unordered_map<size_t, std::pair<double,double>> uv_of;
    if (!planar) {
        if (!tri.IsNull() && tri->HasUVNodes()) {
            double umin=1e300, umax=-1e300, vmin=1e300, vmax=-1e300;
            for (const auto& kv : vmap) {
                gp_Pnt2d p = tri->UVNode(kv.first);
                uv_of[kv.second] = {p.X(), p.Y()};
                umin = std::min(umin, p.X()); umax = std::max(umax, p.X());
                vmin = std::min(vmin, p.Y()); vmax = std::max(vmax, p.Y());
            }
            double du = umax - umin, dvv = vmax - vmin;
            // Reject seam-straddling sub-faces: across a periodic seam the UV nodes split to
            // both ends of the period, so the UV box spans ~a full period and a UV-domain
            // triangulation maps to overlapping/crossed triangles in 3D.  A non-straddling
            // sub-face spans only its (sub-)sector, well under half a period.
            Handle(Geom_Surface) psurf = BRep_Tool::Surface(face);
            bool seam_straddle = false;
            if (!psurf.IsNull()) {
                if (psurf->IsUPeriodic() && du > 0.5 * psurf->UPeriod()) seam_straddle = true;
                if (psurf->IsVPeriodic() && dvv > 0.5 * psurf->VPeriod()) seam_straddle = true;
            }
            if (seam_straddle) {
                spdlog::debug("    seg {} retri bail: seam-straddling sub-face (UV span {:.4g} x {:.4g})",
                              seg_idx, du, dvv);
                reason = CrBail::SelfIntersect;
                return {};
            }
            if (du > 1e-300 && dvv > 1e-300) {
                for (auto& kv : uv_of) {
                    kv.second.first  = (kv.second.first  - umin) / du;
                    kv.second.second = (kv.second.second - vmin) / dvv;
                }
                use_uv = true;
            }
        }
        if (!use_uv) {
            spdlog::debug("    seg {} retri bail: not planar (dev={:.4g} ext={:.4g}) and no UV",
                          seg_idx, max_dev, ext);
            reason = CrBail::NotPlanar;
            return {};
        }
    }

    auto to2d = [&](size_t vi) -> std::pair<double,double> {
        if (use_uv) {
            auto it = uv_of.find(vi);
            return it != uv_of.end() ? it->second : std::make_pair(0.0, 0.0);
        }
        return plane2d(vi);
    };

    CRCDT cdt;
    std::unordered_map<size_t, CRCDT::Vertex_handle> vh;
    for (size_t vi : verts) {
        auto uv = to2d(vi);
        CRCDT::Vertex_handle h = cdt.insert(CRCDT::Point(uv.first, uv.second));
        h->info() = vi;
        vh[vi] = h;
    }
    size_t n_pts = cdt.number_of_vertices();  // baseline before constraints

    for (const auto& s : segs) {
        auto a = vh.find(s.first);
        auto b = vh.find(s.second);
        if (a != vh.end() && b != vh.end() && a->second != b->second)
            cdt.insert_constraint(a->second, b->second);
    }

    // If the projected boundary self-intersected, CGAL inserts crossing vertices (with
    // no info) — the count grows.  Compare against the post-point baseline (same method),
    // so this is independent of whether number_of_vertices() counts the infinite vertex.
    if (cdt.number_of_vertices() != n_pts) {
        spdlog::debug("    seg {} retri bail: verts {} -> {} after constraints (n_verts_inserted={})",
                      seg_idx, n_pts, (size_t)cdt.number_of_vertices(), verts.size());
        reason = CrBail::SelfIntersect;
        return {};
    }

    cr_mark_domains(cdt);

    std::vector<std::array<size_t,3>> out;
    for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
        if (f->info() < 1 || (f->info() % 2) == 0) continue;  // outside the domain
        size_t a = f->vertex(0)->info(), b = f->vertex(1)->info(), cc = f->vertex(2)->info();
        if (a >= vb.size() || b >= vb.size() || cc >= vb.size()) {  // defensive
            reason = CrBail::SelfIntersect;
            return {};
        }
        std::array<double,3> tn = cr_tri_normal(vb, a, b, cc);
        if (tn[0]*n_out[0] + tn[1]*n_out[1] + tn[2]*n_out[2] < 0.0) std::swap(b, cc);
        out.push_back({a, b, cc});
    }
    return out;
}

std::vector<std::pair<Mesh, TopoDS_Face>> tessellate_shape(const TopoDS_Shape& shape,
                                                            double max_surface_error,
                                                            bool repair_inverted_faces)
{
    IMeshTools_Parameters meshing_params;
    meshing_params.Angle = 90.0;
    meshing_params.AngleInterior = 90.0;
    meshing_params.Deflection = max_surface_error;
    meshing_params.DeflectionInterior = max_surface_error;
    meshing_params.InParallel = true;

    BRepTools::Clean(shape);

    BRepMesh_IncrementalMesh mesher(shape, meshing_params);
    std::vector<Point> vertex_buffer;
    std::vector<std::vector<size_t>> face_buffer;

    std::vector<std::pair<Mesh, TopoDS_Face>> tessellation;
    int seg_idx = 0;

    // Inverted-face repair stats (global picture).
    size_t rep_flagged = 0, rep_fixed = 0;
    size_t rep_bail_segs = 0, rep_bail_planar = 0, rep_bail_selfint = 0, rep_bail_incomplete = 0;

    for (TopExp_Explorer iter(shape, TopAbs_FACE); iter.More(); iter.Next()) {
        Mesh mesh;
        vertex_buffer.clear();
        face_buffer.clear();

        TopLoc_Location loc;
        TopoDS_Face face = TopoDS::Face(iter.Current());

        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);
        if (triangulation.IsNull()) {
            continue;
        }

        auto num_verts = triangulation->NbNodes();
        auto num_triangles_total = triangulation->NbTriangles();
        spdlog::debug("  seg {}: orient={} verts={} tris={}",
                      seg_idx, (int)iter.Current().Orientation(), num_verts, num_triangles_total);
        gp_Trsf transform = loc.IsIdentity() ? gp_Trsf() : loc.Transformation();
        std::unordered_map<int, size_t> vertex_map;
        for (int i = 1; i <= num_verts; i++) {
            auto point = triangulation->Node(i);
            point.Transform(transform);
            vertex_map[i] = vertex_buffer.size();
            vertex_buffer.push_back({point.X(), point.Y(), point.Z()});
        }

        const TopAbs_Orientation orientation = iter.Current().Orientation();
        for (int i = 1; i <= num_triangles_total; i++) {
            auto triangle = triangulation->Triangle(i);

            Standard_Integer anId[3];
            triangle.Get(anId[0], anId[1], anId[2]);
            if (orientation == TopAbs_REVERSED) {
                Standard_Integer tmp = anId[1];
                anId[1] = anId[2];
                anId[2] = tmp;
            }

            if (anId[0] < 1 || anId[0] > num_verts ||
                anId[1] < 1 || anId[1] > num_verts ||
                anId[2] < 1 || anId[2] > num_verts)
            {
                spdlog::warn("  Invalid vertex ids: {}, {}, {} of {}", anId[0] - 1, anId[1] - 1,
                             anId[2] - 1, num_verts);
                continue;
            }

            face_buffer.push_back({vertex_map[anId[0]], vertex_map[anId[1]], vertex_map[anId[2]]});
        }

        // Detect bad BRepMesh output and rebuild the face's connectivity from its nodes:
        //  - folded:  a triangle whose winding opposes the face's outward normal;
        //  - holed:   an isolated vertex (used by no triangle) — a missing triangle;
        //  - dropped: a wire (corner) vertex absent from the triangulation — a quad
        //             tessellated as a triangle, leaving a T-junction with its neighbour.
        if (repair_inverted_faces && !face_buffer.empty()) {
            std::array<double,3> n_out = cr_face_outward_normal(face);

            // Fold test: compare each triangle's winding against the LOCAL surface normal
            // at that triangle (from its UV-node centroid), not one face-center normal —
            // otherwise a correctly-tessellated curved face false-positives (its edge
            // triangles legitimately face >90° from the centre).  Fall back to the
            // face-center normal only when UV nodes are unavailable.
            double orient_sign = (orientation == TopAbs_REVERSED) ? -1.0 : 1.0;
            Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
            bool has_uv = triangulation->HasUVNodes() && !surf.IsNull();
            std::unordered_map<size_t, gp_Pnt2d> uv_of;
            if (has_uv)
                for (const auto& kv : vertex_map)
                    uv_of[kv.second] = triangulation->UVNode(kv.first);

            // A triangle is folded if its winding opposes the local surface normal (UV-node
            // centroid → surface D1), falling back to the face-center normal when UV is
            // unavailable or a vertex has no UV (e.g. a re-created corner).
            auto is_folded = [&](size_t a, size_t b, size_t c) -> bool {
                std::array<double,3> tn = cr_tri_normal(vertex_buffer, a, b, c);
                std::array<double,3> ref = n_out;
                if (has_uv) {
                    auto i0 = uv_of.find(a), i1 = uv_of.find(b), i2 = uv_of.find(c);
                    if (i0 != uv_of.end() && i1 != uv_of.end() && i2 != uv_of.end()) {
                        double u = (i0->second.X()+i1->second.X()+i2->second.X())/3.0;
                        double v = (i0->second.Y()+i1->second.Y()+i2->second.Y())/3.0;
                        gp_Pnt pp; gp_Vec du, dv;
                        surf->D1(u, v, pp, du, dv);
                        gp_Vec nn = du.Crossed(dv);
                        ref = {nn.X()*orient_sign, nn.Y()*orient_sign, nn.Z()*orient_sign};
                    }
                }
                return tn[0]*ref[0] + tn[1]*ref[1] + tn[2]*ref[2] < 0.0;
            };

            int nfold = 0;
            for (const auto& t : face_buffer)
                if (is_folded(t[0], t[1], t[2])) nfold++;
            std::vector<char> used(vertex_buffer.size(), 0);
            for (const auto& t : face_buffer)
                for (int j = 0; j < 3; j++) if (t[j] < used.size()) used[t[j]] = 1;
            int unused = 0;
            for (const auto& kv : vertex_map)
                if (kv.second < used.size() && !used[kv.second]) unused++;

            // Dropped wire vertices: a face boundary vertex not present among the
            // triangulation's nodes (within a tolerance scaled to the face extent).
            int dropped = 0;
            {
                double flo[3]={1e300,1e300,1e300}, fhi[3]={-1e300,-1e300,-1e300};
                for (const auto& kv : vertex_map)
                    for (int j=0;j<3;j++){ double cc=(double)vertex_buffer[kv.second][j];
                        flo[j]=std::min(flo[j],cc); fhi[j]=std::max(fhi[j],cc); }
                double fdiag = std::sqrt((fhi[0]-flo[0])*(fhi[0]-flo[0])+(fhi[1]-flo[1])*(fhi[1]-flo[1])+(fhi[2]-flo[2])*(fhi[2]-flo[2]));
                double ftol2 = std::pow(std::max(1e-9, 1e-6*fdiag), 2);
                for (TopExp_Explorer ve(face, TopAbs_VERTEX); ve.More(); ve.Next()) {
                    gp_Pnt wp = BRep_Tool::Pnt(TopoDS::Vertex(ve.Current()));
                    bool found = false;
                    for (const auto& kv : vertex_map) {
                        double dx=(double)vertex_buffer[kv.second][0]-wp.X();
                        double dy=(double)vertex_buffer[kv.second][1]-wp.Y();
                        double dz=(double)vertex_buffer[kv.second][2]-wp.Z();
                        if (dx*dx+dy*dy+dz*dz <= ftol2) { found = true; break; }
                    }
                    if (!found) { dropped++; break; }
                }
            }

            if (nfold > 0 || unused > 0 || dropped > 0) {
                rep_flagged++;
                spdlog::debug("  seg {}: {} folded, {} unused, {} dropped-corner verts, n_out=[{:.3f},{:.3f},{:.3f}], repairing",
                              seg_idx, nfold, unused, dropped, n_out[0], n_out[1], n_out[2]);
                CrBail reason = CrBail::Ok;
                auto fixed = cr_retriangulate_face(face, triangulation, vertex_map,
                                                   vertex_buffer, n_out, seg_idx, reason);
                // Do no harm: only accept a repair that is demonstrably valid —
                //  (1) it doesn't drop coverage (count >= original; a hole would reduce it);
                //  (2) it has no folded triangles of its own (catches overlapping/crossed
                //      triangulations the count alone can't see, e.g. seam-straddling
                //      sub-faces whose UV is discontinuous);
                //  (3) it doesn't introduce a much thinner triangle than the original — a
                //      sliver, typically from incorporating near-duplicate vertices left by
                //      inconsistent subdivision/sewing, which the repair can't cleanly fix.
                auto min_edge_sq = [&](const auto& tris) {
                    double m = std::numeric_limits<double>::max();
                    for (const auto& t : tris) {
                        for (int e = 0; e < 3; e++) {
                            size_t i = t[e], j = t[(e+1)%3];
                            double dx=(double)vertex_buffer[i][0]-(double)vertex_buffer[j][0];
                            double dy=(double)vertex_buffer[i][1]-(double)vertex_buffer[j][1];
                            double dz=(double)vertex_buffer[i][2]-(double)vertex_buffer[j][2];
                            m = std::min(m, dx*dx+dy*dy+dz*dz);
                        }
                    }
                    return m;
                };
                bool fixed_folded = false;
                for (const auto& t : fixed)
                    if (is_folded(t[0], t[1], t[2])) { fixed_folded = true; break; }
                // 0.25 in squared length = 0.5 in length: reject if the smallest edge shrank
                // by more than half relative to the original tessellation.
                bool introduced_sliver = !fixed.empty() &&
                                         min_edge_sq(fixed) < 0.25 * min_edge_sq(face_buffer);

                if (!fixed.empty() && fixed.size() >= face_buffer.size() &&
                    !fixed_folded && !introduced_sliver) {
                    rep_fixed++;
                    spdlog::debug("  repaired tessellation on seg {}: {} -> {} tris",
                                  seg_idx, face_buffer.size(), fixed.size());
                    face_buffer.clear();
                    for (const auto& t : fixed) face_buffer.push_back({t[0], t[1], t[2]});
                } else if (!fixed.empty()) {
                    rep_bail_incomplete++;
                    spdlog::warn("  seg {}: re-triangulation rejected ({} -> {} tris, folded={}, sliver={}); "
                                 "kept original", seg_idx, face_buffer.size(), fixed.size(),
                                 fixed_folded, introduced_sliver);
                } else {
                    if (reason == CrBail::FewSegments)        rep_bail_segs++;
                    else if (reason == CrBail::NotPlanar)     rep_bail_planar++;
                    else                                      rep_bail_selfint++;
                    spdlog::warn("  seg {}: could not safely re-triangulate flagged face; "
                                 "kept original", seg_idx);
                }
            }
        }

        PMP::repair_polygon_soup(vertex_buffer, face_buffer,
                                 CGAL::parameters::geom_traits(PointArray_traits()));

        std::vector<Mesh::Vertex_index> vertex_indices;
        for (const auto& vertex: vertex_buffer) {
            vertex_indices.push_back(
                mesh.add_vertex(Mesh::Point(vertex[0], vertex[1], vertex[2])));
        }
        for (const auto& f: face_buffer) {
            mesh.add_face(vertex_indices[f[0]], vertex_indices[f[1]], vertex_indices[f[2]]);
        }

        if (mesh.number_of_faces() == 0) {
            continue;
        }

        tessellation.push_back(std::make_pair(mesh, face));
        seg_idx++;
    }

    if (repair_inverted_faces && rep_flagged > 0) {
        size_t bailed = rep_bail_segs + rep_bail_planar + rep_bail_selfint + rep_bail_incomplete;
        spdlog::info("  inverted-face repair: {}/{} segments flagged, {} repaired, {} kept original "
                     "(bail: {} few-seg, {} non-planar, {} self-intersect, {} incomplete)",
                     rep_flagged, seg_idx, rep_fixed, bailed,
                     rep_bail_segs, rep_bail_planar, rep_bail_selfint, rep_bail_incomplete);
    }

    return tessellation;
}

size_t get_short_edge_count(const std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                            double min_edge_length)
{
    size_t short_edges = 0;
    double min_sq = min_edge_length * min_edge_length;

    for (const auto& [mesh, face]: tessellation) {
        if (!Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(face)).IsNull()) {
            continue;
        }

        for (auto e: mesh.edges()) {
            if (!mesh.is_border(e)) {
                continue;
            }
            auto h = mesh.halfedge(e);
            auto p_start = mesh.point(mesh.source(h));
            auto p_end   = mesh.point(mesh.target(h));
            if (CGAL::squared_distance(p_start, p_end) < min_sq) {
                short_edges++;
            }
        }
    }

    return short_edges;
}

double find_surface_error_param(const TopoDS_Shape& shape, double min_edge_length,
                                int max_iterations, double conversion_scale)
{
    spdlog::info("finding max surface error param ...");
    max_iterations = std::max(2, std::min(100, max_iterations));

    double threshold_max = min_edge_length;
    double threshold_min = 0;
    // Deflection tuning only counts short border edges — skip the inverted-face repair.
    const bool kNoRepair = false;
    auto tessellation = tessellate_shape(shape, threshold_max, kNoRepair);

    size_t max_short_edges = get_short_edge_count(tessellation, min_edge_length);
    spdlog::info("  initial short edges: {}", max_short_edges);

    for (int i = 0; i < max_iterations; i++) {
        double threshold = (threshold_min + threshold_max) / 2.0;
        tessellation = tessellate_shape(shape, threshold, kNoRepair);
        size_t num_short_edges = get_short_edge_count(tessellation, min_edge_length);
        if (num_short_edges > max_short_edges * 1.01) {
            threshold_min = threshold;
        }
        else {
            threshold_max = threshold;
        }
        spdlog::info("    {} short edges at {}", num_short_edges, threshold * conversion_scale);
    }

    spdlog::info("  found max surface error param = {}", threshold_max * conversion_scale);
    return threshold_max;
}

std::vector<std::pair<Mesh, TopoDS_Face>> boundary_meshes(const TopoDS_Shape& shape)
{
    std::vector<std::pair<Mesh, TopoDS_Face>> result;
    for (TopExp_Explorer face_iter(shape, TopAbs_FACE); face_iter.More(); face_iter.Next()) {
        TopoDS_Face face = TopoDS::Face(face_iter.Current());

        // Walk the outer wire, sampling one point per edge vertex.
        TopExp_Explorer wire_iter(face, TopAbs_WIRE);
        if (!wire_iter.More()) continue;

        std::vector<gp_Pnt> pts;
        BRepTools_WireExplorer wexp(TopoDS::Wire(wire_iter.Current()), face);
        for (; wexp.More(); wexp.Next()) {
            gp_Pnt p = BRep_Tool::Pnt(wexp.CurrentVertex());
            pts.push_back(p);
        }
        if (pts.size() < 3) continue;

        Mesh mesh;
        std::vector<Mesh::Vertex_index> vis;
        vis.reserve(pts.size());
        for (const auto& p : pts)
            vis.push_back(mesh.add_vertex(Mesh::Point(p.X(), p.Y(), p.Z())));

        // Fan triangulation from vis[0].
        bool any = false;
        for (size_t i = 1; i + 1 < vis.size(); i++) {
            auto f = mesh.add_face(vis[0], vis[i], vis[i + 1]);
            if (f != Mesh::null_face()) any = true;
        }
        if (!any) continue;
        result.emplace_back(std::move(mesh), face);
    }
    return result;
}


bool save_shape_as_step(const std::string& path, const TopoDS_Shape& shape) {
    STEPControl_Writer writer;
    auto status = writer.Transfer(shape, STEPControl_AsIs);
    if (status != IFSelect_RetDone) {
        spdlog::error("Error transferring shape to STEP writer.");
        return false;
    }

    status = writer.Write(path.c_str());
    if (status != IFSelect_RetDone) {
        spdlog::error("Error writing STEP file.");
    }

    return true;
}
