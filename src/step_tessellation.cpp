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
#include <filesystem>
#include <fstream>
#include <limits>
#include <list>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <spdlog/spdlog.h>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom_Plane.hxx>
#include <GProp_GProps.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>

namespace PMP = CGAL::Polygon_mesh_processing;


// -- Inverted-face repair (A2-detect-fix) ------------------------------------
// Some valid sub-faces with extreme UV-parameter anisotropy are mis-triangulated by
// BRepMesh (folded / crossed triangles).  We detect triangles whose winding opposes
// the face's outward normal and rebuild the face's connectivity from its EXISTING
// boundary + interior nodes via a constrained Delaunay triangulation on the best-fit
// plane.  Vertex positions are never changed -- only connectivity -- so all geometry
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
// TopAbs_REVERSED faces -- matching the convention BRepMesh's triangles already use.
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
// (crossed seams) or drop wire corners entirely (T-junctions) -- both corruptions we are
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
        // dropped the corner (otherwise the quad loses a node -> T-junction / hole).
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

    // Vertex set = interior nodes (vmap) union boundary nodes from segments (which may
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

    // Near-planar faces triangulate on the best-fit PLANE -- robust even when the
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
    // no info) -- the count grows.  Compare against the post-point baseline (same method),
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

        // Capture UV parameters for each BRepMesh node while Poly_Triangulation data is live.
        // These seed the local-Newton projection in step_projection.cpp.  Dropped-corner
        // vertices created by cr_boundary_segments get NaN and are filled in remeshing.cpp.
        std::unordered_map<size_t, std::pair<double,double>> vertex_uv_buffer;
        if (triangulation->HasUVNodes()) {
            for (const auto& kv : vertex_map) {
                gp_Pnt2d uv = triangulation->UVNode(kv.first);
                vertex_uv_buffer[kv.second] = {uv.X(), uv.Y()};
            }
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
        //  - holed:   an isolated vertex (used by no triangle) -- a missing triangle;
        //  - dropped: a wire (corner) vertex absent from the triangulation -- a quad
        //             tessellated as a triangle, leaving a T-junction with its neighbour.
        if (repair_inverted_faces && !face_buffer.empty()) {
            std::array<double,3> n_out = cr_face_outward_normal(face);

            // Fold test: compare each triangle's winding against the LOCAL surface normal
            // at that triangle (from its UV-node centroid), not one face-center normal --
            // otherwise a correctly-tessellated curved face false-positives (its edge
            // triangles legitimately face >90 deg from the centre).  Fall back to the
            // face-center normal only when UV nodes are unavailable.
            double orient_sign = (orientation == TopAbs_REVERSED) ? -1.0 : 1.0;
            Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
            bool has_uv = triangulation->HasUVNodes() && !surf.IsNull();
            std::unordered_map<size_t, gp_Pnt2d> uv_of;
            if (has_uv)
                for (const auto& kv : vertex_map)
                    uv_of[kv.second] = triangulation->UVNode(kv.first);

            // A triangle is folded if its winding CLEARLY opposes the local surface normal
            // (UV-node centroid -> surface D1), falling back to the face-center normal when UV
            // is unavailable or a vertex has no UV (e.g. a re-created corner).
            //
            // We require the (normalized) angle to be well past 90 deg, not merely past it: a
            // genuine fold is a full winding inversion (cos ~= -1), whereas a near-degenerate
            // sliver's normal direction is dominated by sub-micron node noise and points
            // essentially at random -- it lands near cos ~= 0 and would false-positive on a
            // bare `dot < 0` test (these slivers are harmless: negligible area, ambiguous
            // winding).  -0.25 (~104 deg) cleanly separates real folds from sliver noise.
            constexpr double FOLD_COS_THRESH = -0.25;
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
                double tnlen = std::sqrt(tn[0]*tn[0]+tn[1]*tn[1]+tn[2]*tn[2]);
                double reflen = std::sqrt(ref[0]*ref[0]+ref[1]*ref[1]+ref[2]*ref[2]);
                if (tnlen < 1e-20 || reflen < 1e-20) return false;   // can't orient a degenerate tri
                double cosang = (tn[0]*ref[0]+tn[1]*ref[1]+tn[2]*ref[2]) / (tnlen*reflen);
                return cosang < FOLD_COS_THRESH;
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
                // Do no harm: only accept a repair that is demonstrably valid --
                //  (1) it doesn't drop coverage (count >= original; a hole would reduce it);
                //  (2) it has no folded triangles of its own (catches overlapping/crossed
                //      triangulations the count alone can't see, e.g. seam-straddling
                //      sub-faces whose UV is discontinuous);
                //  (3) it doesn't introduce a much thinner triangle than the original -- a
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

                // Over-cover guard: a valid repair re-triangulates the SAME region, so its
                // total area can't exceed the BREP face area.  A CDT that doesn't respect a
                // concave/looped boundary fills the notch (or a hole) with triangles spanning
                // empty space outside the face -- those are in-plane (not folded) and not
                // slivers, so only an area check catches them.  Curved faces chord *inside*
                // the surface (area <= face area), so the 1% margin is just for planar noise.
                double face_area = 0.0;
                {
                    GProp_GProps fgp;
                    BRepGProp::SurfaceProperties(face, fgp);
                    face_area = fgp.Mass();
                }
                double fixed_area = 0.0;
                for (const auto& t : fixed) {
                    std::array<double,3> nrm = cr_tri_normal(vertex_buffer, t[0], t[1], t[2]);
                    fixed_area += 0.5 * std::sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2]);
                }
                bool overcovers = !fixed.empty() && face_area > 0.0 &&
                                  fixed_area > face_area * 1.01;

                if (!fixed.empty() && fixed.size() >= face_buffer.size() &&
                    !fixed_folded && !introduced_sliver && !overcovers) {
                    rep_fixed++;
                    spdlog::debug("  repaired tessellation on seg {}: {} -> {} tris",
                                  seg_idx, face_buffer.size(), fixed.size());
                    face_buffer.clear();
                    for (const auto& t : fixed) face_buffer.push_back({t[0], t[1], t[2]});
                } else if (!fixed.empty()) {
                    rep_bail_incomplete++;
                    spdlog::warn("  seg {}: re-triangulation rejected ({} -> {} tris, folded={}, sliver={}, "
                                 "overcover={}); kept original", seg_idx, face_buffer.size(), fixed.size(),
                                 fixed_folded, introduced_sliver, overcovers);
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

        // Populate "v:uv" property map for UV-seeded projection.  vertex_indices[i] is the
        // Mesh::Vertex_index for vertex_buffer[i]; vertex_uv_buffer holds UV for the original
        // BRepMesh nodes (by vertex_buffer index).  New vertices appended by cr_boundary_segments
        // (dropped corners) are absent from vertex_uv_buffer and keep the NaN default.
        {
            const double kNaN = std::numeric_limits<double>::quiet_NaN();
            auto uv_prop = mesh.add_property_map<Mesh::Vertex_index,
                                                  std::pair<double,double>>(
                "v:uv", {kNaN, kNaN}).first;
            for (size_t i = 0; i < vertex_indices.size(); i++) {
                auto uit = vertex_uv_buffer.find(i);
                if (uit != vertex_uv_buffer.end())
                    uv_prop[vertex_indices[i]] = uit->second;
            }
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
    // Deflection tuning only counts short border edges -- skip the inverted-face repair.
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


// -- UV-grid CDT tessellation (Phase C2/C3) -----------------------------------
// For analytic curved surfaces (cylinder, cone, torus, ...) the interior of each
// face mesh is replaced by a uniform UV-parameter grid triangulated with CGAL CDT,
// seeded by the fixed border vertices already established by BRepMesh + split/snap.

bool uv_grid_retessellate(Mesh& mesh, const TopoDS_Face& face, int u_steps, int v_steps,
                          double min_edge_length, size_t face_idx, const std::string& dump_dir)
{
    const double kNaNuv = std::numeric_limits<double>::quiet_NaN();
    auto uv_map = mesh.add_property_map<Mesh::Vertex_index, std::pair<double,double>>(
        "v:uv", {kNaNuv, kNaNuv}).first;
    // Verify at least one border vertex has a valid (non-NaN) UV — proxy for map existence.
    bool uv_ok = false;
    for (auto v : mesh.vertices()) {
        if (mesh.is_border(v) && !std::isnan(uv_map[v].first)) { uv_ok = true; break; }
    }
    if (!uv_ok) {
        spdlog::warn("uv_grid_retessellate: mesh has no valid v:uv values on border vertices");
        return false;
    }

    // Re-evaluate UV values for all border vertices using ShapeAnalysis_Surface
    // to fill NaNs and correct any wrapped/corrupted V coordinates produced by BRepMesh.
    {
        Handle(Geom_Surface) geom_surf = BRep_Tool::Surface(face);
        ShapeAnalysis_Surface sa(geom_surf);
        for (auto v : mesh.vertices()) {
            if (!mesh.is_border(v)) continue;
            auto& uv = uv_map[v];
            auto p = mesh.point(v);
            gp_Pnt p3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
            gp_Pnt2d uv2d = sa.ValueOfUV(p3d, 1e-7);

            if (std::isnan(uv.first) || std::isnan(uv.second)) {
                uv = {uv2d.X(), uv2d.Y()};
            } else {
                // CRITICAL CORRECTION: Overwrite V with the true geometric projection.
                // BRepMesh occasionally assigns boundary nodes to the wrong side of the periodic
                // V domain (e.g., giving V=2.9657 to a 3D point lying on the V=0.0 cap).
                // We keep the original U to preserve unrolled periodic seams, but strictly correct V.
                uv.second = uv2d.Y();
            }
        }
    }

    BRepAdaptor_Surface surf(face);

    // Jacobian-corrected UV step sizes: uv_du_step/uv_dv_step give ~target arc length in 3D.
    double uv_u_scale, uv_v_scale, uv_du_step, uv_dv_step;
    {
        double u_mid = (surf.FirstUParameter() + surf.LastUParameter()) / 2.0;
        double v_mid = (surf.FirstVParameter() + surf.LastVParameter()) / 2.0;
        gp_Pnt cp; gp_Vec d_du, d_dv;
        surf.D1(u_mid, v_mid, cp, d_du, d_dv);
        uv_u_scale = std::max(d_du.Magnitude(), 1e-10);
        uv_v_scale = std::max(d_dv.Magnitude(), 1e-10);
        double u_arc = (surf.LastUParameter() - surf.FirstUParameter()) * uv_u_scale;
        double v_arc = (surf.LastVParameter() - surf.FirstVParameter()) * uv_v_scale;
        double target_u = u_steps > 1 ? u_arc / u_steps : u_arc;
        double target_v = v_steps > 1 ? v_arc / v_steps : v_arc;
        double target   = std::max(std::min(target_u, target_v), min_edge_length);
        uv_du_step = target / uv_u_scale;
        uv_dv_step = target / uv_v_scale;
    }

    const bool u_per = surf.IsUPeriodic();
    const double u_period = u_per ? surf.UPeriod() : 0.0;
    const bool v_per = surf.IsVPeriodic();
    const double v_period = v_per ? surf.VPeriod() : 0.0;

    // If the surface is periodic, we must adjust the step sizes so that the period
    // is an EXACT integer multiple of the step size. This guarantees that the grid
    // of interior points is perfectly symmetric and centered.
    if (u_per && u_period > 1e-9) {
        int n_u_steps = std::max(1, (int)std::round(u_period / uv_du_step));
        uv_du_step = u_period / n_u_steps;
    }
    if (v_per && v_period > 1e-9) {
        int n_v_steps = std::max(1, (int)std::round(v_period / uv_dv_step));
        uv_dv_step = v_period / n_v_steps;
    }

    // 1. Extract ordered border loops by walking border halfedges.
    std::unordered_set<size_t> visited_hes;
    std::vector<std::vector<Mesh::Vertex_index>> loops;
    for (auto h : mesh.halfedges()) {
        if (!mesh.is_border(h)) continue;
        if (visited_hes.count(h.idx())) continue;
        std::vector<Mesh::Vertex_index> loop;
        auto cur = h;
        do {
            visited_hes.insert(cur.idx());
            loop.push_back(mesh.target(cur));
            cur = mesh.next(cur);
        } while (cur != h);
        if (loop.size() >= 3) loops.push_back(loop);
    }
    if (loops.empty()) {
        spdlog::warn("uv_grid_retessellate: no border loops found");
        return false;
    }

    // Bail out on a self-touching wire: two DISTINCT border vertices at the same 3D point
    // (a sliver pinching to a model vertex where the boundary touches itself).  There the
    // CDT's convex-hull fill over-covers the thin outside lens between the two tangent
    // boundary curves — the lens opens into the interior across a non-constraint edge, so no
    // domain/centroid test cleanly excludes it, producing border over-coverage (open edges
    // against the adjacent faces).  BRepMesh respects the trim exactly there, so fall back to
    // it for this face.  (The periodic seam is the SAME mesh vertex appearing twice, not two
    // distinct vertices, so it does not trip this check.)
    {
        std::map<std::tuple<double,double,double>, size_t> pinch_pos;
        for (const auto& lp : loops)
            for (auto v : lp) {
                auto p = mesh.point(v);
                auto key = std::make_tuple(std::round(CGAL::to_double(p.x()) * 1e5),
                                           std::round(CGAL::to_double(p.y()) * 1e5),
                                           std::round(CGAL::to_double(p.z()) * 1e5));
                auto it = pinch_pos.find(key);
                if (it != pinch_pos.end() && it->second != v.idx()) {
                    spdlog::info("uv_grid_retessellate [face {}]: self-touching wire at a pinch "
                                 "point — falling back to BRepMesh for this face", face_idx);
                    return false;
                }
                pinch_pos[key] = v.idx();
            }
    }

    // 2. Build UV arrays for each loop using pure geometric projection.
    // We anchor each loop to the U=0 periodic seam to guarantee they unroll
    // into a perfect vertical rectangle, completely eliminating parallelogram skew.
    Handle(Geom_Surface) geom_surf = BRep_Tool::Surface(face);
    ShapeAnalysis_Surface sa(geom_surf);

    std::vector<std::vector<std::pair<double,double>>> loop_uvs(loops.size());
    for (size_t li = 0; li < loops.size(); li++) {
        size_t N = loops[li].size();

        // --- ARRAY ROTATION PASS ---
        // Find the vertex geometrically closest to the U=0 seam and rotate the array
        // so it becomes Index 0. This forces the unroller to anchor both the top
        // and bottom cylinder caps to the exact same starting column.
        if (u_per && u_period > 1e-9) {
            double min_u = std::numeric_limits<double>::max();
            size_t min_idx = 0;
            for (size_t i = 0; i < N; i++) {
                auto p = mesh.point(loops[li][i]);
                gp_Pnt p3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
                gp_Pnt2d uv2d = sa.ValueOfUV(p3d, 1e-7);

                double canon_u = std::fmod(uv2d.X(), u_period);
                if (canon_u < 0) canon_u += u_period;
                // Snap vertices extremely close to 2pi back to 0
                if (u_period - canon_u < 1e-5) canon_u = 0.0;

                if (canon_u < min_u) {
                    min_u = canon_u;
                    min_idx = i;
                }
            }
            std::rotate(loops[li].begin(), loops[li].begin() + min_idx, loops[li].end());
        }

        loop_uvs[li].resize(N + 1); // Resize to hold the duplicate seam vertex

        // Unroll the original N vertices
        for (size_t i = 0; i < N; i++) {
            auto p = mesh.point(loops[li][i]);
            gp_Pnt p3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
            gp_Pnt2d uv2d = sa.ValueOfUV(p3d, 1e-7);

            double u_curr = uv2d.X();
            double v_curr = uv2d.Y();

            if (i > 0) {
                double u_prev = loop_uvs[li][i-1].first;
                double du = u_curr - u_prev;
                while (u_period > 1e-9 && du > u_period * 0.5)  { u_curr -= u_period; du = u_curr - u_prev; }
                while (u_period > 1e-9 && du < -u_period * 0.5) { u_curr += u_period; du = u_curr - u_prev; }

                double v_prev = loop_uvs[li][i-1].second;
                double dv = v_curr - v_prev;
                while (v_period > 1e-9 && dv > v_period * 0.5)  { v_curr -= v_period; dv = v_curr - v_prev; }
                while (v_period > 1e-9 && dv < -v_period * 0.5) { v_curr += v_period; dv = v_curr - v_prev; }
            }
            loop_uvs[li][i] = {u_curr, v_curr};
        }

        // --- PERIODIC SEAM DUPLICATION ---
        if (u_per && u_period > 1e-9) {
            double u_curr = loop_uvs[li][0].first;
            double u_prev = loop_uvs[li][N-1].first;
            double du = u_curr - u_prev;
            while (du > u_period * 0.5)  { u_curr -= u_period; du = u_curr - u_prev; }
            while (du < -u_period * 0.5) { u_curr += u_period; du = u_curr - u_prev; }

            loop_uvs[li][N] = {u_curr, loop_uvs[li][0].second};
            loops[li].push_back(loops[li][0]); // Append physical vertex index
        } else {
            loop_uvs[li].pop_back(); // Remove extra slot if not periodic
        }

        // ==============================================================
        // Alignment Pass: Shift the entire loop so the minimum U is canonical [0, u_period)
        // ==============================================================
        if (u_per && u_period > 1e-9) {
            double min_u = std::numeric_limits<double>::max();
            for (const auto& uv : loop_uvs[li]) {
                min_u = std::min(min_u, uv.first);
            }
            double target_min = std::fmod(min_u, u_period);
            if (target_min < 0.0) target_min += u_period;
            double shift = target_min - min_u;
            if (std::abs(shift) > 1e-9) {
                for (auto& uv : loop_uvs[li]) {
                    uv.first += shift;
                }
            }
        }
    }

    spdlog::debug("uv_grid_retessellate [face {}]: u_scale={:.4f} v_scale={:.4f} "
                  "du_step={:.4f} dv_step={:.4f} u_per={} loops={}",
                  face_idx, uv_u_scale, uv_v_scale, uv_du_step, uv_dv_step,
                  u_per, (int)loops.size());

    // 3. Build CDT. vertex handle → index in unified vertex table.
    //    Border vertices:  index 0 .. N_border-1  → border_mesh_verts[i]
    //    Interior Steiner: index N_border .. N-1  → interior_info[i-N_border]
    std::map<CRCDT::Vertex_handle, size_t> vh_to_idx;
    std::vector<Mesh::Vertex_index> border_mesh_verts;
    std::vector<std::pair<gp_Pnt, std::pair<double,double>>> interior_info;
    int n_interior = 0;
    int n_mark_outside = 0, n_mark_degen = 0;

    {
        CRCDT cdt;

        bool cdt_ok = true;
        // True if any periodic seam-wrap closing edge was left unconstrained.  When set, the
        // border loop is topologically open at the seam, so nesting-level domain marking would
        // leak — we fall back to centroid classification for those (full-cylinder) faces.
        bool seam_skipped = false;
        std::vector<std::vector<CRCDT::Vertex_handle>> loop_vhs(loops.size());
        try {
            for (size_t li = 0; li < loops.size(); li++) {
                loop_vhs[li].reserve(loops[li].size());
                for (size_t vi = 0; vi < loops[li].size(); vi++) {
                    double u = loop_uvs[li][vi].first, v = loop_uvs[li][vi].second;
                    auto vh = cdt.insert(CRCDT::Point(u * uv_u_scale, v * uv_v_scale));
                    if (!vh_to_idx.count(vh)) {
                        vh_to_idx[vh] = border_mesh_verts.size();
                        border_mesh_verts.push_back(loops[li][vi]);
                    }
                    loop_vhs[li].push_back(vh);
                }
            }

            if (cdt.dimension() < 2) {
                spdlog::warn("uv_grid_retessellate: border UV points are degenerate (dim={})",
                             cdt.dimension());
                cdt_ok = false;
            }

            // Insert constrained edges for each border loop.
            // Interior grid edges are NOT constrained — Steiner points are unconstrained,
            // so the Delaunay criterion governs the interior triangulation.  Constrained
            // interior edges caused CDT failures when wire-cutout boundaries crossed them.
            size_t n_after_border = cdt.number_of_vertices();
            for (size_t li = 0; li < loops.size() && cdt_ok; li++) {
                const auto& vhs = loop_vhs[li];
                for (size_t i = 0; i < vhs.size() && cdt_ok; i++) {
                    auto va = vhs[i], vb = vhs[(i + 1) % vhs.size()];
                    if (va == vb) continue;
                    auto pa = va->point(), pb = vb->point();
                    double du = pa.x() - pb.x(), dv = pa.y() - pb.y();
                    if (du*du + dv*dv < 1e-20) continue;
                    // For periodic faces the last loop edge is the seam-wrap closing edge
                    // (Δu ≈ ±u_period).  Skip it — the seam is implicitly enforced by the
                    // convex hull of the border vertices.
                    if (u_per && i == vhs.size() - 1) {
                        double u_period_scaled = u_period * uv_u_scale;
                        if (std::abs(std::abs(du) - u_period_scaled) < u_period_scaled * 0.01) {
                            spdlog::debug("  loop {}: skip seam-wrap closing constraint |du|={:.4f}",
                                          li, std::abs(du));
                            seam_skipped = true;
                            continue;
                        }
                    }
                    cdt.insert_constraint(va, vb);
                }
            }
            if (cdt.number_of_vertices() != n_after_border) {
                spdlog::warn("uv_grid_retessellate: border constraints created extra CDT vertices"
                             " (self-intersecting UV?)");
                cdt_ok = false;
            }

            // Generate interior grid points with Jacobian-corrected UV spacing.
            // BRepClass_FaceClassifier pre-filters points outside the trimmed face
            // (holes, wire cutouts) before inserting them into the CDT.
            if (cdt_ok) {
                // Grid bounds come from the ACTUAL (unrolled) border-loop UV extent, NOT
                // surf.First/LastUParameter().  For a face trimmed into a non-base period
                // (e.g. a 180° hole at u∈[2π,3π]) the adaptor reports the base surface range
                // [0,2π], which would generate the grid in the wrong period — disjoint from
                // the border verts.  The loop UVs are already correctly unrolled, so their
                // min/max give the true trimmed range for both full and partial cylinders.
                double u_min_g = std::numeric_limits<double>::infinity();
                double u_max_g = -std::numeric_limits<double>::infinity();
                double v_min_g = std::numeric_limits<double>::infinity();
                double v_max_g = -std::numeric_limits<double>::infinity();
                for (const auto& lp : loop_uvs)
                    for (const auto& uv : lp) {
                        u_min_g = std::min(u_min_g, uv.first);
                        u_max_g = std::max(u_max_g, uv.first);
                        v_min_g = std::min(v_min_g, uv.second);
                        v_max_g = std::max(v_max_g, uv.second);
                    }

                size_t N_border_now = border_mesh_verts.size();
                int n_inserted = 0, n_filtered = 0;
                // Wrapping base for surf.Value / classifier = the surface's canonical U origin
                // ([0,2π) for a cylinder).  BRepClass_FaceClassifier rejects params outside the
                // canonical range, and a periodic surf.Value is correct at the canonical angle.
                // (Grid *bounds* above stay in the loop's own period to match the border verts.)
                double u_first = surf.FirstUParameter();

                // STRCITLY GENERATE INTERIOR POINTS ONLY
                // We start 1 step inside (uv_du_step) and end at least 0.5 steps before the far boundary.
                // This guarantees that we NEVER place grid vertices on the left or right boundaries,
                // allowing CDT to bridge directly to the original BRep border vertices without duplicates!
                double start_u = u_min_g + uv_du_step;
                double end_u   = u_max_g - uv_du_step * 0.5;

                int i_u = 1; // Start index at 1 since we are strictly in the interior
                for (double u = start_u; u < end_u; u += uv_du_step, ++i_u) {
                    int i_v = 1;
                    for (double v = v_min_g + uv_dv_step; v < v_max_g - uv_dv_step * 0.5; v += uv_dv_step, ++i_v) {

                        // Apply the uniform checkerboard micro-perturbation to prevent co-circular/collinear degeneracy.
                        // Since we are strictly in the interior, we can safely perturb all generated points.
                        double pert_u = u + (((i_u + i_v) % 2 == 0) ? 1e-5 : -1e-5) * uv_du_step;
                        double pert_v = v + ((i_u % 2 == 0) ? 1e-5 : -1e-5) * uv_dv_step;

                        // Map the perturbed coordinate to canonical UV range
                        double u_cls = pert_u;
                        if (u_per && u_period > 1e-9) {
                            double offset = u_cls - u_first;
                            offset = std::fmod(offset, u_period);
                            if (offset < 0.0) offset += u_period;
                            u_cls = u_first + offset;
                        }

                        BRepClass_FaceClassifier clf(face, gp_Pnt2d(u_cls, pert_v), 1e-6);
                        if (clf.State() != TopAbs_IN && clf.State() != TopAbs_ON) {
                            ++n_filtered;
                            continue;
                        }

                        auto vh = cdt.insert(CRCDT::Point(pert_u * uv_u_scale, pert_v * uv_v_scale));

                        if (!vh_to_idx.count(vh)) {
                            gp_Pnt p3d = surf.Value(u_cls, pert_v); // Unscaled UV
                            vh_to_idx[vh] = N_border_now + interior_info.size();
                            interior_info.push_back({p3d, {u_cls, pert_v}});
                            ++n_inserted;
                        }
                    }
                }
                spdlog::debug("  grid: {} inserted, {} filtered (BRepClass)", n_inserted, n_filtered);

                // Degenerate-area threshold (scaled UV²): CGAL can emit zero-area slivers
                // along exactly-collinear boundary chains (cap edges, seam corners).  These
                // overlap the real triangles and corrupt the rebuild, so we drop any triangle
                // whose 2D area is a tiny fraction of a grid cell.
                double cell_area_scaled = (uv_du_step * uv_u_scale) * (uv_dv_step * uv_v_scale);
                double min_tri_2area = cell_area_scaled * 1e-4;
                auto tri_2area = [](CRCDT::Face_handle f) {
                    auto p0 = f->vertex(0)->point(), p1 = f->vertex(1)->point(),
                         p2 = f->vertex(2)->point();
                    return std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                                    (p2.x() - p0.x()) * (p1.y() - p0.y()));
                };

                // Classify a CDT point (scaled UV) against the trimmed face, unscaling and
                // canonicalising U for the periodic surface.
                auto uv_inside = [&](double sx, double sy) -> bool {
                    double cx = sx / uv_u_scale, cy = sy / uv_v_scale;
                    double cx_cls = cx;
                    if (u_per && u_period > 1e-9) {
                        double offset = cx_cls - u_first;
                        offset = std::fmod(offset, u_period);
                        if (offset < 0.0) offset += u_period;
                        cx_cls = u_first + offset;
                    }
                    BRepClass_FaceClassifier clf(face, gp_Pnt2d(cx_cls, cy), 1e-6);
                    return clf.State() == TopAbs_IN || clf.State() == TopAbs_ON;
                };

                // ---- Triangle selection ---------------------------------------------
                //  * Fully-constrained loops (partial cylinders, holes): CDT nesting-level
                //    domain marking — topological and tolerance-free, keeps thin valid slivers
                //    and drops concave bays / bay-spanning fans correctly.
                //  * Seam-open loops (full periodic cylinders, wrap-closing edge skipped): the
                //    open seam makes nesting marking leak, so classify by the triangle centroid.
                if (!seam_skipped) cr_mark_domains(cdt);
                for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
                    if (tri_2area(f) < min_tri_2area) { f->info() = 0; ++n_mark_degen; continue; }
                    bool inside;
                    if (!seam_skipped) {
                        inside = (f->info() % 2) == 1;   // odd nesting level == inside
                    } else {
                        double cx = (f->vertex(0)->point().x() + f->vertex(1)->point().x() +
                                     f->vertex(2)->point().x()) / 3.0;
                        double cy = (f->vertex(0)->point().y() + f->vertex(1)->point().y() +
                                     f->vertex(2)->point().y()) / 3.0;
                        inside = uv_inside(cx, cy);
                    }
                    f->info() = inside ? 1 : 0;
                    if (inside) ++n_interior; else ++n_mark_outside;
                }
                spdlog::debug("  cdt: {} finite faces, {} inside ({} mode)",
                              (int)cdt.number_of_faces(), n_interior,
                              seam_skipped ? "centroid" : "domain");
            }
        } catch (const std::exception& e) {
            spdlog::warn("uv_grid_retessellate: CDT exception: {}", e.what());
            cdt_ok = false;
        } catch (...) {
            spdlog::warn("uv_grid_retessellate: CDT unknown exception");
            cdt_ok = false;
        }

        if (!cdt_ok || n_interior == 0) {
            if (n_interior == 0 && cdt_ok)
                spdlog::warn("uv_grid_retessellate: CDT produced 0 interior triangles");
            return false;
        }

        // Extract the triangles selected by the marking pass above (domain marking for
        // constrained loops, centroid classification for seam-open full-period faces).
        struct CdtTri { size_t a, b, c; };
        std::vector<CdtTri> cdt_tris;
        cdt_tris.reserve(n_interior);
        int n_missing = 0;
        {
            // Selection verdict was recorded in f->info() (1 = keep) by the marking pass.
            for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
                if (f->info() != 1) continue;
                auto it0 = vh_to_idx.find(f->vertex(0));
                auto it1 = vh_to_idx.find(f->vertex(1));
                auto it2 = vh_to_idx.find(f->vertex(2));
                if (it0 == vh_to_idx.end() || it1 == vh_to_idx.end() || it2 == vh_to_idx.end()) {
                    ++n_missing; continue;
                }
                cdt_tris.push_back({it0->second, it1->second, it2->second});
            }
        }
        if (n_missing > 0)
            spdlog::warn("uv_grid_retessellate [face {}]: {} triangles skipped (CDT vertex not in map)",
                         face_idx, n_missing);
        spdlog::debug("  extraction: {} kept, {} outside, {} degenerate, {} missing",
                      cdt_tris.size(), n_mark_outside, n_mark_degen, n_missing);
        if (cdt_tris.empty()) {
            spdlog::warn("uv_grid_retessellate: no mappable interior triangles");
            return false;
        }

        // Diagnostic export of the 2D CDT (--dump-cdt-obj <dir>).
        if (!dump_dir.empty()) {
            std::string obj_path = (std::filesystem::path(dump_dir) /
                                    ("Face_" + std::to_string(face_idx) + "_CDT.obj")).string();
            std::ofstream out_obj(obj_path);
            if (!out_obj.is_open()) {
                spdlog::warn("uv_grid_retessellate [face {}]: failed to open '{}' for CDT dump",
                             face_idx, obj_path);
            } else {
                std::map<CRCDT::Vertex_handle, int> obj_v_idx;
                int idx = 1;
                for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
                    out_obj << "v " << v->point().x() << " " << v->point().y() << " 0\n";
                    obj_v_idx[v] = idx++;
                }
                for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); ++f) {
                    if (f->info() == 1) {   // selection verdict from the marking pass
                        out_obj << "f " << obj_v_idx[f->vertex(0)] << " "
                                        << obj_v_idx[f->vertex(1)] << " "
                                        << obj_v_idx[f->vertex(2)] << "\n";
                    }
                }
                spdlog::debug("uv_grid_retessellate [face {}]: wrote debug CDT to '{}'",
                              face_idx, obj_path);
            }
        }

        size_t N_border = border_mesh_verts.size();

        // 5. Build new mesh.
        auto n_out = cr_face_outward_normal(face);

        Mesh new_mesh;
        auto new_uv = new_mesh.add_property_map<Mesh::Vertex_index, std::pair<double,double>>(
            "v:uv", {kNaNuv, kNaNuv}).first;
        auto new_locked = new_mesh.add_property_map<Mesh::Vertex_index, bool>(
            "v:border_locked", false).first;

        std::vector<Mesh::Vertex_index> vert_new(N_border + interior_info.size(),
                                                  Mesh::null_vertex());
        std::unordered_map<size_t, Mesh::Vertex_index> old_idx_to_new;

        // Deduplicates vertices that evaluate to the exact same 3D location (periodic seams).
        std::map<std::tuple<double,double,double>, Mesh::Vertex_index> pos_to_new_v;

        for (size_t i = 0; i < N_border; i++) {
            auto old_v = border_mesh_verts[i];
            auto it = old_idx_to_new.find(old_v.idx());
            Mesh::Vertex_index nv;
            if (it == old_idx_to_new.end()) {
                nv = new_mesh.add_vertex(mesh.point(old_v));
                auto p = mesh.point(old_v);
                gp_Pnt p3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
                gp_Pnt2d uv2d = sa.ValueOfUV(p3d, 1e-7);
                new_uv[nv] = {uv2d.X(), uv2d.Y()};
                new_locked[nv] = true;
                old_idx_to_new[old_v.idx()] = nv;
            } else {
                nv = it->second;
            }
            vert_new[i] = nv;

            auto p = new_mesh.point(nv);
            double px = std::round(CGAL::to_double(p.x()) * 1e5);
            double py = std::round(CGAL::to_double(p.y()) * 1e5);
            double pz = std::round(CGAL::to_double(p.z()) * 1e5);
            pos_to_new_v[{px, py, pz}] = nv;
        }

        int seam_merge_count = 0;
        for (size_t i = 0; i < interior_info.size(); i++) {
            const auto& p3d = interior_info[i].first;
            // Round to 5 decimals (approx 1 micron) to cleanly merge periodic seams
            double px = std::round(p3d.X() * 1e5);
            double py = std::round(p3d.Y() * 1e5);
            double pz = std::round(p3d.Z() * 1e5);
            auto key = std::make_tuple(px, py, pz);

            auto it = pos_to_new_v.find(key);
            if (it != pos_to_new_v.end()) {
                // If it hits here, it successfully merged the left and right periodic seam!
                vert_new[N_border + i] = it->second;
                seam_merge_count++;
            } else {
                auto nv = new_mesh.add_vertex(K::Point_3(p3d.X(), p3d.Y(), p3d.Z()));
                new_uv[nv] = interior_info[i].second;
                new_locked[nv] = false;
                vert_new[N_border + i] = nv;
                pos_to_new_v[key] = nv;
            }
        }
        spdlog::debug("uv_grid_retessellate [face {}]: successfully merged {} seam/duplicate vertices",
                      face_idx, seam_merge_count);

        // Seam vertices: for periodic (u_closed) faces the seam physical vertex appears
        // at both u=FirstU and u=LastU as two distinct CDT handles, both mapping to the
        // same new-mesh vertex.  CDT triangles adjacent to both handles become identical
        // mesh faces → pre-filter with a vertex-triple deduplication set.
        std::set<std::array<size_t, 3>> seen_faces;
        int added = 0, rejected = 0;
        for (const auto& tri : cdt_tris) {
            if (tri.a >= vert_new.size() || tri.b >= vert_new.size() || tri.c >= vert_new.size()) continue;
            auto va = vert_new[tri.a], vb = vert_new[tri.b], vc = vert_new[tri.c];
            if (va == Mesh::null_vertex() || vb == Mesh::null_vertex() || vc == Mesh::null_vertex()) continue;
            if (va == vb || vb == vc || va == vc) continue;
            std::array<size_t, 3> fkey = {va.idx(), vb.idx(), vc.idx()};
            std::sort(fkey.begin(), fkey.end());
            if (!seen_faces.insert(fkey).second) continue;  // seam duplicate

            // Consistently orient the triangles based on the face's orientation.
            // Since CDT finite faces are already guaranteed to have a consistent CCW winding
            // in 2D space, we simply maintain this consistency globally. Calculating a local
            // 3D dot product against a static face-center normal is incorrect for curved
            // surfaces (like cylinders wrapping 360 degrees) and caused winding inversions
            // that triggered non-manifold rejections.
            if (face.Orientation() == TopAbs_REVERSED) {
                std::swap(vb, vc);
            }
            if (new_mesh.add_face(va, vb, vc) != Mesh::null_face()) added++;
            else rejected++;
        }
        if (rejected > 0)
            spdlog::warn("uv_grid_retessellate [face {}]: add_face rejected {} triangles"
                         " (non-manifold)", face_idx, rejected);
        if (added == 0) {
            spdlog::warn("uv_grid_retessellate: 0 faces added");
            return false;
        }
        spdlog::debug("uv_grid_retessellate [face {}]: {} border verts, {} interior Steiner, {} triangles",
                      face_idx, N_border, interior_info.size(), added);
        mesh = std::move(new_mesh);

        // Border halfedge audit: actual count must equal the sum of border loop sizes.
        // A mismatch means the CDT left gaps (open edges that aren't part of the BREP wire).
        {
            size_t expected_he = 0;
            for (const auto& lp : loops) expected_he += lp.size();
            size_t actual_he = 0;
            for (auto h : mesh.halfedges())
                if (mesh.is_border(h)) actual_he++;
            if (actual_he != expected_he) {
                spdlog::debug("uv_grid_retessellate [face {}]: border halfedge mismatch — "
                              "expected {} got {}",
                              face_idx, expected_he, actual_he);
                // Dump UV of every border halfedge so we can locate the holes.
                auto dump_uv = mesh.property_map<Mesh::Vertex_index,
                                                  std::pair<double,double>>("v:uv").value();
                for (auto h : mesh.halfedges()) {
                    if (!mesh.is_border(h)) continue;
                    auto s = mesh.source(h), t = mesh.target(h);
                    auto [su, sv] = dump_uv[s];
                    auto [tu, tv] = dump_uv[t];
                    spdlog::debug("  bhe: ({:.4f},{:.4f})→({:.4f},{:.4f}) dU={:.4f} dV={:.4f}",
                                  su, sv, tu, tv, tu - su, tv - sv);
                }
            } else
                spdlog::debug("uv_grid_retessellate [face {}]: border halfedges OK ({})",
                              face_idx, actual_he);
        }

        return true;
    }
}

void repair_open_boundary_loops(std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                                 double min_edge_length,
                                 double collapse_area_ratio)
{
    // -- Phase 1: merged position map + per-position vertex reverse lookup ------
    std::map<std::tuple<double,double,double>, int> pid;
    std::vector<std::array<double,3>> pos;
    std::vector<std::vector<std::pair<size_t, Mesh::Vertex_index>>> id_to_verts;

    for (size_t si = 0; si < tessellation.size(); si++) {
        Mesh& mesh = tessellation[si].first;
        for (auto v : mesh.vertices()) {
            auto p = mesh.point(v);
            double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y()), z = CGAL::to_double(p.z());
            auto key = std::make_tuple(x, y, z);
            auto it = pid.find(key);
            int id;
            if (it == pid.end()) {
                id = (int)pos.size();
                pid[key] = id;
                pos.push_back({x, y, z});
                id_to_verts.emplace_back();
            } else {
                id = it->second;
            }
            id_to_verts[id].push_back({si, v});
        }
    }

    // -- Phase 2: edge incidence count + average triangle area ----------------
    std::map<std::pair<int,int>, int> edge_count;
    double total_tri_area = 0.0;
    size_t total_tri_count = 0;

    auto lookup_id = [&](const Mesh::Point& p) -> int {
        double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y()), z = CGAL::to_double(p.z());
        return pid.at(std::make_tuple(x, y, z));
    };

    for (const auto& [mesh, face] : tessellation) {
        for (auto f : mesh.faces()) {
            int ids[3]; int k = 0;
            K::Point_3 pts[3];
            for (auto v : mesh.vertices_around_face(mesh.halfedge(f))) {
                if (k < 3) { ids[k] = lookup_id(mesh.point(v)); pts[k] = mesh.point(v); }
                k++;
            }
            if (k != 3) continue;
            for (int e = 0; e < 3; e++) {
                int a = ids[e], b = ids[(e+1)%3];
                if (a > b) std::swap(a, b);
                edge_count[{a, b}]++;
            }
            auto cr = CGAL::cross_product(pts[1] - pts[0], pts[2] - pts[0]);
            double area = 0.5 * std::sqrt(CGAL::to_double(cr.squared_length()));
            total_tri_area += area;
            total_tri_count++;
        }
    }
    double avg_tri_area = (total_tri_count > 0) ? total_tri_area / total_tri_count : 1.0;
    // Cap the reference by the smallest possible triangle at the minimum edge length
    // so a large average (from big flat faces) doesn't make the collapse threshold too coarse.
    double ref_area = std::min(avg_tri_area, 0.5 * min_edge_length * min_edge_length);

    // -- Phase 3: collect open edges and trace into closed loops --------------
    std::map<int, std::vector<int>> open_adj;
    for (const auto& [e, cnt] : edge_count) {
        if (cnt != 1) continue;
        open_adj[e.first].push_back(e.second);
        open_adj[e.second].push_back(e.first);
    }

    if (open_adj.empty()) return;

    std::unordered_set<int> visited;
    std::vector<std::vector<int>> loops;

    for (const auto& [start, _] : open_adj) {
        if (visited.count(start)) continue;
        std::vector<int> loop;
        int cur = start, prev = -1;
        while (!visited.count(cur)) {
            visited.insert(cur);
            loop.push_back(cur);
            int next = -1;
            for (int nb : open_adj[cur])
                if (nb != prev && !visited.count(nb)) { next = nb; break; }
            if (next == -1) break;
            prev = cur; cur = next;
        }
        if ((int)loop.size() >= 3) loops.push_back(std::move(loop));
    }

    spdlog::info("  [loop repair] {} open boundary loop(s), avg_tri_area={:.3e} ({} tris), "
                 "ref_area={:.3e} (min(avg, 0.5*min_edge^2))",
                 loops.size(), avg_tri_area, total_tri_count, ref_area);

    // -- Phase 4: classify and fix each loop ----------------------------------
    std::unordered_set<size_t> meshes_to_rebuild;

    for (auto& loop : loops) {
        int N = (int)loop.size();

        // Polygon area via fan from centroid
        double cx = 0, cy = 0, cz = 0;
        for (int id : loop) { cx += pos[id][0]; cy += pos[id][1]; cz += pos[id][2]; }
        cx /= N; cy /= N; cz /= N;

        double loop_area = 0.0;
        for (int i = 0; i < N; i++) {
            const auto& a = pos[loop[i]], &b = pos[loop[(i+1)%N]];
            double ax=a[0]-cx, ay=a[1]-cy, az=a[2]-cz;
            double bx=b[0]-cx, by=b[1]-cy, bz=b[2]-cz;
            double xc=ay*bz-az*by, yc=az*bx-ax*bz, zc=ax*by-ay*bx;
            loop_area += 0.5 * std::sqrt(xc*xc + yc*yc + zc*zc);
        }

        // Shortest and longest edge in the loop.
        // min/max ratio characterises degeneracy within the loop itself -- no external
        // reference needed.  A loop with one very short edge (near-zero-span BREP artifact)
        // has ratio << 1; a loop representing genuine missing geometry has roughly
        // uniform edges and ratio near 1.
        double min_edge_len = std::numeric_limits<double>::max();
        double max_edge_len = 0.0;
        int snap_from = -1, snap_to = -1;
        for (int i = 0; i < N; i++) {
            int a = loop[i], b = loop[(i+1)%N];
            const auto& pa = pos[a]; const auto& pb = pos[b];
            double dx=pb[0]-pa[0], dy=pb[1]-pa[1], dz=pb[2]-pa[2];
            double len = std::sqrt(dx*dx+dy*dy+dz*dz);
            if (len < min_edge_len) { min_edge_len = len; snap_from = a; snap_to = b; }
            if (len > max_edge_len) max_edge_len = len;
        }
        spdlog::info("  [loop repair] {} vertices, area={:.3e}, min_edge={:.3e}, "
                     "max_edge={:.3e}, area/ref={:.3e}",
                     N, loop_area, min_edge_len, max_edge_len,
                     loop_area / ref_area);

        if (loop_area < ref_area * collapse_area_ratio) {
            // Collapse: snap snap_from -> snap_to across all meshes
            spdlog::info("    -> collapse: ({:.5f},{:.5f},{:.5f}) -> ({:.5f},{:.5f},{:.5f})",
                         pos[snap_from][0], pos[snap_from][1], pos[snap_from][2],
                         pos[snap_to][0],   pos[snap_to][1],   pos[snap_to][2]);
            K::Point_3 target(pos[snap_to][0], pos[snap_to][1], pos[snap_to][2]);
            int snaps = 0;
            for (auto& [si, v] : id_to_verts[snap_from]) {
                tessellation[si].first.point(v) = target;
                meshes_to_rebuild.insert(si);
                snaps++;
            }
            spdlog::info("    snapped {} vertices", snaps);
        } else {
            // Fill: fan-triangulate from loop[0], add to the segment that owns most loop verts
            spdlog::info("    -> fill: fan-triangulating {} vertices", N);
            std::map<size_t, int> seg_votes;
            for (int id : loop)
                for (auto& [si, v] : id_to_verts[id])
                    seg_votes[si]++;
            size_t best_si = 0;
            int best_votes = -1;
            for (auto& [si, votes] : seg_votes)
                if (votes > best_votes) { best_votes = votes; best_si = si; }

            Mesh& tmesh = tessellation[best_si].first;
            std::vector<Mesh::Vertex_index> hverts;
            for (int id : loop) {
                Mesh::Vertex_index hv = Mesh::null_vertex();
                for (auto& [si, v] : id_to_verts[id])
                    if (si == best_si) { hv = v; break; }
                if (hv == Mesh::null_vertex())
                    hv = tmesh.add_vertex(K::Point_3(pos[id][0], pos[id][1], pos[id][2]));
                hverts.push_back(hv);
            }
            // The loop tracing gives an arbitrary vertex ordering.  add_face requires
            // that the directed halfedges v0->v1, v1->v2, v2->v0 are all free boundary
            // halfedges in the mesh.  If the natural order fails, the reverse winding
            // will succeed because one of the two orientations must match the free side.
            int added = 0;
            for (int i = 1; i < N - 1; i++) {
                auto fh = tmesh.add_face(hverts[0], hverts[i], hverts[i+1]);
                if (fh == Mesh::null_face())
                    fh = tmesh.add_face(hverts[0], hverts[i+1], hverts[i]);
                if (fh != Mesh::null_face())
                    added++;
            }
            spdlog::info("    added {} fill triangles to seg {}", added, best_si);
            meshes_to_rebuild.insert(best_si);
        }
    }

    // -- Phase 5: rebuild affected meshes via repair_polygon_soup --------------
    for (size_t si : meshes_to_rebuild) {
        Mesh& mesh = tessellation[si].first;
        std::vector<Point> vbuf;
        std::vector<std::vector<size_t>> fbuf;
        std::unordered_map<size_t, size_t> vmap;

        for (auto v : mesh.vertices()) {
            vmap[v.idx()] = vbuf.size();
            auto p = mesh.point(v);
            vbuf.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())});
        }
        for (auto f : mesh.faces()) {
            std::vector<size_t> face;
            int k = 0;
            for (auto v : mesh.vertices_around_face(mesh.halfedge(f))) {
                face.push_back(vmap.at(v.idx())); k++;
            }
            if (k == 3) fbuf.push_back(face);
        }

        PMP::repair_polygon_soup(vbuf, fbuf, CGAL::parameters::geom_traits(PointArray_traits()));

        mesh = Mesh{};
        std::vector<Mesh::Vertex_index> new_verts;
        for (const auto& v : vbuf)
            new_verts.push_back(mesh.add_vertex(K::Point_3(v[0], v[1], v[2])));
        for (const auto& f : fbuf)
            if (f.size() == 3) mesh.add_face(new_verts[f[0]], new_verts[f[1]], new_verts[f[2]]);

        // Rebuilt mesh has no UV data; add NaN-filled "v:uv" so the UV fill pass in
        // remeshing.cpp will populate all vertices via ShapeAnalysis_Surface::ValueOfUV.
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        mesh.add_property_map<Mesh::Vertex_index, std::pair<double,double>>(
            "v:uv", {kNaN, kNaN});
    }
}