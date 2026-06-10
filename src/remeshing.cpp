// CGAL must be included before any OCCT header, because OCCT defines a
// conflicting preprocessor macro: #define Handle(Class) opencascade::handle<Class>
// CGAL has an internal Handle class that the macro corrupts.
#define CGAL_NO_PRECONDITIONS
#define CGAL_NO_ASSERTIONS
#define CGAL_NO_WARNINGS
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <precision_mesh/remeshing.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <limits>
#include <mutex>
#include <unordered_map>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepTools.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>

namespace PMP = CGAL::Polygon_mesh_processing;

void split_border_edges(std::vector<Mesh>& meshes, double max_edge_length) {
    spdlog::info("  splitting long border edges ...");

    size_t faces_before = 0;
    for (const auto& mesh: meshes) {
        faces_before += mesh.number_of_faces();
    }

    size_t border_num_before = 0;
    size_t border_num_after = 0;

    for (auto& mesh: meshes) {
        std::vector<EdgeDescriptor> border_edges;
        CGAL::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));
        border_num_before += border_edges.size();
        PMP::split_long_edges(border_edges, max_edge_length, mesh);
        border_edges.clear();
        CGAL::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));
        border_num_after += border_edges.size();
    }

    size_t faces_after = 0;
    for (const auto& mesh: meshes) {
        faces_after += mesh.number_of_faces();
    }

    spdlog::info("    border edges: {} -> {}", border_num_before, border_num_after);
    spdlog::info("    faces: {} -> {}", faces_before, faces_after);
}

void split_crease_edges(std::vector<Mesh>& meshes, double crease_angle, double max_edge_length) {
    spdlog::info("  splitting long crease edges ...");

    size_t faces_before = 0;
    for (const auto& mesh: meshes) {
        faces_before += mesh.number_of_faces();
    }

    size_t crease_num_before = 0;
    size_t crease_num_after = 0;
    std::mutex mutex;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
            size_t num_before = 0;
            size_t num_after = 0;
            for (size_t i = r.begin(); i != r.end(); ++i) {
                Mesh& mesh = meshes[i];

                auto crease_features =
                    mesh.add_property_map<Mesh::Edge_index, bool>("crease", false).first;
                PMP::detect_sharp_edges(mesh, crease_angle, crease_features);

                std::vector<EdgeDescriptor> crease_edges;
                for (const auto& edge: edges(mesh)) {
                    if (get(crease_features, edge)) {
                        crease_edges.push_back(edge);
                    }
                }
                num_before += crease_edges.size();
                PMP::split_long_edges(crease_edges, max_edge_length, mesh);
                mesh.remove_property_map(crease_features);
                crease_edges.clear();

                auto crease_features2 =
                    mesh.add_property_map<Mesh::Edge_index, bool>("crease", false).first;
                PMP::detect_sharp_edges(mesh, crease_angle, crease_features2);
                for (const auto& edge: edges(mesh)) {
                    if (get(crease_features2, edge)) {
                        crease_edges.push_back(edge);
                    }
                }
                num_after += crease_edges.size();
            }
            std::scoped_lock<std::mutex> lock(mutex);
            crease_num_before += num_before;
            crease_num_after += num_after;
        });

    size_t faces_after = 0;
    for (const auto& mesh: meshes) {
        faces_after += mesh.number_of_faces();
    }

    spdlog::info("    crease edges: {} -> {}", crease_num_before, crease_num_after);
    spdlog::info("    faces: {} -> {}", faces_before, faces_after);
}

void snap_border_midpoints_to_brep(std::vector<Mesh>& meshes,
                                   const std::vector<TopoDS_Face>& segments,
                                   double max_edge_length)
{
    if (segments.empty()) return;

    // Build a global position map: exact double-tuple → all (seg, vertex) instances.
    // Border vertices that are bit-identical across adjacent meshes collapse to the same
    // key, so every copy of a chord midpoint receives the same canonical snapped position.
    using PosKey = std::tuple<double,double,double>;
    struct Ref { size_t seg; Mesh::Vertex_index v; };
    std::map<PosKey, std::vector<Ref>> by_pos;

    for (size_t m = 0; m < meshes.size(); m++)
        for (auto v : meshes[m].vertices()) {
            if (!meshes[m].is_border(v)) continue;
            auto p = meshes[m].point(v);
            by_pos[{CGAL::to_double(p.x()),
                    CGAL::to_double(p.y()),
                    CGAL::to_double(p.z())}].push_back({m, v});
        }

    // Per-segment list of non-seam BREP edges with their curves and parameter ranges.
    struct CurveInfo { TopoDS_Edge edge; Standard_Real f, l; Handle(Geom_Curve) crv; };
    std::vector<std::vector<CurveInfo>> seg_curves(meshes.size());
    for (size_t m = 0; m < segments.size() && m < meshes.size(); m++)
        for (TopExp_Explorer ee(segments[m], TopAbs_EDGE); ee.More(); ee.Next()) {
            TopoDS_Edge e = TopoDS::Edge(ee.Current());
            if (BRep_Tool::IsClosed(e, segments[m])) continue;
            Standard_Real f, l;
            Handle(Geom_Curve) crv = BRep_Tool::Curve(e, f, l);
            if (crv.IsNull()) continue;
            if (f > l) std::swap(f, l);
            seg_curves[m].push_back({e, f, l, crv});
        }

    // Classification by a border vertex's own distance to BREP edges in ref0's segment.
    //
    // A border vertex lies, by construction, on the face's wire.  PMP::split_long_edges
    // places split vertices on the straight CHORD between an edge's endpoints, so:
    //
    //   * Straight edges (generator line, U-cut): chord == edge, every split vertex stays
    //     exactly ON the line → distance ≈ 0 ≤ edge_tol → already correct, skip.
    //
    //   * Curved edges (arcs, splines): every split vertex sits off the curve by its local
    //     sagitta → distance > edge_tol → project it back onto the curve.  This must work
    //     for ALL split vertices, not just a single midpoint: an edge several times longer
    //     than max_edge_length is bisected repeatedly, leaving many collinear chord points
    //     whose sagitta (relative to the full arc) far exceeds max_edge_length.  Hence there
    //     is NO upper distance bound and NO neighbour check — both wrongly rejected the
    //     interior chord points and left them on the chord (linear-interpolation artifact).
    const double edge_tol  = 1e-3;

    // Cache (seg, PosKey) → (nearest_edge_idx, distance).
    std::map<std::pair<size_t,PosKey>, std::pair<int,double>> dist_cache;
    auto nearest_edge_for = [&](size_t seg, const PosKey& key) -> std::pair<int,double> {
        auto cache_key = std::make_pair(seg, key);
        auto it = dist_cache.find(cache_key);
        if (it != dist_cache.end()) return it->second;
        if (seg >= seg_curves.size()) { dist_cache[cache_key]={-1,1e300}; return {-1,1e300}; }

        TopoDS_Vertex vs = BRepBuilderAPI_MakeVertex(
            gp_Pnt(std::get<0>(key), std::get<1>(key), std::get<2>(key)));
        int best = -1; double best_d = 1e300;
        for (int i = 0; i < (int)seg_curves[seg].size(); i++) {
            BRepExtrema_DistShapeShape ext;
            ext.LoadS1(seg_curves[seg][i].edge);
            ext.LoadS2(vs); ext.Perform();
            if (ext.IsDone() && ext.Value() < best_d) { best_d = ext.Value(); best = i; }
        }
        dist_cache[cache_key] = {best, best_d};
        return {best, best_d};
    };

    int total_snapped = 0;

    for (auto& [pos, refs] : by_pos) {
        if (refs.empty()) continue;
        const auto& ref0 = refs[0];

        // Find the nearest BREP edge to this border vertex.  If it already lies on that
        // edge (straight-edge split point, or an original wire vertex) leave it alone;
        // otherwise it is an off-curve chord point and must be projected onto the curve.
        auto [e_idx, e_dist] = nearest_edge_for(ref0.seg, pos);
        if (e_idx < 0 || e_dist <= edge_tol) continue;

        // Project onto the identified curve.  NearestPoint() lands on the correct arc point
        // for any chord position; for the rare corner-adjacent point whose nearest edge is
        // ambiguous the projection distance is still minimal, so this is the best snap.
        const auto& ci = seg_curves[ref0.seg][e_idx];
        GeomAPI_ProjectPointOnCurve proj(
            gp_Pnt(std::get<0>(pos), std::get<1>(pos), std::get<2>(pos)),
            ci.crv, ci.f, ci.l);
        if (proj.NbPoints() == 0) continue;

        gp_Pnt q = proj.NearestPoint();
        Mesh::Point snapped(q.X(), q.Y(), q.Z());

        // Apply the same canonical position to every copy of this chord midpoint.
        for (auto& ref : refs)
            meshes[ref.seg].point(ref.v) = snapped;
        total_snapped++;
    }

    if (total_snapped > 0)
        spdlog::info("  snapped {} border midpoints onto BREP curves (chord→arc correction)",
                     total_snapped);
}

void remesh_and_project(
    std::vector<Mesh>& meshes,
    const std::vector<TopoDS_Face>& segments,
    WireProjectorCachePtr wire_projectors,
    std::vector<std::unique_ptr<StepProjector>>& surface_projectors,
    std::vector<std::unique_ptr<StepBorderProjector>>& border_projectors,
    const RemeshParams& params,
    const std::vector<bool>& skip_mask)
{
    double max_remeshing_surface_error =
        std::min(params.max_surface_error, params.min_edge_length * 0.1);

    // Segments flagged in skip_mask (UV-grid CDT tessellated) are already final: their
    // interior is exactly on the analytic surface and their border is shared with the
    // neighbours, so they bypass remeshing, projection, and the repair passes entirely.
    assert(skip_mask.empty() || skip_mask.size() == meshes.size());
    auto skip = [&skip_mask](size_t m) { return !skip_mask.empty() && skip_mask[m]; };

    // Lock all border vertices before remeshing begins.  vertex_is_constrained_map prevents
    // border vertices from being moved by tangential relaxation.  Built once from the
    // pre-remeshing border (which includes any midpoints from split_border_edges).
    // The edge constraint map is rebuilt per-call inside the iteration loop because CGAL
    // writes to it during remeshing (sub-edges of split constrained edges are marked
    // constrained), leaving stale values that would corrupt subsequent iterations.
    if (params.is_step) {
        for (auto& mesh : meshes) {
            auto vcmap = mesh.add_property_map<Mesh::Vertex_index, bool>(
                "v:border_locked", false).first;
            for (auto v : mesh.vertices())
                if (mesh.is_border(v))
                    vcmap[v] = true;
        }
    }

    // Count total border halfedges across all meshes by scanning halfedge flags directly.
    // Does not depend on vertex back-pointer state — reliable at any point in the pipeline.
    auto count_border_halfedges = [&]() -> size_t {
        size_t n = 0;
        for (const auto& mesh : meshes)
            for (auto h : mesh.halfedges())
                if (mesh.is_border(h)) n++;
        return n;
    };

    for (int i = 0; i < params.iterations; i++) {
        spdlog::info("    iteration {}", i + 1);
        spdlog::info("      remeshing ...");

        // Per-segment border halfedge counts — populated unconditionally so the post-remesh
        // border topology restore can use them regardless of log level.
        std::vector<size_t> border_before(meshes.size(), 0);
        for (size_t m = 0; m < meshes.size(); m++)
            for (auto h : meshes[m].halfedges())
                if (meshes[m].is_border(h)) border_before[m]++;
        spdlog::debug("      [border] before remesh: {} border halfedges",
                      count_border_halfedges());

        // Save border vertex positions per segment so that any border vertices incorrectly
        // collapsed by isotropic_remeshing can be re-inserted afterward.
        std::vector<std::vector<Mesh::Point>> saved_bverts(meshes.size());
        if (params.is_step) {
            for (size_t m = 0; m < meshes.size(); m++) {
                if (skip(m)) continue;
                for (auto h : meshes[m].halfedges())
                    if (meshes[m].is_border(h))
                        saved_bverts[m].push_back(meshes[m].point(meshes[m].target(h)));
            }
        }

        // Snapshot border vertex positions before remeshing.  isotropic_remeshing honours
        // protect_constraints and the constraint maps, but we restore positions explicitly
        // as a belt-and-suspenders guarantee of exact position agreement between adjacent
        // segments on their shared boundary edges.
        std::vector<std::unordered_map<size_t, Mesh::Point>> saved_border(meshes.size());
        if (params.is_step) {
            for (size_t m = 0; m < meshes.size(); m++) {
                if (skip(m)) continue;
                const Mesh& mesh = meshes[m];
                for (auto v : mesh.vertices())
                    if (mesh.is_border(v))
                        saved_border[m][v.idx()] = mesh.point(v);
            }
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t m = r.begin(); m != r.end(); ++m) {
                    if (skip(m)) continue;
                    Mesh& mesh = meshes[m];
                    const std::pair edge_min_max{params.min_edge_length, params.max_edge_length};
                    PMP::Adaptive_sizing_field<Mesh> sizing_field(max_remeshing_surface_error,
                                                                  edge_min_max, faces(mesh), mesh);
                    try {
                        if (params.is_step) {
                            auto vcmap = lookup_property_map<Mesh::Vertex_index, bool>(
                                mesh, "v:border_locked");
                            // Rebuild the edge constraint map immediately before each call so
                            // it reflects the current border state.  CGAL writes to this map
                            // during remeshing (marking sub-edges of split constrained edges),
                            // so a map built in a previous iteration would have stale values.
                            auto ecmap = mesh.add_property_map<Mesh::Edge_index, bool>(
                                "e:border_locked", false).first;
                            for (auto e : mesh.edges())
                                ecmap[e] = mesh.is_border(e);
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .protect_constraints(true)
                                                 .edge_is_constrained_map(ecmap)
                                                 .vertex_is_constrained_map(vcmap));
                        }
                        else {
                            auto crease_map =
                                lookup_property_map<Mesh::Edge_index, bool>(mesh, "crease");
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .edge_is_constrained_map(crease_map)
                                                 .protect_constraints(true));
                        }
                    }
                    catch (const std::out_of_range&) {
                        spdlog::warn("Out of range exception when remeshing segment {}.", m);
                    }
                }});

        if (spdlog::should_log(spdlog::level::debug)) {
            size_t total_after = 0;
            for (size_t m = 0; m < meshes.size(); m++) {
                const Mesh& mesh = meshes[m];
                size_t n = 0;
                double min_len = std::numeric_limits<double>::max();
                for (auto h : mesh.halfedges()) {
                    if (!mesh.is_border(h)) continue;
                    n++;
                    auto a = mesh.point(mesh.source(h));
                    auto b = mesh.point(mesh.target(h));
                    double dx = CGAL::to_double(b.x()-a.x()),
                           dy = CGAL::to_double(b.y()-a.y()),
                           dz = CGAL::to_double(b.z()-a.z());
                    min_len = std::min(min_len, std::sqrt(dx*dx+dy*dy+dz*dz));
                }
                total_after += n;
                if (n != border_before[m])
                    spdlog::debug("      [border-change] seg {} {} -> {} halfedges, min_border_edge={:.6f}",
                                  m, border_before[m], n, min_len);
            }
            spdlog::debug("      [border] after remesh:  {} border halfedges", total_after);
        }

        // Restore border vertex positions to their pre-remeshing values.
        if (params.is_step) {
            std::atomic<size_t> total_restored{0};
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m = r.begin(); m != r.end(); ++m) {
                        if (skip(m)) continue;
                        Mesh& mesh = meshes[m];
                        size_t restored = 0;
                        for (auto v : mesh.vertices()) {
                            auto it = saved_border[m].find(v.idx());
                            if (it != saved_border[m].end()) {
                                if (mesh.point(v) != it->second) {
                                    mesh.point(v) = it->second;
                                    restored++;
                                }
                            }
                        }
                        total_restored += restored;
                    }
                });
            if (total_restored)
                spdlog::debug("      restored {} border vertex positions after remeshing",
                              total_restored.load());
        }
        spdlog::debug("      [border] after pos-restore: {} border halfedges", count_border_halfedges());

        // Restore border vertices collapsed by isotropic_remeshing despite the constraint maps.
        // For each segment whose border halfedge count dropped, scan for saved border vertex
        // positions that are no longer on the border, find the spanning border edge, and re-split
        // it.  Position-restore has already run, so the two remaining endpoints (A, C) are at
        // their exact pre-remesh positions — split_long_edges places the new vertex at their
        // midpoint, which equals the saved position for the midpoint-insertion case.  For the
        // general (non-midpoint) case we relocate the new vertex explicitly afterward.
        if (params.is_step) {
            size_t total_restored_bv = 0;
            for (size_t m = 0; m < meshes.size(); m++) {
                if (skip(m)) continue;
                Mesh& mesh = meshes[m];
                size_t cur_bhe = 0;
                for (auto h : mesh.halfedges()) if (mesh.is_border(h)) cur_bhe++;
                if (cur_bhe >= border_before[m]) continue;

                auto vcmap = lookup_property_map<Mesh::Vertex_index, bool>(mesh, "v:border_locked");
                size_t restored = 0;

                for (const auto& saved_pt : saved_bverts[m]) {
                    double px = CGAL::to_double(saved_pt.x()),
                           py = CGAL::to_double(saved_pt.y()),
                           pz = CGAL::to_double(saved_pt.z());
                    // Check whether this position still exists on the border.
                    bool found = false;
                    for (auto h : mesh.halfedges()) {
                        if (!mesh.is_border(h)) continue;
                        auto p = mesh.point(mesh.target(h));
                        if (CGAL::to_double(p.x()) == px &&
                            CGAL::to_double(p.y()) == py &&
                            CGAL::to_double(p.z()) == pz) { found = true; break; }
                    }
                    if (found) continue;

                    // Find the border halfedge whose source/target bracket saved_pt.
                    HalfEdgeDescriptor span_he = mesh.null_halfedge();
                    for (auto h : mesh.halfedges()) {
                        if (!mesh.is_border(h)) continue;
                        auto a = mesh.point(mesh.source(h));
                        auto c = mesh.point(mesh.target(h));
                        double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y()), az = CGAL::to_double(a.z());
                        double cx = CGAL::to_double(c.x()), cy = CGAL::to_double(c.y()), cz = CGAL::to_double(c.z());
                        double dx = cx-ax, dy = cy-ay, dz = cz-az;
                        double len2 = dx*dx+dy*dy+dz*dz;
                        if (len2 < 1e-20) continue;
                        double t = ((px-ax)*dx + (py-ay)*dy + (pz-az)*dz) / len2;
                        if (t < 1e-6 || t > 1.0-1e-6) continue;
                        double ex = ax+t*dx-px, ey = ay+t*dy-py, ez = az+t*dz-pz;
                        if (ex*ex+ey*ey+ez*ez > 1e-10) continue;
                        span_he = h;
                        break;
                    }
                    if (span_he == mesh.null_halfedge()) {
                        spdlog::warn("      [restore-border] seg {} no span found for ({:.6f},{:.6f},{:.6f})",
                                     m, px, py, pz);
                        continue;
                    }

                    // Split the spanning border edge. split_long_edges inserts the new vertex at
                    // the midpoint; we then relocate it to the exact saved position if needed.
                    auto a_pt = mesh.point(mesh.source(span_he));
                    auto c_pt = mesh.point(mesh.target(span_he));
                    double elen = std::sqrt(
                        CGAL::to_double((c_pt.x()-a_pt.x())*(c_pt.x()-a_pt.x()) +
                                        (c_pt.y()-a_pt.y())*(c_pt.y()-a_pt.y()) +
                                        (c_pt.z()-a_pt.z())*(c_pt.z()-a_pt.z())));
                    std::vector<EdgeDescriptor> to_split = {mesh.edge(span_he)};
                    PMP::split_long_edges(to_split, elen * 0.5, mesh);
                    // The new vertex has vcmap=false (default).  Find it and lock it.
                    for (auto v : mesh.vertices()) {
                        if (!mesh.is_border(v) || vcmap[v]) continue;
                        mesh.point(v) = saved_pt;  // relocate to exact saved position
                        vcmap[v] = true;
                        break;
                    }
                    restored++;
                }
                if (restored) {
                    spdlog::info("      [restore-border] seg {} restored {} collapsed border vertices",
                                 m, restored);
                    total_restored_bv += restored;
                }
            }
            if (total_restored_bv)
                spdlog::debug("      [border] after border-restore: {} border halfedges", count_border_halfedges());
        }

        // Repair vertices whose halfedge back-pointer was nulled by the remesher while
        // some halfedges still reference them as a target.  Use two passes so that border
        // halfedges are preferred: a border vertex assigned an interior back-pointer will
        // be misclassified by mesh.is_border(v) in the next iteration's snapshot and by
        // CGAL's internal border checks.  Vertices with null back-pointers that have NO
        // valid halfedge referencing them (truly isolated) are left in place — they are
        // harmless and safe to skip during projection.  We do NOT call collect_garbage()
        // here because the remesher may leave removed-but-not-yet-compacted halfedges that
        // still reference these vertices; triggering compaction would remap those halfedge
        // targets to the wrong vertices.
        std::atomic<size_t> total_repaired{0};
        std::atomic<size_t> total_border_repaired{0};
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t m = r.begin(); m != r.end(); ++m) {
                    if (skip(m)) continue;
                    Mesh& mesh = meshes[m];
                    size_t repaired = 0, border_repaired = 0;
                    // Pass 1: for each null-halfedge vertex, collect the best halfedge
                    // candidate, upgrading to a border halfedge whenever one is found.
                    std::unordered_map<size_t, HalfEdgeDescriptor> best;
                    for (auto h : mesh.halfedges()) {
                        auto v = mesh.target(h);
                        if (mesh.halfedge(v) != mesh.null_halfedge()) continue;
                        auto [it, inserted] = best.emplace(v.idx(), h);
                        if (!inserted && mesh.is_border(h))
                            it->second = h;  // upgrade to border halfedge
                    }
                    // Pass 2: apply back-pointers.
                    for (auto& [vidx, h] : best) {
                        mesh.set_halfedge(Mesh::Vertex_index(vidx), h);
                        repaired++;
                        if (mesh.is_border(h)) border_repaired++;
                    }
                    total_repaired += repaired;
                    total_border_repaired += border_repaired;
                }
            });
        if (total_repaired) {
            spdlog::info("      repaired {} null-halfedge vertices ({} restored as border)",
                         total_repaired.load(), total_border_repaired.load());
        }
        spdlog::debug("      [border] after null-he repair: {} border halfedges", count_border_halfedges());

        // Fill UV for vertices inserted by the remesher (they have NaN from the "v:uv" default).
        // ShapeAnalysis_Surface::ValueOfUV is used rather than BRepExtrema global search because
        // it handles periodicity correctly and is immune to the trimming-boundary false-minimum.
        if (params.is_step) {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m = r.begin(); m != r.end(); ++m) {
                        if (skip(m)) continue;
                        Mesh& mesh = meshes[m];
                        const double kNaN = std::numeric_limits<double>::quiet_NaN();
                        auto uv_map = mesh.add_property_map<Mesh::Vertex_index,
                                                             std::pair<double,double>>(
                            "v:uv", {kNaN, kNaN}).first;
                        Handle(Geom_Surface) surf = BRep_Tool::Surface(segments[m]);
                        if (surf.IsNull()) continue;
                        ShapeAnalysis_Surface sa(surf);

                        // Period-normalization for the UV seeds we're about to fill.
                        // ValueOfUV returns UV in the surface's natural range (e.g. [0,2π])
                        // which may differ from the face's actual parameter range (e.g. [π,3π]).
                        // Normalize so project_to_step's Newton search starts on the right sheet.
                        const bool u_per = surf->IsUPeriodic();
                        const bool v_per = surf->IsVPeriodic();
                        const double u_T  = u_per ? surf->UPeriod() : 0.0;
                        const double v_T  = v_per ? surf->VPeriod() : 0.0;
                        Standard_Real umin_f=0, umax_f=0, vmin_f=0, vmax_f=0;
                        if (u_per || v_per)
                            BRepTools::UVBounds(segments[m], umin_f, umax_f, vmin_f, vmax_f);

                        for (auto v : mesh.vertices()) {
                            if (mesh.halfedge(v) == mesh.null_halfedge()) continue;
                            auto& uv = uv_map[v];
                            if (!std::isnan(uv.first)) continue;
                            auto pt = mesh.point(v);
                            gp_Pnt p(CGAL::to_double(pt.x()),
                                     CGAL::to_double(pt.y()),
                                     CGAL::to_double(pt.z()));
                            gp_Pnt2d uv2 = sa.ValueOfUV(p, 1e-7);
                            double uu = uv2.X(), vv = uv2.Y();
                            if (u_per && u_T > 0.0) {
                                while (uu < umin_f) uu += u_T;
                                while (uu > umin_f + u_T) uu -= u_T;
                            }
                            if (v_per && v_T > 0.0) {
                                while (vv < vmin_f) vv += v_T;
                                while (vv > vmin_f + v_T) vv -= v_T;
                            }
                            uv = {uu, vv};
                        }
                    }
                });
        }

        if (params.is_step && !params.no_projection) {
            spdlog::info("      projecting ...");
            // Weight increases each iteration (1/N, 1/(N-1), ..., 1/1), reaching full
            // projection (weight=1) on the final iteration.  The gradual schedule prevents
            // degenerate faces from aggressively snapping vertices before the mesh topology
            // has had time to adapt.
            double weight = 1.0 / (params.iterations - i);
            std::vector<ProjectionStats> seg_stats(meshes.size());
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m = r.begin(); m != r.end(); ++m) {
                        if (skip(m)) continue;
                        project_to_step(segments[m], meshes[m], wire_projectors,
                                        *surface_projectors[m], *border_projectors[m],
                                        weight, m, &seg_stats[m],
                                        params.max_edge_length);
                    }});

            // Merge per-segment stats and print a single summary for this iteration.
            ProjectionStats global;
            for (const auto& s : seg_stats) global.merge(s);

            if (global.min_delta == std::numeric_limits<double>::max()) global.min_delta = 0.0;
            size_t n_projected = global.n_total - global.n_skipped_null;
            double mean_delta  = n_projected > 0 ? global.sum_delta / n_projected : 0.0;

            spdlog::info("      {} verts ({} border-skip, {} surface, {} occt-skip, {} null-he skipped) "
                         "weight={:.3f} delta min={:.2e} mean={:.2e} max={:.2e} "
                         ">1e-3:{} >1e-2:{} >1e-1:{}",
                         n_projected,
                         global.n_border_skipped, global.n_surface_projected,
                         global.n_occt_skip, global.n_skipped_null,
                         weight,
                         global.min_delta, mean_delta, global.max_delta,
                         global.n_above_1em3, global.n_above_1em2, global.n_above_1em1);

            for (size_t k = 0; k < global.top_records.size(); k++) {
                const auto& r = global.top_records[k];
                if (r.delta < 1e-9) break;
                spdlog::warn("        top{} delta={:.4f} seg={} face={} vert={} "
                             "border={} null_he={} on_edge={} edge_dist={:.4f} "
                             "in=({:.4f},{:.4f},{:.4f}) proj=({:.4f},{:.4f},{:.4f})",
                             k + 1, r.delta, r.segment_index, r.face_code, r.vert_idx,
                             r.is_border, r.null_he, r.on_real_edge, r.edge_dist,
                             r.px, r.py, r.pz, r.qx, r.qy, r.qz);
            }
        }

        spdlog::debug("      [border] after projection: {} border halfedges", count_border_halfedges());

        if (params.on_iteration_done)
            params.on_iteration_done(i + 1);
    }
}
