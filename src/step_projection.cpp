// CGAL must be included before any OCCT header (OCCT defines a conflicting Handle macro).
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <precision_mesh/step_projection.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <set>
#include <unordered_set>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Extrema_GenLocateExtPS.hxx>
#include <Extrema_POnSurf.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>

namespace PMP = CGAL::Polygon_mesh_processing;

StepProjector::StepProjector() {}

StepProjector::StepProjector(const TopoDS_Shape& shape) {
    extrema.LoadS1(shape);
}

void StepProjector::setShape(const TopoDS_Shape& shape) {
    extrema.LoadS1(shape);
}

Mesh::Point StepProjector::operator()(const Mesh::Point& p) {
    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2]));
    extrema.LoadS2(vertex);
    extrema.Perform();
    auto nearest = extrema.PointOnShape1(1);
    return Mesh::Point(nearest.X(), nearest.Y(), nearest.Z());
}

StepBorderProjector::StepBorderProjector() {}

StepBorderProjector::StepBorderProjector(const TopoDS_Face& face) {
    BRep_Builder builder;
    builder.MakeCompound(border);
    for (TopExp_Explorer wire_exp(face, TopAbs_WIRE); wire_exp.More(); wire_exp.Next()) {
        builder.Add(border, TopoDS::Wire(wire_exp.Current()));
    }
    extrema.LoadS1(border);
}

void StepBorderProjector::setFace(const TopoDS_Shape& face) {
    BRep_Builder builder;
    builder.MakeCompound(border);
    for (TopExp_Explorer wire_exp(face, TopAbs_WIRE); wire_exp.More(); wire_exp.Next()) {
        builder.Add(border, TopoDS::Wire(wire_exp.Current()));
    }
    extrema.LoadS1(border);
}

Mesh::Point StepBorderProjector::operator()(const Mesh::Point& p) {
    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2]));
    extrema.LoadS2(vertex);
    extrema.Perform();
    auto nearest = extrema.PointOnShape1(1);
    return Mesh::Point(nearest.X(), nearest.Y(), nearest.Z());
}

TopoDS_Wire get_border_loop_wire(
    const TopoDS_Vertex& vertex,
    const TopTools_IndexedDataMapOfShapeListOfShape& vertex_to_edge_map,
    const TopTools_IndexedDataMapOfShapeListOfShape& edge_to_face_map)
{
    size_t wire_size = 0;
    BRepBuilderAPI_MakeWire wire_maker;
    std::unordered_set<int> evaluated;

    std::deque<TopoDS_Vertex> vertex_queue = { vertex };

    while (!vertex_queue.empty()) {
        auto v = vertex_queue.front();
        vertex_queue.pop_front();

        auto edges = vertex_to_edge_map.FindFromKey(v);
        for (const auto& e: edges) {
            int edge_code = shapeHashCode(e);
            if (!evaluated.insert(edge_code).second) {
                continue;
            }
            auto faces = edge_to_face_map.FindFromKey(e);
            if (faces.Size() == 1) {
                auto single_face = TopoDS::Face(faces.First());
                // Seam edges belong to only one face but are parametric artifacts with no
                // physical boundary -- skip them.
                if (BRep_Tool::IsClosed(TopoDS::Edge(e), single_face)) {
                    continue;
                }
                // Degenerate edges (e.g. cone apex) have zero length -- skip them.
                if (BRep_Tool::Degenerated(TopoDS::Edge(e))) {
                    continue;
                }
            }

            wire_maker.Add(TopoDS::Edge(e));
            wire_size++;

            TopoDS_Vertex v1, v2;
            TopExp::Vertices(TopoDS::Edge(e), v1, v2);
            vertex_queue.push_back(v1);
            vertex_queue.push_back(v2);
        }
    }

    if (wire_size == 0) {
        return TopoDS_Wire();
    }

    spdlog::debug("making wire of size: {}", wire_size);
    return wire_maker.Wire();
}

WireProjectorCachePtr get_edge_vertex_wire_projectors(const TopoDS_Shape& shape) {
    auto wire_projectors =
        std::make_shared<std::unordered_map<int, std::unordered_map<int, StepProjector>>>();

    TopTools_IndexedDataMapOfShapeListOfShape edge_to_face_map;
    TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_to_face_map);

    std::vector<TopoDS_Face> faces;
    for (TopExp_Explorer face_exp(shape, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
        auto face = TopoDS::Face(face_exp.Current());
        int face_code = shapeHashCode(face);
        faces.push_back(face);
        (*wire_projectors)[face_code] = {};
    }

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, faces.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t f = r.begin(); f != r.end(); ++f) {
                auto& face = faces[f];
                int face_code = shapeHashCode(face);

                TopTools_IndexedDataMapOfShapeListOfShape vertex_to_edge_map;
                TopExp::MapShapesAndUniqueAncestors(face, TopAbs_VERTEX, TopAbs_EDGE,
                                                    vertex_to_edge_map);

                auto face_wire_projectors_it = wire_projectors->find(face_code);
                if (face_wire_projectors_it == wire_projectors->end()) {
                    continue;
                }
                auto& face_wire_projectors = face_wire_projectors_it->second;

                for (TopExp_Explorer vertex_exp(face, TopAbs_VERTEX); vertex_exp.More();
                     vertex_exp.Next())
                {
                    auto v = TopoDS::Vertex(vertex_exp.Current());
                    int vertex_code = shapeHashCode(v);

                    if (face_wire_projectors.count(vertex_code)) {
                        continue;
                    }

                    auto wire = get_border_loop_wire(v, vertex_to_edge_map, edge_to_face_map);
                    if (wire.IsNull()) {
                        continue;
                    }

                    StepProjector projector(wire);

                    for (TopExp_Explorer wire_exp(wire, TopAbs_VERTEX); wire_exp.More();
                         wire_exp.Next())
                    {
                        auto wire_v = TopoDS::Vertex(wire_exp.Current());
                        face_wire_projectors[shapeHashCode(wire_v)] = projector;
                    }
                }
            }});

    spdlog::debug("done making wire projectors.");
    return wire_projectors;
}

std::unordered_map<Mesh::Vertex_index, StepProjector>
get_border_vertex_projector_map(const TopoDS_Face& face, Mesh& mesh,
                                WireProjectorCachePtr wire_projectors, double tolerance)
{
    int face_code = shapeHashCode(face);
    const auto& face_wire_projectors = (*wire_projectors)[face_code];

    std::unordered_map<Mesh::Vertex_index, StepProjector> border_vertex_projector_map;

    std::vector<TopoDS_Edge> edges;
    std::vector<BRepExtrema_DistShapeShape> edge_extremas;
    for (TopExp_Explorer edge_exp(face, TopAbs_EDGE); edge_exp.More(); edge_exp.Next()) {
        edge_extremas.push_back(BRepExtrema_DistShapeShape());
        edge_extremas.back().LoadS1(TopoDS::Edge(edge_exp.Current()));
        edges.push_back(TopoDS::Edge(edge_exp.Current()));
    }

    // First pass: assign projectors to border vertices already on a face edge.
    std::deque<Mesh::Vertex_index> vertex_queue;
    for (auto v: mesh.vertices()) {
        if (!mesh.is_border(v)) {
            continue;
        }
        // Null-halfedge vertices are isolated (repaired or removed upstream); skip them
        // here so they are not queued as unresolvable border vertices.
        if (mesh.halfedge(v) == mesh.null_halfedge()) {
            continue;
        }

        auto p = mesh.point(v);
        auto point = gp_Pnt(p[0], p[1], p[2]);
        TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(point);

        size_t nearest_edge_index = 0;
        double nearest_edge_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < edge_extremas.size(); i++) {
            edge_extremas[i].LoadS2(vertex);
            edge_extremas[i].Perform();
            double dist = edge_extremas[i].Value();
            if (dist < nearest_edge_dist) {
                nearest_edge_index = i;
                nearest_edge_dist = dist;
            }
        }

        if (nearest_edge_dist <= tolerance) {
            const auto& nearest_edge = edges[nearest_edge_index];
            TopoDS_Vertex v1, v2;
            TopExp::Vertices(TopoDS::Edge(nearest_edge), v1, v2);
            double d1 = point.Distance(BRep_Tool::Pnt(v1));
            double d2 = point.Distance(BRep_Tool::Pnt(v2));
            TopoDS_Vertex nearest_vertex = (d2 < d1) ? v2 : v1;

            int vertex_code = shapeHashCode(nearest_vertex);
            auto projector_it = face_wire_projectors.find(vertex_code);
            if (projector_it != face_wire_projectors.end()) {
                border_vertex_projector_map[v] = projector_it->second;
            }
        }
        else {
            vertex_queue.push_back(v);
        }
    }

    if (vertex_queue.empty()) {
        return border_vertex_projector_map;
    }

    spdlog::debug("border map size: {}", border_vertex_projector_map.size());
    spdlog::debug("vertex queue size: {}", vertex_queue.size());
    spdlog::debug("num edges: {}", edges.size());

    // Second pass: propagate assignments along the border to unresolved vertices.
    bool made_progress = true;
    while (made_progress) {
        size_t start_size = vertex_queue.size();
        std::deque<Mesh::Vertex_index> remaining;
        while (!vertex_queue.empty()) {
            auto v = vertex_queue.front();
            vertex_queue.pop_front();

            bool found_neighbor = false;
            for (auto halfedge : CGAL::halfedges_around_target(mesh.halfedge(v), mesh)) {
                if (!mesh.is_border(halfedge)) {
                    continue;
                }
                auto neighbor_v = mesh.source(halfedge);
                auto v_it = border_vertex_projector_map.find(neighbor_v);
                if (v_it != border_vertex_projector_map.end()) {
                    border_vertex_projector_map[v] = v_it->second;
                    found_neighbor = true;
                    made_progress = true;
                    break;
                }
            }

            if (!found_neighbor) {
                remaining.push_back(v);
            }
        }

        made_progress = remaining.size() < start_size;
        std::swap(remaining, vertex_queue);

        if (!made_progress && !vertex_queue.empty()) {
            spdlog::warn("Failed to associate {} border vertices with face boundary. face_code={}",
                         vertex_queue.size(), face_code);
            for (auto unresolved_v : vertex_queue) {
                auto p = mesh.point(unresolved_v);
                double px = CGAL::to_double(p.x());
                double py = CGAL::to_double(p.y());
                double pz = CGAL::to_double(p.z());
                TopoDS_Vertex vshape = BRepBuilderAPI_MakeVertex(gp_Pnt(px, py, pz));

                double min_dist = std::numeric_limits<double>::max();
                size_t nearest_idx = 0;
                for (size_t i = 0; i < edge_extremas.size(); i++) {
                    edge_extremas[i].LoadS2(vshape);
                    edge_extremas[i].Perform();
                    if (edge_extremas[i].IsDone() && edge_extremas[i].Value() < min_dist) {
                        min_dist = edge_extremas[i].Value();
                        nearest_idx = i;
                    }
                }
                bool nearest_is_seam = !edges.empty() &&
                    BRep_Tool::IsClosed(edges[nearest_idx], face);

                int border_he = 0, total_he = 0;
                for (auto h : CGAL::halfedges_around_target(
                        mesh.halfedge(unresolved_v), mesh)) {
                    total_he++;
                    if (mesh.is_border(h)) border_he++;
                }

                spdlog::warn("  unresolved vertex ({:.6f},{:.6f},{:.6f}): "
                             "nearest_edge_dist={:.6f} edge={} ({}), "
                             "border_halfedges={}/{}",
                             px, py, pz,
                             min_dist, nearest_idx,
                             nearest_is_seam ? "seam" : "real",
                             border_he, total_he);
            }
        }
    }

    return border_vertex_projector_map;
}

void ProjectionStats::merge(const ProjectionStats& other)
{
    n_total             += other.n_total;
    n_skipped_null      += other.n_skipped_null;
    n_border_skipped    += other.n_border_skipped;
    n_surface_projected += other.n_surface_projected;
    sum_delta           += other.sum_delta;
    min_delta            = std::min(min_delta, other.min_delta);
    max_delta            = std::max(max_delta, other.max_delta);
    n_above_1em3        += other.n_above_1em3;
    n_above_1em2        += other.n_above_1em2;
    n_above_1em1        += other.n_above_1em1;
    n_occt_skip         += other.n_occt_skip;
    for (const auto& r : other.top_records) {
        auto pos = std::lower_bound(top_records.begin(), top_records.end(), r,
            [](const Record& a, const Record& b) { return a.delta > b.delta; });
        top_records.insert(pos, r);
        if (top_records.size() > TOP_N) top_records.pop_back();
    }
}

void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr wire_projectors,
                     StepProjector& surface_projector,
                     StepBorderProjector& border_projector,
                     double weight,
                     size_t segment_index,
                     ProjectionStats* stats_out,
                     double max_projection_delta)
{
    double w1 = std::max(0.0, std::min(1.0, weight));
    double w2 = 1.0 - w1;

    int face_code = shapeHashCode(face);

    auto border_vertex_projector_map =
        get_border_vertex_projector_map(face, mesh, wire_projectors);

    // Build a compound of real (non-seam) BREP edges for geometry-based vertex
    // classification — same pattern as validate_tessellation.  mesh.is_border(v)
    // conflates two classes: true BREP boundary / subdivision cut (→ wire) and
    // periodic seam (interior to the surface → surface).  Seam vertices routed
    // to the wire fallback drift to the wrong edge over iterations.
    // Subdivision cut edges are intentional boundaries and are deliberately included.
    BRep_Builder ec_builder;
    TopoDS_Compound edge_comp;
    ec_builder.MakeCompound(edge_comp);
    bool has_real_edges = false;
    for (TopExp_Explorer ee(face, TopAbs_EDGE); ee.More(); ee.Next()) {
        const TopoDS_Edge& e = TopoDS::Edge(ee.Current());
        if (!BRep_Tool::IsClosed(e, face)) {
            ec_builder.Add(edge_comp, e);
            has_real_edges = true;
        }
    }
    BRepExtrema_DistShapeShape edge_classifier;
    if (has_real_edges) edge_classifier.LoadS1(edge_comp);
    const double edge_tol = 1e-3;  // matches get_border_vertex_projector_map default

    // UV-seeded local surface projection.  Constructed once per face; Perform() is called
    // per interior vertex using the stored "v:uv" as initial guess.  This eliminates the
    // wrong-local-minimum failure of the global BRepExtrema search on complex B-spline faces.
    // Falls back to surface_projector (global) if UV is unavailable or Newton diverges.
    const double kNaN_uv = std::numeric_limits<double>::quiet_NaN();
    auto uv_map = mesh.add_property_map<Mesh::Vertex_index, std::pair<double,double>>(
        "v:uv", {kNaN_uv, kNaN_uv}).first;
    BRepAdaptor_Surface adaptor(face);
    Extrema_GenLocateExtPS local_ext(adaptor, 1e-7, 1e-7);

    ProjectionStats stats;

    for (auto v: mesh.vertices()) {
        stats.n_total++;

        // Vertices repaired upstream (remesh_and_project restores null halfedge back-
        // pointers after each remesh iteration).  Truly isolated ones that could not be
        // repaired (no valid halfedge references them) are skipped here: they have no
        // incident faces so projection has no visible effect, and routing them through the
        // wire projector path would be wrong.
        if (mesh.halfedge(v) == mesh.null_halfedge()) {
            stats.n_skipped_null++;
            continue;
        }

        auto input = mesh.point(v);
        double px = CGAL::to_double(input.x()), py = CGAL::to_double(input.y()),
               pz = CGAL::to_double(input.z());

        // Use mesh.is_border as a cheap pre-filter (interior vertices are never on a
        // real boundary), then confirm geometrically so seam vertices are classified
        // as surface rather than edge.
        bool on_real_edge = false;
        double edge_check_dist = -1.0;
        if (has_real_edges && mesh.is_border(v)) {
            TopoDS_Vertex vshape = BRepBuilderAPI_MakeVertex(
                gp_Pnt(input.x(), input.y(), input.z()));
            edge_classifier.LoadS2(vshape);
            edge_classifier.Perform();
            if (edge_classifier.IsDone()) {
                edge_check_dist = edge_classifier.Value();
                on_real_edge = edge_check_dist <= edge_tol;
            }
        }

        // Border vertices are static: they come from BRepMesh (exact BREP placement) or
        // split_border_edges (chord midpoint, close enough).  Projecting them independently
        // per segment causes adjacent segments to diverge by small floating-point amounts,
        // breaking the soup-based watertightness check.  Skip them entirely.
        if (on_real_edge) {
            stats.n_border_skipped++;
            continue;
        }

        Mesh::Point projected;
        bool used_local_uv = false;
        {
            // Try UV-seeded local Newton search first (avoids wrong-local-minimum on B-splines).
            const auto& stored_uv = uv_map[v];
            if (!std::isnan(stored_uv.first)) {
                local_ext.Perform(gp_Pnt(px, py, pz), stored_uv.first, stored_uv.second);
                if (local_ext.IsDone()) {
                    const gp_Pnt& lp = local_ext.Point().Value();
                    projected = Mesh::Point(lp.X(), lp.Y(), lp.Z());
                    double ru, rv;
                    local_ext.Point().Parameter(ru, rv);
                    uv_map[v] = {ru, rv};
                    used_local_uv = true;
                }
            }
            if (!used_local_uv) {
                projected = surface_projector(input);
            }
        }

        double qx = CGAL::to_double(projected.x()), qy = CGAL::to_double(projected.y()),
               qz = CGAL::to_double(projected.z());
        double delta = gp_Pnt(px, py, pz).Distance(gp_Pnt(qx, qy, qz));

        // Guard against OCCT extrema converging to the wrong local minimum: only applies to
        // the global-search fallback path (used_local_uv=false).  With stored UV, local Newton
        // covers the vast majority of interior vertices and never triggers this false-minimum.
        if (!used_local_uv &&
            max_projection_delta > 0.0 &&
            delta > max_projection_delta &&
            surface_projector.extrema.IsDone() &&
            surface_projector.extrema.NbSolution() > 0 &&
            surface_projector.extrema.SupportTypeShape1(1) != BRepExtrema_IsInFace) {
            stats.n_occt_skip++;
            spdlog::debug("      occt-guard skip (global fallback): seg={} face={} vert={} delta={:.4f}",
                          segment_index, face_code, v.idx(), delta);
            continue;
        }

        stats.n_surface_projected++;

        stats.sum_delta += delta;
        if (delta < stats.min_delta) stats.min_delta = delta;
        if (delta > stats.max_delta) stats.max_delta = delta;
        if (delta > 1e-3) stats.n_above_1em3++;
        if (delta > 1e-2) stats.n_above_1em2++;
        if (delta > 1e-1) stats.n_above_1em1++;

        if (stats_out &&
            (stats.top_records.size() < ProjectionStats::TOP_N ||
             delta > stats.top_records.back().delta)) {
            ProjectionStats::Record rec{delta, px, py, pz, qx, qy, qz,
                                        mesh.is_border(v),
                                        mesh.halfedge(v) == mesh.null_halfedge(),
                                        on_real_edge, edge_check_dist,
                                        v.idx(), segment_index, face_code};
            auto pos = std::lower_bound(stats.top_records.begin(), stats.top_records.end(),
                rec, [](const ProjectionStats::Record& a, const ProjectionStats::Record& b) {
                    return a.delta > b.delta;
                });
            stats.top_records.insert(pos, rec);
            if (stats.top_records.size() > ProjectionStats::TOP_N)
                stats.top_records.pop_back();
        }

        mesh.point(v) = Mesh::Point(
            w1 * projected.x() + w2 * input.x(),
            w1 * projected.y() + w2 * input.y(),
            w1 * projected.z() + w2 * input.z());
    }

    if (stats_out) *stats_out = std::move(stats);
}

double get_distance_to_face(const TopoDS_Face& face, double x, double y, double z) {
    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(x, y, z));
    BRepExtrema_DistShapeShape dist(vertex, face);
    dist.Perform();
    if (dist.IsDone() && dist.Value() >= 0) {
        return dist.Value();
    }
    spdlog::error("Failed to calculate distance.");
    return -1;
}

TessellationValidation validate_tessellation(const std::vector<Mesh>& meshes,
                                             const std::vector<TopoDS_Face>& segments,
                                             const std::vector<TopoDS_Face>& edge_faces,
                                             double tolerance, int samples_per_tri)
{
    auto dist_to = [](BRepExtrema_DistShapeShape& ext, const Mesh::Point& p) -> double {
        TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(gp_Pnt(p.x(), p.y(), p.z()));
        ext.LoadS2(v);
        ext.Perform();
        return (ext.IsDone() && ext.NbSolution() > 0) ? ext.Value() : -1.0;
    };

    const size_t n = std::min(meshes.size(), segments.size());

    // Per-segment partials (one slot each -> no races); mean_surface_error holds the sum.
    std::vector<TessellationValidation> partials(n);

    // Each segment is independent and builds its own BRepExtrema instances (OCCT is safe
    // per-instance), so segments validate in parallel -- same pattern as remesh_and_project.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t>& rng) {
        for (size_t m = rng.begin(); m != rng.end(); ++m) {
            const Mesh& mesh = meshes[m];
            const TopoDS_Face& face = segments[m];
            TessellationValidation& p = partials[m];
            double err_sum = 0.0;

            // Validate against the ORIGINAL face's REAL edges, excluding the periodic seam
            // (a seam is interior to the continuous surface, not a part boundary) -- and NOT
            // against the segment's mesh border, which also includes subdivision cuts.  When
            // there was no subdivision, edge_faces[m] == segments[m].
            const TopoDS_Face& edge_face = (m < edge_faces.size()) ? edge_faces[m] : face;
            BRep_Builder builder;
            TopoDS_Compound edge_comp;
            builder.MakeCompound(edge_comp);
            for (TopExp_Explorer ee(edge_face, TopAbs_EDGE); ee.More(); ee.Next()) {
                const TopoDS_Edge& e = TopoDS::Edge(ee.Current());
                if (BRep_Tool::IsClosed(e, edge_face)) continue;   // seam edge -- not a boundary
                builder.Add(edge_comp, e);
            }
            BRepExtrema_DistShapeShape edge_ext;
            edge_ext.LoadS1(edge_comp);
            BRepExtrema_DistShapeShape face_ext;
            face_ext.LoadS1(face);

            for (auto v : mesh.vertices()) {
                const auto& pt = mesh.point(v);
                // Classify by proximity to a REAL original edge, not by mesh.is_border:
                // a vertex on a real edge is a boundary vertex; everything else (true
                // interior, seam, subdivision cut) must lie on the surface.
                double de = dist_to(edge_ext, pt);
                if (de >= 0 && de <= tolerance) {
                    p.border_verts++;
                    p.max_border_edge_dist = std::max(p.max_border_edge_dist, de);
                } else {
                    p.interior_verts++;
                    double d = dist_to(face_ext, pt);
                    if (d >= 0) {
                        p.max_interior_face_dist = std::max(p.max_interior_face_dist, d);
                        if (d > tolerance) {
                            p.misclassified_interior++;
                            spdlog::debug("    [offsurface] seg {} vert ({:.4f},{:.4f},{:.4f}) is {:.6f} off the original surface",
                                          m, CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()), d);
                        }
                    }
                }
            }

            // Surface error: sample triangle interiors (centroid, + edge midpoints) and
            // measure to the BREP face -- captures chord deviation between vertices.
            for (auto f : mesh.faces()) {
                p.total_tris++;
                auto h0 = mesh.halfedge(f);
                auto h1 = mesh.next(h0);
                auto h2 = mesh.next(h1);
                const auto& a = mesh.point(mesh.source(h0));
                const auto& b = mesh.point(mesh.source(h1));
                const auto& c = mesh.point(mesh.source(h2));

                std::vector<Mesh::Point> samples;
                samples.emplace_back((a.x()+b.x()+c.x())/3.0, (a.y()+b.y()+c.y())/3.0, (a.z()+b.z()+c.z())/3.0);
                if (samples_per_tri >= 4) {
                    samples.emplace_back((a.x()+b.x())/2.0, (a.y()+b.y())/2.0, (a.z()+b.z())/2.0);
                    samples.emplace_back((b.x()+c.x())/2.0, (b.y()+c.y())/2.0, (b.z()+c.z())/2.0);
                    samples.emplace_back((c.x()+a.x())/2.0, (c.y()+a.y())/2.0, (c.z()+a.z())/2.0);
                }
                for (const auto& s : samples) {
                    double d = dist_to(face_ext, s);
                    if (d >= 0) {
                        p.max_surface_error = std::max(p.max_surface_error, d);
                        err_sum += d;
                        p.surface_samples++;
                    }
                }
            }
            p.mean_surface_error = err_sum;  // sum for now; combined below
        }
    });

    // Reduce per-segment partials.
    TessellationValidation r;
    r.segments = n;
    double err_sum = 0.0;
    for (const auto& p : partials) {
        r.total_tris            += p.total_tris;
        r.border_verts          += p.border_verts;
        r.interior_verts        += p.interior_verts;
        r.misclassified_border  += p.misclassified_border;
        r.misclassified_interior+= p.misclassified_interior;
        r.surface_samples       += p.surface_samples;
        r.max_border_edge_dist   = std::max(r.max_border_edge_dist, p.max_border_edge_dist);
        r.max_interior_face_dist = std::max(r.max_interior_face_dist, p.max_interior_face_dist);
        r.max_surface_error      = std::max(r.max_surface_error, p.max_surface_error);
        err_sum += p.mean_surface_error;  // holds this segment's sum
    }
    r.mean_surface_error = r.surface_samples ? err_sum / r.surface_samples : 0.0;

    // Watertightness, measured orientation-independently from the triangle soup (vertices
    // merged by exact position).  We deliberately do NOT use merge_meshes() / is_border here:
    // its incremental add_face() silently drops faces on orientation/non-manifold conflicts,
    // fabricating border edges that aren't real holes.  An undirected edge incident to exactly
    // one triangle is a real open boundary (crack); incident to >2 is non-manifold.
    std::map<std::tuple<double,double,double>, int> pid;
    std::vector<std::array<double,3>> pos;
    auto id_of = [&](const Mesh::Point& p) {
        double x=CGAL::to_double(p.x()), y=CGAL::to_double(p.y()), z=CGAL::to_double(p.z());
        auto key = std::make_tuple(x,y,z);
        auto it = pid.find(key);
        if (it != pid.end()) return it->second;
        int id = (int)pos.size(); pid[key]=id; pos.push_back({x,y,z}); return id;
    };
    std::map<std::pair<int,int>, int> edge_count;
    for (const auto& mesh : meshes)
        for (auto f : mesh.faces()) {
            int ids[3]; int k = 0;
            for (auto v : mesh.vertices_around_face(mesh.halfedge(f))) { if (k < 3) ids[k] = id_of(mesh.point(v)); k++; }
            if (k != 3) continue;
            for (int e = 0; e < 3; e++) { int a = ids[e], b = ids[(e+1)%3]; if (a > b) std::swap(a,b); edge_count[{a,b}]++; }
        }

    auto seg_vertex = [&](const Mesh& mesh, double X, double Y, double Z) {
        for (auto v : mesh.vertices()) {
            auto pt = mesh.point(v);
            double dx=CGAL::to_double(pt.x())-X, dy=CGAL::to_double(pt.y())-Y, dz=CGAL::to_double(pt.z())-Z;
            if (dx*dx+dy*dy+dz*dz <= tolerance*tolerance) return v;
        }
        return Mesh::null_vertex();
    };
    int non_manifold_edges = 0;
    for (const auto& kv : edge_count) {
        if (kv.second > 2) { non_manifold_edges++; continue; }
        if (kv.second != 1) continue;
        r.open_boundary_edges++;
        const auto& A = pos[kv.first.first]; const auto& B = pos[kv.first.second];
        double len = std::sqrt((B[0]-A[0])*(B[0]-A[0])+(B[1]-A[1])*(B[1]-A[1])+(B[2]-A[2])*(B[2]-A[2]));
        std::string owners;
        for (size_t m = 0; m < n; m++) {
            auto v0 = seg_vertex(meshes[m], A[0],A[1],A[2]);
            auto v1 = seg_vertex(meshes[m], B[0],B[1],B[2]);
            bool h0 = (v0 != Mesh::null_vertex()), h1 = (v1 != Mesh::null_vertex());
            if (h0 && h1)        owners += " " + std::to_string(m) + "(both," + std::to_string(meshes[m].number_of_faces()) + "t)";
            else if (h0 || h1)   owners += " " + std::to_string(m) + (h0 ? "(p0)" : "(p1)");
        }
        spdlog::debug("    [openedge] ({:.4f},{:.4f},{:.4f}) -> ({:.4f},{:.4f},{:.4f}) len={:.5f} owners:{}",
                      A[0],A[1],A[2], B[0],B[1],B[2], len, owners);
    }
    if (non_manifold_edges)
        spdlog::warn("    non-manifold edges (incident>2): {}", non_manifold_edges);

    // Exact local-topology dump for sliver/needle open edges (len < tolerance).  For each such
    // corner, list EVERY triangle (across all segment meshes) that touches it, with its TRUE
    // owning segment id (read from the mesh, not position-matched) and each vertex's distance to
    // its segment's nearest BREP corner (d2c).  Resolves whether near-duplicate corner reps sit
    // within one face (a within-face needle) or split across faces.
    std::vector<std::array<double,3>> sliver_corners;
    for (const auto& kv : edge_count) {
        if (kv.second != 1) continue;
        const auto& A = pos[kv.first.first]; const auto& B = pos[kv.first.second];
        double dx=B[0]-A[0], dy=B[1]-A[1], dz=B[2]-A[2];
        if (dx*dx+dy*dy+dz*dz >= tolerance*tolerance) continue;     // not a needle
        std::array<double,3> P = {(A[0]+B[0])*0.5, (A[1]+B[1])*0.5, (A[2]+B[2])*0.5};
        bool dup = false;
        for (const auto& q : sliver_corners) {
            double qx=q[0]-P[0], qy=q[1]-P[1], qz=q[2]-P[2];
            if (qx*qx+qy*qy+qz*qz <= tolerance*tolerance) { dup = true; break; }
        }
        if (!dup) sliver_corners.push_back(P);
    }
    for (const auto& P : sliver_corners) {
        spdlog::info("  [slivercheck] corner ({:.7f},{:.7f},{:.7f}) -- triangles touching it:",
                     P[0], P[1], P[2]);
        std::set<size_t> sliver_segs;
        for (size_t m = 0; m < n; m++) {
            const Mesh& mesh = meshes[m];
            std::vector<std::array<double,3>> corners;
            const TopoDS_Face& ef = (m < edge_faces.size()) ? edge_faces[m] : segments[m];
            for (TopExp_Explorer ve(ef, TopAbs_VERTEX); ve.More(); ve.Next()) {
                gp_Pnt cp = BRep_Tool::Pnt(TopoDS::Vertex(ve.Current()));
                corners.push_back({cp.X(), cp.Y(), cp.Z()});
            }
            auto d2c = [&](const std::array<double,3>& v) -> double {
                double best = -1.0;
                for (const auto& c : corners) {
                    double d = std::sqrt((c[0]-v[0])*(c[0]-v[0])+(c[1]-v[1])*(c[1]-v[1])+(c[2]-v[2])*(c[2]-v[2]));
                    if (best < 0 || d < best) best = d;
                }
                return best;
            };
            for (auto f : mesh.faces()) {
                std::array<std::array<double,3>,3> vc;
                int k = 0; bool touches = false;
                for (auto v : mesh.vertices_around_face(mesh.halfedge(f))) {
                    if (k >= 3) { k++; break; }
                    auto pt = mesh.point(v);
                    vc[k] = {CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z())};
                    double tx=vc[k][0]-P[0], ty=vc[k][1]-P[1], tz=vc[k][2]-P[2];
                    if (tx*tx+ty*ty+tz*tz <= tolerance*tolerance) touches = true;
                    k++;
                }
                if (k != 3 || !touches) continue;
                sliver_segs.insert(m);
                spdlog::info("    seg {} tri: ({:.7f},{:.7f},{:.7f})[d2c={:.2e}] "
                             "({:.7f},{:.7f},{:.7f})[d2c={:.2e}] ({:.7f},{:.7f},{:.7f})[d2c={:.2e}]",
                             m, vc[0][0],vc[0][1],vc[0][2], d2c(vc[0]),
                             vc[1][0],vc[1][1],vc[1][2], d2c(vc[1]),
                             vc[2][0],vc[2][1],vc[2][2], d2c(vc[2]));
            }
        }

        // BREP topology dump for the faces at this sliver corner.
        // Build a vertex identity table: assign each unique BREP vertex object (by IsSame)
        // a local id, then report which segments share each object and its 3D position.
        // This directly reveals whether positions A and B at the hole are the same BREP
        // vertex (pcurve evaluation issue) or distinct vertex objects (BREP defect).
        std::vector<TopoDS_Vertex> uniq;
        auto vid = [&](const TopoDS_Vertex& v) -> int {
            for (int i = 0; i < (int)uniq.size(); i++)
                if (uniq[i].IsSame(v)) return i;
            uniq.push_back(v); return (int)uniq.size() - 1;
        };
        std::map<int, std::vector<size_t>> vid_segs;
        for (size_t m : sliver_segs) {
            for (TopExp_Explorer ve(segments[m], TopAbs_VERTEX); ve.More(); ve.Next()) {
                int id = vid(TopoDS::Vertex(ve.Current()));
                if (vid_segs[id].empty() || vid_segs[id].back() != m)
                    vid_segs[id].push_back(m);
            }
        }
        spdlog::info("  [brep-topology] vertex identity across sliver segs:");
        for (auto& [id, segs] : vid_segs) {
            gp_Pnt p = BRep_Tool::Pnt(uniq[id]);
            std::string ss; for (auto s : segs) ss += " " + std::to_string(s);
            spdlog::info("    V{}: ({:.7f},{:.7f},{:.7f}) in segs{}", id, p.X(), p.Y(), p.Z(), ss);
        }

        spdlog::info("  [brep-topology] wire edges per sliver seg:");
        for (size_t m : sliver_segs) {
            const TopoDS_Face& face = segments[m];
            spdlog::info("    seg {}:", m);
            for (TopExp_Explorer we(face, TopAbs_WIRE); we.More(); we.Next()) {
                for (BRepTools_WireExplorer ee(TopoDS::Wire(we.Current()), face); ee.More(); ee.Next()) {
                    const TopoDS_Edge& e = ee.Current();
                    TopoDS_Vertex ev1 = ee.CurrentVertex();
                    TopoDS_Vertex ev2 = TopExp::LastVertex(e, Standard_True);
                    Standard_Real f = 0, l = 0;
                    bool has_curve = !BRep_Tool::Curve(e, f, l).IsNull();
                    bool degen = BRep_Tool::Degenerated(e);
                    std::string ctype = "null";
                    if (has_curve) {
                        BRepAdaptor_Curve ac(e);
                        switch (ac.GetType()) {
                            case GeomAbs_Line:         ctype = "line";    break;
                            case GeomAbs_Circle:       ctype = "circle";  break;
                            case GeomAbs_BSplineCurve: ctype = "bspline"; break;
                            default:                   ctype = "other";   break;
                        }
                    }
                    spdlog::info("      {} V{}->V{} {} degen={} [{:.17g},{:.17g}] span={:.3e}",
                                 ee.Orientation() == TopAbs_FORWARD ? "FWD" : "REV",
                                 vid(ev1), vid(ev2), ctype, degen, f, l, l - f);
                }
            }
        }
    }

    return r;
}
