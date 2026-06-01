// CGAL must be included before any OCCT header (OCCT defines a conflicting Handle macro).
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <precision_mesh/step_projection.h>

#include <algorithm>
#include <deque>
#include <limits>
#include <unordered_set>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
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
                // physical boundary — skip them.
                if (BRep_Tool::IsClosed(TopoDS::Edge(e), single_face)) {
                    continue;
                }
                // Degenerate edges (e.g. cone apex) have zero length — skip them.
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
            spdlog::warn("Failed to associate {} border vertices with face boundary.",
                         vertex_queue.size());
        }
    }

    return border_vertex_projector_map;
}

void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr wire_projectors,
                     StepProjector& surface_projector,
                     StepBorderProjector& border_projector,
                     double weight)
{
    double w1 = std::max(0.0, std::min(1.0, weight));
    double w2 = 1.0 - w1;

    auto border_vertex_projector_map =
        get_border_vertex_projector_map(face, mesh, wire_projectors);

    for (auto v: mesh.vertices()) {
        auto input = mesh.point(v);

        Mesh::Point projected;
        if (mesh.is_border(v)) {
            auto wire_projector_it = border_vertex_projector_map.find(v);
            if (wire_projector_it == border_vertex_projector_map.end()) {
                spdlog::warn("Failed to find projector corresponding to mesh border vertex.");
                projected = border_projector(input);
            }
            else {
                projected = wire_projector_it->second(input);
            }
        }
        else {
            projected = surface_projector(input);
        }

        mesh.point(v) = Mesh::Point(
            w1 * projected.x() + w2 * input.x(),
            w1 * projected.y() + w2 * input.y(),
            w1 * projected.z() + w2 * input.z());
    }
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
                                             double tolerance, int samples_per_tri)
{
    auto dist_to = [](BRepExtrema_DistShapeShape& ext, const Mesh::Point& p) -> double {
        TopoDS_Vertex v = BRepBuilderAPI_MakeVertex(gp_Pnt(p.x(), p.y(), p.z()));
        ext.LoadS2(v);
        ext.Perform();
        return (ext.IsDone() && ext.NbSolution() > 0) ? ext.Value() : -1.0;
    };

    const size_t n = std::min(meshes.size(), segments.size());

    // Per-segment partials (one slot each → no races); mean_surface_error holds the sum.
    std::vector<TessellationValidation> partials(n);

    // Each segment is independent and builds its own BRepExtrema instances (OCCT is safe
    // per-instance), so segments validate in parallel — same pattern as remesh_and_project.
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t>& rng) {
        for (size_t m = rng.begin(); m != rng.end(); ++m) {
            const Mesh& mesh = meshes[m];
            const TopoDS_Face& face = segments[m];
            TessellationValidation& p = partials[m];
            double err_sum = 0.0;

            BRep_Builder builder;
            TopoDS_Compound edge_comp;
            builder.MakeCompound(edge_comp);
            for (TopExp_Explorer ee(face, TopAbs_EDGE); ee.More(); ee.Next())
                builder.Add(edge_comp, ee.Current());
            BRepExtrema_DistShapeShape edge_ext;
            edge_ext.LoadS1(edge_comp);
            BRepExtrema_DistShapeShape face_ext;
            face_ext.LoadS1(face);

            for (auto v : mesh.vertices()) {
                const auto& pt = mesh.point(v);
                if (mesh.is_border(v)) {
                    p.border_verts++;
                    double d = dist_to(edge_ext, pt);
                    if (d >= 0) {
                        p.max_border_edge_dist = std::max(p.max_border_edge_dist, d);
                        if (d > tolerance) p.misclassified_border++;
                    }
                } else {
                    p.interior_verts++;
                    double d = dist_to(face_ext, pt);
                    if (d >= 0) {
                        p.max_interior_face_dist = std::max(p.max_interior_face_dist, d);
                        if (d > tolerance) p.misclassified_interior++;
                    }
                }
            }

            // Surface error: sample triangle interiors (centroid, + edge midpoints) and
            // measure to the BREP face — captures chord deviation between vertices.
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

    return r;
}
