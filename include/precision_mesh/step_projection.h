#pragma once

#include <algorithm>
#include <deque>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <CGAL/boost/graph/iterator.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>

#include <precision_mesh/mesh_util.h>
#include <precision_mesh/step_types.h>

// Projects mesh vertices onto the nearest point of a STEP shape.
template<class Mesh>
struct StepProjector {

    StepProjector() {}
    StepProjector(const TopoDS_Shape& shape) {
        extrema.LoadS1(shape);
    }

    void setShape(const TopoDS_Shape& shape) {
        extrema.LoadS1(shape);
    }

    typename Mesh::Point operator()(const typename Mesh::Point& p) {
        TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2]));
        extrema.LoadS2(vertex);
        extrema.Perform();
        auto nearest = extrema.PointOnShape1(1);
        return typename Mesh::Point(nearest.X(), nearest.Y(), nearest.Z());
    }

    BRepExtrema_DistShapeShape extrema;
};

// Projects mesh border vertices onto the wire boundary of a STEP face.
template<class Mesh>
struct StepBorderProjector {
    StepBorderProjector() {}
    StepBorderProjector(const TopoDS_Face& face) {
        BRep_Builder builder;
        builder.MakeCompound(border);
        for (TopExp_Explorer wire_exp(face, TopAbs_WIRE); wire_exp.More(); wire_exp.Next()) {
            builder.Add(border, TopoDS::Wire(wire_exp.Current()));
        }
        extrema.LoadS1(border);
    }

    void setFace(const TopoDS_Shape& face) {
        BRep_Builder builder;
        builder.MakeCompound(border);
        for (TopExp_Explorer wire_exp(face, TopAbs_WIRE); wire_exp.More(); wire_exp.Next()) {
            builder.Add(border, TopoDS::Wire(wire_exp.Current()));
        }
        extrema.LoadS1(border);
    }

    typename Mesh::Point operator()(const typename Mesh::Point& p) {
        TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2]));
        extrema.LoadS2(vertex);
        extrema.Perform();
        auto nearest = extrema.PointOnShape1(1);
        return typename Mesh::Point(nearest.X(), nearest.Y(), nearest.Z());
    }

    TopoDS_Compound border;
    BRepExtrema_DistShapeShape extrema;
};

// Per-face map from STEP vertex hash → wire projector, shared across all segments.
template<class Mesh>
using WireProjectorCachePtr =
    std::shared_ptr<std::unordered_map<int, std::unordered_map<int, StepProjector<Mesh>>>>;

// Build a wire from all boundary edges reachable from a given STEP vertex, skipping
// seam and degenerate edges which are not physical boundaries.
inline TopoDS_Wire get_border_loop_wire(
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
                // Seam edges (e.g. the U=0/2π closure of a full 360° cylinder)
                // belong to only one face but appear twice in its wire with
                // opposing orientations.  They are parametric artifacts with no
                // corresponding physical surface boundary, so skip them.
                if (BRep_Tool::IsClosed(TopoDS::Edge(e), single_face)) {
                    continue;
                }
                // Degenerate edges (e.g. the apex of a cone or sphere) have zero
                // length; adding them to BRepBuilderAPI_MakeWire would corrupt the
                // wire, so skip them.
                if (BRep_Tool::Degenerated(TopoDS::Edge(e))) {
                    continue;
                }
                // Free outer-boundary edges of open shells are also Size()==1
                // but represent real surface boundaries, so fall through.
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

// Build a per-face cache of wire projectors, one per STEP vertex, in parallel.
template<class Mesh>
WireProjectorCachePtr<Mesh> get_edge_vertex_wire_projectors(const TopoDS_Shape& shape) {
    auto wire_projectors =
        std::make_shared<std::unordered_map<int, std::unordered_map<int, StepProjector<Mesh>>>>();

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

                    auto projector = StepProjector<Mesh>(wire);

                    for (TopExp_Explorer wire_exp(wire, TopAbs_VERTEX); wire_exp.More();
                         wire_exp.Next())
                    {
                        auto wire_v = TopoDS::Vertex(wire_exp.Current());
                        int wire_vertex_code = shapeHashCode(wire_v);
                        face_wire_projectors[wire_vertex_code] = projector;
                    }
                }
            }});

    spdlog::debug("done making wire projectors.");
    return wire_projectors;
}

// Associate each border vertex of mesh with the wire projector for its nearest
// STEP edge boundary, propagating assignments to vertices not directly on an edge
// via their already-assigned neighbors.
template<class Mesh>
std::unordered_map<typename Mesh::Vertex_index, StepProjector<Mesh>>
get_border_vertex_projector_map(const TopoDS_Face& face, Mesh& mesh,
                                WireProjectorCachePtr<Mesh> wire_projectors,
                                double tolerance = 1e-3)
{
    int face_code = shapeHashCode(face);
    const auto& face_wire_projectors = (*wire_projectors)[face_code];

    std::unordered_map<typename Mesh::Vertex_index, StepProjector<Mesh>> border_vertex_projector_map;

    std::vector<TopoDS_Edge> edges;
    std::vector<BRepExtrema_DistShapeShape> edge_extremas;
    for (TopExp_Explorer edge_exp(face, TopAbs_EDGE); edge_exp.More(); edge_exp.Next()) {
        edge_extremas.push_back(BRepExtrema_DistShapeShape());
        edge_extremas.back().LoadS1(TopoDS::Edge(edge_exp.Current()));
        edges.push_back(TopoDS::Edge(edge_exp.Current()));
    }

    // First pass: assign projectors to border vertices that already lie on a face edge.
    std::deque<typename Mesh::Vertex_index> vertex_queue;
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

    // Second pass: propagate assignments along the border to vertices not on an edge.
    bool made_progress = true;
    while (made_progress) {
        size_t start_size = vertex_queue.size();
        std::deque<typename Mesh::Vertex_index> remaining;
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

// Project all mesh vertices toward the STEP surface, blending by weight [0,1].
// Border vertices are projected onto their associated wire boundary; interior
// vertices are projected onto the surface.  Weight ramps up each iteration so
// vertices reach the surface exactly on the final pass.
template<class Mesh>
void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr<Mesh> wire_projectors,
                     StepProjector<Mesh>& surface_projector,
                     StepBorderProjector<Mesh>& border_projector,
                     double weight = 1.0)
{
    double w1 = std::max(0.0, std::min(1.0, weight));
    double w2 = 1.0 - w1;

    auto border_vertex_projector_map =
        get_border_vertex_projector_map<Mesh>(face, mesh, wire_projectors);

    for (auto v: mesh.vertices()) {
        auto input = mesh.point(v);

        typename Mesh::Point projected;
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

        mesh.point(v) = typename Mesh::Point(
            w1 * projected.x() + w2 * input.x(),
            w1 * projected.y() + w2 * input.y(),
            w1 * projected.z() + w2 * input.z());
    }
}

// Returns the distance from point (x,y,z) to the nearest point on face, or -1 on failure.
inline double get_distance_to_face(const TopoDS_Face& face, double x, double y, double z) {
    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(x, y, z));
    BRepExtrema_DistShapeShape dist(vertex, face);
    dist.Perform();
    if (dist.IsDone() && dist.Value() >= 0) {
        return dist.Value();
    }
    spdlog::error("Failed to calculate distance.");
    return -1;
}
