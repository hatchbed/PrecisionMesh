#pragma once

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <gp_Pnt.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

#include <precision_mesh/step_subdivision.h>

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

template<class Mesh>
using WireProjectorCachePtr = std::shared_ptr<std::unordered_map<int, std::unordered_map<int, StepProjector<Mesh>>>>;

TopoDS_Wire get_border_loop_wire(
    const TopoDS_Vertex& vertex,
    const TopTools_IndexedDataMapOfShapeListOfShape& vertex_to_edge_map,
    const TopTools_IndexedDataMapOfShapeListOfShape& edge_to_face_map)
{
    size_t wire_size = 0;
    BRepBuilderAPI_MakeWire wire_maker;
    std::unordered_set<int> evaluated;

    // Get the vertices of the target edge
    std::deque<TopoDS_Vertex> vertex_queue = { vertex };

    while (!vertex_queue.empty()) {
        auto v = vertex_queue.front();
        vertex_queue.pop_front();

        auto edges = vertex_to_edge_map.FindFromKey(v);
        for (const auto& e: edges) {
            int edge_code = shapeHashCode(e);
            if (!evaluated.insert(edge_code).second) {
                // ignore already evaluated edge
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
                // Degenerate edges (e.g. the apex of a cone or sphere) are also
                // single-face and have zero length; adding them to
                // BRepBuilderAPI_MakeWire would corrupt the wire, so skip them.
                if (BRep_Tool::Degenerated(TopoDS::Edge(e))) {
                    continue;
                }
                // Free outer-boundary edges of open shells are also Size()==1
                // but represent real surface boundaries that mesh border vertices
                // must project onto, so fall through and include them in the wire.
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
    return wire_maker.Wire();;
}

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
                TopExp::MapShapesAndUniqueAncestors(face, TopAbs_VERTEX, TopAbs_EDGE, vertex_to_edge_map);

                auto face_wire_projectors_it = wire_projectors->find(face_code);
                if (face_wire_projectors_it == wire_projectors->end()) {
                    continue;
                }
                auto& face_wire_projectors = face_wire_projectors_it->second;

                for (TopExp_Explorer vertex_exp(face, TopAbs_VERTEX); vertex_exp.More(); vertex_exp.Next()) {
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

                    for (TopExp_Explorer wire_exp(wire, TopAbs_VERTEX); wire_exp.More(); wire_exp.Next()) {
                        auto wire_v = TopoDS::Vertex(wire_exp.Current());
                        int wire_vertex_code = shapeHashCode(wire_v);
                        face_wire_projectors[wire_vertex_code] = projector;
                    }
                }
            }});

    spdlog::debug("done making wire projectors.");

    return wire_projectors;
}

template<class Mesh>
std::unordered_map<typename Mesh::Vertex_index, StepProjector<Mesh>> get_border_vertex_projector_map(
    const TopoDS_Face& face, Mesh& mesh,
    WireProjectorCachePtr<Mesh> wire_projectors, double tolerance=1e-3)
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

    // first pass to associate mesh border points already on the face border with projectors
    std::deque<typename Mesh::Vertex_index> vertex_queue;
    for (auto v: mesh.vertices()) {
        if (mesh.is_border(v)) {
            auto p = mesh.point(v);
            auto point = gp_Pnt(p[0], p[1], p[2]);
            TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(point);

            // check if the vertex lies on a border edge
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
                // if so, find the nearest edge vertex and the associate projector

                const auto& nearest_edge = edges[nearest_edge_index];
                TopoDS_Vertex v1, v2;
                TopExp::Vertices(TopoDS::Edge(nearest_edge), v1, v2);
                double d1 = point.Distance(BRep_Tool::Pnt(v1));
                double d2 = point.Distance(BRep_Tool::Pnt(v2));
                TopoDS_Vertex nearest_vertex = v1;
                if (d2 < d1) {
                    nearest_vertex = v2;
                }

                int vertex_code = shapeHashCode(nearest_vertex);
                auto projector_it = face_wire_projectors.find(vertex_code);
                if (projector_it == face_wire_projectors.end()) {
                    continue;
                }

                border_vertex_projector_map[v] = projector_it->second;
            }
            else {
                // otherwise queue the vertex for the second pass
                vertex_queue.push_back(v);
            }
        }
    }

    if (vertex_queue.empty()) {
        return border_vertex_projector_map;
    }

    spdlog::debug("border map size: {}", border_vertex_projector_map.size());
    spdlog::debug("vertex queue size: {}", vertex_queue.size());
    spdlog::debug("num edges: {}", edges.size());

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

                    // TODO(malban): track wire parameter to limit search space when projecting
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
            spdlog::warn("Failed to associate {} border vertices with face boundary.", vertex_queue.size());
        }
    }

    return border_vertex_projector_map;
}


template<class Mesh>
void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr<Mesh> wire_projectors,
                     StepProjector<Mesh>& surface_projector,
                     StepBorderProjector<Mesh>& border_projector,
                     double weight = 1.0)
{
    double w1 = std::max(0.0, std::min(1.0, weight));
    double w2 = 1.0 - w1;

    auto border_vertex_projector_map = get_border_vertex_projector_map<Mesh>(face, mesh,
                                                                             wire_projectors);

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

        double nx = w1 * projected.x() + w2 * input.x();
        double ny = w1 * projected.y() + w2 * input.y();
        double nz = w1 * projected.z() + w2 * input.z();

        mesh.point(v) = typename Mesh::Point(nx, ny, nz);
    }
}

double get_distance_to_face(const TopoDS_Face& face, double x, double y, double z) {

    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(x, y, z));
    BRepExtrema_DistShapeShape distance_calculator(vertex, face);
    distance_calculator.Perform();
    if (distance_calculator.IsDone() && distance_calculator.Value() >= 0) {
        return distance_calculator.Value();
    }

    spdlog::error("Failed to calculate distance.");
    return -1;
}


std::vector<TopoDS_Face> get_faces(TopoDS_Shape& shape) {
    std::vector<TopoDS_Face> faces;
    for (TopExp_Explorer iter(shape, TopAbs_FACE); iter.More(); iter.Next()) {
        faces.push_back(TopoDS::Face(iter.Current()));
    }
    return faces;
}

int get_face_type(const TopoDS_Face& face) {
    auto surface = BRep_Tool::Surface(face);

    auto plane = Handle(Geom_Plane)::DownCast(surface);
    if (!plane.IsNull()) {
        return 1;
    }

    auto cylinder = Handle(Geom_CylindricalSurface)::DownCast(surface);
    if (!cylinder.IsNull()) {
        return 2;
    }

    auto sphere = Handle(Geom_SphericalSurface)::DownCast(surface);
    if (!sphere.IsNull()) {
        return 3;
    }

    auto torous = Handle(Geom_ToroidalSurface)::DownCast(surface);
    if (!torous.IsNull()) {
        return 4;
    }

    auto cone = Handle(Geom_ConicalSurface)::DownCast(surface);
    if (!cone.IsNull()) {
        return 5;
    }

    GeomAdaptor_Surface adaptor(surface);
    GeomAbs_SurfaceType surface_type = adaptor.GetType();
    if (surface_type == GeomAbs_Plane) {
        return 1;
    }
    else if (surface_type == GeomAbs_Cylinder) {
        return 2;
    }
    else if (surface_type == GeomAbs_Sphere) {
        return 3;
    }
    else if (surface_type == GeomAbs_Torus) {
        return 4;
    }
    else if (surface_type == GeomAbs_Cone) {
        return 5;
    }

    return 0;   
}

float get_face_area(const TopoDS_Face& face) {
    GProp_GProps surface_props;
    BRepGProp::SurfaceProperties(face, surface_props);
    return surface_props.Mass();
}

std::array<float, 3> get_face_centroid(const TopoDS_Face& face) {
    GProp_GProps surface_props;
    BRepGProp::SurfaceProperties(face, surface_props);
    auto centroid = surface_props.CentreOfMass();
    return { static_cast<float>(centroid.X()), static_cast<float>(centroid.Y()), static_cast<float>(centroid.Z()) };
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

template<class Mesh>
std::vector<std::pair<Mesh, TopoDS_Face>> tessellate_shape(const TopoDS_Shape& shape, double max_surface_error) {
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

        // copy vertices
        auto num_verts = triangulation->NbNodes();
        gp_Trsf transform = loc.Transformation();
        std::unordered_map<int, size_t> vertex_map;
        for (int i = 1; i <= num_verts; i++) {
            auto point = triangulation->Node(i);
            point.Transform(transform);
            vertex_map[i] = vertex_buffer.size();
            vertex_buffer.push_back({point.X(), point.Y(), point.Z()});
        }

        auto num_triangles = triangulation->NbTriangles();
        const TopAbs_Orientation orientation = iter.Current().Orientation();
        for (int i = 1; i <= num_triangles; i++) {
            auto triangle = triangulation->Triangle(i);

            Standard_Integer anId[3];
            triangle.Get(anId[0], anId[1], anId[2]);
            if (orientation == TopAbs_REVERSED) {
                // Swap 1, 2.
                Standard_Integer aTmpIdx = anId[1];
                anId[1] = anId[2];
                anId[2] = aTmpIdx;
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

        PMP::repair_polygon_soup(vertex_buffer, face_buffer,
                                 CGAL::parameters::geom_traits(PointArray_traits()));

        // create mesh from vertex and face buffers
        std::vector<typename Mesh::Vertex_index> vertex_indices;
        for (const auto& vertex: vertex_buffer) {
            vertex_indices.push_back(mesh.add_vertex(typename Mesh::Point(vertex[0], vertex[1], vertex[2])));
        }
        for (size_t i = 0; i < face_buffer.size(); i++) {
            const auto& face = face_buffer[i];
            mesh.add_face(vertex_indices[face[0]], vertex_indices[face[1]], vertex_indices[face[2]]);
        }

        if (mesh.number_of_faces() == 0) {
            continue;
        }

        tessellation.push_back(std::make_pair(mesh, face));
    }

    return tessellation;
}

size_t get_short_edge_count(const std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation, double min_edge_length) {
    size_t short_edges = 0;
    double min_edge_length_sq = min_edge_length * min_edge_length;
    for (const auto& [mesh, face]: tessellation) {
        auto plane = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(face));
        if (!plane.IsNull()) {
            // ignore planes
            continue;
        }

        for (auto e: mesh.edges()) {
            if (!mesh.is_border(e)) {
                continue;
            }

            auto h = mesh.halfedge(e);
            auto v_start = mesh.source(h);
            auto v_end = mesh.target(h);
            auto p_start = mesh.point(v_start);
            auto p_end = mesh.point(v_end);
            if (CGAL::squared_distance(p_start, p_end) < min_edge_length_sq) {
                short_edges++;
            }
        }
    }

    return short_edges;
}

// Auto-tune the OCCT tessellation deflection parameter for a given shape when
// the user has not specified an explicit --max-surface-error.
//
// OCCT's BRepMesh_IncrementalMesh accepts a deflection (max surface error) that
// controls tessellation density.  PrecisionMesh is concerned with uniform triangle
// sizes bounded by min/max edge length.  Edges that are too long can be split
// during remeshing, but there is no mechanism to merge edges that are already too
// short.  Some short border edges are unavoidable regardless of deflection — for
// example a face whose width is less than min_edge_length will always produce short
// edges at its boundary.  These are the "inherent" short edges.
//
// The goal is to find the finest deflection (most accurate surface representation)
// that does not introduce MORE short border edges than this inherent baseline.
// Going finer improves surface accuracy; the short-edge count acts as the brake
// that prevents the deflection from producing a mesh that cannot be remeshed
// cleanly.  For flat surfaces, deflection has no effect on triangle count, so the
// search converges to a very small value harmlessly.  For curved surfaces (e.g. a
// sphere), the search converges to the deflection where arc segments on shared
// edges just reach min_edge_length — naturally sizing the tessellation to the
// requested edge-length budget.
//
// The search is a binary search over [0, min_edge_length].  The baseline short-edge
// count is measured at the coarsest deflection (min_edge_length).  A 1% tolerance
// on that count allows for noise while still preventing a meaningful increase.
template<class Mesh>
double find_surface_error_param(const TopoDS_Shape& shape, double min_edge_length, int max_iterations=10,
                                double conversion_scale=1)
{
    spdlog::info("finding max surface error param ...");
    max_iterations = std::max(2, std::min(100, max_iterations));

    double threshold_max = min_edge_length;
    double threshold_min = 0;
    auto tessellation = tessellate_shape<Mesh>(shape, threshold_max);

    size_t max_short_edges = get_short_edge_count(tessellation, min_edge_length);
    spdlog::info("  initial short edges: {}", max_short_edges);

    for (int i = 0; i < max_iterations; i++) {
        double threshold = (threshold_min + threshold_max) / 2.0;
        tessellation = tessellate_shape<Mesh>(shape, threshold);
        size_t num_short_edges = get_short_edge_count(tessellation, min_edge_length);
        if (num_short_edges > max_short_edges * 1.01) {
            // This deflection is too fine — it introduces new short edges beyond
            // the inherent baseline.  Raise the lower bound to exclude it.
            threshold_min = threshold;
        }
        else {
            // This deflection is acceptable.  Lower the upper bound and try finer.
            threshold_max = threshold;
        }
        spdlog::info("    {} short edges at {}", num_short_edges, threshold * conversion_scale);
    }

    spdlog::info("  found max surface error param = {}", threshold_max * conversion_scale);

    return threshold_max;
}
