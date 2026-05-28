#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <BRep_Tool.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>

#include <precision_mesh/mesh_util.h>

template<class Mesh>
std::vector<std::pair<Mesh, TopoDS_Face>> tessellate_shape(const TopoDS_Shape& shape,
                                                            double max_surface_error)
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
        gp_Trsf transform = loc.IsIdentity() ? gp_Trsf() : loc.Transformation();
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

        PMP::repair_polygon_soup(vertex_buffer, face_buffer,
                                 CGAL::parameters::geom_traits(PointArray_traits()));

        std::vector<typename Mesh::Vertex_index> vertex_indices;
        for (const auto& vertex: vertex_buffer) {
            vertex_indices.push_back(mesh.add_vertex(
                typename Mesh::Point(vertex[0], vertex[1], vertex[2])));
        }
        for (const auto& f: face_buffer) {
            mesh.add_face(vertex_indices[f[0]], vertex_indices[f[1]], vertex_indices[f[2]]);
        }

        if (mesh.number_of_faces() == 0) {
            continue;
        }

        tessellation.push_back(std::make_pair(mesh, face));
    }

    return tessellation;
}

// Count border edges shorter than min_edge_length, ignoring planar faces.
size_t get_short_edge_count(const std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                            double min_edge_length);

template<class Mesh>
double find_surface_error_param(const TopoDS_Shape& shape, double min_edge_length,
                                int max_iterations = 10, double conversion_scale = 1)
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

bool save_shape_as_step(const std::string& path, const TopoDS_Shape& shape);
