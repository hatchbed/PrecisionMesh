#pragma once

#include <string>
#include <vector>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <precision_mesh/mesh_util.h>

template <typename Traits>
bool saveComponentsToStl(const std::string& path, const std::vector<Mesh>& components, float scale=1) {
    Traits traits;

    std::vector<typename Mesh::Point> vertices;
    std::vector<std::vector<size_t>> faces;

    typedef typename Traits::Less_xyz_3 Less_xyz_3;
    typedef std::map<Mesh::Point, size_t, Less_xyz_3> VertexIndex;
    VertexIndex vertex_lookup(traits.less_xyz_3_object());

    for (size_t i = 0; i < components.size(); i++) {
        const auto& mesh = components[i];
        for (const auto& f : mesh.faces()) {
            std::vector<size_t> face;
            for (const auto& v : mesh.vertices_around_face(mesh.halfedge(f))) {
                const auto& point = mesh.point(v);
                auto result = vertex_lookup.find(point);
                if (result == vertex_lookup.end()) {
                    size_t index = vertices.size();
                    vertex_lookup[point] = index;
                    face.push_back(index);
                    vertices.push_back({ 
                        point[0] * scale,
                        point[1] * scale,
                        point[2] * scale
                    });
                }
                else {
                  face.push_back(result->second);
                }
            }
            faces.push_back(face);
        }
    }

    CGAL::IO::write_STL(path, vertices, faces, CGAL::parameters::stream_precision(17));

    return true;
}
