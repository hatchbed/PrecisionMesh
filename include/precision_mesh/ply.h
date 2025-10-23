#pragma once

#include <fstream>
#include <iostream>
#include <unordered_map>

#include <tinyply.h>
#include <TopoDS_Face.hxx>

#include <precision_mesh/mesh_util.h>
#include <precision_mesh/step_util.h>

template <typename Traits>
bool saveComponentsToPly(const std::string& path, const std::vector<Mesh>& components, 
                         const std::vector<TopoDS_Face>& step_faces,
                         const std::unordered_map<size_t, int>& component_map) 
{
    Traits traits;

    std::vector<std::array<float, 3>> vertices;
    std::vector<std::array<uint32_t, 3>> faces;
    std::vector<int32_t> component_ids;
    std::vector<float> component_areas(step_faces.size(), 0.0f);
    std::vector<std::array<float, 3>> component_centers(step_faces.size(), {0.0f, 0.0f, 0.0f});
    std::vector<int32_t> component_types(step_faces.size(), 0);

    typedef typename Traits::Less_xyz_3 Less_xyz_3;
    typedef std::map<Mesh::Point, uint32_t, Less_xyz_3> VertexIndex;
    VertexIndex vertex_lookup(traits.less_xyz_3_object());

    for (size_t i = 0; i < components.size(); i++) {
        int component_id = i;
        auto component_id_it = component_map.find(i);
        if (component_id_it != component_map.end()) {
            component_id = component_id_it->second;
        }

        printf("component %zu of %zu, id = %d\n", i, components.size(), component_id);

        const auto& mesh = components[i];
        for (const auto& f : mesh.faces()) {
            size_t face_idx = 0;
            std::array<uint32_t, 3> face;
            for (const auto& v : mesh.vertices_around_face(mesh.halfedge(f))) {
                if (face_idx >= 3) {
                    break;
                }
                const auto& point = mesh.point(v);
                auto result = vertex_lookup.find(point);
                if (result == vertex_lookup.end()) {
                    uint32_t index = static_cast<uint32_t>(vertices.size());
                    vertices.push_back({static_cast<float>(point[0]),
                                        static_cast<float>(point[1]),
                                        static_cast<float>(point[2])});
                    vertex_lookup[point] = index;
                    face[face_idx++] = index;
                }
                else {
                  face[face_idx++] = result->second;
                }
            }

            if (face_idx < 3) {
                continue;
            }

            component_ids.push_back(component_id);
            faces.push_back(face);
        }
    }

    printf("writing metadata for %zu faces ...\n", step_faces.size());
    for (size_t i = 0; i < step_faces.size(); i++) {
        component_types[i] = get_face_type(step_faces[i]);
        component_areas[i] = get_face_area(step_faces[i]);
        component_centers[i] = get_face_centroid(step_faces[i]);
    }

    std::filebuf fb;
    //fb.open(path, std::ios::out | std::ios::binary);
    fb.open(path, std::ios::out);
    std::ostream output_stream(&fb);

    tinyply::PlyFile ply_file;
    ply_file.add_properties_to_element("vertex", { "x", "y", "z" },
        tinyply::Type::FLOAT32, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    ply_file.add_properties_to_element("face", { "vertex_index" },
        tinyply::Type::UINT32, faces.size(), reinterpret_cast<uint8_t*>(faces.data()), tinyply::Type::UINT8, 3);
    ply_file.add_properties_to_element("face", { "component" },
        tinyply::Type::INT32, component_ids.size(), reinterpret_cast<uint8_t*>(component_ids.data()), tinyply::Type::INVALID, 0);
    ply_file.add_properties_to_element("component", { "type" },
        tinyply::Type::INT32, component_types.size(), reinterpret_cast<uint8_t*>(component_types.data()), tinyply::Type::INVALID, 0);
    ply_file.add_properties_to_element("component", { "area" },
        tinyply::Type::FLOAT32, component_areas.size(), reinterpret_cast<uint8_t*>(component_areas.data()), tinyply::Type::INVALID, 0);
    ply_file.add_properties_to_element("component", { "cx", "cy", "cz" },
        tinyply::Type::FLOAT32, component_centers.size(), reinterpret_cast<uint8_t*>(component_centers.data()), tinyply::Type::INVALID, 0);
    //ply_file.write(output_stream, true);
    ply_file.write(output_stream, false);

    return true;
}
