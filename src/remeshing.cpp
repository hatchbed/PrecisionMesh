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
#include <mutex>

#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <boost/iterator/function_output_iterator.hpp>

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
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));
        border_num_before += border_edges.size();
        PMP::split_long_edges(border_edges, max_edge_length, mesh);
        border_edges.clear();
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
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

void remesh_and_project(
    std::vector<Mesh>& meshes,
    const std::vector<TopoDS_Face>& segments,
    WireProjectorCachePtr wire_projectors,
    std::vector<std::unique_ptr<StepProjector>>& surface_projectors,
    std::vector<std::unique_ptr<StepBorderProjector>>& border_projectors,
    const RemeshParams& params)
{
    double max_remeshing_surface_error =
        std::min(params.max_surface_error, params.min_edge_length * 0.1);

    for (int i = 0; i < params.iterations; i++) {
        spdlog::info("    iteration {}", i + 1);
        spdlog::info("      remeshing ...");
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t m = r.begin(); m != r.end(); ++m) {
                    Mesh& mesh = meshes[m];
                    const std::pair edge_min_max{params.min_edge_length, params.max_edge_length};
                    PMP::Adaptive_sizing_field<Mesh> sizing_field(max_remeshing_surface_error,
                                                                  edge_min_max, faces(mesh), mesh);
                    try {
                        if (params.is_step) {
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .protect_constraints(true));
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

        if (params.is_step && !params.no_projection) {
            spdlog::info("      projecting ...");
            // Weight increases each iteration (1/N, 1/(N-1), ..., 1/1), reaching full
            // projection (weight=1) on the final iteration.  The gradual schedule prevents
            // degenerate faces from aggressively snapping vertices before the mesh topology
            // has had time to adapt.
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m = r.begin(); m != r.end(); ++m) {
                        double weight = 1.0 / (params.iterations - i);
                        project_to_step(segments[m], meshes[m], wire_projectors,
                                              *surface_projectors[m], *border_projectors[m], weight);
                    }});
        }
    }
}
