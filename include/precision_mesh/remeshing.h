#pragma once

#include <memory>
#include <vector>

#include <TopoDS_Face.hxx>

#include <precision_mesh/mesh_util.h>
#include <precision_mesh/step_projection.h>

struct RemeshParams {
    double min_edge_length;
    double max_edge_length;
    double max_surface_error;
    int iterations;
    bool is_step;
    bool no_projection;
};

// Split border edges longer than max_edge_length on each mesh in the collection.
void split_border_edges(std::vector<Mesh>& meshes, double max_edge_length);

// Detect sharp crease edges and split those longer than max_edge_length.
// Only meaningful for generic (non-STEP) meshes; the crease property map is
// added, used, and removed inside this function.
void split_crease_edges(std::vector<Mesh>& meshes, double crease_angle, double max_edge_length);

// Run adaptive isotropic remeshing for the requested number of iterations,
// optionally projecting vertices back onto the STEP surface after each pass.
// For non-STEP meshes, pass empty/null projector containers and is_step=false.
void remesh_and_project(
    std::vector<Mesh>& meshes,
    const std::vector<TopoDS_Face>& segments,
    WireProjectorCachePtr wire_projectors,
    std::vector<std::unique_ptr<StepProjector>>& surface_projectors,
    std::vector<std::unique_ptr<StepBorderProjector>>& border_projectors,
    const RemeshParams& params);
