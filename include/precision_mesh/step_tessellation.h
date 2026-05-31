#pragma once

#include <string>
#include <utility>
#include <vector>

// mesh_util.h (CGAL) must come before any OCCT header that defines the Handle macro.
#include <precision_mesh/mesh_util.h>

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

// repair_inverted_faces: re-triangulate faces where BRepMesh produced folded / holed /
// dropped-corner output (see step_tessellation.cpp).
std::vector<std::pair<Mesh, TopoDS_Face>> tessellate_shape(const TopoDS_Shape& shape,
                                                            double max_surface_error,
                                                            bool repair_inverted_faces = true);

// Build minimal fan-polygon meshes from each face's boundary wire vertices,
// bypassing BRepMesh entirely.  Useful for visualising raw subdivision geometry.
std::vector<std::pair<Mesh, TopoDS_Face>> boundary_meshes(const TopoDS_Shape& shape);

// Count border edges shorter than min_edge_length, ignoring planar faces.
size_t get_short_edge_count(const std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                            double min_edge_length);

// Auto-tune the OCCT tessellation deflection so the resulting mesh does not
// introduce short border edges beyond the inherent baseline.
double find_surface_error_param(const TopoDS_Shape& shape, double min_edge_length,
                                int max_iterations = 10, double conversion_scale = 1);

bool save_shape_as_step(const std::string& path, const TopoDS_Shape& shape);
