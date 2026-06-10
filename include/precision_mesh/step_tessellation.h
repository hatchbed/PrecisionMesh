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

// Replace the interior triangulation of `mesh` with a UV-grid CDT tessellation.
// Border vertices (mesh.is_border(v)) are fixed constraints; interior is filled
// by Jacobian-corrected Steiner points evaluated on `face` via BRepClass_FaceClassifier.
// Requires "v:uv" property map to exist on `mesh` (set by tessellate_shape).
// Returns false (mesh unchanged) if CDT fails or produces no triangles.
// dump_dir: when non-empty, write the 2D CDT to <dump_dir>/Face_<face_idx>_CDT.obj
// for diagnosis (the directory must already exist).
bool uv_grid_retessellate(Mesh& mesh, const TopoDS_Face& face, int u_steps, int v_steps,
                          double min_edge_length, size_t face_idx = 0,
                          const std::string& dump_dir = {});

// Post-tessellation open boundary loop repair.  Merges the soup by exact vertex
// position, traces open boundary edges into closed loops, and classifies each by
// comparing the loop's polygon area to ref_area = min(avg_tri_area, 0.5*min_edge^2):
//   loop_area < ref_area * collapse_area_ratio -> collapse (snap the shortest edge's
//       vertex pair together and rebuild the affected segment meshes)
//   otherwise                                  -> fill (fan-triangulate the loop)
// Capping by 0.5*min_edge_length^2 ensures the threshold stays tight even when
// the average triangle area is large (e.g. a part with big flat faces).
void repair_open_boundary_loops(std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                                 double min_edge_length,
                                 double collapse_area_ratio = 0.01);
