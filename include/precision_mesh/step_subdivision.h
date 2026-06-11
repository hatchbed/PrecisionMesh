#pragma once

#include <tuple>
#include <utility>
#include <vector>

#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <TopoDS_Face.hxx>

#include <precision_mesh/step_types.h>

enum class CurvedFaceType { None, Cylinder, Cone, Torus, Revolution, Extrusion };

struct FaceTessellationSteps {
    CurvedFaceType type = CurvedFaceType::None;
    int u_steps = 1;
    int v_steps = 1;
};

// Classify a face and return the UV subdivision steps that would be applied.
// Returns type=None and u=v=1 for non-regularly-curved faces (planes, splines, etc.).
// CDT-eligible face types use this to determine the interior grid density for
// uv_grid_retessellate().
FaceTessellationSteps get_face_tessellation_steps(const TopoDS_Face& face,
                                                   double min_edge_length,
                                                   double max_edge_length,
                                                   double max_surface_error);

// Single source of truth for which face types use the UV-grid CDT path
// (uv_grid_retessellate) instead of BRepAlgoAPI_Splitter subdivision + isotropic
// remeshing.  Such faces skip subdivision entirely; faces whose steps are u=v=1
// additionally need no CDT (the BRepMesh interior is already at target density).
// Extend here as cone/torus/revolution/extrusion CDT support lands.
inline bool cdt_eligible(const FaceTessellationSteps& steps)
{
    return steps.type == CurvedFaceType::Cylinder ||
           steps.type == CurvedFaceType::Cone;
}

// Compute the number of U (angular) and V (axial) parametric subdivisions for a
// cylindrical face such that the resulting sub-faces satisfy edge length and surface
// error constraints.  Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_cylinder_steps(const Geom_CylindricalSurface& surface,
                                          const TopoDS_Face& face,
                                          double min_edge_length,
                                          double max_edge_length,
                                          double max_surface_error);

// Compute subdivision steps for a conical face.  Uses the radius one ring segment
// from the narrow end as the characteristic radius to avoid over-tessellating the
// apex.  Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_cone_steps(const Geom_ConicalSurface& surface,
                                      const TopoDS_Face& face,
                                      double min_edge_length,
                                      double max_edge_length,
                                      double max_surface_error);

// Compute U subdivision steps for a general surface of revolution by sampling the
// profile curve to find the minimum radius from the axis.  V subdivision is left to
// isotropic remeshing.  Returns {u_steps, 1}.
std::pair<int,int> compute_revolution_steps(const Geom_SurfaceOfRevolution& surface,
                                            const TopoDS_Face& face,
                                            double min_edge_length,
                                            double max_edge_length,
                                            double max_surface_error);

// Compute V (extrusion direction) subdivision steps for a linear extrusion face.
// The extrusion direction is always straight so only the edge length constraint
// applies.  Returns {1, v_steps}.
std::pair<int,int> compute_extrusion_steps(const Geom_SurfaceOfLinearExtrusion& surface,
                                           const TopoDS_Face& face,
                                           double max_edge_length);

// Compute subdivision steps for a toroidal face using the outer equator radius for
// U and the tube cross-section radius for V, with a shared diagonal budget.
// Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_torus_steps(const Geom_ToroidalSurface& surface,
                                       const TopoDS_Face& face,
                                       double min_edge_length,
                                       double max_edge_length,
                                       double max_surface_error);

// Split a single STEP face into a grid of sub-faces along parametric iso-lines.
// Returns the original face unchanged if both u_steps and v_steps are < 2.
std::vector<TopoDS_Face> subdivide_face(const TopoDS_Face& face, int u_steps, int v_steps);

// Pre-subdivide all faces in a STEP shape according to curvature and edge length
// constraints, then sew the results back into a compound.  Returns the new shape
// and a map from each sub-face back to its original face index.
std::tuple<TopoDS_Shape, FaceMap> subdivide_step_shape(TopoDS_Shape& shape,
                                                        double min_edge_length,
                                                        double max_edge_length,
                                                        double max_surface_error);
