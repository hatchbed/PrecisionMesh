#pragma once

#include <utility>

#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <TopoDS_Face.hxx>

enum class CurvedFaceType { None, Cylinder, Cone, Torus, Revolution, Extrusion };

struct FaceTessellationSteps {
    CurvedFaceType type = CurvedFaceType::None;
    int u_steps = 1;
    int v_steps = 1;
};

// Classify a face and return the UV step counts that size its CDT interior grid.
// Returns type=None and u=v=1 for non-regularly-curved faces (planes, splines, etc.).
// CDT-eligible face types use this to determine the interior grid density for
// uv_grid_retessellate().
FaceTessellationSteps get_face_tessellation_steps(const TopoDS_Face& face,
                                                   double min_edge_length,
                                                   double max_edge_length,
                                                   double max_surface_error);

// Single source of truth for which analytic face types use the UV-grid CDT path
// (uv_grid_retessellate).  Every other face (plane, sphere, B-spline/Bezier patch,
// offset, ...) keeps its BRepMesh interior and is refined by isotropic remeshing.
inline bool cdt_eligible(const FaceTessellationSteps& steps)
{
    return steps.type == CurvedFaceType::Cylinder ||
           steps.type == CurvedFaceType::Cone ||
           steps.type == CurvedFaceType::Torus ||
           steps.type == CurvedFaceType::Revolution ||
           steps.type == CurvedFaceType::Extrusion;
}

// Compute the number of U (angular) and V (axial) parametric step counts for a
// cylindrical face such that the CDT tessellation satisfies edge length and surface
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

// Compute subdivision steps for a general surface of revolution.  U (angular) is sized
// from the minimum profile radius to the axis; V (profile) from the profile arc length
// and its tightest radius of curvature (sagitta).  Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_revolution_steps(const Geom_SurfaceOfRevolution& surface,
                                            const TopoDS_Face& face,
                                            double min_edge_length,
                                            double max_edge_length,
                                            double max_surface_error);

// Compute subdivision steps for a surface of linear extrusion S(u,v)=C(u)+v*Dir.
// U is the (generally curved) profile curve and V is the straight extrusion ruling.
// V is sized purely by the edge-length budget (zero surface error along a straight
// line); U is sized from the profile arc length AND its tightest radius of curvature
// (sagitta), mirroring the revolution profile.  Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_extrusion_steps(const Geom_SurfaceOfLinearExtrusion& surface,
                                           const TopoDS_Face& face,
                                           double min_edge_length,
                                           double max_edge_length,
                                           double max_surface_error);

// Compute subdivision steps for a toroidal face using the outer equator radius for
// U and the tube cross-section radius for V, with a shared diagonal budget.
// Returns {u_steps, v_steps}, both >= 1.
std::pair<int,int> compute_torus_steps(const Geom_ToroidalSurface& surface,
                                       const TopoDS_Face& face,
                                       double min_edge_length,
                                       double max_edge_length,
                                       double max_surface_error);
