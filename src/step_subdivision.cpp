#include <precision_mesh/step_subdivision.h>
#include <precision_mesh/step_tessellation.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

#include <spdlog/spdlog.h>

#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <Geom_Curve.hxx>
#include <GeomLProp_CLProps.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

std::pair<int,int> compute_cylinder_steps(const Geom_CylindricalSurface& cylinder,
                                          const TopoDS_Face& face,
                                          double min_edge_length,
                                          double max_edge_length,
                                          double max_surface_error)
{
    int u_steps = 1;
    int v_steps = 1;

    auto dir = cylinder.Axis().Direction();
    double radius = cylinder.Radius();

    spdlog::debug("cylinder:");
    spdlog::debug("  radius: {}", radius);
    spdlog::debug("  axis: {}, {}, {}", dir.X(), dir.Y(), dir.Z());

    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("  U (angle): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (height): {} -> {} ({})", v1, v2, v2 - v1);
    spdlog::debug("  max edge length: {}", max_edge_length);

    // Adjust max edge length down to account for the hypotenuse of the right
    // triangle formed by the U and V edges.  Edge lengths cannot exceed the
    // chord diameter.
    double max_u_edge_length = max_edge_length / std::sqrt(2.0);
    max_u_edge_length = std::min(max_u_edge_length, 2 * radius);
    double min_u_edge_length = std::min(min_edge_length, 2 * radius);

    spdlog::debug("  max u edge length: {}", max_u_edge_length);
    spdlog::debug("  min u edge length: {}", min_u_edge_length);

    double angle = 0;

    if (max_surface_error / radius <= 2) {
        double max_angle = 2 * std::acos(1 - max_surface_error / radius);
        spdlog::debug("  max surface error angle: {}", max_angle);
        angle = max_angle;
    }

    double min_angle = 2 * std::asin(min_u_edge_length / (2 * radius));
    spdlog::debug("  min edge length angle: {}", min_angle);
    if (min_angle > angle) {
        angle = min_angle;
    }

    double max_angle = 2 * std::asin(max_u_edge_length / (2 * radius));
    spdlog::debug("  max edge length angle: {}", max_angle);
    if (max_angle < angle) {
        angle = max_angle;
    }

    if (angle > 0) {
        spdlog::debug("  step size: {}", angle);
        int steps = std::ceil((u2 - u1) / angle);
        spdlog::debug("  steps: {}", steps);
        angle = (u2 - u1) / steps;
        spdlog::debug("  actual step size: {}", angle);
        if (steps > 1) {
            u_steps = steps;
        }
    }

    double u_length = 2 * radius * std::sin(angle / 2);
    spdlog::debug("  u length: {}", u_length);

    double max_v_edge_length =
        std::sqrt(std::max(0.0, max_edge_length * max_edge_length - u_length * u_length));
    spdlog::debug("  max v length: {}", max_v_edge_length);

    if (v2 - v1 > max_v_edge_length) {
        v_steps = static_cast<int>(std::ceil((v2 - v1) / max_v_edge_length));
    }

    spdlog::debug("  v length: {}", (v2 - v1) / v_steps);

    return {u_steps, v_steps};
}

std::pair<int,int> compute_cone_steps(const Geom_ConicalSurface& cone,
                                      const TopoDS_Face& face,
                                      double min_edge_length,
                                      double max_edge_length,
                                      double max_surface_error)
{
    int u_steps = 1;
    int v_steps = 1;

    double ref_radius = cone.RefRadius();
    double semi_angle = cone.SemiAngle();

    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("cone:");
    spdlog::debug("  ref radius: {}", ref_radius);
    spdlog::debug("  semi angle: {}", semi_angle);
    spdlog::debug("  U (angle): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (slant height): {} -> {} ({})", v1, v2, v2 - v1);

    // r(v) = RefRadius + v * sin(SemiAngle)
    double r_v1 = ref_radius + v1 * std::sin(semi_angle);
    double r_v2 = ref_radius + v2 * std::sin(semi_angle);
    double r_min = std::min(r_v1, r_v2);
    spdlog::debug("  r at v1: {}, r at v2: {}", r_v1, r_v2);

    // Use the radius one ring segment from the narrow end as the characteristic
    // radius.  This is always positive (even for a true apex where r_min==0) and
    // prevents over-tessellating the narrow tip.
    double r_ring = r_min + max_edge_length * std::sin(std::abs(semi_angle));
    spdlog::debug("  r_ring (one step from narrow end): {}", r_ring);

    double max_u_edge_length = std::min(max_edge_length / std::sqrt(2.0), 2 * r_ring);
    double min_u_edge_length = std::min(min_edge_length, 2 * r_ring);
    spdlog::debug("  max u edge length: {}", max_u_edge_length);
    spdlog::debug("  min u edge length: {}", min_u_edge_length);

    double angle = 0;

    if (max_surface_error / r_ring <= 2) {
        angle = 2 * std::acos(1 - max_surface_error / r_ring);
        spdlog::debug("  max surface error angle: {}", angle);
    }

    double min_angle = 2 * std::asin(min_u_edge_length / (2 * r_ring));
    if (min_angle > angle) {
        angle = min_angle;
    }

    double max_u_angle = 2 * std::asin(max_u_edge_length / (2 * r_ring));
    if (max_u_angle < angle) {
        angle = max_u_angle;
    }

    if (angle > 0) {
        int steps = static_cast<int>(std::ceil((u2 - u1) / angle));
        spdlog::debug("  U steps: {}", steps);
        angle = (u2 - u1) / steps;
        spdlog::debug("  actual U angle: {}", angle);
        if (steps > 1) {
            u_steps = steps;
        }
    }

    // Cone generator lines are straight, so V edges lie exactly on the surface
    // regardless of step size.  Only the diagonal budget applies.
    double u_length = 2 * r_ring * std::sin(angle / 2);
    spdlog::debug("  u chord at r_ring: {}", u_length);

    double max_v_edge_length =
        std::sqrt(std::max(0.0, max_edge_length * max_edge_length - u_length * u_length));
    spdlog::debug("  max v edge length: {}", max_v_edge_length);

    if (max_v_edge_length > 0 && (v2 - v1) > max_v_edge_length) {
        v_steps = static_cast<int>(std::ceil((v2 - v1) / max_v_edge_length));
        spdlog::debug("  V steps: {}", v_steps);
    }

    return {u_steps, v_steps};
}

std::pair<int,int> compute_revolution_steps(const Geom_SurfaceOfRevolution& revolution,
                                            const TopoDS_Face& face,
                                            double min_edge_length,
                                            double max_edge_length,
                                            double max_surface_error)
{
    int u_steps = 1;

    gp_Pnt axis_origin = revolution.Axis().Location();
    gp_Vec axis_dir(revolution.Axis().Direction());
    Handle(Geom_Curve) basis_curve = revolution.BasisCurve();

    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("surface of revolution:");
    spdlog::debug("  U (angle): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (profile): {} -> {} ({})", v1, v2, v2 - v1);

    // Sample the profile curve to find the minimum radius from the revolution axis
    // across the face's V range.  This characterizes the tightest U-direction curvature.
    // Analytical revolved surfaces (cylinder, cone, sphere, torus) are already caught
    // by earlier DownCasts and never reach here.
    const int num_samples = 20;
    double r_min = std::numeric_limits<double>::max();
    for (int i = 0; i <= num_samples; i++) {
        double v = v1 + (v2 - v1) * i / num_samples;
        gp_Pnt pt;
        basis_curve->D0(v, pt);
        gp_Vec OP(axis_origin, pt);
        double along = OP.Dot(axis_dir);
        gp_Vec perp = OP - axis_dir * along;
        r_min = std::min(r_min, perp.Magnitude());
    }
    spdlog::debug("  min radius from axis: {}", r_min);

    // Skip U subdivision if the profile passes through or very near the axis — the
    // apex region is degenerate and isotropic remeshing handles it better than
    // uniform angular cuts.
    if (r_min > min_edge_length) {
        double max_u_edge_length = std::min(max_edge_length / std::sqrt(2.0), 2 * r_min);
        double min_u_edge_length = std::min(min_edge_length, 2 * r_min);

        double angle = 0;
        if (max_surface_error / r_min <= 2) {
            angle = 2 * std::acos(1 - max_surface_error / r_min);
            spdlog::debug("  max surface error angle: {}", angle);
        }

        double min_angle = 2 * std::asin(min_u_edge_length / (2 * r_min));
        if (min_angle > angle) {
            angle = min_angle;
        }

        double max_u_angle = 2 * std::asin(max_u_edge_length / (2 * r_min));
        if (max_u_angle < angle) {
            angle = max_u_angle;
        }

        if (angle > 0) {
            int steps = static_cast<int>(std::ceil((u2 - u1) / angle));
            spdlog::debug("  U steps: {}", steps);
            if (steps > 1) {
                u_steps = steps;
            }
        }
    }

    // V (profile) sizing.  Unlike a cylinder/cone whose V is a straight ruling, the
    // revolution profile is a general curve: its arc length is NOT proportional to the
    // parameter, and it carries its own curvature.  The UV-grid CDT path owns the
    // interior (it does not lean on isotropic remeshing for V), so size V rows from
    // BOTH the profile arc length (edge-length budget) and the profile's tightest
    // radius of curvature (sagitta) — mirroring compute_torus_steps for the tube.
    int v_steps = 1;
    {
        const int num_v_samples = 40;
        double profile_len = 0.0;
        double rho_min = std::numeric_limits<double>::max();
        gp_Vec d_prev;
        GeomLProp_CLProps clp(basis_curve, 2, 1e-7);
        for (int i = 0; i <= num_v_samples; i++) {
            double v = v1 + (v2 - v1) * i / num_v_samples;
            gp_Pnt p; gp_Vec d1;
            basis_curve->D1(v, p, d1);
            if (i > 0) profile_len += 0.5 * (d1.Magnitude() + d_prev.Magnitude())
                                          * (v2 - v1) / num_v_samples;
            d_prev = d1;
            clp.SetParameter(v);
            if (clp.IsTangentDefined()) {
                double k = clp.Curvature();
                if (k > 1e-12) rho_min = std::min(rho_min, 1.0 / k);
            }
        }
        spdlog::debug("  profile arc length: {}", profile_len);
        spdlog::debug("  profile min radius of curvature: {}", rho_min);

        double v_chord = max_edge_length;
        if (rho_min < std::numeric_limits<double>::max() &&
            max_surface_error / rho_min <= 2) {
            double v_angle = 2 * std::acos(1 - max_surface_error / rho_min);
            v_chord = std::min(v_chord, 2 * rho_min * std::sin(v_angle / 2));
        }
        if (v_chord > 0 && profile_len > v_chord) {
            v_steps = static_cast<int>(std::ceil(profile_len / v_chord));
            spdlog::debug("  V steps: {}", v_steps);
        }
    }

    return {u_steps, v_steps};
}

std::pair<int,int> compute_extrusion_steps(const Geom_SurfaceOfLinearExtrusion& extrusion,
                                           const TopoDS_Face& face,
                                           double min_edge_length,
                                           double max_edge_length,
                                           double max_surface_error)
{
    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("surface of linear extrusion:");
    spdlog::debug("  U (profile): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (extrusion): {} -> {} ({})", v1, v2, v2 - v1);

    // V is the extrusion direction — a straight ruling with zero surface error.
    // Subdivide purely by max_edge_length (|Dir| is unit, so v IS the 3D distance).
    int v_steps = 1;
    if ((v2 - v1) > max_edge_length) {
        v_steps = static_cast<int>(std::ceil((v2 - v1) / max_edge_length));
        spdlog::debug("  V steps: {}", v_steps);
    }

    // U is the profile curve C(u): a general curve like the revolution profile, so size
    // U rows from BOTH the profile arc length (edge-length budget) and the profile's
    // tightest radius of curvature (sagitta).  The UV-grid CDT path owns the interior
    // for extrusion faces, so U is no longer left to isotropic remeshing.  (Mirrors the
    // V-profile branch of compute_revolution_steps.)
    int u_steps = 1;
    Handle(Geom_Curve) basis_curve = extrusion.BasisCurve();
    if (!basis_curve.IsNull()) {
        const int num_u_samples = 40;
        double profile_len = 0.0;
        double rho_min = std::numeric_limits<double>::max();
        gp_Vec d_prev;
        GeomLProp_CLProps clp(basis_curve, 2, 1e-7);
        for (int i = 0; i <= num_u_samples; i++) {
            double u = u1 + (u2 - u1) * i / num_u_samples;
            gp_Pnt p; gp_Vec d1;
            basis_curve->D1(u, p, d1);
            if (i > 0) profile_len += 0.5 * (d1.Magnitude() + d_prev.Magnitude())
                                          * (u2 - u1) / num_u_samples;
            d_prev = d1;
            clp.SetParameter(u);
            if (clp.IsTangentDefined()) {
                double k = clp.Curvature();
                if (k > 1e-12) rho_min = std::min(rho_min, 1.0 / k);
            }
        }
        spdlog::debug("  profile arc length: {}", profile_len);
        spdlog::debug("  profile min radius of curvature: {}", rho_min);

        double u_chord = max_edge_length;
        if (rho_min < std::numeric_limits<double>::max() &&
            max_surface_error / rho_min <= 2) {
            double u_angle = 2 * std::acos(1 - max_surface_error / rho_min);
            u_chord = std::min(u_chord, 2 * rho_min * std::sin(u_angle / 2));
        }
        u_chord = std::max(u_chord, min_edge_length);
        if (u_chord > 0 && profile_len > u_chord) {
            u_steps = static_cast<int>(std::ceil(profile_len / u_chord));
            spdlog::debug("  U steps: {}", u_steps);
        }
    }

    return {u_steps, v_steps};
}

std::pair<int,int> compute_torus_steps(const Geom_ToroidalSurface& torus,
                                       const TopoDS_Face& face,
                                       double min_edge_length,
                                       double max_edge_length,
                                       double max_surface_error)
{
    int u_steps = 1;
    int v_steps = 1;

    double major_radius = torus.MajorRadius();
    double minor_radius = torus.MinorRadius();
    double outer_radius = major_radius + minor_radius;

    spdlog::debug("torus:");
    spdlog::debug("  major radius: {}", major_radius);
    spdlog::debug("  minor radius: {}", minor_radius);

    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("  U (ring angle): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (tube angle): {} -> {} ({})", v1, v2, v2 - v1);
    spdlog::debug("  max edge length: {}", max_edge_length);

    // Surface error and chord length are worst-case at the outer equator (v=0)
    // where the ring radius is R+r.  Using outer_radius is conservative.
    double max_u_edge_length = max_edge_length / std::sqrt(2.0);
    max_u_edge_length = std::min(max_u_edge_length, 2 * outer_radius);
    double min_u_edge_length = std::min(min_edge_length, 2 * outer_radius);
    spdlog::debug("  max u edge length: {}", max_u_edge_length);
    spdlog::debug("  min u edge length: {}", min_u_edge_length);

    double u_angle = 0;
    if (max_surface_error / outer_radius <= 2) {
        double max_angle = 2 * std::acos(1 - max_surface_error / outer_radius);
        spdlog::debug("  max surface error angle (U): {}", max_angle);
        u_angle = max_angle;
    }

    double min_u_angle = 2 * std::asin(min_u_edge_length / (2 * outer_radius));
    spdlog::debug("  min edge length angle (U): {}", min_u_angle);
    if (min_u_angle > u_angle) {
        u_angle = min_u_angle;
    }

    double max_u_angle = 2 * std::asin(max_u_edge_length / (2 * outer_radius));
    spdlog::debug("  max edge length angle (U): {}", max_u_angle);
    if (max_u_angle < u_angle) {
        u_angle = max_u_angle;
    }

    if (u_angle > 0) {
        int steps = std::ceil((u2 - u1) / u_angle);
        spdlog::debug("  U steps: {}", steps);
        u_angle = (u2 - u1) / steps;
        spdlog::debug("  actual U angle: {}", u_angle);
        if (steps > 1) {
            u_steps = steps;
        }
    }

    double u_chord = 2 * outer_radius * std::sin(u_angle / 2);
    spdlog::debug("  u chord: {}", u_chord);

    // Two constraints for V: diagonal budget and tube cross-section curvature.
    double max_v_chord = std::sqrt(
        std::max(0.0, max_edge_length * max_edge_length - u_chord * u_chord));
    spdlog::debug("  max v chord from diagonal: {}", max_v_chord);

    if (max_surface_error / minor_radius <= 2) {
        double v_angle_from_error = 2 * std::acos(1 - max_surface_error / minor_radius);
        double v_chord_from_error = 2 * minor_radius * std::sin(v_angle_from_error / 2);
        spdlog::debug("  max v chord from surface error: {}", v_chord_from_error);
        max_v_chord = std::min(max_v_chord, v_chord_from_error);
    }

    max_v_chord = std::min(max_v_chord, 2 * minor_radius);

    if (max_v_chord > 0) {
        double v_angle_max = 2 * std::asin(max_v_chord / (2 * minor_radius));
        spdlog::debug("  max V angle: {}", v_angle_max);
        if (v_angle_max > 0 && (v2 - v1) > v_angle_max) {
            v_steps = static_cast<int>(std::ceil((v2 - v1) / v_angle_max));
            spdlog::debug("  V steps: {}", v_steps);
        }
    }

    return {u_steps, v_steps};
}

FaceTessellationSteps get_face_tessellation_steps(const TopoDS_Face& face,
                                                   double min_edge_length,
                                                   double max_edge_length,
                                                   double max_surface_error)
{
    auto surface = BRep_Tool::Surface(face);
    FaceTessellationSteps result;

    if (auto h = Handle(Geom_CylindricalSurface)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Cylinder;
        std::tie(result.u_steps, result.v_steps) = compute_cylinder_steps(*h, face, min_edge_length,
                                                                           max_edge_length, max_surface_error);
    } else if (auto h = Handle(Geom_ConicalSurface)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Cone;
        std::tie(result.u_steps, result.v_steps) = compute_cone_steps(*h, face, min_edge_length,
                                                                       max_edge_length, max_surface_error);
    } else if (auto h = Handle(Geom_SurfaceOfRevolution)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Revolution;
        std::tie(result.u_steps, result.v_steps) = compute_revolution_steps(*h, face, min_edge_length,
                                                                             max_edge_length, max_surface_error);
    } else if (auto h = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Extrusion;
        std::tie(result.u_steps, result.v_steps) = compute_extrusion_steps(*h, face, min_edge_length,
                                                                            max_edge_length, max_surface_error);
    } else if (auto h = Handle(Geom_ToroidalSurface)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Torus;
        std::tie(result.u_steps, result.v_steps) = compute_torus_steps(*h, face, min_edge_length,
                                                                        max_edge_length, max_surface_error);
    }

    return result;
}
