#include <precision_mesh/step_subdivision.h>
#include <precision_mesh/step_tessellation.h>

#include <cmath>
#include <vector>

#include <spdlog/spdlog.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepTools.hxx>
#include <Geom_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <ShapeFix_Edge.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>

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

    return {u_steps, 1};
}

std::pair<int,int> compute_extrusion_steps(const Geom_SurfaceOfLinearExtrusion& /*extrusion*/,
                                           const TopoDS_Face& face,
                                           double max_edge_length)
{
    int v_steps = 1;

    Standard_Real u1, u2, v1, v2;
    BRepTools::UVBounds(face, u1, u2, v1, v2);

    spdlog::debug("surface of linear extrusion:");
    spdlog::debug("  U (profile): {} -> {} ({})", u1, u2, u2 - u1);
    spdlog::debug("  V (extrusion): {} -> {} ({})", v1, v2, v2 - v1);

    // V is the extrusion direction — always a straight line with zero surface error.
    // Subdivide purely by max_edge_length.
    if ((v2 - v1) > max_edge_length) {
        v_steps = static_cast<int>(std::ceil((v2 - v1) / max_edge_length));
        spdlog::debug("  V steps: {}", v_steps);
    }

    return {1, v_steps};
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

std::vector<TopoDS_Face> subdivide_face(const TopoDS_Face& face, int u_steps, int v_steps)
{
    if (u_steps < 2 && v_steps < 2) {
        return { face };
    }

    Standard_Real u_first, u_last, v_first, v_last;
    BRepTools::UVBounds(face, u_first, u_last, v_first, v_last);

    auto surface = BRep_Tool::Surface(face);

    double u_range = u_last - u_first;
    double v_range = v_last - v_first;

    double u_step_size = u_range / u_steps;
    double v_step_size = v_range / v_steps;

    spdlog::debug("subdividing face ({}-{}[{}], {}-{}[{}])", u_first, u_last, u_steps,
                  v_first, v_last, v_steps);

    ShapeFix_Edge edge_fix;

    // Cut extent in each direction: extend a little past the face bounds so the cut fully
    // crosses non-periodic boundaries.  In a PERIODIC direction, however, the cut span must
    // never reach a full period: on a face that wraps (or nearly wraps) the whole period, an
    // extended cut runs past the seam and overlaps itself, and the splitter then leaves the
    // seam-adjacent cells un-split (collapsing a whole band into one giant cell).  Cap the
    // extension by the gap to the period.  This subsumes every case with no special-casing:
    //  - full-period face (gap 0)        -> extension 0 -> span [first,last], closes cleanly;
    //  - almost-full face (tiny gap)     -> tiny extension, span stays just under a period;
    //  - small arc (large gap)           -> full 0.01 extension as before.
    // It also handles faces that span ~2pi-eps with NO seam edge (a STEP 359.98-degree
    // cylinder), which neither a UPeriod-tolerance test nor seam detection catches.
    double u_ext = 0.01, v_ext = 0.01;
    if (!surface.IsNull() && surface->IsUPeriodic())
        u_ext = std::min(0.01, std::max(0.0, (surface->UPeriod() - u_range) * 0.49));
    if (!surface.IsNull() && surface->IsVPeriodic())
        v_ext = std::min(0.01, std::max(0.0, (surface->VPeriod() - v_range) * 0.49));
    // U-cut lines (constant U) span V; V-cut lines (constant V) span U.
    double v_span_lo = v_first - v_ext;
    double v_span_len = v_range + 2 * v_ext;
    double u_span_lo = u_first - u_ext;
    double u_span_len = u_range + 2 * u_ext;

    TopTools_ListOfShape cut_tools;
    for (int u_step = 1; u_step < u_steps; u_step++) {
        double u_val = u_first + u_step * u_step_size;
        auto v_line = new Geom2d_Line(gp_Pnt2d(u_val, v_span_lo), gp_Dir2d(0, 1));
        TopoDS_Edge v_edge = BRepBuilderAPI_MakeEdge(v_line, surface, 0, v_span_len);
        edge_fix.FixAddCurve3d(v_edge);
        cut_tools.Append(v_edge);
    }
    for (int v_step = 1; v_step < v_steps; v_step++) {
        double v_val = v_first + v_step * v_step_size;
        auto u_line = new Geom2d_Line(gp_Pnt2d(u_span_lo, v_val), gp_Dir2d(1, 0));
        TopoDS_Edge u_edge = BRepBuilderAPI_MakeEdge(u_line, surface, 0, u_span_len);
        edge_fix.FixAddCurve3d(u_edge);
        cut_tools.Append(u_edge);
    }

    TopTools_ListOfShape cut_args;
    cut_args.Append(face);

    BRepAlgoAPI_Splitter splitter;
    splitter.SetNonDestructive(true);
    splitter.SetRunParallel(true);
    splitter.SetArguments(cut_args);
    splitter.SetTools(cut_tools);
    splitter.Build();

    const TopAbs_Orientation orig_orient = face.Orientation();
    std::vector<TopoDS_Face> subdivs;
    auto modified = splitter.Modified(face);
    for (const auto& shape: modified) {
        if (shape.ShapeType() == TopAbs_FACE) {
            TopoDS_Face subface = TopoDS::Face(shape);
            if (subface.Orientation() != orig_orient) {
                subface = TopoDS::Face(subface.Reversed());
            }
            subdivs.push_back(subface);
        }
        else {
            spdlog::warn("Modified shape is not a face.");
        }
    }

    // Note: subdivs.size() < u_steps*v_steps is normal for trimmed faces — grid cells
    // outside the trimmed UV region produce no sub-face.  Logged at debug only.
    spdlog::debug("subdivs created: {} (full-grid would be {}, u_steps={} v_steps={})",
                  subdivs.size(), u_steps * v_steps, u_steps, v_steps);

    return subdivs;
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
        std::tie(result.u_steps, result.v_steps) = compute_extrusion_steps(*h, face, max_edge_length);
    } else if (auto h = Handle(Geom_ToroidalSurface)::DownCast(surface); !h.IsNull()) {
        result.type = CurvedFaceType::Torus;
        std::tie(result.u_steps, result.v_steps) = compute_torus_steps(*h, face, min_edge_length,
                                                                        max_edge_length, max_surface_error);
    }

    return result;
}

std::tuple<TopoDS_Shape, FaceMap> subdivide_step_shape(TopoDS_Shape& shape,
                                                        double min_edge_length,
                                                        double max_edge_length,
                                                        double max_surface_error)
{
    BRep_Builder builder;
    TopoDS_Compound new_shape;
    builder.MakeCompound(new_shape);

    FaceMap pre_sewn_face_map;

    int original_face_id = 0;
    for (TopExp_Explorer iter(shape, TopAbs_FACE); iter.More(); iter.Next()) {
        TopoDS_Face face = TopoDS::Face(iter.Current());

        auto steps = get_face_tessellation_steps(face, min_edge_length, max_edge_length, max_surface_error);
        int u_steps = steps.u_steps;
        int v_steps = steps.v_steps;

        // CDT-eligible face types bypass the BRepAlgoAPI_Splitter subdivision path.
        // Their interior is filled by uv_grid_retessellate() after tessellation.
        if (cdt_eligible(steps)) {
            u_steps = 1;
            v_steps = 1;
        }

        auto subdivs = subdivide_face(face, u_steps, v_steps);
        spdlog::debug("face {}: {} sub-faces, u_steps={} v_steps={} face_orient={}",
                      original_face_id, subdivs.size(), u_steps, v_steps, (int)face.Orientation());
        for (const auto& subdiv: subdivs) {
            builder.Add(new_shape, subdiv);
            pre_sewn_face_map[subdiv] = original_face_id;
        }

        original_face_id++;
    }

    // Sewing tolerance must scale with the model (geometry is in native STEP units here).
    // A hardcoded 0.1 was unit-agnostic (0.1 m for a metre-unit file) and, even in mm, was
    // coarse enough to mis-merge tiny tip-region edges and break cross-strip adjacency,
    // which fed BRepMesh inconsistent shared-edge nodes and produced twisted tessellations.
    double sew_tolerance = min_edge_length * 0.01;
    spdlog::debug("sewing tolerance: {} (native units)", sew_tolerance);
    BRepBuilderAPI_Sewing sewing;
    sewing.SetTolerance(sew_tolerance);
    sewing.Add(new_shape);
    sewing.Perform();

    auto sewed = sewing.SewedShape();

    FaceMap face_map;
    for (const auto& [face, id] : pre_sewn_face_map) {
        TopoDS_Shape modified = sewing.ModifiedSubShape(face);
        if (!modified.IsNull()) {
            for (TopExp_Explorer face_exp(modified, TopAbs_FACE); face_exp.More(); face_exp.Next()) {
                face_map[TopoDS::Face(face_exp.Current())] = id;
            }
        }
    }

    return std::make_tuple(sewed, face_map);
}
