#include <precision_mesh/step_face_util.h>

#include <cstdio>
#include <cmath>
#include <vector>

#include <spdlog/spdlog.h>

#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <BRepGProp.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

// -- sample_face_wire --------------------------------------------------------

std::vector<float> sample_face_wire(const TopoDS_Face& face, double deflection) {
    std::vector<float> xyz;
    for (TopExp_Explorer e_exp(face, TopAbs_EDGE); e_exp.More(); e_exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(e_exp.Current());

        // Skip degenerate edges (e.g., cone apex)
        Standard_Real first, last;
        if (BRep_Tool::Curve(edge, first, last).IsNull()) continue;
        if (last - first < 1e-15) continue;

        BRepAdaptor_Curve curve(edge);
        GCPnts_TangentialDeflection sampler(curve, 0.2, deflection, 2, 1.0e-9);

        int n = sampler.NbPoints();
        if (n < 2) continue;

        gp_Pnt prev = sampler.Value(1);
        for (int i = 2; i <= n; i++) {
            gp_Pnt curr = sampler.Value(i);
            xyz.push_back((float)prev.X()); xyz.push_back((float)prev.Y()); xyz.push_back((float)prev.Z());
            xyz.push_back((float)curr.X()); xyz.push_back((float)curr.Y()); xyz.push_back((float)curr.Z());
            prev = curr;
        }
    }
    return xyz;
}

// -- sample_free_edges --------------------------------------------------------

std::vector<float> sample_free_edges(const TopoDS_Shape& shape, double deflection) {
    std::vector<float> xyz;
    TopTools_IndexedDataMapOfShapeListOfShape edge_faces;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);

    int n_null_curve = 0, n_degenerate = 0, n_sampled = 0;
    for (int i = 1; i <= edge_faces.Extent(); i++) {
        const TopTools_ListOfShape& faces = edge_faces.FindFromIndex(i);
        if (faces.Extent() >= 2) continue;                       // shared -- not a free edge
        TopoDS_Edge edge = TopoDS::Edge(edge_faces.FindKey(i));
        if (faces.Extent() >= 1 &&
            BRep_Tool::IsClosed(edge, TopoDS::Face(faces.First()))) continue;   // periodic seam

        Standard_Real first, last;
        if (BRep_Tool::Curve(edge, first, last).IsNull()) { n_null_curve++; continue; }
        if (last - first < 1e-15) { n_degenerate++; continue; }

        BRepAdaptor_Curve curve(edge);
        GCPnts_TangentialDeflection sampler(curve, 0.2, deflection, 2, 1.0e-9);
        int n = sampler.NbPoints();
        if (n < 2) { n_degenerate++; continue; }

        n_sampled++;
        gp_Pnt prev = sampler.Value(1);
        for (int j = 2; j <= n; j++) {
            gp_Pnt curr = sampler.Value(j);
            xyz.push_back((float)prev.X()); xyz.push_back((float)prev.Y()); xyz.push_back((float)prev.Z());
            xyz.push_back((float)curr.X()); xyz.push_back((float)curr.Y()); xyz.push_back((float)curr.Z());
            prev = curr;
        }
    }
    if (n_null_curve || n_degenerate)
        spdlog::debug("    free-edge sampling: {} rendered, {} skipped (null curve), {} skipped (degenerate/short)",
                      n_sampled, n_null_curve, n_degenerate);
    return xyz;
}

// -- Helpers ------------------------------------------------------------------

static std::string fmt_g(double v, int sig = 4) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*g", sig, v);
    return buf;
}

static std::string curve_type_name(const Handle(Geom_Curve)& c) {
    if (c.IsNull()) return "none";
    GeomAdaptor_Curve a(c);
    switch (a.GetType()) {
        case GeomAbs_Line:         return "line";
        case GeomAbs_Circle:       return "circle";
        case GeomAbs_Ellipse:      return "ellipse";
        case GeomAbs_Hyperbola:    return "hyperbola";
        case GeomAbs_Parabola:     return "parabola";
        case GeomAbs_BezierCurve:  return "bezier";
        case GeomAbs_BSplineCurve: return "bspline";
        case GeomAbs_OffsetCurve:  return "offset";
        default:                   return "curve";
    }
}

// -- Public API ----------------------------------------------------------------

std::vector<TopoDS_Face> get_faces(const TopoDS_Shape& shape) {
    std::vector<TopoDS_Face> faces;
    for (TopExp_Explorer iter(shape, TopAbs_FACE); iter.More(); iter.Next()) {
        faces.push_back(TopoDS::Face(iter.Current()));
    }
    return faces;
}

int get_face_type(const TopoDS_Face& face) {
    auto surface = BRep_Tool::Surface(face);

    if (!Handle(Geom_Plane)::DownCast(surface).IsNull())             return 1;
    if (!Handle(Geom_CylindricalSurface)::DownCast(surface).IsNull()) return 2;
    if (!Handle(Geom_SphericalSurface)::DownCast(surface).IsNull())   return 3;
    if (!Handle(Geom_ToroidalSurface)::DownCast(surface).IsNull())    return 4;
    if (!Handle(Geom_ConicalSurface)::DownCast(surface).IsNull())     return 5;

    GeomAdaptor_Surface adaptor(surface);
    switch (adaptor.GetType()) {
        case GeomAbs_Plane:    return 1;
        case GeomAbs_Cylinder: return 2;
        case GeomAbs_Sphere:   return 3;
        case GeomAbs_Torus:    return 4;
        case GeomAbs_Cone:     return 5;
        default:               return 0;
    }
}

float get_face_area(const TopoDS_Face& face) {
    GProp_GProps props;
    BRepGProp::SurfaceProperties(face, props);
    return props.Mass();
}

std::array<float, 3> get_face_centroid(const TopoDS_Face& face) {
    GProp_GProps props;
    BRepGProp::SurfaceProperties(face, props);
    auto c = props.CentreOfMass();
    return { static_cast<float>(c.X()), static_cast<float>(c.Y()), static_cast<float>(c.Z()) };
}

std::string get_face_description(const TopoDS_Face& face) {
    auto surface = BRep_Tool::Surface(face);
    if (surface.IsNull()) return "None";

    GeomAdaptor_Surface a(surface);
    switch (a.GetType()) {
        case GeomAbs_Plane:
            return "Plane";

        case GeomAbs_Cylinder:
            return "Cylinder  r=" + fmt_g(a.Cylinder().Radius());

        case GeomAbs_Cone: {
            auto cone = a.Cone();
            return "Cone  half-angle=" + fmt_g(cone.SemiAngle() * 180.0 / M_PI, 3)
                   + " deg  r=" + fmt_g(cone.RefRadius());
        }

        case GeomAbs_Sphere:
            return "Sphere  r=" + fmt_g(a.Sphere().Radius());

        case GeomAbs_Torus: {
            auto t = a.Torus();
            return "Torus  R=" + fmt_g(t.MajorRadius())
                   + "  r=" + fmt_g(t.MinorRadius());
        }

        case GeomAbs_BezierSurface: {
            auto b = a.Bezier();
            return "Bezier  deg=" + std::to_string(b->UDegree())
                   + "x" + std::to_string(b->VDegree())
                   + "  " + std::to_string(b->NbUPoles())
                   + "x" + std::to_string(b->NbVPoles()) + " poles";
        }

        case GeomAbs_BSplineSurface: {
            auto b = a.BSpline();
            std::string s = "BSpline  deg=" + std::to_string(b->UDegree())
                            + "x" + std::to_string(b->VDegree())
                            + "  " + std::to_string(b->NbUPoles())
                            + "x" + std::to_string(b->NbVPoles()) + " poles";
            if (b->IsURational() || b->IsVRational()) s += "  rational";
            return s;
        }

        case GeomAbs_SurfaceOfRevolution: {
            auto rev = Handle(Geom_SurfaceOfRevolution)::DownCast(surface);
            if (!rev.IsNull())
                return "Revolution  (" + curve_type_name(rev->BasisCurve()) + " profile)";
            return "Revolution";
        }

        case GeomAbs_SurfaceOfExtrusion: {
            auto ext = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface);
            if (!ext.IsNull())
                return "Extrusion  (" + curve_type_name(ext->BasisCurve()) + " profile)";
            return "Extrusion";
        }

        case GeomAbs_OffsetSurface: {
            auto off = Handle(Geom_OffsetSurface)::DownCast(surface);
            if (!off.IsNull())
                return "Offset  d=" + fmt_g(off->Offset());
            return "Offset";
        }

        default:
            return "Other";
    }
}
