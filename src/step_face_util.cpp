#include <precision_mesh/step_face_util.h>

#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>

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
