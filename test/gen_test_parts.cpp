// Generator for synthetic STEP test parts exercising specific tessellation paths.
// Usage: gen_test_parts <shape> [output.step]
// Shapes:
//   party_hat   Hollow cone shell: outer cone with apex (30 deg semi-angle), identical
//               inner cone offset down the axis, planar annular ring at the base.
//               Exercises apex-containing cone faces in the UV-grid CDT path.

#include <cmath>
#include <cstdio>
#include <map>
#include <string>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>

namespace {

TopoDS_Shape make_party_hat()
{
    const double semi_angle = 30.0 * M_PI / 180.0;
    const double base_radius = 30.0;                            // mm
    const double height = base_radius / std::tan(semi_angle);   // apex at z = height
    const double axial_offset = 2.0;                            // wall: inner cone down by this

    gp_Ax2 axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    TopoDS_Shape outer = BRepPrimAPI_MakeCone(axis, base_radius, 0.0, height).Shape();

    gp_Trsf down;
    down.SetTranslation(gp_Vec(0, 0, -axial_offset));
    TopoDS_Shape inner = BRepBuilderAPI_Transform(outer, down, /*copy=*/true).Shape();

    // outer minus inner: outer cone face (apex at z=height), inner cone face (apex at
    // z=height-offset), and the base annulus ring between them at z=0.
    return BRepAlgoAPI_Cut(outer, inner).Shape();
}

}  // namespace

int main(int argc, char** argv)
{
    std::map<std::string, TopoDS_Shape (*)()> shapes = {
        {"party_hat", make_party_hat},
    };

    if (argc < 2 || !shapes.count(argv[1])) {
        std::fprintf(stderr, "usage: %s <shape> [output.step]\nshapes:\n", argv[0]);
        for (const auto& [name, fn] : shapes) std::fprintf(stderr, "  %s\n", name.c_str());
        return 1;
    }
    std::string name = argv[1];
    std::string output = argc > 2 ? argv[2] : name + ".step";

    TopoDS_Shape shape = shapes[name]();
    if (shape.IsNull()) {
        std::fprintf(stderr, "failed to build shape '%s'\n", name.c_str());
        return 1;
    }

    STEPControl_Writer writer;
    if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone ||
        writer.Write(output.c_str()) != IFSelect_RetDone) {
        std::fprintf(stderr, "failed to write '%s'\n", output.c_str());
        return 1;
    }
    std::printf("wrote %s: %s\n", name.c_str(), output.c_str());
    return 0;
}
