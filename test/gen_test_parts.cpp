// Generator for synthetic STEP test parts exercising specific tessellation paths.
// Usage: gen_test_parts <shape> [output.step]
// Shapes:
//   party_hat     Hollow cone shell: outer cone with apex (30 deg semi-angle), identical
//                 inner cone offset down the axis, planar annular ring at the base.
//                 Exercises apex-containing cone faces in the UV-grid CDT path.
//   torus_quarter Quarter-tube torus (R=20mm, r=5mm, v in [0, pi/2]), full azimuthal
//                 revolution.  Exercises torus CDT (varying row_radius, u-seam).
//   torus_half    Half-tube torus (v in [-pi/2, pi/2]), symmetric about outer equator,
//                 full azimuthal revolution.  Wider v span than torus_quarter.
//   donut         Complete torus (both u and v full-period).  Fully-closed face —
//                 CDT skips (no border loops), falls back to BRepMesh.
//   donut_c       C-shaped torus: half azimuthal revolution (u in [0, pi]), full tube
//                 cross-section (v full-period).  CDT skips (v-seam not yet implemented),
//                 falls back to BRepMesh.
//   cube          20mm box, all planar faces.  Exercises the all-planar early-return
//                 path in find_surface_error_param (no arc-forced nodes to calibrate).


#include <cmath>
#include <cstdio>
#include <map>
#include <string>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
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

TopoDS_Shape make_torus_quarter()
{
    // Quarter-tube torus: R=20mm (tube center to axis), r=5mm (tube radius).
    // v in [0, π/2] — outer top quarter of the tube cross-section.
    // Full azimuthal revolution (u in [0, 2π]) exercises the u-seam CDT path.
    // row_radius varies from R+r (v=0, outer equator) to R (v=π/2, tube top).
    return BRepPrimAPI_MakeTorus(20.0, 5.0, 0.0, M_PI / 2.0).Shape();
}

TopoDS_Shape make_torus_half()
{
    // Half-tube torus: R=20mm, r=5mm, v in [-π/2, π/2].
    // Symmetric about the outer equator (v=0): row_radius peaks at R+r on the
    // equator and falls to R at both edges — exercises the full row_radius variation.
    // Full azimuthal revolution (u in [0, 2π]) — u-seam CDT path.
    return BRepPrimAPI_MakeTorus(20.0, 5.0, -M_PI / 2.0, M_PI / 2.0).Shape();
}

TopoDS_Shape make_donut()
{
    // Complete torus: R=20mm, r=5mm.  Both u and v are full-period — the BREP
    // face has no border loops.  CDT skips (fully-closed face), BRepMesh interior.
    return BRepPrimAPI_MakeTorus(20.0, 5.0).Shape();
}

TopoDS_Shape make_donut_c()
{
    // C-shaped torus: R=20mm, r=5mm, u in [0, π] (half azimuthal revolution).
    // The tube cross-section is fully closed (v full-period), so CDT skips the
    // lateral face (v-seam not yet implemented) and falls back to BRepMesh.
    // The two planar semi-disk caps at u=0 and u=π are non-CDT faces.
    return BRepPrimAPI_MakeTorus(20.0, 5.0, M_PI).Shape();
}

TopoDS_Shape make_cube()
{
    // 20mm box, all planar faces.  No arc-forced nodes — exercises the early-return
    // path in find_surface_error_param (no non-planar faces to calibrate against).
    return BRepPrimAPI_MakeBox(20.0, 20.0, 20.0).Shape();
}

}  // namespace

int main(int argc, char** argv)
{
    std::map<std::string, TopoDS_Shape (*)()> shapes = {
        {"party_hat",     make_party_hat},
        {"torus_quarter", make_torus_quarter},
        {"torus_half",    make_torus_half},
        {"donut",         make_donut},
        {"donut_c",       make_donut_c},
        {"cube",          make_cube},
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
