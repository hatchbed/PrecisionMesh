#pragma once

#include <array>
#include <string>
#include <vector>

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

std::vector<TopoDS_Face> get_faces(const TopoDS_Shape& shape);
int get_face_type(const TopoDS_Face& face);
float get_face_area(const TopoDS_Face& face);
std::array<float, 3> get_face_centroid(const TopoDS_Face& face);

// Returns a compact human-readable description of the surface geometry,
// e.g. "Cylinder  r=1.234", "BSpline  deg=3×3  8×6 poles  rational".
std::string get_face_description(const TopoDS_Face& face);

// Sample the boundary wires of a TopoDS_Face into 3D line segments for GL_LINES
// rendering (interleaved x,y,z endpoint pairs).  deflection is the maximum chord
// height used for curved-edge approximation (same units as the model).
std::vector<float> sample_face_wire(const TopoDS_Face& face, double deflection);

// Sample the FREE boundary edges of a shape (edges adjacent to <2 faces, excluding
// periodic seam edges) into 3D line segments for GL_LINES rendering.  These are the
// real open boundaries / un-stitched edges of the BREP topology.
std::vector<float> sample_free_edges(const TopoDS_Shape& shape, double deflection);
