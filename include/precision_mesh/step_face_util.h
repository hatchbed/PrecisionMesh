#pragma once

#include <array>
#include <vector>

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

std::vector<TopoDS_Face> get_faces(const TopoDS_Shape& shape);
int get_face_type(const TopoDS_Face& face);
float get_face_area(const TopoDS_Face& face);
std::array<float, 3> get_face_centroid(const TopoDS_Face& face);
