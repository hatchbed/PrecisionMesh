#pragma once

#include <climits>

#include <Standard_Version.hxx>
#include <TopoDS_Shape.hxx>

// HashCode(upper) was removed in OCCT 7.8; use std::hash instead.
#if OCC_VERSION_HEX >= 0x070800
inline int shapeHashCode(const TopoDS_Shape& s) { return static_cast<int>(std::hash<TopoDS_Shape>{}(s)); }
#else
inline int shapeHashCode(const TopoDS_Shape& s) { return s.HashCode(INT_MAX); }
#endif
