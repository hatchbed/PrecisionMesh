#pragma once

#include <climits>
#include <unordered_map>

#include <Standard_Version.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

// HashCode(upper) was removed in OCCT 7.8; use std::hash instead.
#if OCC_VERSION_HEX >= 0x070800
inline int shapeHashCode(const TopoDS_Shape& s) { return static_cast<int>(std::hash<TopoDS_Shape>{}(s)); }
#else
inline int shapeHashCode(const TopoDS_Shape& s) { return s.HashCode(INT_MAX); }
#endif

struct FaceHasher {
    std::size_t operator()(const TopoDS_Face& face) const {
        return static_cast<std::size_t>(shapeHashCode(face));
    }
};

struct FaceEqual {
    bool operator()(const TopoDS_Face& a, const TopoDS_Face& b) const {
        return a.IsSame(b);
    }
};

typedef std::unordered_map<TopoDS_Face, int, FaceHasher, FaceEqual> FaceMap;
