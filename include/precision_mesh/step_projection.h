#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

// mesh_util.h (CGAL) must come before any OCCT header that defines the Handle macro.
#include <precision_mesh/mesh_util.h>

#include <BRepExtrema_DistShapeShape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

#include <precision_mesh/step_types.h>

// Projects mesh vertices onto the nearest point of a STEP shape.
struct StepProjector {
    StepProjector();
    StepProjector(const TopoDS_Shape& shape);
    void setShape(const TopoDS_Shape& shape);
    Mesh::Point operator()(const Mesh::Point& p);
    BRepExtrema_DistShapeShape extrema;
};

// Projects mesh border vertices onto the wire boundary of a STEP face.
struct StepBorderProjector {
    StepBorderProjector();
    StepBorderProjector(const TopoDS_Face& face);
    void setFace(const TopoDS_Shape& face);
    Mesh::Point operator()(const Mesh::Point& p);
    TopoDS_Compound border;
    BRepExtrema_DistShapeShape extrema;
};

// Per-face map from STEP vertex hash → wire projector, shared across all segments.
using WireProjectorCachePtr =
    std::shared_ptr<std::unordered_map<int, std::unordered_map<int, StepProjector>>>;

// Build a wire from all boundary edges reachable from a given STEP vertex, skipping
// seam and degenerate edges which are not physical boundaries.
TopoDS_Wire get_border_loop_wire(
    const TopoDS_Vertex& vertex,
    const TopTools_IndexedDataMapOfShapeListOfShape& vertex_to_edge_map,
    const TopTools_IndexedDataMapOfShapeListOfShape& edge_to_face_map);

// Build a per-face cache of wire projectors, one per STEP vertex, in parallel.
WireProjectorCachePtr get_edge_vertex_wire_projectors(const TopoDS_Shape& shape);

// Associate each border vertex of mesh with the wire projector for its nearest
// STEP edge boundary, propagating assignments via already-assigned neighbors.
std::unordered_map<Mesh::Vertex_index, StepProjector>
get_border_vertex_projector_map(const TopoDS_Face& face, Mesh& mesh,
                                WireProjectorCachePtr wire_projectors,
                                double tolerance = 1e-3);

// Project all mesh vertices toward the STEP surface, blending by weight [0,1].
void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr wire_projectors,
                     StepProjector& surface_projector,
                     StepBorderProjector& border_projector,
                     double weight = 1.0);

// Returns the distance from point (x,y,z) to the nearest point on face, or -1 on failure.
double get_distance_to_face(const TopoDS_Face& face, double x, double y, double z);
