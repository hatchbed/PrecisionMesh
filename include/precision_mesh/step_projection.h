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

// Result of validate_tessellation (all distances in model units).
struct TessellationValidation {
    size_t segments        = 0;
    size_t total_tris      = 0;
    size_t border_verts    = 0;   // mesh-border vertices (edge-origin)
    size_t interior_verts  = 0;   // mesh-interior vertices (surface-origin)
    double max_border_edge_dist    = 0.0;  // worst border vertex → nearest STEP edge
    double max_interior_face_dist  = 0.0;  // worst interior vertex → its STEP face
    size_t misclassified_border    = 0;    // border verts farther than tol from any edge
    size_t misclassified_interior  = 0;    // interior verts farther than tol from face
    double max_surface_error  = 0.0;  // worst triangle-interior sample → BREP face
    double mean_surface_error = 0.0;  // mean over all samples
    size_t surface_samples    = 0;
    size_t open_boundary_edges = 0;   // border edges of the merged mesh (cracks / true boundary)
};

// Validate the tessellation against the BREP: per vertex, check that border vertices lie
// on a STEP edge and interior vertices on the STEP face (invariants 2 & 3); sample each
// triangle's interior (centroid; + edge midpoints if samples_per_tri >= 4) and measure the
// surface error to the BREP face; and count open boundary edges of the merged mesh
// (watertightness).  Intended as an opt-in validation pass — BRep distance queries make it
// slow on large meshes.
// `edge_faces` is parallel to `segments`: the ORIGINAL (pre-subdivision) face each segment
// came from.  Vertices are classified against that face's REAL edges (seam edges excluded) —
// not the segment mesh's border — so periodic seams and subdivision cuts (which are interior
// to the continuous surface, not true part boundaries) are validated against the surface, not
// against a boundary.  Pass `segments` itself when there was no subdivision.
TessellationValidation validate_tessellation(const std::vector<Mesh>& meshes,
                                             const std::vector<TopoDS_Face>& segments,
                                             const std::vector<TopoDS_Face>& edge_faces,
                                             double tolerance,
                                             int samples_per_tri = 4);
