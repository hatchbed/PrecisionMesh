#pragma once

#include <cstddef>
#include <limits>
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

// Per-face map from STEP vertex hash -> wire projector, shared across all segments.
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

// Per-call statistics collected by project_to_step; merge across segments with merge().
struct ProjectionStats {
    struct Record {
        double delta        = 0.0;
        double px = 0, py = 0, pz = 0;   // input position
        double qx = 0, qy = 0, qz = 0;   // projected position
        bool   is_border    = false;
        bool   null_he      = false;
        bool   on_real_edge = false;
        double edge_dist    = -1.0;
        size_t vert_idx     = 0;
        size_t segment_index = 0;
        int    face_code    = 0;
    };

    static constexpr size_t TOP_N = 20;

    size_t n_total             = 0;
    size_t n_skipped_null      = 0;
    size_t n_border_skipped    = 0;  // border vertices: static, not projected
    size_t n_surface_projected = 0;
    double sum_delta           = 0.0;
    double min_delta           = std::numeric_limits<double>::max();
    double max_delta           = 0.0;
    size_t n_above_1em3        = 0;  // delta > 1e-3
    size_t n_above_1em2        = 0;  // delta > 1e-2
    size_t n_above_1em1        = 0;  // delta > 1e-1
    size_t n_occt_skip         = 0;  // surface projector returned a face-boundary point (skipped)
    std::vector<Record> top_records;  // sorted descending by delta, capped at TOP_N

    void merge(const ProjectionStats& other);
};

// Project all mesh vertices toward the STEP surface, blending by weight [0,1].
// If stats_out is non-null, per-vertex statistics are written there (no printing).
// segment_index is stored in records for identification after merging.
// max_projection_delta: if > 0, surface-projected vertices whose nearest point on the
// face is a trimming-boundary edge/vertex AND whose displacement exceeds this threshold
// are skipped (OCCT extrema wrong-local-minimum guard).  Pass params.max_edge_length.
void project_to_step(const TopoDS_Face& face, Mesh& mesh,
                     WireProjectorCachePtr wire_projectors,
                     StepProjector& surface_projector,
                     StepBorderProjector& border_projector,
                     double weight = 1.0,
                     size_t segment_index = 0,
                     ProjectionStats* stats_out = nullptr,
                     double max_projection_delta = 0.0);

// Returns the distance from point (x,y,z) to the nearest point on face, or -1 on failure.
double get_distance_to_face(const TopoDS_Face& face, double x, double y, double z);

// Result of validate_tessellation (all distances in model units).
struct TessellationValidation {
    size_t segments        = 0;
    size_t total_tris      = 0;
    size_t border_verts    = 0;   // mesh-border vertices (edge-origin)
    size_t interior_verts  = 0;   // mesh-interior vertices (surface-origin)
    double max_border_edge_dist    = 0.0;  // worst border vertex -> nearest STEP edge
    double max_interior_face_dist  = 0.0;  // worst interior vertex -> its STEP face
    size_t misclassified_border    = 0;    // border verts farther than tol from any edge
    size_t misclassified_interior  = 0;    // interior verts farther than tol from face
    double max_surface_error  = 0.0;  // worst triangle-interior sample -> BREP face
    double mean_surface_error = 0.0;  // mean over all samples
    size_t surface_samples    = 0;    // 0 when triangle sampling was disabled
    size_t open_boundary_edges = 0;   // border edges of the merged mesh (cracks / true boundary)
    size_t non_manifold_edges  = 0;   // undirected edges incident to >2 triangles
    size_t flipped_triangles   = 0;   // triangles whose winding opposes the BREP face normal
};

// Validate the tessellation against the BREP: per vertex, check that border vertices lie
// on a STEP edge and interior vertices on the STEP face (invariants 2 & 3); sample each
// triangle's interior (centroid if samples_per_tri >= 1; + edge midpoints if >= 4) and
// measure the surface error to the BREP face; count open boundary edges of the merged
// mesh (watertightness); and check each triangle's winding against the BREP face normal
// at its centroid (flipped_triangles -- orientation bugs are invisible to the
// position/watertightness metrics but break shading and inside/outside tests).
// samples_per_tri == 0 skips the per-triangle surface-error sampling entirely (it
// dominates the cost on large meshes); vertex placement, watertightness, and winding
// are always checked.  Intended as an opt-in validation pass -- BRep distance queries
// make it slow on large meshes.
// `edge_faces` is parallel to `segments`: the ORIGINAL (pre-subdivision) face each segment
// came from.  Vertices are classified against that face's REAL edges (seam edges excluded) --
// not the segment mesh's border -- so periodic seams and subdivision cuts (which are interior
// to the continuous surface, not true part boundaries) are validated against the surface, not
// against a boundary.  Pass `segments` itself when there was no subdivision.
TessellationValidation validate_tessellation(const std::vector<Mesh>& meshes,
                                             const std::vector<TopoDS_Face>& segments,
                                             const std::vector<TopoDS_Face>& edge_faces,
                                             double tolerance,
                                             int samples_per_tri = 4);
