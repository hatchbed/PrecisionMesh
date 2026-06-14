#pragma once

#include <memory>
#include <string>
#include <vector>
#include <cstddef>

// NOTE: No GL/GLFW/ImGui/OCCT headers here — pimpl pattern keeps them
// isolated from headers that include OCCT, which defines Handle() as a macro
// that corrupts third-party template names.
#include <precision_mesh/mesh_util.h>

struct ViewerSegment {
    const Mesh* mesh;
    int         face_type;       // 0=Other 1=Plane 2=Cylinder 3=Sphere 4=Torus 5=Cone
    int         segment_index;   // index of the CAD face (one segment per face)
    float       area;            // surface area in model units
    std::string surface_desc;    // human-readable surface description
    std::string tess_approach;   // tessellation method used for this segment
    // Pre-tessellation BRep wire of this face's boundary for the hover overlay
    // (GL_LINES xyz pairs, sampled from the STEP geometry, not the mesh).
    std::vector<float> face_wire;
};

class Viewer {
public:
    Viewer(int width, int height, const std::string& title);
    ~Viewer();

    Viewer(const Viewer&) = delete;
    Viewer& operator=(const Viewer&) = delete;

    // Called from the worker thread — non-blocking.  Converts meshes to flat
    // GPU data and posts to a pending slot.  The render loop uploads it on the
    // next frame.
    // `free_brep_edges` are GL_LINES xyz pairs for the input shape's free boundary edges
    // (drawn bold yellow).  Open edges of the merged triangle soup (the real cracks) are
    // computed from the segment meshes and drawn bold magenta.  `healed_edges` are
    // GL_LINES xyz pairs for edges that were open before loop repair and are now closed
    // (drawn bold green).
    void update(const std::vector<ViewerSegment>& segments,
                const std::string& stage_label,
                const std::vector<float>& free_brep_edges = {},
                const std::vector<float>& healed_edges = {});

    // Called by the worker thread when all processing is complete.  After this,
    // the next window-close event exits the event loop normally.
    void notify_done();

    // Main-thread event loop.  Blocks until processing is done AND the window
    // has been closed by the user.  Must be called from the main thread.
    void run_event_loop();

    // Show the window if it was hidden.  Thread-safe — sets an atomic flag
    // that the render loop checks on the next iteration.
    void show();

    // Public forward declaration so GLFW callbacks (non-member functions) can
    // name the type.  The full definition is in viewer.cpp.
    struct Impl;

private:
    std::unique_ptr<Impl> impl_;
};
