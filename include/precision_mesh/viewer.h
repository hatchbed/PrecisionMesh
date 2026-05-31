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
    int         face_type;           // 0=Other 1=Plane 2=Cylinder 3=Sphere 4=Torus 5=Cone
    int         segment_index;       // index within the subdivided tessellation
    int         original_face_index; // index of the original un-subdivided CAD face
    int         subdiv_count;        // total number of tessellation segments for this CAD face
    float       area;                // surface area in model units
    std::string surface_desc;        // human-readable surface description
    // Pre-tessellation BRep wire samples for hover overlays (GL_LINES xyz pairs).
    // Sampled directly from the STEP geometry, not derived from the mesh.
    std::vector<float> subseg_wire;    // boundary of this sub-face (red hover)
    std::vector<float> orig_face_wire; // boundary of the original un-subdivided face (green hover)
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
    void update(const std::vector<ViewerSegment>& segments,
                const std::string& stage_label);

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
