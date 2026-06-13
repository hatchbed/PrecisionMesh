// GLEW must come before any other GL headers.
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <precision_mesh/viewer.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef __unix__
#include <sys/select.h>
#include <unistd.h>
#endif

#include <Eigen/Dense>

#include <spdlog/spdlog.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ── Vertex layout (9 floats = 36 bytes) ─────────────────────────────────────
struct GpuVertex {
    float x, y, z;
    float nx, ny, nz;
    float face_type;   // 0–5 (CAD primitive type)
    float orig_mod;    // original_face_index % 12  (for "By Segment" = original CAD face)
    float subdiv_mod;  // segment_index % 12         (for "By Subdivision" = tessellation sub-face)
};
static_assert(sizeof(GpuVertex) == 36, "GpuVertex layout changed");

// ── Color palettes ───────────────────────────────────────────────────────────
static const float kFaceTypeColors[6][3] = {
    {0.60f, 0.60f, 0.60f},  // 0 other    — gray
    {0.40f, 0.60f, 0.85f},  // 1 plane    — steel blue
    {0.30f, 0.75f, 0.45f},  // 2 cylinder — forest green
    {1.00f, 0.65f, 0.20f},  // 3 sphere   — warm orange
    {0.75f, 0.40f, 0.90f},  // 4 torus    — purple
    {0.90f, 0.40f, 0.40f},  // 5 cone     — coral red
};
static const float kSegColors[12][3] = {
    {0.92f, 0.26f, 0.21f}, {0.91f, 0.46f, 0.13f}, {1.00f, 0.76f, 0.03f},
    {0.30f, 0.69f, 0.31f}, {0.13f, 0.59f, 0.95f}, {0.61f, 0.15f, 0.69f},
    {0.00f, 0.74f, 0.83f}, {0.55f, 0.76f, 0.29f}, {1.00f, 0.60f, 0.00f},
    {0.40f, 0.23f, 0.72f}, {0.90f, 0.18f, 0.49f}, {0.47f, 0.33f, 0.28f},
};

// ── Shader sources ───────────────────────────────────────────────────────────
static const char* kTriVertSrc = R"glsl(
#version 330 core
layout(location=0) in vec3  a_pos;
layout(location=1) in vec3  a_normal;
layout(location=2) in float a_face_type;
layout(location=3) in float a_orig_mod;
layout(location=4) in float a_subdiv_mod;

uniform mat4 u_mvp;
uniform mat3 u_normal_mat;

out vec3  v_normal_view;
out float v_face_type;
out float v_orig_mod;
out float v_subdiv_mod;

void main() {
    gl_Position   = u_mvp * vec4(a_pos, 1.0);
    v_normal_view = normalize(u_normal_mat * a_normal);
    v_face_type   = a_face_type;
    v_orig_mod    = a_orig_mod;
    v_subdiv_mod  = a_subdiv_mod;
}
)glsl";

static const char* kTriFragSrc = R"glsl(
#version 330 core
in vec3  v_normal_view;
in float v_face_type;
in float v_orig_mod;
in float v_subdiv_mod;

// color_mode: 0=uniform  1=face_type  2=by_segment(orig CAD face)  3=by_subdivision
uniform int   u_color_mode;
uniform vec3  u_light_dir;
uniform float u_ambient;
uniform vec4  u_override_color; // alpha>0.5 → use this color (wireframe edge pass)
uniform vec3  u_face_type_colors[6];
uniform vec3  u_palette[12];
uniform vec3  u_uniform_color;

out vec4 frag_color;

void main() {
    vec3 base;
    if (u_color_mode == 1)
        base = u_face_type_colors[clamp(int(v_face_type), 0, 5)];
    else if (u_color_mode == 2)
        base = u_palette[clamp(int(v_orig_mod),   0, 11)];
    else if (u_color_mode == 3)
        base = u_palette[clamp(int(v_subdiv_mod), 0, 11)];
    else
        base = u_uniform_color;

    float diff  = max(dot(v_normal_view, u_light_dir), 0.0);
    float light = u_ambient + (1.0 - u_ambient) * diff;
    vec3  col   = base * light;

    if (u_override_color.a > 0.5)
        col = u_override_color.rgb;

    frag_color = vec4(col, 1.0);
}
)glsl";

static const char* kLineVertSrc = R"glsl(
#version 330 core
layout(location=0) in vec3 a_pos;
uniform mat4 u_mvp;
void main() { gl_Position = u_mvp * vec4(a_pos, 1.0); }
)glsl";

static const char* kLineFragSrc = R"glsl(
#version 330 core
uniform vec4 u_color;
out vec4 frag_color;
void main() { frag_color = u_color; }
)glsl";

// ── Shader compilation helpers ───────────────────────────────────────────────
static GLuint compile_shader(GLenum type, const char* src) {
    GLuint sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    GLint ok = 0;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[512];
        glGetShaderInfoLog(sh, sizeof(log), nullptr, log);
        spdlog::error("Shader compile error: {}", log);
    }
    return sh;
}

static GLuint link_program(const char* vert, const char* frag) {
    GLuint vs = compile_shader(GL_VERTEX_SHADER,   vert);
    GLuint fs = compile_shader(GL_FRAGMENT_SHADER, frag);
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[512];
        glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        spdlog::error("Program link error: {}", log);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

// ── Camera math (Eigen, no glm) ───────────────────────────────────────────────
static Eigen::Matrix4f make_perspective(float fov_y, float aspect, float near, float far) {
    float f = 1.0f / std::tan(fov_y * 0.5f);
    Eigen::Matrix4f m = Eigen::Matrix4f::Zero();
    m(0,0) =  f / aspect;
    m(1,1) =  f;
    m(2,2) = (far + near) / (near - far);
    m(2,3) = (2.0f * far * near) / (near - far);
    m(3,2) = -1.0f;
    return m;
}

static Eigen::Matrix4f make_look_at(Eigen::Vector3f eye, Eigen::Vector3f center,
                                     Eigen::Vector3f up) {
    Eigen::Vector3f f = (center - eye).normalized();
    Eigen::Vector3f r = f.cross(up).normalized();
    Eigen::Vector3f u = r.cross(f);
    Eigen::Matrix4f m = Eigen::Matrix4f::Identity();
    m(0,0)= r.x(); m(0,1)= r.y(); m(0,2)= r.z(); m(0,3)= -r.dot(eye);
    m(1,0)= u.x(); m(1,1)= u.y(); m(1,2)= u.z(); m(1,3)= -u.dot(eye);
    m(2,0)=-f.x(); m(2,1)=-f.y(); m(2,2)=-f.z(); m(2,3)=  f.dot(eye);
    m(3,0)= 0;     m(3,1)= 0;     m(3,2)= 0;     m(3,3)=  1;
    return m;
}

// ── Impl definition ───────────────────────────────────────────────────────────
struct Viewer::Impl {
    GLFWwindow* window = nullptr;

    // Shader programs
    GLuint tri_prog  = 0;
    GLuint line_prog = 0;

    // Triangle mesh GL objects
    GLuint tri_vao = 0, tri_vbo = 0;
    GLsizei tri_vert_count = 0;

    // Border edge GL objects
    GLuint border_vao = 0, border_vbo = 0;
    GLsizei border_vert_count = 0;

    // Free BREP boundary edges (bold yellow), open triangle-soup edges (bold magenta),
    // non-manifold edges (bold orange), and healed edges (bold green).
    GLuint free_edge_vao = 0, free_edge_vbo = 0;
    GLsizei free_edge_vert_count = 0;
    GLuint open_edge_vao = 0, open_edge_vbo = 0;
    GLsizei open_edge_vert_count = 0;
    GLuint nm_edge_vao = 0, nm_edge_vbo = 0;
    GLsizei nm_edge_vert_count = 0;
    GLuint healed_edge_vao = 0, healed_edge_vbo = 0;
    GLsizei healed_edge_vert_count = 0;

    // ── Per-triangle metadata (one entry per face) ───────────────────────────
    struct TriangleMeta {
        int         face_type;
        int         original_face_index;
        int         segment_index;
        int         subdiv_count;
        float       area;
        std::string surface_desc;
        std::string tess_approach;
    };

    // ── Pending snapshot (worker → render) ───────────────────────────────────
    struct Snapshot {
        std::vector<GpuVertex>    verts;
        std::vector<TriangleMeta> meta;      // parallel to verts/3 (one per triangle)
        std::vector<float>        border_xyz;  // x,y,z triples
        std::vector<float>        free_edge_xyz;   // free BREP boundary edges (yellow)
        std::vector<float>        open_edge_xyz;   // open triangle-soup edges (magenta)
        std::vector<float>        nm_edge_xyz;     // non-manifold edges (orange, count > 2)
        std::vector<float>        healed_edge_xyz; // loop-repaired edges (green)
        // segment_index → all border edges (for hover red wire)
        std::unordered_map<int, std::vector<float>> seg_borders;
        // original_face_index → CAD-boundary border edges (for hover green wire)
        std::unordered_map<int, std::vector<float>> orig_face_borders;
        std::string               label;
        size_t                    tri_count  = 0;
        size_t                    seg_count  = 0;
        float                     bbox_min[3] = {0,0,0};
        float                     bbox_max[3] = {0,0,0};
        bool                      valid = false;
    };
    Snapshot   pending;
    bool       has_pending = false;
    std::mutex pending_mtx;

    // ── CPU-side vertex positions for ray picking ─────────────────────────────
    // Kept in sync with the GPU verts after each upload.
    std::vector<GpuVertex>    pick_verts;
    std::vector<TriangleMeta> pick_meta;

    // ── Hover state ───────────────────────────────────────────────────────────
    int         hover_face_type    = -1;
    int         hover_orig_face    = -1;
    int         hover_subdiv       = -1;
    int         hover_subdiv_count = 0;
    float       hover_area         = 0.0f;
    std::string hover_surface_desc;
    std::string hover_tess_approach;
    double      last_pick_mx = -1e9, last_pick_my = -1e9;

    // ── Hover wire overlay GL objects ─────────────────────────────────────────
    GLuint  hover_subseg_vao = 0, hover_subseg_vbo = 0;   // red: hovered subsegment border
    GLsizei hover_subseg_vert_count = 0;
    GLuint  hover_face_vao = 0, hover_face_vbo = 0;       // green: hovered original face border
    GLsizei hover_face_vert_count = 0;

    // CPU-side border data (segment_index / orig_face_index → xyz triples)
    std::unordered_map<int, std::vector<float>> seg_border_lines;
    std::unordered_map<int, std::vector<float>> face_border_lines;

    // Track last-uploaded hover state to avoid redundant GL uploads
    int last_uploaded_hover_subdiv    = -1;
    int last_uploaded_hover_orig_face = -1;

    // Current stats (updated after upload)
    std::string current_label = "Waiting...";
    size_t      current_tris  = 0;
    size_t      current_segs  = 0;

    // ── Camera ───────────────────────────────────────────────────────────────
    Eigen::Vector3f cam_target   = {0,0,0};
    float           cam_radius   = 1.0f;
    float           cam_azimuth  = 0.3f;
    float           cam_elevation= 0.4f;
    float           cam_fov_y    = (float)(M_PI / 4.0);
    float           cam_near     = 0.001f;
    float           cam_far      = 1e6f;

    // ── Mouse state ──────────────────────────────────────────────────────────
    bool   lmb_down = false, mmb_down = false;
    double last_mx  = 0,     last_my  = 0;

    // ── UI state ─────────────────────────────────────────────────────────────
    int  display_mode  = 2;     // 0=solid 1=wireframe 2=solid+edges
    int  color_mode    = 1;     // 0=uniform 1=face_type 2=segment
    bool show_borders    = true;
    bool show_free_edges   = true;  // yellow: input BREP free boundary edges
    bool show_open_edges   = true;  // magenta: open triangle-soup edges (cracks)
    bool show_nm_edges     = true;  // orange: non-manifold edges (incident to 3+ triangles)
    bool show_healed_edges = true;  // green: edges closed by loop repair

    // ── Camera state ─────────────────────────────────────────────────────────
    bool camera_fitted = false;  // auto-fit only on the first upload

    // ── Threading ────────────────────────────────────────────────────────────
    std::atomic<bool> processing_done{false};
    std::atomic<bool> reopen_requested{false};
    std::atomic<bool> stop_stdin{false};
    std::thread       stdin_thread;

    // ── Helpers ──────────────────────────────────────────────────────────────
    Eigen::Vector3f cam_eye() const {
        return cam_target + cam_radius * Eigen::Vector3f(
            std::cos(cam_elevation) * std::sin(cam_azimuth),
            std::sin(cam_elevation),
            std::cos(cam_elevation) * std::cos(cam_azimuth));
    }

    Eigen::Matrix4f view_matrix() const {
        return make_look_at(cam_eye(), cam_target, Eigen::Vector3f::UnitY());
    }

    Eigen::Matrix4f proj_matrix(float aspect) const {
        return make_perspective(cam_fov_y, aspect, cam_near, cam_far);
    }

    void fit_camera(const Snapshot& s) {
        float cx = (s.bbox_min[0] + s.bbox_max[0]) * 0.5f;
        float cy = (s.bbox_min[1] + s.bbox_max[1]) * 0.5f;
        float cz = (s.bbox_min[2] + s.bbox_max[2]) * 0.5f;
        cam_target = {cx, cy, cz};
        float dx = s.bbox_max[0] - s.bbox_min[0];
        float dy = s.bbox_max[1] - s.bbox_min[1];
        float dz = s.bbox_max[2] - s.bbox_min[2];
        float diag = std::sqrt(dx*dx + dy*dy + dz*dz);
        cam_radius = std::max(diag * 0.75f, 1e-6f);
        cam_near   = cam_radius * 0.001f;
        cam_far    = cam_radius * 100.0f;
    }
};

// ── GLFW callbacks ────────────────────────────────────────────────────────────
static void cb_close(GLFWwindow* win) {
    auto* impl = static_cast<Viewer::Impl*>(glfwGetWindowUserPointer(win));
    if (!impl->processing_done.load()) {
        glfwSetWindowShouldClose(win, GLFW_FALSE);
        glfwHideWindow(win);
        spdlog::info("Viewer hidden — type 'v' + Enter in the terminal to reopen.");
    }
    // else: allow normal close (glfwWindowShouldClose stays true)
}

static void cb_scroll(GLFWwindow* win, double /*dx*/, double dy) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    auto* impl = static_cast<Viewer::Impl*>(glfwGetWindowUserPointer(win));
    impl->cam_radius *= std::exp((float)(-dy * 0.12));
    impl->cam_radius  = std::max(impl->cam_radius, impl->cam_near * 2.0f);
}

static void cb_mouse_button(GLFWwindow* win, int button, int action, int /*mods*/) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    auto* impl = static_cast<Viewer::Impl*>(glfwGetWindowUserPointer(win));
    double mx, my;
    glfwGetCursorPos(win, &mx, &my);
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        impl->lmb_down = (action == GLFW_PRESS);
        impl->last_mx  = mx;
        impl->last_my  = my;
    }
    if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
        impl->mmb_down = (action == GLFW_PRESS);
        impl->last_mx  = mx;
        impl->last_my  = my;
    }
}

static void cb_cursor(GLFWwindow* win, double mx, double my) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    auto* impl = static_cast<Viewer::Impl*>(glfwGetWindowUserPointer(win));
    double dx = mx - impl->last_mx;
    double dy = my - impl->last_my;
    impl->last_mx = mx;
    impl->last_my = my;

    int w, h;
    glfwGetWindowSize(win, &w, &h);
    if (w == 0 || h == 0) return;

    bool shift = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT)  == GLFW_PRESS ||
                 glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;

    bool pan = impl->mmb_down || (impl->lmb_down && shift);
    bool orbit = impl->lmb_down && !shift;

    if (orbit) {
        impl->cam_azimuth   -= (float)(dx / w * 2.0 * M_PI);
        impl->cam_elevation += (float)(dy / h * M_PI);
        const float elev_max = (float)(M_PI * 0.499);
        impl->cam_elevation = std::max(-elev_max, std::min(elev_max, impl->cam_elevation));
    }
    if (pan) {
        Eigen::Vector3f eye  = impl->cam_eye();
        Eigen::Vector3f fwd  = (impl->cam_target - eye).normalized();
        Eigen::Vector3f right = fwd.cross(Eigen::Vector3f::UnitY()).normalized();
        Eigen::Vector3f up    = right.cross(fwd).normalized();
        float scale = impl->cam_radius * std::tan(impl->cam_fov_y * 0.5f) * 2.0f / h;
        impl->cam_target -= right * (float)(dx * scale);
        impl->cam_target += up    * (float)(dy * scale);
    }
}

static void cb_key(GLFWwindow* win, int key, int /*scan*/, int action, int /*mods*/) {
    if (action != GLFW_PRESS) return;
    auto* impl = static_cast<Viewer::Impl*>(glfwGetWindowUserPointer(win));
    if (key == GLFW_KEY_F) {
        std::lock_guard<std::mutex> lock(impl->pending_mtx);
        // Re-fit to whatever is currently displayed — re-read last snapshot bbox.
        // We track last bbox in the impl for this purpose.
        // (handled inside run_event_loop after upload)
    }
}

// ── Viewer construction / destruction ────────────────────────────────────────
Viewer::Viewer(int width, int height, const std::string& title)
    : impl_(std::make_unique<Impl>())
{
    if (!glfwInit()) {
        spdlog::error("glfwInit() failed");
        return;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    impl_->window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (!impl_->window) {
        spdlog::error("glfwCreateWindow() failed");
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(impl_->window);
    glfwSwapInterval(1);  // vsync

    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        spdlog::error("glewInit() failed: {}", (const char*)glewGetErrorString(err));
        return;
    }
    // glewInit sometimes generates a spurious GL_INVALID_ENUM on core profile; clear it.
    glGetError();

    glfwSetWindowUserPointer(impl_->window, impl_.get());
    glfwSetWindowCloseCallback(impl_->window,     cb_close);
    glfwSetScrollCallback(impl_->window,          cb_scroll);
    glfwSetMouseButtonCallback(impl_->window,     cb_mouse_button);
    glfwSetCursorPosCallback(impl_->window,       cb_cursor);
    glfwSetKeyCallback(impl_->window,             cb_key);

    // ImGui setup
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(impl_->window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // Compile shaders
    impl_->tri_prog  = link_program(kTriVertSrc,  kTriFragSrc);
    impl_->line_prog = link_program(kLineVertSrc, kLineFragSrc);

    // Allocate VAO/VBO for triangles
    glGenVertexArrays(1, &impl_->tri_vao);
    glGenBuffers(1, &impl_->tri_vbo);
    glBindVertexArray(impl_->tri_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->tri_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          (void*)offsetof(GpuVertex, x));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          (void*)offsetof(GpuVertex, nx));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          (void*)offsetof(GpuVertex, face_type));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          (void*)offsetof(GpuVertex, orig_mod));
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          (void*)offsetof(GpuVertex, subdiv_mod));
    glBindVertexArray(0);

    // Allocate VAO/VBO for border lines (vec3 positions only)
    glGenVertexArrays(1, &impl_->border_vao);
    glGenBuffers(1, &impl_->border_vbo);
    glBindVertexArray(impl_->border_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->border_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    // Free BREP boundary edges (yellow) and open soup edges (magenta)
    glGenVertexArrays(1, &impl_->free_edge_vao);
    glGenBuffers(1, &impl_->free_edge_vbo);
    glBindVertexArray(impl_->free_edge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->free_edge_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    glGenVertexArrays(1, &impl_->open_edge_vao);
    glGenBuffers(1, &impl_->open_edge_vbo);
    glBindVertexArray(impl_->open_edge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->open_edge_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    glGenVertexArrays(1, &impl_->healed_edge_vao);
    glGenBuffers(1, &impl_->healed_edge_vbo);
    glBindVertexArray(impl_->healed_edge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->healed_edge_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    glGenVertexArrays(1, &impl_->nm_edge_vao);
    glGenBuffers(1, &impl_->nm_edge_vbo);
    glBindVertexArray(impl_->nm_edge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->nm_edge_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    // Hover subsegment border (red)
    glGenVertexArrays(1, &impl_->hover_subseg_vao);
    glGenBuffers(1, &impl_->hover_subseg_vbo);
    glBindVertexArray(impl_->hover_subseg_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->hover_subseg_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    // Hover original-face border (green)
    glGenVertexArrays(1, &impl_->hover_face_vao);
    glGenBuffers(1, &impl_->hover_face_vbo);
    glBindVertexArray(impl_->hover_face_vao);
    glBindBuffer(GL_ARRAY_BUFFER, impl_->hover_face_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
}

Viewer::~Viewer() {
    if (!impl_->window) return;

    // Stop stdin thread
    impl_->stop_stdin.store(true);
    if (impl_->stdin_thread.joinable())
        impl_->stdin_thread.join();

    // Cleanup GL
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteVertexArrays(1, &impl_->tri_vao);
    glDeleteBuffers(1, &impl_->tri_vbo);
    glDeleteVertexArrays(1, &impl_->border_vao);
    glDeleteBuffers(1, &impl_->border_vbo);
    glDeleteVertexArrays(1, &impl_->free_edge_vao);
    glDeleteBuffers(1, &impl_->free_edge_vbo);
    glDeleteVertexArrays(1, &impl_->open_edge_vao);
    glDeleteBuffers(1, &impl_->open_edge_vbo);
    glDeleteVertexArrays(1, &impl_->healed_edge_vao);
    glDeleteBuffers(1, &impl_->healed_edge_vbo);
    glDeleteVertexArrays(1, &impl_->nm_edge_vao);
    glDeleteBuffers(1, &impl_->nm_edge_vbo);
    glDeleteVertexArrays(1, &impl_->hover_subseg_vao);
    glDeleteBuffers(1, &impl_->hover_subseg_vbo);
    glDeleteVertexArrays(1, &impl_->hover_face_vao);
    glDeleteBuffers(1, &impl_->hover_face_vbo);
    glDeleteProgram(impl_->tri_prog);
    glDeleteProgram(impl_->line_prog);

    glfwDestroyWindow(impl_->window);
    glfwTerminate();
}

// ── Viewer::update (worker thread) ───────────────────────────────────────────
void Viewer::update(const std::vector<ViewerSegment>& segments,
                    const std::string& label,
                    const std::vector<float>& free_brep_edges,
                    const std::vector<float>& healed_edges)
{
    if (!impl_->window) return;

    Impl::Snapshot snap;
    snap.label     = label;
    snap.seg_count = segments.size();
    snap.valid     = true;
    snap.free_edge_xyz   = free_brep_edges; // bold yellow: input BREP free boundary edges
    snap.healed_edge_xyz = healed_edges;   // bold green: edges closed by loop repair

    // Open triangle-soup edges (bold magenta): undirected edges incident to exactly one
    // triangle across ALL segment meshes, merged by exact position — the real cracks (this
    // matches the validator's watertightness measure, not per-mesh is_border).
    {
        std::map<std::tuple<double,double,double>, int> pid;
        std::vector<std::array<double,3>> pos;
        auto id_of = [&](const Mesh::Point& p) {
            double x=(double)p.x(), y=(double)p.y(), z=(double)p.z();
            auto key = std::make_tuple(x,y,z);
            auto it = pid.find(key);
            if (it != pid.end()) return it->second;
            int id = (int)pos.size(); pid[key]=id; pos.push_back({x,y,z}); return id;
        };
        std::map<std::pair<int,int>, int> ecount;
        for (const auto& seg : segments) {
            if (!seg.mesh) continue;
            const Mesh& m = *seg.mesh;
            for (auto f : m.faces()) {
                int ids[3]; int k = 0;
                for (auto v : m.vertices_around_face(m.halfedge(f))) { if (k<3) ids[k]=id_of(m.point(v)); k++; }
                if (k != 3) continue;
                for (int e=0;e<3;e++){ int a=ids[e],b=ids[(e+1)%3]; if(a>b) std::swap(a,b); ecount[{a,b}]++; }
            }
        }
        for (const auto& kv : ecount) {
            const auto& A = pos[kv.first.first]; const auto& B = pos[kv.first.second];
            auto emit = [&](std::vector<float>& dst) {
                dst.push_back((float)A[0]); dst.push_back((float)A[1]); dst.push_back((float)A[2]);
                dst.push_back((float)B[0]); dst.push_back((float)B[1]); dst.push_back((float)B[2]);
            };
            if (kv.second == 1) emit(snap.open_edge_xyz);
            if (kv.second  > 2) emit(snap.nm_edge_xyz);
        }
    }

    float bmin[3] = {1e30f, 1e30f, 1e30f};
    float bmax[3] = {-1e30f,-1e30f,-1e30f};

    // ── Detect subdivision boundary edges ──────────────────────────────────
    // A border edge shared between two segments of the SAME original CAD face
    // is a subdivision boundary.  All other border edges are real CAD boundaries
    // and are the only ones rendered in the border line overlay.
    //
    // Key: canonicalized (sorted) pair of float-precision endpoint coordinates.
    // Border edges shared within a face group appear twice (once per segment),
    // so a count > 1 identifies them as subdivision boundaries.

    using Vec3f   = std::array<float, 3>;
    using EdgeKey = std::pair<Vec3f, Vec3f>;

    auto make_edge_key = [](const Mesh::Point& a, const Mesh::Point& b) -> EdgeKey {
        Vec3f pa = {(float)a.x(), (float)a.y(), (float)a.z()};
        Vec3f pb = {(float)b.x(), (float)b.y(), (float)b.z()};
        if (pa > pb) std::swap(pa, pb);
        return {pa, pb};
    };

    std::set<EdgeKey> subdiv_edges;
    {
        // Group segment indices by original_face_index
        std::unordered_map<int, std::vector<size_t>> face_groups;
        for (size_t i = 0; i < segments.size(); i++)
            face_groups[segments[i].original_face_index].push_back(i);

        for (auto& [orig_idx, group] : face_groups) {
            if (group.size() <= 1) continue;
            std::map<EdgeKey, int> edge_count;
            for (size_t si : group) {
                if (!segments[si].mesh) continue;
                const Mesh& m = *segments[si].mesh;
                for (auto e : m.edges()) {
                    if (!m.is_border(e)) continue;
                    auto h = m.halfedge(e, 0);
                    edge_count[make_edge_key(m.point(m.source(h)), m.point(m.target(h)))]++;
                }
            }
            for (auto& [key, cnt] : edge_count)
                if (cnt > 1) subdiv_edges.insert(key);
        }
    }

    // ── Build vertex and border-line data ──────────────────────────────────
    for (const auto& seg : segments) {
        if (!seg.mesh) continue;
        const Mesh& mesh = *seg.mesh;

        float ft       = (float)std::max(0, std::min(5, seg.face_type));
        float orig_m   = (float)(((seg.original_face_index % 12) + 12) % 12);
        float subdiv_m = (float)(((seg.segment_index       % 12) + 12) % 12);

        // Hover wire overlays come from pre-sampled BRep geometry, not mesh border edges.
        snap.seg_borders[seg.segment_index] = seg.subseg_wire;
        if (!snap.orig_face_borders.count(seg.original_face_index))
            snap.orig_face_borders[seg.original_face_index] = seg.orig_face_wire;

        // Global CAD-boundary overlay (cyan): mesh border edges that are not subdivision seams.
        for (auto e : mesh.edges()) {
            if (!mesh.is_border(e)) continue;
            auto h  = mesh.halfedge(e, 0);
            auto p0 = mesh.point(mesh.source(h));
            auto p1 = mesh.point(mesh.target(h));
            if (subdiv_edges.count(make_edge_key(p0, p1))) continue;
            snap.border_xyz.push_back((float)p0.x());
            snap.border_xyz.push_back((float)p0.y());
            snap.border_xyz.push_back((float)p0.z());
            snap.border_xyz.push_back((float)p1.x());
            snap.border_xyz.push_back((float)p1.y());
            snap.border_xyz.push_back((float)p1.z());
        }

        // Flat triangle soup
        for (auto face : mesh.faces()) {
            auto h0 = mesh.halfedge(face);
            auto h1 = mesh.next(h0);
            auto h2 = mesh.next(h1);
            auto p0 = mesh.point(mesh.source(h0));
            auto p1 = mesh.point(mesh.source(h1));
            auto p2 = mesh.point(mesh.source(h2));

            auto e1 = p1 - p0;
            auto e2 = p2 - p0;
            auto n  = CGAL::cross_product(e1, e2);
            double nlen = std::sqrt((double)n.squared_length());
            float nx = 0, ny = 1, nz = 0;
            if (nlen > 1e-15) {
                nx = (float)(n.x() / nlen);
                ny = (float)(n.y() / nlen);
                nz = (float)(n.z() / nlen);
            }

            for (const auto& p : {p0, p1, p2}) {
                float fx = (float)p.x(), fy = (float)p.y(), fz = (float)p.z();
                bmin[0] = std::min(bmin[0], fx); bmin[1] = std::min(bmin[1], fy); bmin[2] = std::min(bmin[2], fz);
                bmax[0] = std::max(bmax[0], fx); bmax[1] = std::max(bmax[1], fy); bmax[2] = std::max(bmax[2], fz);
                snap.verts.push_back({fx, fy, fz, nx, ny, nz, ft, orig_m, subdiv_m});
            }
            snap.meta.push_back({seg.face_type, seg.original_face_index, seg.segment_index,
                                 seg.subdiv_count, seg.area, seg.surface_desc,
                                 seg.tess_approach});
        }
        snap.tri_count += mesh.number_of_faces();
    }

    if (!snap.verts.empty()) {
        std::memcpy(snap.bbox_min, bmin, 12);
        std::memcpy(snap.bbox_max, bmax, 12);
    }

    {
        std::lock_guard<std::mutex> lock(impl_->pending_mtx);
        impl_->pending     = std::move(snap);
        impl_->has_pending = true;
    }
}

// ── Upload snapshot to GPU (render thread) ────────────────────────────────────
static void upload_snapshot(Viewer::Impl& impl, const Viewer::Impl::Snapshot& s) {
    // Triangle VBO
    glBindBuffer(GL_ARRAY_BUFFER, impl.tri_vbo);
    if (!s.verts.empty()) {
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(s.verts.size() * sizeof(GpuVertex)),
                     s.verts.data(), GL_DYNAMIC_DRAW);
    } else {
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.tri_vert_count = (GLsizei)s.verts.size();

    // Border line VBO
    glBindBuffer(GL_ARRAY_BUFFER, impl.border_vbo);
    if (!s.border_xyz.empty()) {
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(s.border_xyz.size() * sizeof(float)),
                     s.border_xyz.data(), GL_DYNAMIC_DRAW);
    } else {
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.border_vert_count = (GLsizei)(s.border_xyz.size() / 3);

    // Free BREP boundary edges (yellow)
    glBindBuffer(GL_ARRAY_BUFFER, impl.free_edge_vbo);
    if (!s.free_edge_xyz.empty())
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(s.free_edge_xyz.size() * sizeof(float)),
                     s.free_edge_xyz.data(), GL_DYNAMIC_DRAW);
    else
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.free_edge_vert_count = (GLsizei)(s.free_edge_xyz.size() / 3);

    // Open triangle-soup edges (magenta)
    glBindBuffer(GL_ARRAY_BUFFER, impl.open_edge_vbo);
    if (!s.open_edge_xyz.empty())
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(s.open_edge_xyz.size() * sizeof(float)),
                     s.open_edge_xyz.data(), GL_DYNAMIC_DRAW);
    else
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.open_edge_vert_count = (GLsizei)(s.open_edge_xyz.size() / 3);

    // Non-manifold edges (orange)
    glBindBuffer(GL_ARRAY_BUFFER, impl.nm_edge_vbo);
    if (!s.nm_edge_xyz.empty())
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(s.nm_edge_xyz.size() * sizeof(float)),
                     s.nm_edge_xyz.data(), GL_DYNAMIC_DRAW);
    else
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.nm_edge_vert_count = (GLsizei)(s.nm_edge_xyz.size() / 3);

    // Healed edges (green)
    glBindBuffer(GL_ARRAY_BUFFER, impl.healed_edge_vbo);
    if (!s.healed_edge_xyz.empty())
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(s.healed_edge_xyz.size() * sizeof(float)),
                     s.healed_edge_xyz.data(), GL_DYNAMIC_DRAW);
    else
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    impl.healed_edge_vert_count = (GLsizei)(s.healed_edge_xyz.size() / 3);

    impl.current_label = s.label;
    impl.current_tris  = s.tri_count;
    impl.current_segs  = s.seg_count;
    if (!impl.camera_fitted) {
        impl.fit_camera(s);
        impl.camera_fitted = true;
    }

    // Copy for CPU ray picking
    impl.pick_verts = s.verts;
    impl.pick_meta  = s.meta;
    impl.hover_face_type    = impl.hover_orig_face = impl.hover_subdiv = -1;
    impl.hover_subdiv_count = 0;
    impl.hover_area         = 0.0f;
    impl.hover_surface_desc.clear();
    impl.hover_tess_approach.clear();
    impl.last_pick_mx = impl.last_pick_my = -1e9;

    // Copy border maps for hover queries and clear hover GL buffers
    impl.seg_border_lines  = s.seg_borders;
    impl.face_border_lines = s.orig_face_borders;
    impl.last_uploaded_hover_subdiv    = -1;
    impl.last_uploaded_hover_orig_face = -1;
    glBindBuffer(GL_ARRAY_BUFFER, impl.hover_subseg_vbo);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    impl.hover_subseg_vert_count = 0;
    glBindBuffer(GL_ARRAY_BUFFER, impl.hover_face_vbo);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    impl.hover_face_vert_count = 0;
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// ── Ray–triangle intersection (Möller–Trumbore) ───────────────────────────────
static bool ray_tri(const Eigen::Vector3f& orig, const Eigen::Vector3f& dir,
                    const Eigen::Vector3f& v0,   const Eigen::Vector3f& v1,
                    const Eigen::Vector3f& v2,   float& t)
{
    constexpr float EPS = 1e-8f;
    Eigen::Vector3f e1  = v1 - v0;
    Eigen::Vector3f e2  = v2 - v0;
    Eigen::Vector3f h   = dir.cross(e2);
    float a = e1.dot(h);
    if (std::abs(a) < EPS) return false;
    float f = 1.0f / a;
    Eigen::Vector3f s = orig - v0;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) return false;
    Eigen::Vector3f q = s.cross(e1);
    float v = f * dir.dot(q);
    if (v < 0.0f || u + v > 1.0f) return false;
    t = f * e2.dot(q);
    return t > EPS;
}

// ── Hover picking ─────────────────────────────────────────────────────────────
static const char* kFaceTypeNames[] = {
    "Other", "Plane", "Cylinder", "Sphere", "Torus", "Cone"
};
static const float kPanelW = 220.0f;

static void pick_under_cursor(Viewer::Impl& impl, int fb_w, int fb_h,
                               double mx, double my)
{
    // No geometry or same position as last pick
    if (impl.pick_meta.empty()) return;
    if (mx == impl.last_pick_mx && my == impl.last_pick_my) return;
    impl.last_pick_mx = mx;
    impl.last_pick_my = my;

    // Mouse is over the ImGui panel — clear hover
    if (mx >= fb_w - kPanelW) {
        impl.hover_face_type = impl.hover_orig_face = impl.hover_subdiv = -1;
        return;
    }

    // Unproject cursor to world-space ray
    float ndc_x = (float)(2.0 * mx / fb_w - 1.0);
    float ndc_y = (float)(1.0 - 2.0 * my / fb_h);
    float aspect = (fb_h > 0) ? (float)fb_w / fb_h : 1.0f;

    Eigen::Matrix4f inv_vp = (impl.proj_matrix(aspect) * impl.view_matrix()).inverse();
    auto unproj = [&](float z) -> Eigen::Vector3f {
        Eigen::Vector4f p = inv_vp * Eigen::Vector4f(ndc_x, ndc_y, z, 1.0f);
        return p.head<3>() / p.w();
    };
    Eigen::Vector3f ray_orig = unproj(-1.0f);
    Eigen::Vector3f ray_dir  = (unproj(1.0f) - ray_orig).normalized();

    // Linear scan — Möller–Trumbore over flat triangle soup
    float t_min   = 1e30f;
    int   hit_tri = -1;
    size_t ntri   = impl.pick_meta.size();

    for (size_t i = 0; i < ntri; i++) {
        size_t vi = i * 3;
        Eigen::Vector3f v0(impl.pick_verts[vi+0].x, impl.pick_verts[vi+0].y, impl.pick_verts[vi+0].z);
        Eigen::Vector3f v1(impl.pick_verts[vi+1].x, impl.pick_verts[vi+1].y, impl.pick_verts[vi+1].z);
        Eigen::Vector3f v2(impl.pick_verts[vi+2].x, impl.pick_verts[vi+2].y, impl.pick_verts[vi+2].z);
        float t;
        if (ray_tri(ray_orig, ray_dir, v0, v1, v2, t) && t < t_min) {
            t_min   = t;
            hit_tri = (int)i;
        }
    }

    if (hit_tri >= 0) {
        const auto& m           = impl.pick_meta[hit_tri];
        impl.hover_face_type     = m.face_type;
        impl.hover_orig_face     = m.original_face_index;
        impl.hover_subdiv        = m.segment_index;
        impl.hover_subdiv_count  = m.subdiv_count;
        impl.hover_area          = m.area;
        impl.hover_surface_desc  = m.surface_desc;
        impl.hover_tess_approach = m.tess_approach;
    } else {
        impl.hover_face_type    = impl.hover_orig_face = impl.hover_subdiv = -1;
        impl.hover_subdiv_count = 0;
        impl.hover_area         = 0.0f;
        impl.hover_surface_desc.clear();
        impl.hover_tess_approach.clear();
    }
}

// ── Upload hover border overlays (render thread, called after pick) ───────────
static void upload_hover_overlays(Viewer::Impl& impl) {
    bool subdiv_changed = (impl.hover_subdiv != impl.last_uploaded_hover_subdiv);
    bool face_changed   = (impl.hover_orig_face != impl.last_uploaded_hover_orig_face);
    if (!subdiv_changed && !face_changed) return;

    if (subdiv_changed) {
        impl.last_uploaded_hover_subdiv = impl.hover_subdiv;
        glBindBuffer(GL_ARRAY_BUFFER, impl.hover_subseg_vbo);
        auto it = (impl.hover_subdiv >= 0)
                  ? impl.seg_border_lines.find(impl.hover_subdiv)
                  : impl.seg_border_lines.end();
        if (it != impl.seg_border_lines.end() && !it->second.empty()) {
            glBufferData(GL_ARRAY_BUFFER,
                         (GLsizeiptr)(it->second.size() * sizeof(float)),
                         it->second.data(), GL_DYNAMIC_DRAW);
            impl.hover_subseg_vert_count = (GLsizei)(it->second.size() / 3);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
            impl.hover_subseg_vert_count = 0;
        }
    }

    if (face_changed) {
        impl.last_uploaded_hover_orig_face = impl.hover_orig_face;
        glBindBuffer(GL_ARRAY_BUFFER, impl.hover_face_vbo);
        auto it = (impl.hover_orig_face >= 0 && impl.hover_subdiv_count > 1)
                  ? impl.face_border_lines.find(impl.hover_orig_face)
                  : impl.face_border_lines.end();
        if (it != impl.face_border_lines.end() && !it->second.empty()) {
            glBufferData(GL_ARRAY_BUFFER,
                         (GLsizeiptr)(it->second.size() * sizeof(float)),
                         it->second.data(), GL_DYNAMIC_DRAW);
            impl.hover_face_vert_count = (GLsizei)(it->second.size() / 3);
        } else {
            glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
            impl.hover_face_vert_count = 0;
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// ── render_frame (render thread) ─────────────────────────────────────────────
static void render_frame(Viewer::Impl& impl) {
    int fb_w, fb_h;
    glfwGetFramebufferSize(impl.window, &fb_w, &fb_h);
    glViewport(0, 0, fb_w, fb_h);
    glClearColor(0.18f, 0.18f, 0.18f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // ── 3D mesh rendering ────────────────────────────────────────────────────
    if (impl.tri_vert_count > 0) {
        float aspect = (fb_w > 0 && fb_h > 0) ? (float)fb_w / fb_h : 1.0f;
        Eigen::Matrix4f view = impl.view_matrix();
        Eigen::Matrix4f proj = impl.proj_matrix(aspect);
        Eigen::Matrix4f mvp  = proj * view;
        Eigen::Matrix3f nm   = view.topLeftCorner<3,3>().inverse().transpose();

        glUseProgram(impl.tri_prog);

        auto setM4 = [&](const char* n, const Eigen::Matrix4f& m) {
            glUniformMatrix4fv(glGetUniformLocation(impl.tri_prog, n), 1, GL_FALSE, m.data());
        };
        auto setM3 = [&](const char* n, const Eigen::Matrix3f& m) {
            glUniformMatrix3fv(glGetUniformLocation(impl.tri_prog, n), 1, GL_FALSE, m.data());
        };
        auto setI  = [&](const char* n, int v)   { glUniform1i(glGetUniformLocation(impl.tri_prog, n), v); };
        auto setF  = [&](const char* n, float v)  { glUniform1f(glGetUniformLocation(impl.tri_prog, n), v); };
        auto set3  = [&](const char* n, float a, float b, float c) {
            glUniform3f(glGetUniformLocation(impl.tri_prog, n), a, b, c);
        };
        auto set4  = [&](const char* n, float a, float b, float c, float d) {
            glUniform4f(glGetUniformLocation(impl.tri_prog, n), a, b, c, d);
        };

        setM4("u_mvp", mvp);
        setM3("u_normal_mat", nm);
        setI("u_color_mode",  impl.color_mode);
        setF("u_ambient",     0.25f);
        set3("u_light_dir", 0.577f, 0.577f, 0.577f);
        set3("u_uniform_color", 0.7f, 0.7f, 0.7f);
        set4("u_override_color", 0, 0, 0, 0);  // no override

        for (int i = 0; i < 6; i++) {
            std::string name = std::string("u_face_type_colors[") + std::to_string(i) + "]";
            glUniform3f(glGetUniformLocation(impl.tri_prog, name.c_str()),
                        kFaceTypeColors[i][0], kFaceTypeColors[i][1], kFaceTypeColors[i][2]);
        }
        for (int i = 0; i < 12; i++) {
            std::string name = std::string("u_palette[") + std::to_string(i) + "]";
            glUniform3f(glGetUniformLocation(impl.tri_prog, name.c_str()),
                        kSegColors[i][0], kSegColors[i][1], kSegColors[i][2]);
        }

        glBindVertexArray(impl.tri_vao);

        bool solid_pass = (impl.display_mode == 0 || impl.display_mode == 2);
        bool edge_pass  = (impl.display_mode == 1 || impl.display_mode == 2);

        if (solid_pass) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            if (edge_pass) {
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(1.0f, 1.0f);
            }
            glDrawArrays(GL_TRIANGLES, 0, impl.tri_vert_count);
            if (edge_pass) glDisable(GL_POLYGON_OFFSET_FILL);
        }
        if (edge_pass) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            set4("u_override_color", 0.12f, 0.12f, 0.12f, 1.0f);
            glDrawArrays(GL_TRIANGLES, 0, impl.tri_vert_count);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        glBindVertexArray(0);

        // Border edge overlay
        if (impl.show_borders && impl.border_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"),
                               1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"),
                        0.0f, 1.0f, 1.0f, 1.0f);
            glLineWidth(2.0f);
            glBindVertexArray(impl.border_vao);
            glDrawArrays(GL_LINES, 0, impl.border_vert_count);
            glBindVertexArray(0);
            glLineWidth(1.0f);
        }

        // Free BREP boundary edges (bold yellow) — drawn on top, depth-test off so they're
        // always visible.
        if (impl.show_free_edges && impl.free_edge_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"), 1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"), 1.0f, 0.92f, 0.0f, 1.0f);
            glLineWidth(4.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.free_edge_vao);
            glDrawArrays(GL_LINES, 0, impl.free_edge_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }

        // Healed edges (bold green) -- previously-open edges closed by loop repair.
        if (impl.show_healed_edges && impl.healed_edge_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"), 1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"), 0.0f, 0.9f, 0.2f, 1.0f);
            glLineWidth(4.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.healed_edge_vao);
            glDrawArrays(GL_LINES, 0, impl.healed_edge_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }

        // Open triangle-soup edges (bold magenta) -- remaining cracks.
        if (impl.show_open_edges && impl.open_edge_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"), 1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"), 1.0f, 0.0f, 1.0f, 1.0f);
            glLineWidth(4.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.open_edge_vao);
            glDrawArrays(GL_LINES, 0, impl.open_edge_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }

        // Non-manifold edges (bold orange) -- edges incident to 3+ triangles; drawn last.
        if (impl.show_nm_edges && impl.nm_edge_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"), 1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"), 1.0f, 0.5f, 0.0f, 1.0f);
            glLineWidth(4.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.nm_edge_vao);
            glDrawArrays(GL_LINES, 0, impl.nm_edge_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }

        // Hover original-face border (green) — drawn before subseg so red is on top
        if (impl.hover_face_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"),
                               1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"),
                        0.15f, 1.0f, 0.15f, 1.0f);
            glLineWidth(3.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.hover_face_vao);
            glDrawArrays(GL_LINES, 0, impl.hover_face_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }

        // Hover subsegment border (red)
        if (impl.hover_subseg_vert_count > 0) {
            glUseProgram(impl.line_prog);
            glUniformMatrix4fv(glGetUniformLocation(impl.line_prog, "u_mvp"),
                               1, GL_FALSE, mvp.data());
            glUniform4f(glGetUniformLocation(impl.line_prog, "u_color"),
                        1.0f, 0.15f, 0.15f, 1.0f);
            glLineWidth(3.0f);
            glDisable(GL_DEPTH_TEST);
            glBindVertexArray(impl.hover_subseg_vao);
            glDrawArrays(GL_LINES, 0, impl.hover_subseg_vert_count);
            glBindVertexArray(0);
            glEnable(GL_DEPTH_TEST);
            glLineWidth(1.0f);
        }
    }

    // ── ImGui panel ──────────────────────────────────────────────────────────
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    const float panel_w = 220.0f;
    ImGui::SetNextWindowPos({(float)fb_w - panel_w, 0.0f}, ImGuiCond_Always);
    ImGui::SetNextWindowSize({panel_w, (float)fb_h}, ImGuiCond_Always);
    ImGui::Begin("##controls", nullptr,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                 ImGuiWindowFlags_NoMove     | ImGuiWindowFlags_NoCollapse);

    ImGui::TextUnformatted("Stage:");
    ImGui::SameLine();
    ImGui::TextWrapped("%s", impl.current_label.c_str());
    if (!impl.processing_done.load()) {
        ImGui::TextDisabled("Processing...");
    }
    ImGui::Separator();

    ImGui::TextUnformatted("Display");
    ImGui::RadioButton("Solid",        &impl.display_mode, 0);
    ImGui::RadioButton("Wireframe",    &impl.display_mode, 1);
    ImGui::RadioButton("Solid+Edges",  &impl.display_mode, 2);
    ImGui::Separator();

    ImGui::TextUnformatted("Color");
    ImGui::RadioButton("Uniform",        &impl.color_mode, 0);
    ImGui::RadioButton("By Face Type",   &impl.color_mode, 1);
    ImGui::RadioButton("By Segment",     &impl.color_mode, 2);
    ImGui::RadioButton("By Subdivision", &impl.color_mode, 3);
    ImGui::Separator();

    ImGui::Checkbox("CAD Boundaries", &impl.show_borders);
    ImGui::Checkbox("Free Edges (yellow)", &impl.show_free_edges);
    ImGui::Checkbox("Open Edges (magenta)", &impl.show_open_edges);
    ImGui::Checkbox("Non-Manifold (orange)", &impl.show_nm_edges);
    ImGui::Checkbox("Healed Edges (green)", &impl.show_healed_edges);
    ImGui::Separator();

    ImGui::TextDisabled("Segments:  %zu", impl.current_segs);
    ImGui::TextDisabled("Triangles: %zu", impl.current_tris);
    ImGui::Separator();
    ImGui::TextUnformatted("Hover");
    if (impl.hover_face_type >= 0) {
        // Surface description (wraps if long)
        ImGui::TextWrapped("  %s", impl.hover_surface_desc.c_str());
        // Segment = original CAD face
        ImGui::Text("  Segment: %d", impl.hover_orig_face);
        // Subdivisions of that face
        if (impl.hover_subdiv_count > 1) {
            ImGui::Text("  Subdiv:  %d  (%d total)", impl.hover_subdiv, impl.hover_subdiv_count);
        } else {
            ImGui::TextDisabled("  (no subdivisions)");
        }
        // Tessellation approach
        if (!impl.hover_tess_approach.empty())
            ImGui::Text("  Method:  %s", impl.hover_tess_approach.c_str());
        // Area
        if (impl.hover_area > 0.0f)
            ImGui::Text("  Area:    %.4g", impl.hover_area);
    } else {
        ImGui::TextDisabled("  (none)");
    }
    ImGui::Separator();

    ImGui::TextDisabled("F: fit camera");
    ImGui::TextDisabled("LMB: orbit");
    ImGui::TextDisabled("MMB/Shift+LMB: pan");
    ImGui::TextDisabled("Scroll: zoom");

    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// ── Viewer::notify_done ───────────────────────────────────────────────────────
void Viewer::notify_done() {
    impl_->processing_done.store(true);
}

// ── Viewer::show ──────────────────────────────────────────────────────────────
void Viewer::show() {
    impl_->reopen_requested.store(true);
}

// ── Viewer::run_event_loop ────────────────────────────────────────────────────
void Viewer::run_event_loop() {
    if (!impl_->window) return;

    // Stdin reader thread — only when running interactively
#if defined(__unix__) || defined(__linux__)
    if (isatty(STDIN_FILENO)) {
        impl_->stdin_thread = std::thread([this]() {
            while (!impl_->stop_stdin.load()) {
                fd_set fds;
                FD_ZERO(&fds);
                FD_SET(STDIN_FILENO, &fds);
                struct timeval tv{0, 100000};  // 100 ms
                if (select(STDIN_FILENO + 1, &fds, nullptr, nullptr, &tv) > 0) {
                    int c = getchar();
                    if (c == 'v' || c == 'V')
                        impl_->reopen_requested.store(true);
                }
            }
        });
    }
#endif

    while (true) {
        bool done = glfwWindowShouldClose(impl_->window) &&
                    impl_->processing_done.load();
        if (done) break;

        // Reopen request from stdin thread
        if (impl_->reopen_requested.exchange(false)) {
            glfwShowWindow(impl_->window);
            glfwFocusWindow(impl_->window);
        }

        // F key: re-fit camera (check outside ImGui, via GLFW key state)
        if (glfwGetKey(impl_->window, GLFW_KEY_F) == GLFW_PRESS) {
            std::lock_guard<std::mutex> lock(impl_->pending_mtx);
            if (impl_->pending.valid)
                impl_->fit_camera(impl_->pending);
        }

        // Check for pending snapshot
        {
            std::lock_guard<std::mutex> lock(impl_->pending_mtx);
            if (impl_->has_pending) {
                upload_snapshot(*impl_, impl_->pending);
                impl_->has_pending = false;
            }
        }

        bool visible = glfwGetWindowAttrib(impl_->window, GLFW_VISIBLE);
        if (!visible) {
            glfwWaitEventsTimeout(0.05);
            continue;
        }

        glfwPollEvents();

        {
            int fb_w, fb_h;
            glfwGetFramebufferSize(impl_->window, &fb_w, &fb_h);
            double mx, my;
            glfwGetCursorPos(impl_->window, &mx, &my);
            pick_under_cursor(*impl_, fb_w, fb_h, mx, my);
            upload_hover_overlays(*impl_);
        }

        render_frame(*impl_);
        glfwSwapBuffers(impl_->window);
    }

    // Clean up stdin thread
    impl_->stop_stdin.store(true);
    if (impl_->stdin_thread.joinable())
        impl_->stdin_thread.join();
}
