#pragma once

// Face metric analysis for the measured-metric CDT path (Workstream D).
//
// Samples the first fundamental form (E,F,G) and principal curvatures on a grid over a
// face's UV bounds to decide (a) whether CDT is worthwhile at all (cdt_beneficial) and
// (b) whether a face is a coherent sweep we can develop (sweep_coherent), and to provide
// the arc-length tables the numeric development map / ring sampler consume.
//
// Every quantity here only informs *layout* (where Steiner points land); vertices are
// always evaluated on the true surface via surf.Value(u,v), so a metric error never moves
// a vertex off the geometry.

#include <vector>
#include <TopoDS_Face.hxx>

struct FaceMetricAnalysis {
    // sampled on an n_u × n_v grid over the face's (border-derived) UV bounds
    double k_max_abs = 0.0;        // max |principal curvature| over samples
    double aniso_ratio = 1.0;      // max over samples of |k_max|/max(|k_min|, eps)
    double shear_mean = 0.0;       // mean of |F|/sqrt(E*G)  (|cos| of iso-curve angle)
    double shear_max = 0.0;
    double su_cv_along_ring = 0.0; // max over rings of stddev(|S_u|)/mean(|S_u|)
    double sv_cv_along_col = 0.0;  // same for |S_v| along iso-u columns
    double su_min = 0.0, sv_min = 0.0;   // degenerate-edge detection (vs median)
    double su_med = 0.0, sv_med = 0.0;
    std::vector<double> u_samples, v_samples;  // grid coordinates (ascending)
    std::vector<double> arclen_u;  // cumulative 3D arc length along iso-v at v_mid
    std::vector<double> arclen_v;  // cumulative 3D arc length along iso-u at u_mid
    std::vector<double> ring_len;  // total 3D length of each sampled iso-v ring
    int n_undefined = 0;           // samples where curvature was undefined
    bool valid = false;            // false if too many evaluations failed
};

// Sample the metric on an n_u × n_v grid over [u_min,u_max] × [v_min,v_max].
FaceMetricAnalysis analyze_face_metric(const TopoDS_Face& face,
                                       double u_min, double u_max,
                                       double v_min, double v_max,
                                       int n_u = 17, int n_v = 17);

// "Should this face use CDT at all?" — from the sagitta relation error ≈ L²·κ/8:
//   L_req = sqrt(8 * max_surface_error / k_max_abs);
//   CDT is beneficial iff max_edge_length > L_req (curvature tight relative to target).
bool cdt_beneficial(const FaceMetricAnalysis& m, double max_edge_length,
                    double max_surface_error);

// "Can we develop this face?" — a coherent, non-degenerate, low-shear sweep whose
// |S_u| is roughly constant along each ring (the v1 sampler assumption).
bool sweep_coherent(const FaceMetricAnalysis& m);

// Linear interpolation of a monotone (params -> s) table, and its inverse.  Both are
// used by the development map / ring sampler (Phase D4).
double interp_arclen(const std::vector<double>& params,
                     const std::vector<double>& s, double x);
double param_at_arclen(const std::vector<double>& params,
                       const std::vector<double>& s, double arclen);
