#include "precision_mesh/surface_metrics.h"

#include <algorithm>
#include <cmath>

#include <BRepAdaptor_Surface.hxx>
#include <BRepLProp_SLProps.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace {

double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// coefficient of variation stddev/mean of a sample set (0 if mean ~ 0)
double cv_of(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= v.size();
    if (std::abs(mean) < 1e-12) return 0.0;
    double var = 0.0;
    for (double x : v) var += (x - mean) * (x - mean);
    var /= v.size();
    return std::sqrt(var) / std::abs(mean);
}

}  // namespace

FaceMetricAnalysis analyze_face_metric(const TopoDS_Face& face,
                                       double u_min, double u_max,
                                       double v_min, double v_max,
                                       int n_u, int n_v) {
    FaceMetricAnalysis m;
    if (n_u < 3) n_u = 3;
    if (n_v < 3) n_v = 3;
    if (!(u_max > u_min) || !(v_max > v_min)) return m;  // invalid bounds -> valid=false

    BRepAdaptor_Surface surf(face);

    m.u_samples.resize(n_u);
    m.v_samples.resize(n_v);
    for (int i = 0; i < n_u; i++) m.u_samples[i] = u_min + (u_max - u_min) * i / (n_u - 1);
    for (int j = 0; j < n_v; j++) m.v_samples[j] = v_min + (v_max - v_min) * j / (n_v - 1);

    // su[i][j] = |S_u|, sv[i][j] = |S_v| at sample (u_i, v_j).
    std::vector<std::vector<double>> su(n_u, std::vector<double>(n_v, 0.0));
    std::vector<std::vector<double>> sv(n_u, std::vector<double>(n_v, 0.0));
    std::vector<double> su_all, sv_all, shear_all;
    su_all.reserve(n_u * n_v);
    sv_all.reserve(n_u * n_v);

    for (int i = 0; i < n_u; i++) {
        for (int j = 0; j < n_v; j++) {
            double u = m.u_samples[i], v = m.v_samples[j];
            gp_Pnt p;
            gp_Vec d1u, d1v;
            surf.D1(u, v, p, d1u, d1v);
            double E = d1u.Dot(d1u), F = d1u.Dot(d1v), G = d1v.Dot(d1v);
            double su_ij = std::sqrt(std::max(E, 0.0));
            double sv_ij = std::sqrt(std::max(G, 0.0));
            su[i][j] = su_ij;
            sv[i][j] = sv_ij;
            su_all.push_back(su_ij);
            sv_all.push_back(sv_ij);
            double denom = std::sqrt(std::max(E * G, 0.0));
            shear_all.push_back(denom > 1e-12 ? std::abs(F) / denom : 1.0);

            BRepLProp_SLProps props(surf, u, v, 2, 1e-9);
            if (props.IsCurvatureDefined()) {
                double k1 = std::abs(props.MaxCurvature());
                double k2 = std::abs(props.MinCurvature());
                double kmaj = std::max(k1, k2), kmin = std::min(k1, k2);
                m.k_max_abs = std::max(m.k_max_abs, kmaj);
                double ar = kmaj / std::max(kmin, 1e-9);
                m.aniso_ratio = std::max(m.aniso_ratio, ar);
            } else {
                m.n_undefined++;
            }
        }
    }

    // Shear stats.
    m.shear_max = 0.0;
    double shear_sum = 0.0;
    for (double s : shear_all) { m.shear_max = std::max(m.shear_max, s); shear_sum += s; }
    m.shear_mean = shear_all.empty() ? 0.0 : shear_sum / shear_all.size();

    // |S_u| / |S_v| extremes.
    m.su_min = *std::min_element(su_all.begin(), su_all.end());
    m.sv_min = *std::min_element(sv_all.begin(), sv_all.end());
    m.su_med = median_of(su_all);
    m.sv_med = median_of(sv_all);

    // |S_u| variation along each iso-v ring (row over u); |S_v| variation along each
    // iso-u column.
    m.su_cv_along_ring = 0.0;
    for (int j = 0; j < n_v; j++) {
        std::vector<double> ring(n_u);
        for (int i = 0; i < n_u; i++) ring[i] = su[i][j];
        m.su_cv_along_ring = std::max(m.su_cv_along_ring, cv_of(ring));
    }
    m.sv_cv_along_col = 0.0;
    for (int i = 0; i < n_u; i++) {
        std::vector<double> col(n_v);
        for (int j = 0; j < n_v; j++) col[j] = sv[i][j];
        m.sv_cv_along_col = std::max(m.sv_cv_along_col, cv_of(col));
    }

    // Arc-length tables (3D chord accumulation).
    int j_mid = n_v / 2, i_mid = n_u / 2;
    auto val = [&](double u, double v) { return surf.Value(u, v); };

    m.arclen_u.assign(n_u, 0.0);
    for (int i = 1; i < n_u; i++)
        m.arclen_u[i] = m.arclen_u[i - 1] +
            val(m.u_samples[i], m.v_samples[j_mid]).Distance(
                val(m.u_samples[i - 1], m.v_samples[j_mid]));

    m.arclen_v.assign(n_v, 0.0);
    for (int j = 1; j < n_v; j++)
        m.arclen_v[j] = m.arclen_v[j - 1] +
            val(m.u_samples[i_mid], m.v_samples[j]).Distance(
                val(m.u_samples[i_mid], m.v_samples[j - 1]));

    m.ring_len.assign(n_v, 0.0);
    for (int j = 0; j < n_v; j++) {
        double len = 0.0;
        for (int i = 1; i < n_u; i++)
            len += val(m.u_samples[i], m.v_samples[j]).Distance(
                       val(m.u_samples[i - 1], m.v_samples[j]));
        m.ring_len[j] = len;
    }

    // Valid if the geometry sampled non-degenerately and most curvatures were defined.
    int total = n_u * n_v;
    m.valid = (m.su_med > 1e-9) && (m.sv_med > 1e-9) &&
              (m.n_undefined <= total / 2);
    return m;
}

bool cdt_beneficial(const FaceMetricAnalysis& m, double max_edge_length,
                    double max_surface_error) {
    if (!m.valid) return false;
    if (m.k_max_abs < 1e-9) return false;          // ~flat -> isotropic remeshing suffices
    double L_req = std::sqrt(8.0 * max_surface_error / m.k_max_abs);
    return max_edge_length > L_req;                // target chords tighter than curvature -> CDT
}

bool sweep_coherent(const FaceMetricAnalysis& m) {
    if (!m.valid) return false;
    return m.su_min > 0.05 * m.su_med &&          // no collapsed edge (degenerate patch)
           m.sv_min > 0.05 * m.sv_med &&
           m.shear_max  < 0.5 &&                  // iso-curves never closer than ~60°
           m.shear_mean < 0.25 &&                 //  ... and mostly near-orthogonal
           m.su_cv_along_ring < 0.3 &&            // v1 sampler assumption (relax in D4.1)
           m.sv_cv_along_col  < 0.3;
}

double interp_arclen(const std::vector<double>& params,
                     const std::vector<double>& s, double x) {
    if (params.empty()) return 0.0;
    if (x <= params.front()) return s.front();
    if (x >= params.back())  return s.back();
    auto it = std::lower_bound(params.begin(), params.end(), x);
    size_t k = it - params.begin();
    double d = params[k] - params[k - 1];
    double t = d > 0 ? (x - params[k - 1]) / d : 0.0;
    return s[k - 1] + t * (s[k] - s[k - 1]);
}

double param_at_arclen(const std::vector<double>& params,
                       const std::vector<double>& s, double arclen) {
    if (s.empty()) return 0.0;
    if (arclen <= s.front()) return params.front();
    if (arclen >= s.back())  return params.back();
    auto it = std::lower_bound(s.begin(), s.end(), arclen);
    size_t k = it - s.begin();
    double d = s[k] - s[k - 1];
    double t = d > 0 ? (arclen - s[k - 1]) / d : 0.0;
    return params[k - 1] + t * (params[k] - params[k - 1]);
}
