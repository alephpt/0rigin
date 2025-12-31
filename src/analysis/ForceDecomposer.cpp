/**
 * ForceDecomposer.cpp
 *
 * Implementation of force decomposition analysis for temporal force measurement
 */

#include "analysis/ForceDecomposer.h"
#include <cmath>
#include <numeric>
#include <algorithm>

ForceDecomposer::ForceComponents ForceDecomposer::decompose(
    double dp_x, double dp_y,
    double dt,
    double grad_R_x, double grad_R_y,
    double delta) {

    ForceComponents forces;

    // Total force: F_total = dp/dt
    forces.F_total_x = dp_x / dt;
    forces.F_total_y = dp_y / dt;

    // Spatial force: F_s = -Δ·∇R
    forces.F_spatial_x = -delta * grad_R_x;
    forces.F_spatial_y = -delta * grad_R_y;

    // Temporal force: F_t = F_total - F_s (residual)
    forces.F_temporal_x = forces.F_total_x - forces.F_spatial_x;
    forces.F_temporal_y = forces.F_total_y - forces.F_spatial_y;

    return forces;
}

std::tuple<double, double> ForceDecomposer::computeSpatialGradient(
    const std::vector<float>& R_field,
    int Nx, int Ny,
    double pos_x, double pos_y,
    double dx) {

    // Get grid indices (with periodic boundaries)
    int ix = static_cast<int>(std::round(pos_x));
    int iy = static_cast<int>(std::round(pos_y));

    // Periodic boundary conditions
    ix = (ix % Nx + Nx) % Nx;
    iy = (iy % Ny + Ny) % Ny;

    // Neighboring indices with periodic boundaries
    int ix_p = (ix + 1) % Nx;
    int ix_m = (ix - 1 + Nx) % Nx;
    int iy_p = (iy + 1) % Ny;
    int iy_m = (iy - 1 + Ny) % Ny;

    // Grid indices
    int idx_xp = iy * Nx + ix_p;
    int idx_xm = iy * Nx + ix_m;
    int idx_yp = iy_p * Nx + ix;
    int idx_ym = iy_m * Nx + ix;

    // Centered finite differences
    double grad_R_x = (R_field[idx_xp] - R_field[idx_xm]) / (2.0 * dx);
    double grad_R_y = (R_field[idx_yp] - R_field[idx_ym]) / (2.0 * dx);

    return {grad_R_x, grad_R_y};
}

std::tuple<double, double> ForceDecomposer::computeTemporalGradient(
    const std::vector<float>& dR_dt_field,
    int Nx, int Ny,
    double pos_x, double pos_y,
    double dx) {

    // Same implementation as computeSpatialGradient, but operates on ∂R/∂t field
    // This gives ∇(∂R/∂t), which appears in temporal force term

    int ix = static_cast<int>(std::round(pos_x));
    int iy = static_cast<int>(std::round(pos_y));

    // Periodic boundary conditions
    ix = (ix % Nx + Nx) % Nx;
    iy = (iy % Ny + Ny) % Ny;

    // Neighboring indices
    int ix_p = (ix + 1) % Nx;
    int ix_m = (ix - 1 + Nx) % Nx;
    int iy_p = (iy + 1) % Ny;
    int iy_m = (iy - 1 + Ny) % Ny;

    // Grid indices
    int idx_xp = iy * Nx + ix_p;
    int idx_xm = iy * Nx + ix_m;
    int idx_yp = iy_p * Nx + ix;
    int idx_ym = iy_m * Nx + ix;

    // Centered finite differences for temporal gradient
    double grad_dRdt_x = (dR_dt_field[idx_xp] - dR_dt_field[idx_xm]) / (2.0 * dx);
    double grad_dRdt_y = (dR_dt_field[idx_yp] - dR_dt_field[idx_ym]) / (2.0 * dx);

    return {grad_dRdt_x, grad_dRdt_y};
}

double ForceDecomposer::computeCorrelation(
    const std::vector<double>& series1,
    const std::vector<double>& series2) {

    if (series1.size() != series2.size() || series1.empty()) {
        return 0.0;
    }

    size_t N = series1.size();

    // Compute means
    double mean1 = std::accumulate(series1.begin(), series1.end(), 0.0) / N;
    double mean2 = std::accumulate(series2.begin(), series2.end(), 0.0) / N;

    // Compute covariance and variances
    double cov = 0.0;
    double var1 = 0.0;
    double var2 = 0.0;

    for (size_t i = 0; i < N; i++) {
        double diff1 = series1[i] - mean1;
        double diff2 = series2[i] - mean2;

        cov += diff1 * diff2;
        var1 += diff1 * diff1;
        var2 += diff2 * diff2;
    }

    // Pearson correlation: ρ = cov / √(var1 · var2)
    double denominator = std::sqrt(var1 * var2);

    if (denominator < 1e-10) {
        return 0.0;  // Avoid division by zero
    }

    return cov / denominator;
}
