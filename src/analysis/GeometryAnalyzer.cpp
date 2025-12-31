#include "GeometryAnalyzer.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

// ============================================================================
// Metric Tensor Computation
// ============================================================================

GeometryAnalyzer::Metric GeometryAnalyzer::computeMetric(double R) {
    // Guard against metric singularity (R → 0)
    R = std::max(R, R_MIN);

    Metric metric;

    // Diagonal metric: ds² = -R²dt² + dx² + dy²
    metric.g_tt = -R * R;
    metric.g_xx = 1.0;
    metric.g_yy = 1.0;
    metric.g_tx = 0.0;
    metric.g_ty = 0.0;

    // Inverse metric: g^μν (for raising indices)
    metric.g_inv_tt = -1.0 / (R * R);
    metric.g_inv_xx = 1.0;
    metric.g_inv_yy = 1.0;

    return metric;
}

Eigen::Matrix3d GeometryAnalyzer::getMetricMatrix(double R) {
    R = std::max(R, R_MIN);

    Eigen::Matrix3d g = Eigen::Matrix3d::Zero();
    g(0, 0) = -R * R;  // g_tt
    g(1, 1) = 1.0;     // g_xx
    g(2, 2) = 1.0;     // g_yy

    return g;
}

// ============================================================================
// Gradient and Laplacian Operators
// ============================================================================

Eigen::Vector2d GeometryAnalyzer::computeRGradient(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    int ix, int iy,
    double dx, double dy)
{
    // Boundary check for centered stencil
    if (ix < BOUNDARY_GUARD || ix >= Nx - BOUNDARY_GUARD ||
        iy < BOUNDARY_GUARD || iy >= Ny - BOUNDARY_GUARD) {
        return Eigen::Vector2d::Zero();
    }

    // Centered finite differences
    // ∂R/∂x ≈ (R[i+1,j] - R[i-1,j]) / (2Δx)
    int idx = iy * Nx + ix;
    double dR_dx = (R_field[idx + 1] - R_field[idx - 1]) / (2.0 * dx);
    double dR_dy = (R_field[idx + Nx] - R_field[idx - Nx]) / (2.0 * dy);

    return Eigen::Vector2d(dR_dx, dR_dy);
}

double GeometryAnalyzer::computeLaplacian(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    int ix, int iy,
    double dx, double dy)
{
    // Boundary check
    if (ix < BOUNDARY_GUARD || ix >= Nx - BOUNDARY_GUARD ||
        iy < BOUNDARY_GUARD || iy >= Ny - BOUNDARY_GUARD) {
        return 0.0;
    }

    // 5-point stencil Laplacian
    // ∇²R = (R[i+1] + R[i-1] - 2R[i]) / dx² + (R[j+1] + R[j-1] - 2R[j]) / dy²
    int idx = iy * Nx + ix;
    double R_center = R_field[idx];
    double R_left = R_field[idx - 1];
    double R_right = R_field[idx + 1];
    double R_down = R_field[idx - Nx];
    double R_up = R_field[idx + Nx];

    double d2R_dx2 = (R_right + R_left - 2.0 * R_center) / (dx * dx);
    double d2R_dy2 = (R_up + R_down - 2.0 * R_center) / (dy * dy);

    return d2R_dx2 + d2R_dy2;
}

// ============================================================================
// Christoffel Symbol Computation
// ============================================================================

GeometryAnalyzer::ChristoffelSymbols GeometryAnalyzer::computeChristoffel(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    int ix, int iy,
    double dx, double dy)
{
    ChristoffelSymbols gamma;

    // Get R value and gradient at grid point
    int idx = iy * Nx + ix;
    double R = std::max(R_field[idx], R_MIN);
    Eigen::Vector2d grad_R = computeRGradient(R_field, Nx, Ny, ix, iy, dx, dy);
    double dR_dx = grad_R(0);
    double dR_dy = grad_R(1);

    // For diagonal metric g_μν = diag(-R², 1, 1):
    //
    // Γ^λ_μν = (1/2) g^λσ (∂_μ g_νσ + ∂_ν g_μσ - ∂_σ g_μν)
    //
    // Key derivatives:
    // ∂g_tt/∂x = ∂(-R²)/∂x = -2R ∂R/∂x
    // ∂g_tt/∂y = ∂(-R²)/∂y = -2R ∂R/∂y

    // Time-time components (dominant for non-relativistic particles)
    // Γ^x_00 = (1/2) g^xx ∂_x g_00 = (1/2)(1)(-2R ∂R/∂x) = -R ∂R/∂x
    // Γ^y_00 = (1/2) g^yy ∂_y g_00 = (1/2)(1)(-2R ∂R/∂y) = -R ∂R/∂y
    gamma.Gamma_x_tt = -R * dR_dx;
    gamma.Gamma_y_tt = -R * dR_dy;

    // Mixed time-space components
    // Γ^t_0x = (1/2) g^tt ∂_x g_00 = (1/2)(-1/R²)(-2R ∂R/∂x) = (1/R) ∂R/∂x
    // Γ^t_0y = (1/2) g^tt ∂_y g_00 = (1/2)(-1/R²)(-2R ∂R/∂y) = (1/R) ∂R/∂y
    gamma.Gamma_t_tx = dR_dx / R;
    gamma.Gamma_t_ty = dR_dy / R;

    // Space-space components (second order in ∂R, usually small)
    // Γ^x_xx = -(1/2) g^xx ∂_x g_tt = -(1/2)(1)(-2R ∂R/∂x) = R ∂R/∂x
    // Γ^y_yy = -(1/2) g^yy ∂_y g_tt = -(1/2)(1)(-2R ∂R/∂y) = R ∂R/∂y
    gamma.Gamma_x_xx = R * dR_dx;
    gamma.Gamma_y_yy = R * dR_dy;

    // Off-diagonal spatial terms
    // Γ^x_yy = -(1/2) g^xx ∂_y g_tt = -(1/2)(1)(-2R ∂R/∂y) = R ∂R/∂y
    // Γ^y_xx = -(1/2) g^yy ∂_x g_tt = -(1/2)(1)(-2R ∂R/∂x) = R ∂R/∂x
    gamma.Gamma_x_xy = R * dR_dy;
    gamma.Gamma_y_xy = R * dR_dx;

    return gamma;
}

// ============================================================================
// Interpolation Methods
// ============================================================================

double GeometryAnalyzer::interpolateR(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double x, double y,
    double L_domain)
{
    // Convert physical coordinates to grid coordinates
    double dx_grid = L_domain / static_cast<double>(Nx);
    double dy_grid = L_domain / static_cast<double>(Ny);

    double fx = x / dx_grid;
    double fy = y / dy_grid;

    // Grid indices for bilinear interpolation
    int ix0 = static_cast<int>(std::floor(fx));
    int iy0 = static_cast<int>(std::floor(fy));
    int ix1 = ix0 + 1;
    int iy1 = iy0 + 1;

    // Clamp to grid bounds
    ix0 = std::clamp(ix0, 0, Nx - 1);
    ix1 = std::clamp(ix1, 0, Nx - 1);
    iy0 = std::clamp(iy0, 0, Ny - 1);
    iy1 = std::clamp(iy1, 0, Ny - 1);

    // Interpolation weights
    double wx = fx - static_cast<double>(ix0);
    double wy = fy - static_cast<double>(iy0);
    wx = std::clamp(wx, 0.0, 1.0);
    wy = std::clamp(wy, 0.0, 1.0);

    // Bilinear interpolation
    double R00 = R_field[iy0 * Nx + ix0];
    double R10 = R_field[iy0 * Nx + ix1];
    double R01 = R_field[iy1 * Nx + ix0];
    double R11 = R_field[iy1 * Nx + ix1];

    double R_interp = (1.0 - wx) * (1.0 - wy) * R00 +
                      wx * (1.0 - wy) * R10 +
                      (1.0 - wx) * wy * R01 +
                      wx * wy * R11;

    return std::max(R_interp, R_MIN);
}

Eigen::Vector2d GeometryAnalyzer::interpolateRGradient(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double x, double y,
    double L_domain)
{
    // Convert physical coordinates to grid coordinates
    double dx_grid = L_domain / static_cast<double>(Nx);
    double dy_grid = L_domain / static_cast<double>(Ny);

    double fx = x / dx_grid;
    double fy = y / dy_grid;

    int ix = static_cast<int>(std::round(fx));
    int iy = static_cast<int>(std::round(fy));

    // Clamp to interior (for centered stencil)
    ix = std::clamp(ix, BOUNDARY_GUARD, Nx - BOUNDARY_GUARD - 1);
    iy = std::clamp(iy, BOUNDARY_GUARD, Ny - BOUNDARY_GUARD - 1);

    // Compute gradient at nearest grid point
    return computeRGradient(R_field, Nx, Ny, ix, iy, dx_grid, dy_grid);
}

GeometryAnalyzer::ChristoffelSymbols GeometryAnalyzer::computeChristoffelInterpolated(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double x, double y,
    double L_domain)
{
    // Get R value and gradient via interpolation
    double R = interpolateR(R_field, Nx, Ny, x, y, L_domain);
    Eigen::Vector2d grad_R = interpolateRGradient(R_field, Nx, Ny, x, y, L_domain);
    double dR_dx = grad_R(0);
    double dR_dy = grad_R(1);

    ChristoffelSymbols gamma;

    // Compute Christoffel symbols using interpolated values
    gamma.Gamma_x_tt = -R * dR_dx;
    gamma.Gamma_y_tt = -R * dR_dy;
    gamma.Gamma_t_tx = dR_dx / R;
    gamma.Gamma_t_ty = dR_dy / R;
    gamma.Gamma_x_xx = R * dR_dx;
    gamma.Gamma_y_yy = R * dR_dy;
    gamma.Gamma_x_xy = R * dR_dy;
    gamma.Gamma_y_xy = R * dR_dx;

    return gamma;
}

// ============================================================================
// Geodesic Acceleration
// ============================================================================

Eigen::Vector2d GeometryAnalyzer::computeGeodesicAcceleration(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double x, double y,
    double vx, double vy,
    double L_domain)
{
    // Get Christoffel symbols at particle position
    ChristoffelSymbols gamma = computeChristoffelInterpolated(
        R_field, Nx, Ny, x, y, L_domain);

    // Geodesic equation: d²x^i/dτ² + Γ^i_μν (dx^μ/dτ)(dx^ν/dτ) = 0
    //
    // For coordinate 4-velocity v^μ = (dt/dτ, dx/dτ, dy/dτ):
    // - In non-relativistic limit: dt/dτ ≈ 1, dx/dτ ≈ vx, dy/dτ ≈ vy
    //
    // Acceleration components:
    // a^x = -Γ^x_μν v^μ v^ν
    //     = -Γ^x_00 - 2Γ^x_0x vx - 2Γ^x_0y vy - Γ^x_xx vx² - Γ^x_yy vy² - 2Γ^x_xy vx vy
    //
    // a^y = -Γ^y_μν v^μ v^ν
    //     = -Γ^y_00 - 2Γ^y_0x vx - 2Γ^y_0y vy - Γ^y_xx vx² - Γ^y_yy vy² - 2Γ^y_xy vx vy

    // Note: Γ^x_0x = Γ^x_0y = 0 for diagonal metric (only Γ^t_0x ≠ 0)
    // Simplified: a^x = -Γ^x_00 - Γ^x_xx vx² - Γ^x_yy vy²
    //             a^y = -Γ^y_00 - Γ^y_xx vx² - Γ^y_yy vy²

    double ax = -gamma.Gamma_x_tt - gamma.Gamma_x_xx * vx * vx - gamma.Gamma_x_xy * vy * vy;
    double ay = -gamma.Gamma_y_tt - gamma.Gamma_y_xy * vx * vx - gamma.Gamma_y_yy * vy * vy;

    return Eigen::Vector2d(ax, ay);
}

// ============================================================================
// Curvature Diagnostics
// ============================================================================

double GeometryAnalyzer::computeRiemannScalar(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    int ix, int iy,
    double dx, double dy)
{
    // For diagonal metric g_μν = diag(-R², 1, 1):
    //
    // Ricci scalar R = g^μν R_μν (scalar curvature)
    //
    // Weak field approximation:
    // R ≈ -2∇²(ln R) / R² = -2(∇²R)/R + 2|∇R|²/R²
    //
    // This is an approximate formula valid for slowly varying R

    int idx = iy * Nx + ix;
    double R = std::max(R_field[idx], R_MIN);

    Eigen::Vector2d grad_R = computeRGradient(R_field, Nx, Ny, ix, iy, dx, dy);
    double laplacian_R = computeLaplacian(R_field, Nx, Ny, ix, iy, dx, dy);

    double grad_R_squared = grad_R.squaredNorm();

    // R ≈ -2(∇²R)/R + 2|∇R|²/R²
    double riemann_scalar = -2.0 * laplacian_R / R + 2.0 * grad_R_squared / (R * R);

    return riemann_scalar;
}

double GeometryAnalyzer::computeGaussianCurvature(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    int ix, int iy,
    double dx, double dy)
{
    // Gaussian curvature K for 2D spatial slice (x, y)
    //
    // For pure spatial metric h_ij = diag(1, 1): K = 0 (flat)
    //
    // However, R-field induces effective curvature via time component:
    // K_eff ≈ (∇²R) / R (leading order in weak field)
    //
    // This measures spatial curvature induced by spacetime geometry

    int idx = iy * Nx + ix;
    double R = std::max(R_field[idx], R_MIN);

    double laplacian_R = computeLaplacian(R_field, Nx, Ny, ix, iy, dx, dy);

    // K_eff ≈ (∇²R) / R
    double gaussian_curvature = laplacian_R / R;

    return gaussian_curvature;
}
