/**
 * GeodesicIntegrator.cpp
 *
 * Implementation of geodesic equation verification in TRD curved spacetime
 */

#include "GeodesicIntegrator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

GeodesicIntegrator::GeodesicIntegrator(int Nx, int Ny, double delta)
    : _Nx(Nx), _Ny(Ny), _delta(delta) {
    if (Nx <= 0 || Ny <= 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }
}

GeodesicIntegrator::~GeodesicIntegrator() {
}

GeodesicIntegrator::MetricTensor GeodesicIntegrator::computeMetric(
    double x, double y,
    const std::vector<double>& R_field,
    const std::vector<double>& v_field) const {

    MetricTensor metric;

    // Interpolate R and v at (x, y)
    double R = interpolateField(x, y, R_field);
    double v = interpolateField(x, y, v_field);

    // Clamp R to valid range (must be > 1 for valid metric)
    R = std::max(1.0001, R);
    // Clamp v to [0, 1)
    v = std::max(0.0, std::min(0.9999, v));

    // TRD metric: g_μν = R² × diag(-(1-v²), 1, 1, 0)
    double R2 = R * R;
    double g00_coeff = 1.0 - v * v;

    metric.g00 = -g00_coeff * R2;  // Time component: negative signature
    metric.g11 = R2;               // x spatial component
    metric.g22 = R2;               // y spatial component
    metric.g33 = 0.0;              // No time dimension in 2D spatial

    // Off-diagonal components (all zero for diagonal metric)
    metric.g01 = metric.g02 = metric.g03 = 0.0;
    metric.g12 = metric.g13 = metric.g23 = 0.0;

    return metric;
}

double GeodesicIntegrator::interpolateField(
    double x, double y,
    const std::vector<double>& field) const {

    // Clamp to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(_Nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(_Ny - 1), y));

    int xi = static_cast<int>(x);
    int yi = static_cast<int>(y);

    // Check bounds
    if (xi < 0 || xi >= _Nx - 1 || yi < 0 || yi >= _Ny - 1) {
        // Out of bounds, use nearest neighbor
        int xi_clamp = std::max(0, std::min(_Nx - 1, xi));
        int yi_clamp = std::max(0, std::min(_Ny - 1, yi));
        return field[yi_clamp * _Nx + xi_clamp];
    }

    // Bilinear interpolation
    double fx = x - xi;
    double fy = y - yi;

    double f00 = field[yi * _Nx + xi];
    double f10 = field[yi * _Nx + (xi + 1)];
    double f01 = field[(yi + 1) * _Nx + xi];
    double f11 = field[(yi + 1) * _Nx + (xi + 1)];

    double f0 = f00 * (1.0 - fx) + f10 * fx;
    double f1 = f01 * (1.0 - fx) + f11 * fx;

    return f0 * (1.0 - fy) + f1 * fy;
}

void GeodesicIntegrator::fieldGradient(
    double x, double y,
    const std::vector<double>& field,
    double& grad_x, double& grad_y,
    double h) const {

    // Finite difference: ∂f/∂x ≈ (f(x+h) - f(x-h)) / 2h
    grad_x = (interpolateField(x + h, y, field) - interpolateField(x - h, y, field)) / (2.0 * h);
    grad_y = (interpolateField(x, y + h, field) - interpolateField(x, y - h, field)) / (2.0 * h);
}

void GeodesicIntegrator::fieldSecondGradient(
    double x, double y,
    const std::vector<double>& field,
    double& d2_xx, double& d2_yy,
    double h) const {

    // Finite difference: ∂²f/∂x² ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
    double f_center = interpolateField(x, y, field);
    d2_xx = (interpolateField(x + h, y, field) - 2.0 * f_center + interpolateField(x - h, y, field)) / (h * h);
    d2_yy = (interpolateField(x, y + h, field) - 2.0 * f_center + interpolateField(x, y - h, field)) / (h * h);
}

GeodesicIntegrator::ChristoffelSymbols GeodesicIntegrator::computeChristoffel(
    double x, double y,
    const std::vector<double>& R_field,
    const std::vector<double>& v_field,
    double h) const {

    ChristoffelSymbols christoffel;

    // Get metric at center point
    MetricTensor g_center = computeMetric(x, y, R_field, v_field);

    // Compute metric derivatives using finite differences
    // ∂g_μν/∂x^ρ at (x, y)

    // For diagonal metric g_μν = diag(g00, g11, g22, g33)
    // Compute g_μν and its derivatives

    double R = interpolateField(x, y, R_field);
    double v = interpolateField(x, y, v_field);
    R = std::max(1.0001, R);
    v = std::max(0.0, std::min(0.9999, v));

    // Derivatives of R and v
    double dR_dx, dR_dy;
    fieldGradient(x, y, R_field, dR_dx, dR_dy, h);

    double dv_dx, dv_dy;
    fieldGradient(x, y, v_field, dv_dx, dv_dy, h);

    // Second derivatives
    double d2R_dxx, d2R_dyy;
    fieldSecondGradient(x, y, R_field, d2R_dxx, d2R_dyy, h);

    // Compute g^μν (inverse metric)
    // For diagonal g_μν, g^μν = diag(1/g00, 1/g11, 1/g22, 1/g33)
    double g00 = g_center.g00;
    double g11 = g_center.g11;
    double g22 = g_center.g22;

    // Handle g00 carefully (it's negative)
    double g00_inv = (std::abs(g00) > 1e-10) ? 1.0 / g00 : 0.0;
    double g11_inv = (std::abs(g11) > 1e-10) ? 1.0 / g11 : 0.0;
    double g22_inv = (std::abs(g22) > 1e-10) ? 1.0 / g22 : 0.0;

    double R2 = R * R;
    double g00_coeff = 1.0 - v * v;

    // ∂g_00/∂x = ∂[-(1-v²)R²]/∂x = -[(1-v²)·2R·∂R/∂x + R²·(-2v·∂v/∂x)]
    double dg00_dx = -(g00_coeff * 2.0 * R * dR_dx + R2 * (-2.0 * v * dv_dx));
    double dg00_dy = -(g00_coeff * 2.0 * R * dR_dy + R2 * (-2.0 * v * dv_dy));

    // ∂g_11/∂x = ∂[R²]/∂x = 2R·∂R/∂x
    double dg11_dx = 2.0 * R * dR_dx;
    double dg11_dy = 2.0 * R * dR_dy;

    // ∂g_22/∂x = ∂[R²]/∂x = 2R·∂R/∂x
    double dg22_dx = 2.0 * R * dR_dx;
    double dg22_dy = 2.0 * R * dR_dy;

    // Now compute Christoffel symbols: Γ^μ_νλ = (1/2)g^μρ(∂_ν g_ρλ + ∂_λ g_νρ - ∂_ρ g_νλ)
    // For 2D spatial (ignore time index 0, work with spatial 1,2)

    // Non-zero Christoffel symbols for diagonal metric:

    // Γ^1_11 = (1/2)g^11(∂_1 g_11 + ∂_1 g_11 - ∂_1 g_11) = (1/2)g^11(∂_1 g_11)
    // Actually: Γ^1_11 = (1/2)g^11 ∂_1 g_11 (for spatial coordinates)
    christoffel.value[1][1][1] = 0.5 * g11_inv * dg11_dx;

    // Γ^1_22 = -(1/2)g^11(∂_2 g_22) = -(1/2)g^11 ∂_2 g_22
    christoffel.value[1][2][2] = -0.5 * g11_inv * dg22_dy;

    // Γ^1_12 = Γ^1_21 = (1/2)g^11 ∂_2 g_11
    christoffel.value[1][1][2] = 0.5 * g11_inv * dg11_dy;
    christoffel.value[1][2][1] = 0.5 * g11_inv * dg11_dy;

    // Γ^2_22 = (1/2)g^22 ∂_2 g_22
    christoffel.value[2][2][2] = 0.5 * g22_inv * dg22_dy;

    // Γ^2_11 = -(1/2)g^22 ∂_1 g_11
    christoffel.value[2][1][1] = -0.5 * g22_inv * dg11_dx;

    // Γ^2_12 = Γ^2_21 = (1/2)g^22 ∂_1 g_22
    christoffel.value[2][1][2] = 0.5 * g22_inv * dg22_dx;
    christoffel.value[2][2][1] = 0.5 * g22_inv * dg22_dx;

    return christoffel;
}

void GeodesicIntegrator::computeGeodesicAcceleration(
    const std::array<double, 2>& pos,
    const std::array<double, 2>& vel,
    const ChristoffelSymbols& christoffel,
    double& ax, double& ay) const {

    double vx = vel[0];
    double vy = vel[1];

    // Geodesic equation: a^μ = d²x^μ/dτ² = -Γ^μ_νλ(dx^ν/dτ)(dx^λ/dτ)
    // For spatial coordinates (ignoring time index 0):

    // a^1 = -Γ^1_11 vx·vx - 2·Γ^1_12 vx·vy - Γ^1_22 vy·vy
    // a^2 = -Γ^2_11 vx·vx - 2·Γ^2_12 vx·vy - Γ^2_22 vy·vy

    ax = -(christoffel.value[1][1][1] * vx * vx +
           2.0 * christoffel.value[1][1][2] * vx * vy +
           christoffel.value[1][2][2] * vy * vy);

    ay = -(christoffel.value[2][1][1] * vx * vx +
           2.0 * christoffel.value[2][1][2] * vx * vy +
           christoffel.value[2][2][2] * vy * vy);
}

std::vector<GeodesicIntegrator::GeodesicPoint>
GeodesicIntegrator::integrateGeodesic(
    const std::array<double, 2>& initial_pos,
    const std::array<double, 2>& initial_vel,
    const std::vector<double>& R_field,
    const std::vector<double>& v_field,
    double dt, int num_steps) const {

    std::vector<GeodesicPoint> trajectory;
    trajectory.reserve(num_steps + 1);

    // Initial state
    std::array<double, 2> pos = initial_pos;
    std::array<double, 2> vel = initial_vel;
    double t = 0.0;

    // Store initial point
    GeodesicPoint point;
    point.t = t;
    point.x = pos[0];
    point.y = pos[1];
    point.vx = vel[0];
    point.vy = vel[1];
    point.acceleration_magnitude = 0.0;
    trajectory.push_back(point);

    // RK4 integration
    for (int step = 0; step < num_steps; ++step) {
        // K1: derivatives at current point
        ChristoffelSymbols christ_k1 = computeChristoffel(pos[0], pos[1], R_field, v_field);
        double ax_k1, ay_k1;
        computeGeodesicAcceleration(pos, vel, christ_k1, ax_k1, ay_k1);

        // K2: midpoint derivatives
        std::array<double, 2> pos_k2 = {pos[0] + 0.5 * dt * vel[0], pos[1] + 0.5 * dt * vel[1]};
        std::array<double, 2> vel_k2 = {vel[0] + 0.5 * dt * ax_k1, vel[1] + 0.5 * dt * ay_k1};
        ChristoffelSymbols christ_k2 = computeChristoffel(pos_k2[0], pos_k2[1], R_field, v_field);
        double ax_k2, ay_k2;
        computeGeodesicAcceleration(pos_k2, vel_k2, christ_k2, ax_k2, ay_k2);

        // K3: midpoint derivatives (with k2 velocity)
        std::array<double, 2> pos_k3 = {pos[0] + 0.5 * dt * vel_k2[0], pos[1] + 0.5 * dt * vel_k2[1]};
        std::array<double, 2> vel_k3 = {vel[0] + 0.5 * dt * ax_k2, vel[1] + 0.5 * dt * ay_k2};
        ChristoffelSymbols christ_k3 = computeChristoffel(pos_k3[0], pos_k3[1], R_field, v_field);
        double ax_k3, ay_k3;
        computeGeodesicAcceleration(pos_k3, vel_k3, christ_k3, ax_k3, ay_k3);

        // K4: end point derivatives
        std::array<double, 2> pos_k4 = {pos[0] + dt * vel_k3[0], pos[1] + dt * vel_k3[1]};
        std::array<double, 2> vel_k4 = {vel[0] + dt * ax_k3, vel[1] + dt * ay_k3};
        ChristoffelSymbols christ_k4 = computeChristoffel(pos_k4[0], pos_k4[1], R_field, v_field);
        double ax_k4, ay_k4;
        computeGeodesicAcceleration(pos_k4, vel_k4, christ_k4, ax_k4, ay_k4);

        // Update: x += dt * v, v += dt * a (weighted average)
        pos[0] += dt * (vel[0] + 2.0 * vel_k2[0] + 2.0 * vel_k3[0] + vel_k4[0]) / 6.0;
        pos[1] += dt * (vel[1] + 2.0 * vel_k2[1] + 2.0 * vel_k3[1] + vel_k4[1]) / 6.0;

        vel[0] += dt * (ax_k1 + 2.0 * ax_k2 + 2.0 * ax_k3 + ax_k4) / 6.0;
        vel[1] += dt * (ay_k1 + 2.0 * ay_k2 + 2.0 * ay_k3 + ay_k4) / 6.0;

        // Clamp position to grid
        clampCoordinate(pos[0], pos[1]);

        // Store trajectory point
        t += dt;
        point.t = t;
        point.x = pos[0];
        point.y = pos[1];
        point.vx = vel[0];
        point.vy = vel[1];
        point.acceleration_magnitude = std::sqrt(ax_k4 * ax_k4 + ay_k4 * ay_k4);
        trajectory.push_back(point);
    }

    return trajectory;
}

std::vector<double> GeodesicIntegrator::compareTrajectories(
    const std::vector<std::array<double, 3>>& dirac_trajectory,
    const std::vector<std::array<double, 3>>& geodesic_trajectory) const {

    std::vector<double> deviations;

    // Compare positions at each time step
    size_t min_steps = std::min(dirac_trajectory.size(), geodesic_trajectory.size());

    for (size_t i = 0; i < min_steps; ++i) {
        double x_dirac = dirac_trajectory[i][1];
        double y_dirac = dirac_trajectory[i][2];
        double x_geo = geodesic_trajectory[i][1];
        double y_geo = geodesic_trajectory[i][2];

        // Displacement: |r_dirac - r_geodesic|
        double dx = x_dirac - x_geo;
        double dy = y_dirac - y_geo;
        double displacement = std::sqrt(dx * dx + dy * dy);

        // Reference: |r_geodesic|
        double reference = std::sqrt(x_geo * x_geo + y_geo * y_geo);

        // Relative deviation
        double deviation = (reference > 1e-6) ? (displacement / reference) : displacement;
        deviations.push_back(deviation);
    }

    return deviations;
}
