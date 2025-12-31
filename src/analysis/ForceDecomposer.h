/**
 * ForceDecomposer.h
 *
 * Phase 4 Test 4.2: Temporal Force Measurement
 *
 * Decomposes momentum change into spatial and temporal force components:
 * - Spatial force: F_s = -Δ·∇R (known from existing theory)
 * - Temporal force: F_t = dp/dt - F_s (residual after spatial contribution)
 *
 * Physics Hypothesis:
 * Time gradient creates force: F_temp ∝ -∇(∂R/∂t)
 * Momentum has spatial AND temporal contributions.
 *
 * Prediction:
 * |F_temp| ≠ 0 in time-varying R-fields (especially during vortex motion)
 * Correlation: ρ(F_temp, ∂R/∂t) > 0.5
 */

#pragma once

#include <vector>
#include <tuple>

class ForceDecomposer {
public:
    struct ForceComponents {
        double F_spatial_x;    // Spatial force x-component: -Δ·∂R/∂x
        double F_spatial_y;    // Spatial force y-component: -Δ·∂R/∂y
        double F_temporal_x;   // Temporal force x-component: dp_x/dt - F_spatial_x
        double F_temporal_y;   // Temporal force y-component: dp_y/dt - F_spatial_y
        double F_total_x;      // Total force x: dp_x/dt
        double F_total_y;      // Total force y: dp_y/dt
    };

    /**
     * Decompose momentum change into spatial and temporal forces
     *
     * Physics:
     * - Spatial force: F_s = -Δ·∇R (gradient of synchronization field)
     * - Total force: F_total = dp/dt (rate of momentum change)
     * - Temporal force: F_t = F_total - F_s (residual)
     *
     * @param dp_x Momentum change in x (current - previous)
     * @param dp_y Momentum change in y (current - previous)
     * @param dt Timestep
     * @param grad_R_x R-field gradient x-component at particle center
     * @param grad_R_y R-field gradient y-component at particle center
     * @param delta Vacuum potential Δ
     * @return ForceComponents struct with spatial, temporal, and total forces
     */
    static ForceComponents decompose(
        double dp_x, double dp_y,
        double dt,
        double grad_R_x, double grad_R_y,
        double delta);

    /**
     * Compute spatial R-field gradient at position
     *
     * Uses centered finite differences:
     * ∂R/∂x ≈ (R(x+dx) - R(x-dx)) / (2·dx)
     *
     * @param R_field Synchronization field R(x,y) [size: Nx*Ny]
     * @param Nx Grid width
     * @param Ny Grid height
     * @param pos_x Position x (grid units)
     * @param pos_y Position y (grid units)
     * @param dx Grid spacing (default: 1.0)
     * @return (∇R_x, ∇R_y) tuple
     */
    static std::tuple<double, double> computeSpatialGradient(
        const std::vector<float>& R_field,
        int Nx, int Ny,
        double pos_x, double pos_y,
        double dx = 1.0);

    /**
     * Compute temporal R-field gradient at position
     *
     * Uses time derivative field: ∂R/∂t
     * Then computes spatial gradient: ∇(∂R/∂t)
     *
     * @param dR_dt_field Time derivative field ∂R/∂t [size: Nx*Ny]
     * @param Nx Grid width
     * @param Ny Grid height
     * @param pos_x Position x (grid units)
     * @param pos_y Position y (grid units)
     * @param dx Grid spacing (default: 1.0)
     * @return (∇(∂R/∂t)_x, ∇(∂R/∂t)_y) tuple
     */
    static std::tuple<double, double> computeTemporalGradient(
        const std::vector<float>& dR_dt_field,
        int Nx, int Ny,
        double pos_x, double pos_y,
        double dx = 1.0);

    /**
     * Compute correlation between temporal force and R-field time derivative
     *
     * Physics: F_temp should correlate with ∂R/∂t
     * Strong correlation (ρ > 0.5) validates temporal force hypothesis
     *
     * @param F_temporal_series Time series of temporal force magnitudes
     * @param dR_dt_series Time series of ∂R/∂t at particle center
     * @return Pearson correlation coefficient ρ ∈ [-1, 1]
     */
    static double computeCorrelation(
        const std::vector<double>& F_temporal_series,
        const std::vector<double>& dR_dt_series);
};
