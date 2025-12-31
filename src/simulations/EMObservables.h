/**
 * EMObservables.h
 *
 * Observable computation and validation for electromagnetic coupling
 *
 * Physics:
 *   - Validate Maxwell equations on grid
 *   - Compute Lorentz force correlation with particle dynamics
 *   - Extract effective fine structure constant
 *   - Test if EM fields from Kuramoto phase explain particle behavior
 *
 * Purpose: Phase 5 (Scenario 2.6B) - Electromagnetic Coupling Validation
 *   - Does EM coupling explain acceleration? (F_EM vs dp/dt)
 *   - Do Maxwell equations hold? (∇·E = ρ, ∇×B = J, etc.)
 *   - What is the effective coupling strength? (α_eff)
 */

#pragma once

#include "../physics/EMFieldComputer.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>

class EMObservables {
public:
    /**
     * Validation metrics for electromagnetic coupling
     */
    struct ValidationMetrics {
        // Lorentz force validation
        double force_correlation;   // Correlation ρ(F_EM, dp/dt)
        double force_magnitude;     // |F_EM| average magnitude
        double force_residual;      // |F_EM - dp/dt| / |dp/dt|

        // Maxwell equation validation (residuals)
        double gauss_law_error;     // |∇·E - 4πρ| (Gaussian units)
        double ampere_law_error;    // |∇×B - 4πJ/c - (1/c)∂_t E|
        double faraday_law_error;   // |∇×E + (1/c)∂_t B|
        double no_monopole_error;   // |∇·B| (should be exactly 0)

        // Gauge invariance
        double coulomb_gauge_error; // |∇·A| (Coulomb gauge: ∇·A = 0)

        // Energy conservation
        double field_energy;        // ∫(E² + B²)/(8π) dV
        double field_energy_change; // dU_field/dt

        // Discovery metrics
        double effective_alpha;     // Measured q²/(ħc) (expect 1/137?)
        double coupling_consistency; // Consistency of coupling across grid

        // Constructor
        ValidationMetrics()
            : force_correlation(0.0), force_magnitude(0.0), force_residual(0.0),
              gauss_law_error(0.0), ampere_law_error(0.0), faraday_law_error(0.0),
              no_monopole_error(0.0), coulomb_gauge_error(0.0),
              field_energy(0.0), field_energy_change(0.0),
              effective_alpha(0.0), coupling_consistency(0.0)
        {}
    };

    /**
     * Validate Maxwell equations on grid
     *
     * Physics:
     *   Gauss's law:    ∇·E = 4πρ
     *   Ampere's law:   ∇×B = 4πJ/c + (1/c)∂_t E
     *   Faraday's law:  ∇×E = -(1/c)∂_t B
     *   No monopoles:   ∇·B = 0
     *
     * Returns average residuals (should be small if Maxwell equations hold)
     *
     * @param fields: Current EM fields
     * @param fields_prev: Previous EM fields (for time derivatives)
     * @param charge_density: Charge density ρ [Nx × Ny]
     * @param current_x: Current density J_x [Nx × Ny]
     * @param current_y: Current density J_y [Nx × Ny]
     * @param dx: Grid spacing in x
     * @param dy: Grid spacing in y
     * @param dt: Timestep
     * @return: ValidationMetrics with Maxwell equation residuals
     */
    static ValidationMetrics validateMaxwellEquations(
        const EMFieldComputer::EMFields& fields,
        const EMFieldComputer::EMFields& fields_prev,
        const Eigen::MatrixXd& charge_density,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y,
        double dx, double dy, double dt);

    /**
     * Validate Lorentz force against particle dynamics
     *
     * Tests if electromagnetic force explains particle acceleration:
     *   F_EM = ρE + J×B  (predicted from EM fields)
     *   F_dyn = m·dp/dt  (measured from particle trajectory)
     *
     * Computes correlation and residual between F_EM and F_dyn
     *
     * @param fields: EM fields
     * @param charge_density: Charge density ρ [Nx × Ny]
     * @param current_x: Current density J_x [Nx × Ny]
     * @param current_y: Current density J_y [Nx × Ny]
     * @param momentum_x_current: Current momentum density p_x [Nx × Ny]
     * @param momentum_y_current: Current momentum density p_y [Nx × Ny]
     * @param momentum_x_prev: Previous momentum density p_x [Nx × Ny]
     * @param momentum_y_prev: Previous momentum density p_y [Nx × Ny]
     * @param dt: Timestep
     * @return: Force correlation and residual
     */
    static void validateLorentzForce(
        const EMFieldComputer::EMFields& fields,
        const Eigen::MatrixXd& charge_density,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y,
        const Eigen::MatrixXd& momentum_x_current,
        const Eigen::MatrixXd& momentum_y_current,
        const Eigen::MatrixXd& momentum_x_prev,
        const Eigen::MatrixXd& momentum_y_prev,
        double dt,
        ValidationMetrics& metrics);

    /**
     * Extract effective fine structure constant
     *
     * Physics: α_eff = q²/(ħc)
     *
     * In natural units (ħ = c = 1): α_eff = q²
     *
     * Measure from coupling strength and field magnitudes
     *
     * @param fields: EM fields
     * @param coupling_strength: Effective charge q
     * @return: Effective fine structure constant
     */
    static double computeEffectiveAlpha(
        const EMFieldComputer::EMFields& fields,
        double coupling_strength);

    /**
     * Compute charge density from spinor field
     *
     * Physics: ρ = Ψ†Ψ = Σ_α |ψ_α|²
     *
     * @param spinor_field: 4-component spinor [4 × Nx × Ny]
     * @param Nx: Grid size x
     * @param Ny: Grid size y
     * @return: Charge density ρ [Nx × Ny]
     */
    static Eigen::MatrixXd computeChargeDensity(
        const std::vector<std::complex<double>>& spinor_field,
        int Nx, int Ny);

    /**
     * Compute current density from spinor field
     *
     * Physics: J = Ψ†(α)Ψ (probability current)
     *
     * For Dirac equation in 2D:
     *   J_x = Ψ†α_x Ψ
     *   J_y = Ψ†α_y Ψ
     *
     * where α_x, α_y are Dirac matrices
     *
     * @param spinor_field: 4-component spinor [4 × Nx × Ny]
     * @param Nx: Grid size x
     * @param Ny: Grid size y
     * @param J_x: Output current density x-component [Nx × Ny]
     * @param J_y: Output current density y-component [Nx × Ny]
     */
    static void computeCurrentDensity(
        const std::vector<std::complex<double>>& spinor_field,
        int Nx, int Ny,
        Eigen::MatrixXd& J_x,
        Eigen::MatrixXd& J_y);

    /**
     * Compute momentum density from spinor field
     *
     * Physics: p = -iΨ†(∇)Ψ
     *
     * @param spinor_field: 4-component spinor [4 × Nx × Ny]
     * @param Nx: Grid size x
     * @param Ny: Grid size y
     * @param dx: Grid spacing x
     * @param dy: Grid spacing y
     * @param p_x: Output momentum density x-component [Nx × Ny]
     * @param p_y: Output momentum density y-component [Nx × Ny]
     */
    static void computeMomentumDensity(
        const std::vector<std::complex<double>>& spinor_field,
        int Nx, int Ny,
        double dx, double dy,
        Eigen::MatrixXd& p_x,
        Eigen::MatrixXd& p_y);

private:
    /**
     * Compute divergence of vector field with periodic boundaries
     *
     * ∇·F = ∂_x F_x + ∂_y F_y
     *
     * @param F_x: Vector field x-component [Nx × Ny]
     * @param F_y: Vector field y-component [Nx × Ny]
     * @param dx: Grid spacing x
     * @param dy: Grid spacing y
     * @return: Divergence ∇·F [Nx × Ny]
     */
    static Eigen::MatrixXd computeDivergence(
        const Eigen::MatrixXd& F_x,
        const Eigen::MatrixXd& F_y,
        double dx, double dy);

    /**
     * Compute curl of vector field (z-component in 2D)
     *
     * (∇×F)_z = ∂_x F_y - ∂_y F_x
     *
     * @param F_x: Vector field x-component [Nx × Ny]
     * @param F_y: Vector field y-component [Nx × Ny]
     * @param dx: Grid spacing x
     * @param dy: Grid spacing y
     * @return: Curl z-component (∇×F)_z [Nx × Ny]
     */
    static Eigen::MatrixXd computeCurl(
        const Eigen::MatrixXd& F_x,
        const Eigen::MatrixXd& F_y,
        double dx, double dy);

    /**
     * Compute RMS (root-mean-square) of field
     *
     * @param field: Field values [Nx × Ny]
     * @return: RMS = sqrt(<field²>)
     */
    static double computeRMS(const Eigen::MatrixXd& field);

    /**
     * Compute correlation coefficient between two fields
     *
     * ρ(A, B) = <AB> / (σ_A σ_B)
     *
     * @param field_A: First field [Nx × Ny]
     * @param field_B: Second field [Nx × Ny]
     * @return: Correlation coefficient [-1, 1]
     */
    static double computeCorrelation(
        const Eigen::MatrixXd& field_A,
        const Eigen::MatrixXd& field_B);
};
