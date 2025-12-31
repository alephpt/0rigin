/**
 * LorentzTransform.h
 *
 * Lorentz transformation methods for boosted frame analysis
 *
 * Physics:
 *   - Coordinate transforms: (x', t') = Λ(β)(x, t)
 *   - Field transforms: R'(x', t'), ψ'(x', t')
 *   - Observable transforms: E', p', m'
 *   - Lorentz invariants: E² - p² = m²
 *
 * Purpose: Enable Phase 2.5B analysis
 *   - Measure R-field in boosted frame (particle rest frame)
 *   - Validate m(v) = γ·Δ·R formula in different frames
 *   - Test Lorentz covariance of SMFT mass formula
 */

#pragma once

#include <eigen3/Eigen/Dense>
#include <complex>
#include <vector>
#include <cmath>
#include <cstdint>

class LorentzTransform {
public:
    LorentzTransform() = default;
    ~LorentzTransform() = default;

    // ========================================================================
    // COORDINATE TRANSFORMS
    // ========================================================================

    /**
     * Get Lorentz boost matrix Λ(β) for 2D+1 spacetime
     *
     * Matrix form (ct', x', y') = Λ(β_x, β_y) · (ct, x, y)
     *
     * For pure boost in x-direction (β_y = 0):
     *   [γ      -γβ     0  ]
     *   [-γβ     γ      0  ]
     *   [0       0      1  ]
     *
     * @param beta_x: Boost velocity in x-direction (units of c)
     * @param beta_y: Boost velocity in y-direction (units of c)
     * @return: 3×3 Lorentz transformation matrix
     */
    Eigen::Matrix3d getBoostMatrix(double beta_x, double beta_y) const;

    /**
     * Get Lorentz gamma factor: γ = 1 / √(1 - β²)
     *
     * @param beta: Total boost velocity magnitude (units of c)
     * @return: Lorentz factor γ ≥ 1
     */
    double getGamma(double beta) const;

    /**
     * Transform single spacetime coordinate
     *
     * Physics: (ct', x', y') = Λ(β) · (ct, x, y)
     *
     * @param ct: Time coordinate (in units where c=1)
     * @param x: Spatial x coordinate
     * @param y: Spatial y coordinate
     * @param beta_x: Boost velocity in x
     * @param beta_y: Boost velocity in y
     * @return: Transformed coordinates (ct', x', y')
     */
    Eigen::Vector3d transformCoordinate(double ct, double x, double y,
                                        double beta_x, double beta_y) const;

    /**
     * Transform grid coordinates for field analysis
     *
     * Given a field R(x,y,t) in lab frame, compute grid coordinates
     * for boosted frame R'(x',y',t').
     *
     * @param x_grid: Lab frame x coordinates [Nx × Ny]
     * @param y_grid: Lab frame y coordinates [Nx × Ny]
     * @param t: Current time in lab frame
     * @param beta_x: Boost velocity in x
     * @param beta_y: Boost velocity in y
     * @return: pair of (x'_grid, y'_grid) in boosted frame
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    transformGridCoordinates(const Eigen::MatrixXd& x_grid,
                            const Eigen::MatrixXd& y_grid,
                            double t,
                            double beta_x,
                            double beta_y) const;

    // ========================================================================
    // FIELD TRANSFORMS
    // ========================================================================

    /**
     * Transform scalar field to boosted frame
     *
     * Physics: Scalar fields transform as R'(x') = R(x)
     * (Lorentz scalars are frame-independent)
     *
     * Note: This interpolates R(x,y,t) onto the boosted grid (x',y',t')
     *
     * @param R_field: Scalar field in lab frame [Nx × Ny]
     * @param x_grid: Lab frame x coordinates
     * @param y_grid: Lab frame y coordinates
     * @param x_prime_grid: Boosted frame x' coordinates
     * @param y_prime_grid: Boosted frame y' coordinates
     * @return: R'(x',y') interpolated to boosted frame grid
     */
    Eigen::MatrixXd transformScalarField(const Eigen::MatrixXd& R_field,
                                         const Eigen::MatrixXd& x_grid,
                                         const Eigen::MatrixXd& y_grid,
                                         const Eigen::MatrixXd& x_prime_grid,
                                         const Eigen::MatrixXd& y_prime_grid) const;

    /**
     * Transform spinor field to boosted frame
     *
     * Physics: Dirac spinors transform as ψ'(x') = S(Λ) ψ(x)
     * where S(Λ) is the spinor representation of the Lorentz boost.
     *
     * For 2D Dirac equation with 4-component spinor:
     *   S(Λ) = exp(-½ωᵢⱼ Σⁱⱼ)
     *   where ωᵢⱼ are boost rapidity parameters
     *   and Σⁱⱼ are spin generators
     *
     * @param psi_field: Spinor field in lab frame [4 × Nx × Ny]
     * @param x_grid: Lab frame x coordinates
     * @param y_grid: Lab frame y coordinates
     * @param x_prime_grid: Boosted frame x' coordinates
     * @param y_prime_grid: Boosted frame y' coordinates
     * @param beta_x: Boost velocity in x
     * @param beta_y: Boost velocity in y
     * @return: ψ'(x',y') transformed to boosted frame
     */
    std::vector<Eigen::MatrixXcd>
    transformSpinorField(const std::vector<std::complex<double>>& psi_field,
                        uint32_t Nx, uint32_t Ny,
                        const Eigen::MatrixXd& x_grid,
                        const Eigen::MatrixXd& y_grid,
                        const Eigen::MatrixXd& x_prime_grid,
                        const Eigen::MatrixXd& y_prime_grid,
                        double beta_x,
                        double beta_y) const;

    // ========================================================================
    // OBSERVABLE TRANSFORMS
    // ========================================================================

    /**
     * Transform energy-momentum 4-vector
     *
     * Physics: (E', pₓ', pᵧ') = Λ(β) · (E, pₓ, pᵧ)
     *
     * Lorentz invariant: E'² - p'² = E² - p² = m²
     *
     * @param E: Energy in lab frame
     * @param px: Momentum in x direction
     * @param py: Momentum in y direction
     * @param beta_x: Boost velocity in x
     * @param beta_y: Boost velocity in y
     * @return: (E', p'ₓ, p'ᵧ) in boosted frame
     */
    Eigen::Vector3d transformEnergyMomentum(double E, double px, double py,
                                            double beta_x, double beta_y) const;

    /**
     * Compute Lorentz invariant mass
     *
     * Physics: m² = E² - p²
     *
     * This should be frame-independent.
     *
     * @param E: Energy
     * @param px: Momentum in x
     * @param py: Momentum in y
     * @return: Invariant mass m
     */
    double computeInvariantMass(double E, double px, double py) const;

    /**
     * Verify Lorentz covariance
     *
     * Tests if physical quantity is correctly transformed:
     *   - Compute quantity in lab frame
     *   - Transform to boosted frame
     *   - Compute quantity in boosted frame
     *   - Compare (should be equal for Lorentz scalars)
     *
     * @param m_lab: Mass in lab frame
     * @param m_boosted: Mass in boosted frame
     * @param tolerance: Relative tolerance (default: 5%)
     * @return: true if |m_lab - m_boosted| / m_lab < tolerance
     */
    bool verifyLorentzCovariance(double m_lab, double m_boosted,
                                 double tolerance = 0.05) const;

    // ========================================================================
    // PHYSICS ANALYSIS
    // ========================================================================

    /**
     * Measure effective mass in boosted frame
     *
     * Physics: m_eff = Δ · R
     *
     * In boosted frame: m'_eff = Δ · R'
     *
     * Lorentz covariance test: m'_eff should equal m_eff
     * (mass is a Lorentz scalar)
     *
     * @param R_boosted: Order parameter in boosted frame [Nx × Ny]
     * @param Delta: Mass gap parameter (typically 1.0)
     * @return: Average effective mass <m'_eff> = Δ · <R'>
     */
    double measureEffectiveMass(const Eigen::MatrixXd& R_boosted,
                               double Delta = 1.0) const;

    /**
     * Compute velocity from particle trajectory
     *
     * Physics: β = dx/dt (in units where c=1)
     *
     * @param x1: Position at time t1
     * @param x2: Position at time t2
     * @param t1: Initial time
     * @param t2: Final time
     * @return: Velocity β = Δx / Δt
     */
    double computeVelocity(double x1, double x2, double t1, double t2) const;

    /**
     * Test gamma factor from measured velocity
     *
     * Physics: γ = 1 / √(1 - v²/c²)
     *
     * Validation: Measure v from trajectory, compute γ,
     *             compare with relativistic mass γ·m₀
     *
     * @param velocity: Measured velocity (units of c)
     * @param gamma_measured: Measured gamma factor from mass
     * @param tolerance: Relative tolerance (default: 5%)
     * @return: true if |γ_computed - γ_measured| / γ_computed < tolerance
     */
    bool testGammaFactor(double velocity, double gamma_measured,
                        double tolerance = 0.05) const;

private:
    // Bilinear interpolation helper for field transforms
    double interpolateBilinear(const Eigen::MatrixXd& field,
                              const Eigen::MatrixXd& x_grid,
                              const Eigen::MatrixXd& y_grid,
                              double x_target,
                              double y_target) const;

    // Spinor boost matrix S(Λ) for Dirac field transformation
    Eigen::Matrix4cd getSpinorBoostMatrix(double beta_x, double beta_y) const;
};
