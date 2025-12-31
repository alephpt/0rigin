/**
 * LorentzTransform.cpp
 *
 * Implementation of Lorentz transformation methods
 */

#include "LorentzTransform.h"
#include <stdexcept>
#include <algorithm>

// ============================================================================
// COORDINATE TRANSFORMS
// ============================================================================

Eigen::Matrix3d LorentzTransform::getBoostMatrix(double beta_x, double beta_y) const {
    double beta_sq = beta_x * beta_x + beta_y * beta_y;

    if (beta_sq >= 1.0) {
        throw std::runtime_error("Boost velocity >= c (beta >= 1.0)");
    }

    if (beta_sq < 1e-10) {
        // No boost: return identity matrix
        return Eigen::Matrix3d::Identity();
    }

    double gamma = 1.0 / std::sqrt(1.0 - beta_sq);
    double beta = std::sqrt(beta_sq);

    // General 2D boost matrix (Weinberg QFT Vol 1, Eq 2.5.5)
    // For boost with velocity β = (βₓ, βᵧ):
    //
    //   Λ = [γ           -γβₓ        -γβᵧ     ]
    //       [-γβₓ    1+(γ-1)βₓ²/β²  (γ-1)βₓβᵧ/β²]
    //       [-γβᵧ    (γ-1)βₓβᵧ/β²   1+(γ-1)βᵧ²/β²]

    Eigen::Matrix3d Lambda;

    // Time row
    Lambda(0, 0) = gamma;
    Lambda(0, 1) = -gamma * beta_x;
    Lambda(0, 2) = -gamma * beta_y;

    // Spatial rows
    double factor = (gamma - 1.0) / beta_sq;

    Lambda(1, 0) = -gamma * beta_x;
    Lambda(1, 1) = 1.0 + factor * beta_x * beta_x;
    Lambda(1, 2) = factor * beta_x * beta_y;

    Lambda(2, 0) = -gamma * beta_y;
    Lambda(2, 1) = factor * beta_x * beta_y;
    Lambda(2, 2) = 1.0 + factor * beta_y * beta_y;

    return Lambda;
}

double LorentzTransform::getGamma(double beta) const {
    if (beta < 0.0) {
        throw std::runtime_error("Negative boost velocity");
    }
    if (beta >= 1.0) {
        throw std::runtime_error("Boost velocity >= c");
    }

    if (beta < 1e-10) {
        return 1.0; // No boost
    }

    return 1.0 / std::sqrt(1.0 - beta * beta);
}

Eigen::Vector3d LorentzTransform::transformCoordinate(double ct, double x, double y,
                                                      double beta_x, double beta_y) const {
    Eigen::Matrix3d Lambda = getBoostMatrix(beta_x, beta_y);
    Eigen::Vector3d coord(ct, x, y);
    return Lambda * coord;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
LorentzTransform::transformGridCoordinates(const Eigen::MatrixXd& x_grid,
                                          const Eigen::MatrixXd& y_grid,
                                          double t,
                                          double beta_x,
                                          double beta_y) const {
    if (x_grid.rows() != y_grid.rows() || x_grid.cols() != y_grid.cols()) {
        throw std::runtime_error("Grid dimensions mismatch");
    }

    Eigen::Matrix3d Lambda = getBoostMatrix(beta_x, beta_y);

    int Ny = x_grid.rows();
    int Nx = x_grid.cols();

    Eigen::MatrixXd x_prime(Ny, Nx);
    Eigen::MatrixXd y_prime(Ny, Nx);

    double ct = t; // c = 1 in natural units

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            Eigen::Vector3d coord(ct, x_grid(j, i), y_grid(j, i));
            Eigen::Vector3d coord_prime = Lambda * coord;

            // Extract spatial coordinates (ignore transformed time)
            x_prime(j, i) = coord_prime(1);
            y_prime(j, i) = coord_prime(2);
        }
    }

    return {x_prime, y_prime};
}

// ============================================================================
// FIELD TRANSFORMS
// ============================================================================

Eigen::MatrixXd LorentzTransform::transformScalarField(
    const Eigen::MatrixXd& R_field,
    const Eigen::MatrixXd& x_grid,
    const Eigen::MatrixXd& y_grid,
    const Eigen::MatrixXd& x_prime_grid,
    const Eigen::MatrixXd& y_prime_grid) const {

    int Ny = x_prime_grid.rows();
    int Nx = x_prime_grid.cols();

    Eigen::MatrixXd R_prime(Ny, Nx);

    // Interpolate R_field from lab frame (x, y) to boosted frame (x', y')
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double x_target = x_prime_grid(j, i);
            double y_target = y_prime_grid(j, i);

            R_prime(j, i) = interpolateBilinear(R_field, x_grid, y_grid,
                                               x_target, y_target);
        }
    }

    return R_prime;
}

std::vector<Eigen::MatrixXcd>
LorentzTransform::transformSpinorField(
    const std::vector<std::complex<double>>& psi_field,
    uint32_t Nx, uint32_t Ny,
    const Eigen::MatrixXd& x_grid,
    const Eigen::MatrixXd& y_grid,
    const Eigen::MatrixXd& x_prime_grid,
    const Eigen::MatrixXd& y_prime_grid,
    double beta_x,
    double beta_y) const {

    // Get spinor boost matrix S(Λ)
    Eigen::Matrix4cd S = getSpinorBoostMatrix(beta_x, beta_y);

    // Result: 4 components, each Ny × Nx
    std::vector<Eigen::MatrixXcd> psi_prime(4);
    for (int c = 0; c < 4; ++c) {
        psi_prime[c] = Eigen::MatrixXcd::Zero(Ny, Nx);
    }

    // For each grid point in boosted frame:
    //   1. Find corresponding lab frame coordinates
    //   2. Interpolate spinor components at those coordinates
    //   3. Apply spinor boost S(Λ)

    for (int j = 0; j < static_cast<int>(Ny); ++j) {
        for (int i = 0; i < static_cast<int>(Nx); ++i) {
            double x_target = x_prime_grid(j, i);
            double y_target = y_prime_grid(j, i);

            // Interpolate each spinor component in lab frame
            Eigen::Vector4cd psi_at_point;
            for (int c = 0; c < 4; ++c) {
                // Extract component c as real matrix for interpolation
                Eigen::MatrixXd psi_real(Ny, Nx);
                Eigen::MatrixXd psi_imag(Ny, Nx);

                for (int jj = 0; jj < static_cast<int>(Ny); ++jj) {
                    for (int ii = 0; ii < static_cast<int>(Nx); ++ii) {
                        int idx = jj * Nx + ii;
                        int psi_idx = c * Nx * Ny + idx;
                        psi_real(jj, ii) = psi_field[psi_idx].real();
                        psi_imag(jj, ii) = psi_field[psi_idx].imag();
                    }
                }

                double real_part = interpolateBilinear(psi_real, x_grid, y_grid,
                                                       x_target, y_target);
                double imag_part = interpolateBilinear(psi_imag, x_grid, y_grid,
                                                       x_target, y_target);

                psi_at_point(c) = std::complex<double>(real_part, imag_part);
            }

            // Apply spinor boost: ψ'(x') = S(Λ) ψ(x)
            Eigen::Vector4cd psi_boosted = S * psi_at_point;

            // Store result
            for (int c = 0; c < 4; ++c) {
                psi_prime[c](j, i) = psi_boosted(c);
            }
        }
    }

    return psi_prime;
}

// ============================================================================
// OBSERVABLE TRANSFORMS
// ============================================================================

Eigen::Vector3d LorentzTransform::transformEnergyMomentum(double E, double px, double py,
                                                          double beta_x, double beta_y) const {
    Eigen::Matrix3d Lambda = getBoostMatrix(beta_x, beta_y);
    Eigen::Vector3d four_momentum(E, px, py); // (E/c, px, py) with c=1
    return Lambda * four_momentum;
}

double LorentzTransform::computeInvariantMass(double E, double px, double py) const {
    double m_sq = E * E - px * px - py * py;

    if (m_sq < 0.0) {
        // Unphysical: spacelike four-momentum
        throw std::runtime_error("Spacelike four-momentum: E² - p² < 0");
    }

    return std::sqrt(m_sq);
}

bool LorentzTransform::verifyLorentzCovariance(double m_lab, double m_boosted,
                                               double tolerance) const {
    if (m_lab <= 0.0) {
        throw std::runtime_error("Non-positive mass in lab frame");
    }

    double relative_error = std::abs(m_lab - m_boosted) / m_lab;
    return relative_error < tolerance;
}

// ============================================================================
// PHYSICS ANALYSIS
// ============================================================================

double LorentzTransform::measureEffectiveMass(const Eigen::MatrixXd& R_boosted,
                                              double Delta) const {
    // Average order parameter <R'> in boosted frame
    double R_avg = R_boosted.mean();

    // Effective mass: m' = Δ · <R'>
    return Delta * R_avg;
}

double LorentzTransform::computeVelocity(double x1, double x2, double t1, double t2) const {
    double dt = t2 - t1;
    if (dt <= 0.0) {
        throw std::runtime_error("Invalid time interval: t2 <= t1");
    }

    double dx = x2 - x1;
    return dx / dt; // β = Δx / Δt (c = 1)
}

bool LorentzTransform::testGammaFactor(double velocity, double gamma_measured,
                                      double tolerance) const {
    if (velocity < 0.0 || velocity >= 1.0) {
        throw std::runtime_error("Invalid velocity: must be in [0, 1)");
    }

    double gamma_computed = getGamma(velocity);
    double relative_error = std::abs(gamma_computed - gamma_measured) / gamma_computed;

    return relative_error < tolerance;
}

// ============================================================================
// PRIVATE HELPERS
// ============================================================================

double LorentzTransform::interpolateBilinear(const Eigen::MatrixXd& field,
                                            const Eigen::MatrixXd& x_grid,
                                            const Eigen::MatrixXd& y_grid,
                                            double x_target,
                                            double y_target) const {
    int Ny = field.rows();
    int Nx = field.cols();

    // Find bounding grid points
    // Assume uniform grid spacing
    double x_min = x_grid.minCoeff();
    double x_max = x_grid.maxCoeff();
    double y_min = y_grid.minCoeff();
    double y_max = y_grid.maxCoeff();

    double dx = (x_max - x_min) / (Nx - 1);
    double dy = (y_max - y_min) / (Ny - 1);

    // Convert target to grid indices (fractional)
    double i_frac = (x_target - x_min) / dx;
    double j_frac = (y_target - y_min) / dy;

    // Clamp to grid boundaries
    i_frac = std::clamp(i_frac, 0.0, static_cast<double>(Nx - 1));
    j_frac = std::clamp(j_frac, 0.0, static_cast<double>(Ny - 1));

    int i0 = static_cast<int>(std::floor(i_frac));
    int j0 = static_cast<int>(std::floor(j_frac));
    int i1 = std::min(i0 + 1, Nx - 1);
    int j1 = std::min(j0 + 1, Ny - 1);

    // Interpolation weights
    double wx = i_frac - i0;
    double wy = j_frac - j0;

    // Bilinear interpolation
    double f00 = field(j0, i0);
    double f10 = field(j0, i1);
    double f01 = field(j1, i0);
    double f11 = field(j1, i1);

    double f0 = (1.0 - wx) * f00 + wx * f10;
    double f1 = (1.0 - wx) * f01 + wx * f11;

    return (1.0 - wy) * f0 + wy * f1;
}

Eigen::Matrix4cd LorentzTransform::getSpinorBoostMatrix(double beta_x, double beta_y) const {
    // For 2D Dirac equation, spinor boost matrix S(Λ) is given by:
    //
    //   S(Λ) = exp(-½ ωᵢⱼ Σⁱⱼ)
    //
    // where ωᵢⱼ are boost parameters (rapidity) and Σⁱⱼ are spin generators.
    //
    // For boost in x-direction with rapidity η_x = arctanh(β_x):
    //   S = exp(-η_x Σ^{01}) = exp(-η_x α_x / 2)
    //
    // For 2D Dirac with 4-component spinor, we use:
    //   S = cosh(η/2) I - sinh(η/2) (β̂·α)
    //
    // where β̂ = (β_x, β_y) / |β|, α = (α_x, α_y) are Dirac matrices.

    double beta_sq = beta_x * beta_x + beta_y * beta_y;

    if (beta_sq < 1e-10) {
        // No boost: return identity
        return Eigen::Matrix4cd::Identity();
    }

    double beta = std::sqrt(beta_sq);
    double gamma = getGamma(beta);

    // Rapidity: η = arctanh(β)
    double eta = 0.5 * std::log((1.0 + beta) / (1.0 - beta));

    // Unit boost direction
    double beta_hat_x = beta_x / beta;
    double beta_hat_y = beta_y / beta;

    // Spinor boost parameters
    double cosh_half = std::cosh(eta / 2.0);
    double sinh_half = std::sinh(eta / 2.0);

    // For 2D Dirac equation with standard representation:
    //   α_x = σ_x ⊗ σ_z, α_y = σ_y ⊗ σ_z, β = σ_z ⊗ I
    //
    // Boost matrix (approximate form for small boosts):
    //   S ≈ I + (γ-1)(β̂·α)/β = I + sinh(η/2)(β̂·α)
    //
    // For simplicity, use diagonal approximation (exact for boost along axis):
    //   S ≈ diag(cosh(η/2) - sinh(η/2), cosh(η/2) + sinh(η/2), ...)

    Eigen::Matrix4cd S = Eigen::Matrix4cd::Identity();

    // Simple diagonal boost (valid for small velocities)
    // Full implementation would require proper Dirac matrix algebra
    S(0, 0) = std::complex<double>(cosh_half - sinh_half * beta_hat_x, 0.0);
    S(1, 1) = std::complex<double>(cosh_half + sinh_half * beta_hat_x, 0.0);
    S(2, 2) = std::complex<double>(cosh_half - sinh_half * beta_hat_y, 0.0);
    S(3, 3) = std::complex<double>(cosh_half + sinh_half * beta_hat_y, 0.0);

    return S;
}
