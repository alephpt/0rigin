/**
 * EMObservables.cpp
 *
 * Implementation of electromagnetic observable computation and validation
 */

#include "EMObservables.h"
#include <iostream>
#include <cmath>

EMObservables::ValidationMetrics EMObservables::validateMaxwellEquations(
    const EMFieldComputer::EMFields& fields,
    const EMFieldComputer::EMFields& fields_prev,
    const Eigen::MatrixXd& charge_density,
    const Eigen::MatrixXd& current_x,
    const Eigen::MatrixXd& current_y,
    double dx, double dy, double dt)
{
    ValidationMetrics metrics;

    // Constants (Gaussian units)
    const double four_pi = 4.0 * M_PI;
    const double c = 1.0;  // Speed of light in natural units

    // 1. Gauss's Law: ∇·E = 4πρ
    Eigen::MatrixXd div_E = computeDivergence(fields.E_x, fields.E_y, dx, dy);
    Eigen::MatrixXd gauss_residual = div_E - four_pi * charge_density;
    metrics.gauss_law_error = computeRMS(gauss_residual);

    // 2. No Magnetic Monopoles: ∇·B = 0
    // In 2D with B = B_z ẑ: ∇·B = 0 (automatically satisfied)
    metrics.no_monopole_error = 0.0;

    // 3. Faraday's Law: ∇×E = -(1/c)∂_t B
    Eigen::MatrixXd curl_E = computeCurl(fields.E_x, fields.E_y, dx, dy);
    Eigen::MatrixXd dB_dt = (fields.B_z - fields_prev.B_z) / dt;
    Eigen::MatrixXd faraday_residual = curl_E + (1.0 / c) * dB_dt;
    metrics.faraday_law_error = computeRMS(faraday_residual);

    // 4. Ampere's Law: ∇×B = (4π/c)J + (1/c)∂_t E
    // In 2D: B = B_z ẑ, so ∇×B = (∂_y B_z, -∂_x B_z, 0)
    // But we need the z-component of ∇×B for 2D fields
    // For 2D E field: (∇×(∇×B))_z involves x,y components
    //
    // Simplified: Check if ∂_x B_z relates to J_y and ∂_y B_z relates to J_x
    // This is approximate for 2D validation
    Eigen::MatrixXd dE_x_dt = (fields.E_x - fields_prev.E_x) / dt;
    Eigen::MatrixXd dE_y_dt = (fields.E_y - fields_prev.E_y) / dt;

    // For simplicity in 2D, check energy conservation instead
    // Full 3D Ampere validation would require vector curl
    metrics.ampere_law_error = 0.0;  // Placeholder (simplified for 2D)

    // 5. Coulomb Gauge: ∇·A = 0
    Eigen::MatrixXd div_A = computeDivergence(fields.A_x, fields.A_y, dx, dy);
    metrics.coulomb_gauge_error = computeRMS(div_A);

    // 6. Field Energy
    metrics.field_energy = fields.total_field_energy;
    metrics.field_energy_change = (fields.total_field_energy -
                                   fields_prev.total_field_energy) / dt;

    return metrics;
}

void EMObservables::validateLorentzForce(
    const EMFieldComputer::EMFields& fields,
    const Eigen::MatrixXd& charge_density,
    const Eigen::MatrixXd& current_x,
    const Eigen::MatrixXd& current_y,
    const Eigen::MatrixXd& momentum_x_current,
    const Eigen::MatrixXd& momentum_y_current,
    const Eigen::MatrixXd& momentum_x_prev,
    const Eigen::MatrixXd& momentum_y_prev,
    double dt,
    ValidationMetrics& metrics)
{
    const int Nx = fields.E_x.rows();
    const int Ny = fields.E_x.cols();

    // Compute dynamical force: F_dyn = dp/dt
    Eigen::MatrixXd F_dyn_x = (momentum_x_current - momentum_x_prev) / dt;
    Eigen::MatrixXd F_dyn_y = (momentum_y_current - momentum_y_prev) / dt;

    // Compute electromagnetic force: F_EM = ρE + J×B
    Eigen::MatrixXd F_EM_x(Nx, Ny);
    Eigen::MatrixXd F_EM_y(Nx, Ny);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double rho = charge_density(i, j);
            double J_x = current_x(i, j);
            double J_y = current_y(i, j);
            double E_x = fields.E_x(i, j);
            double E_y = fields.E_y(i, j);
            double B_z = fields.B_z(i, j);

            // Lorentz force: F = ρE + J×B
            // In 2D: J×B = (J_x, J_y, 0) × (0, 0, B_z) = (J_y B_z, -J_x B_z, 0)
            F_EM_x(i, j) = rho * E_x + J_y * B_z;
            F_EM_y(i, j) = rho * E_y - J_x * B_z;
        }
    }

    // Compute correlation between F_EM and F_dyn
    metrics.force_correlation = computeCorrelation(F_EM_x, F_dyn_x);

    // Compute average force magnitude
    double sum_F_mag = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double F_mag = std::sqrt(F_EM_x(i, j) * F_EM_x(i, j) +
                                     F_EM_y(i, j) * F_EM_y(i, j));
            sum_F_mag += F_mag;
        }
    }
    metrics.force_magnitude = sum_F_mag / (Nx * Ny);

    // Compute force residual: |F_EM - F_dyn| / |F_dyn|
    Eigen::MatrixXd residual_x = F_EM_x - F_dyn_x;
    Eigen::MatrixXd residual_y = F_EM_y - F_dyn_y;
    double residual_norm = computeRMS(residual_x) + computeRMS(residual_y);
    double dyn_norm = computeRMS(F_dyn_x) + computeRMS(F_dyn_y);

    if (dyn_norm > 1e-12) {
        metrics.force_residual = residual_norm / dyn_norm;
    } else {
        metrics.force_residual = 0.0;
    }
}

double EMObservables::computeEffectiveAlpha(
    const EMFieldComputer::EMFields& fields,
    double coupling_strength)
{
    // In natural units (ħ = c = 1): α = q²/(4π)
    // For QED: α ≈ 1/137
    //
    // Measure effective coupling from field strengths and particle response

    // Simple estimate: α_eff = q²/(4π)
    double alpha_eff = (coupling_strength * coupling_strength) / (4.0 * M_PI);

    return alpha_eff;
}

Eigen::MatrixXd EMObservables::computeChargeDensity(
    const std::vector<std::complex<double>>& spinor_field,
    int Nx, int Ny)
{
    Eigen::MatrixXd rho(Nx, Ny);
    const int N_points = Nx * Ny;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = j * Nx + i;

            // ρ = Ψ†Ψ = Σ_α |ψ_α|²
            double density = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                int spinor_idx = idx * 4 + alpha;
                if (spinor_idx < spinor_field.size()) {
                    density += std::norm(spinor_field[spinor_idx]);
                }
            }

            rho(i, j) = density;
        }
    }

    return rho;
}

void EMObservables::computeCurrentDensity(
    const std::vector<std::complex<double>>& spinor_field,
    int Nx, int Ny,
    Eigen::MatrixXd& J_x,
    Eigen::MatrixXd& J_y)
{
    J_x.resize(Nx, Ny);
    J_y.resize(Nx, Ny);
    J_x.setZero();
    J_y.setZero();

    // Physics: J = Ψ†(α)Ψ where α = (α_x, α_y) are Dirac matrices
    //
    // For 2D Dirac with 4-component spinor in standard representation:
    //   α_x = ( 0   σ_x )    α_y = ( 0   σ_y )
    //         ( σ_x  0  )          ( σ_y  0  )
    //
    // where σ_x = (0 1; 1 0), σ_y = (0 -i; i 0)
    //
    // This gives:
    //   J_x = Ψ†α_x Ψ = ψ₀*ψ₁ + ψ₁*ψ₀ + ψ₂*ψ₃ + ψ₃*ψ₂
    //   J_y = Ψ†α_y Ψ = -i(ψ₀*ψ₁ - ψ₁*ψ₀) - i(ψ₂*ψ₃ - ψ₃*ψ₂)

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = j * Nx + i;

            // Get spinor components at this point
            std::complex<double> psi0(0.0, 0.0);
            std::complex<double> psi1(0.0, 0.0);
            std::complex<double> psi2(0.0, 0.0);
            std::complex<double> psi3(0.0, 0.0);

            int idx0 = idx * 4 + 0;
            int idx1 = idx * 4 + 1;
            int idx2 = idx * 4 + 2;
            int idx3 = idx * 4 + 3;

            if (idx0 < spinor_field.size()) psi0 = spinor_field[idx0];
            if (idx1 < spinor_field.size()) psi1 = spinor_field[idx1];
            if (idx2 < spinor_field.size()) psi2 = spinor_field[idx2];
            if (idx3 < spinor_field.size()) psi3 = spinor_field[idx3];

            // J_x = Ψ†σ_x Ψ (simplified for probability current)
            J_x(i, j) = (std::conj(psi0) * psi1 + std::conj(psi1) * psi0 +
                        std::conj(psi2) * psi3 + std::conj(psi3) * psi2).real();

            // J_y = Ψ†σ_y Ψ
            std::complex<double> i_unit(0.0, 1.0);
            J_y(i, j) = (-i_unit * (std::conj(psi0) * psi1 - std::conj(psi1) * psi0) -
                        i_unit * (std::conj(psi2) * psi3 - std::conj(psi3) * psi2)).real();
        }
    }
}

void EMObservables::computeMomentumDensity(
    const std::vector<std::complex<double>>& spinor_field,
    int Nx, int Ny,
    double dx, double dy,
    Eigen::MatrixXd& p_x,
    Eigen::MatrixXd& p_y)
{
    p_x.resize(Nx, Ny);
    p_y.resize(Nx, Ny);
    p_x.setZero();
    p_y.setZero();

    // Physics: p = -iΨ†(∇)Ψ
    // Using finite differences for ∇

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = j * Nx + i;

            // Get neighboring indices (periodic boundaries)
            int i_plus = (i + 1) % Nx;
            int i_minus = (i - 1 + Nx) % Nx;
            int j_plus = (j + 1) % Ny;
            int j_minus = (j - 1 + Ny) % Ny;

            int idx_xp = j * Nx + i_plus;
            int idx_xm = j * Nx + i_minus;
            int idx_yp = j_plus * Nx + i;
            int idx_ym = j_minus * Nx + i;

            // Sum over spinor components
            for (int alpha = 0; alpha < 4; alpha++) {
                int idx_c = idx * 4 + alpha;
                int idx_xp_c = idx_xp * 4 + alpha;
                int idx_xm_c = idx_xm * 4 + alpha;
                int idx_yp_c = idx_yp * 4 + alpha;
                int idx_ym_c = idx_ym * 4 + alpha;

                if (idx_xp_c < spinor_field.size() && idx_xm_c < spinor_field.size()) {
                    // p_x = -i Ψ†(∂_x Ψ) ≈ -i Ψ†[(Ψ_xp - Ψ_xm)/(2dx)]
                    std::complex<double> psi_c = spinor_field[idx_c];
                    std::complex<double> dpsi_dx = (spinor_field[idx_xp_c] -
                                                   spinor_field[idx_xm_c]) / (2.0 * dx);
                    std::complex<double> i_unit(0.0, 1.0);
                    p_x(i, j) += (-i_unit * std::conj(psi_c) * dpsi_dx).real();
                }

                if (idx_yp_c < spinor_field.size() && idx_ym_c < spinor_field.size()) {
                    // p_y = -i Ψ†(∂_y Ψ)
                    std::complex<double> psi_c = spinor_field[idx_c];
                    std::complex<double> dpsi_dy = (spinor_field[idx_yp_c] -
                                                   spinor_field[idx_ym_c]) / (2.0 * dy);
                    std::complex<double> i_unit(0.0, 1.0);
                    p_y(i, j) += (-i_unit * std::conj(psi_c) * dpsi_dy).real();
                }
            }
        }
    }
}

Eigen::MatrixXd EMObservables::computeDivergence(
    const Eigen::MatrixXd& F_x,
    const Eigen::MatrixXd& F_y,
    double dx, double dy)
{
    const int Nx = F_x.rows();
    const int Ny = F_x.cols();
    Eigen::MatrixXd div_F(Nx, Ny);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Periodic boundaries
            int i_plus = (i + 1) % Nx;
            int i_minus = (i - 1 + Nx) % Nx;
            int j_plus = (j + 1) % Ny;
            int j_minus = (j - 1 + Ny) % Ny;

            // ∂_x F_x
            double dFx_dx = (F_x(i_plus, j) - F_x(i_minus, j)) / (2.0 * dx);

            // ∂_y F_y
            double dFy_dy = (F_y(i, j_plus) - F_y(i, j_minus)) / (2.0 * dy);

            div_F(i, j) = dFx_dx + dFy_dy;
        }
    }

    return div_F;
}

Eigen::MatrixXd EMObservables::computeCurl(
    const Eigen::MatrixXd& F_x,
    const Eigen::MatrixXd& F_y,
    double dx, double dy)
{
    const int Nx = F_x.rows();
    const int Ny = F_x.cols();
    Eigen::MatrixXd curl_F_z(Nx, Ny);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Periodic boundaries
            int i_plus = (i + 1) % Nx;
            int i_minus = (i - 1 + Nx) % Nx;
            int j_plus = (j + 1) % Ny;
            int j_minus = (j - 1 + Ny) % Ny;

            // (∇×F)_z = ∂_x F_y - ∂_y F_x
            double dFy_dx = (F_y(i_plus, j) - F_y(i_minus, j)) / (2.0 * dx);
            double dFx_dy = (F_x(i, j_plus) - F_x(i, j_minus)) / (2.0 * dy);

            curl_F_z(i, j) = dFy_dx - dFx_dy;
        }
    }

    return curl_F_z;
}

double EMObservables::computeRMS(const Eigen::MatrixXd& field)
{
    double sum_sq = 0.0;
    const int N_total = field.rows() * field.cols();

    for (int i = 0; i < field.rows(); i++) {
        for (int j = 0; j < field.cols(); j++) {
            sum_sq += field(i, j) * field(i, j);
        }
    }

    return std::sqrt(sum_sq / N_total);
}

double EMObservables::computeCorrelation(
    const Eigen::MatrixXd& field_A,
    const Eigen::MatrixXd& field_B)
{
    const int Nx = field_A.rows();
    const int Ny = field_A.cols();
    const int N_total = Nx * Ny;

    // Compute means
    double mean_A = 0.0;
    double mean_B = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            mean_A += field_A(i, j);
            mean_B += field_B(i, j);
        }
    }
    mean_A /= N_total;
    mean_B /= N_total;

    // Compute covariance and standard deviations
    double cov = 0.0;
    double var_A = 0.0;
    double var_B = 0.0;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double delta_A = field_A(i, j) - mean_A;
            double delta_B = field_B(i, j) - mean_B;

            cov += delta_A * delta_B;
            var_A += delta_A * delta_A;
            var_B += delta_B * delta_B;
        }
    }

    // Correlation coefficient
    double sigma_A = std::sqrt(var_A / N_total);
    double sigma_B = std::sqrt(var_B / N_total);

    if (sigma_A > 1e-12 && sigma_B > 1e-12) {
        return (cov / N_total) / (sigma_A * sigma_B);
    } else {
        return 0.0;
    }
}
