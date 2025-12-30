/**
 * EMFieldComputer.cpp
 *
 * Implementation of electromagnetic field extraction from Kuramoto phase
 */

#include "EMFieldComputer.h"
#include <cmath>
#include <algorithm>

EMFieldComputer::EMFields EMFieldComputer::computeFromPhase(
    const Eigen::MatrixXd& theta_current,
    const Eigen::MatrixXd& theta_previous,
    const Eigen::MatrixXd& R_field,
    double dx, double dy, double dt,
    RegularizationType reg_type)
{
    const int Nx = theta_current.rows();
    const int Ny = theta_current.cols();

    EMFields fields(Nx, Ny);

    // Compute temporal derivative: φ = ∂_t θ
    // Use backward difference: ∂_t θ ≈ (θ_current - θ_previous) / dt
    // CRITICAL FIX: Wrap phase difference to [-π, π] to handle branch cuts
    Eigen::MatrixXd phase_diff = theta_current - theta_previous;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double diff = phase_diff(i, j);
            // Wrap to [-π, π]
            while (diff > M_PI) diff -= 2.0 * M_PI;
            while (diff < -M_PI) diff += 2.0 * M_PI;
            phase_diff(i, j) = diff;
        }
    }
    fields.phi = phase_diff / dt;

    // Compute spatial derivatives: A = ∇θ (with regularization)
    // Use conjugate product method for numerical stability at vortex cores
    // A_x = ∂_x θ, A_y = ∂_y θ
    fields.A_x = computePhaseGradientConjugateProduct(R_field, theta_current, dx, 0, reg_type);
    fields.A_y = computePhaseGradientConjugateProduct(R_field, theta_current, dy, 1, reg_type);

    // Compute field strengths E, B
    computeFieldStrengths(fields, dx, dy);

    // Compute diagnostics
    computeDiagnostics(fields);

    // Compute total field energy
    fields.total_field_energy = computeFieldEnergy(fields, dx, dy);

    return fields;
}

EMFieldComputer::EMFields EMFieldComputer::computeFromPhase(
    const std::vector<float>& theta_current,
    const std::vector<float>& theta_previous,
    const std::vector<float>& R_current,
    int Nx, int Ny,
    double dx, double dy, double dt,
    RegularizationType reg_type)
{
    // Convert std::vector<float> to Eigen::MatrixXd (row-major layout)
    // std::vector layout: theta[j*Nx + i] = value at grid point (i,j)
    // Eigen::MatrixXd layout: mat(i,j) = value at grid point (i,j)

    Eigen::MatrixXd theta_curr_mat(Nx, Ny);
    Eigen::MatrixXd theta_prev_mat(Nx, Ny);
    Eigen::MatrixXd R_curr_mat(Nx, Ny);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = j * Nx + i;  // Row-major indexing
            theta_curr_mat(i, j) = static_cast<double>(theta_current[idx]);
            theta_prev_mat(i, j) = static_cast<double>(theta_previous[idx]);
            R_curr_mat(i, j) = static_cast<double>(R_current[idx]);
        }
    }

    // Call Eigen version with regularization
    return computeFromPhase(theta_curr_mat, theta_prev_mat, R_curr_mat, dx, dy, dt, reg_type);
}

void EMFieldComputer::computeFieldStrengths(
    EMFields& fields,
    double dx, double dy)
{
    const int Nx = fields.phi.rows();
    const int Ny = fields.phi.cols();

    // Compute temporal derivatives of vector potential
    // Note: For first timestep, we can't compute ∂_t A accurately
    // For now, we'll use spatial approach: E = -∇φ - ∂_t A
    // We approximate ∂_t A ≈ 0 for first order (can improve later)

    // E = -∇φ (neglecting ∂_t A term for simplicity in first implementation)
    // E_x = -∂_x φ
    fields.E_x = -computeSpatialDerivative(fields.phi, dx, 0);

    // E_y = -∂_y φ
    fields.E_y = -computeSpatialDerivative(fields.phi, dy, 1);

    // B = ∇×A (in 2D: B_z = ∂_x A_y - ∂_y A_x)
    Eigen::MatrixXd dAy_dx = computeSpatialDerivative(fields.A_y, dx, 0);
    Eigen::MatrixXd dAx_dy = computeSpatialDerivative(fields.A_x, dy, 1);
    fields.B_z = dAy_dx - dAx_dy;
}

Eigen::Vector2d EMFieldComputer::computeLorentzForce(
    const EMFields& fields,
    double charge_density,
    const Eigen::Vector2d& current_density,
    int ix, int iy)
{
    // Lorentz force: f = ρE + J×B
    // In 2D with B = B_z ẑ:
    //   J×B = (J_x, J_y, 0) × (0, 0, B_z) = (J_y·B_z, -J_x·B_z, 0)

    double E_x = fields.E_x(ix, iy);
    double E_y = fields.E_y(ix, iy);
    double B_z = fields.B_z(ix, iy);

    double f_x = charge_density * E_x + current_density.y() * B_z;
    double f_y = charge_density * E_y - current_density.x() * B_z;

    return Eigen::Vector2d(f_x, f_y);
}

double EMFieldComputer::computeFieldEnergy(
    const EMFields& fields,
    double dx, double dy)
{
    // Energy density (Gaussian units): u = (E² + B²) / (8π)
    // Total energy: U = ∫ u dV = ∫ (E² + B²) / (8π) dx dy
    const double inv_8pi = 1.0 / (8.0 * M_PI);
    const int Nx = fields.E_x.rows();
    const int Ny = fields.E_x.cols();

    double total_energy = 0.0;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double E_x = fields.E_x(i, j);
            double E_y = fields.E_y(i, j);
            double B_z = fields.B_z(i, j);

            // Check for NaN/Inf in fields
            if (std::isnan(E_x) || std::isnan(E_y) || std::isnan(B_z) ||
                std::isinf(E_x) || std::isinf(E_y) || std::isinf(B_z)) {
                std::cerr << "[ERROR] Invalid EM field at (" << i << "," << j << "): "
                         << "E=(" << E_x << "," << E_y << "), B=" << B_z << std::endl;
                return std::numeric_limits<double>::quiet_NaN();
            }

            double E2 = E_x * E_x + E_y * E_y;
            double B2 = B_z * B_z;

            total_energy += (E2 + B2) * inv_8pi;
        }
    }

    // Multiply by cell volume
    total_energy *= dx * dy;

    // Final validation
    if (std::isnan(total_energy) || std::isinf(total_energy)) {
        std::cerr << "[ERROR] EM field energy is " << total_energy << std::endl;
    }

    return total_energy;
}

EMFieldComputer::PoyntingVector EMFieldComputer::computePoyntingVector(
    const EMFields& fields,
    double dx, double dy)
{
    // Poynting vector: S = (1/(4π)) E × B
    // In 2D with B = B_z ẑ:
    //   S_x = (1/(4π)) E_y · B_z
    //   S_y = (1/(4π)) (-E_x · B_z)
    //   |S| = (1/(4π)) sqrt(S_x² + S_y²)
    //
    // Interpretation: Energy flux (power per unit area) in the EM field

    const double inv_4pi = 1.0 / (4.0 * M_PI);
    const int Nx = fields.E_x.rows();
    const int Ny = fields.E_x.cols();

    PoyntingVector poynting(Nx, Ny);

    double max_flux_magnitude = 0.0;
    double total_flux_magnitude = 0.0;

    // Compute Poynting vector components and statistics
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double E_x = fields.E_x(i, j);
            double E_y = fields.E_y(i, j);
            double B_z = fields.B_z(i, j);

            // Poynting vector components
            poynting.S_x(i, j) = inv_4pi * E_y * B_z;
            poynting.S_y(i, j) = inv_4pi * (-E_x * B_z);

            // Flux magnitude at this point
            double flux_magnitude = std::sqrt(
                poynting.S_x(i, j) * poynting.S_x(i, j) +
                poynting.S_y(i, j) * poynting.S_y(i, j)
            );

            // Update max flux
            max_flux_magnitude = std::max(max_flux_magnitude, flux_magnitude);

            // Accumulate for total flux
            total_flux_magnitude += flux_magnitude;
        }
    }

    // Compute statistics
    poynting.max_flux = max_flux_magnitude;
    poynting.avg_flux = total_flux_magnitude / (Nx * Ny);

    // Total flux: integrate |S| over domain area
    // total_flux = ∫ |S| dA ≈ Σ |S| · dx·dy
    poynting.total_flux = total_flux_magnitude * dx * dy;

    return poynting;
}

void EMFieldComputer::computeDiagnostics(EMFields& fields)
{
    const int Nx = fields.E_x.rows();
    const int Ny = fields.E_x.cols();
    const int N_total = Nx * Ny;

    fields.max_E = 0.0;
    fields.max_B = 0.0;
    double sum_E = 0.0;
    double sum_B = 0.0;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // Electric field magnitude
            double E_mag = std::sqrt(
                fields.E_x(i, j) * fields.E_x(i, j) +
                fields.E_y(i, j) * fields.E_y(i, j)
            );

            // Magnetic field magnitude (just |B_z| in 2D)
            double B_mag = std::abs(fields.B_z(i, j));

            // Update max values
            fields.max_E = std::max(fields.max_E, E_mag);
            fields.max_B = std::max(fields.max_B, B_mag);

            // Accumulate for average
            sum_E += E_mag;
            sum_B += B_mag;
        }
    }

    fields.avg_E = sum_E / N_total;
    fields.avg_B = sum_B / N_total;
}

double EMFieldComputer::getPeriodicValue(
    const Eigen::MatrixXd& field,
    int i, int j,
    int offset_i, int offset_j)
{
    const int Nx = field.rows();
    const int Ny = field.cols();

    // Apply periodic boundary conditions
    int i_wrapped = (i + offset_i + Nx) % Nx;
    int j_wrapped = (j + offset_j + Ny) % Ny;

    return field(i_wrapped, j_wrapped);
}

Eigen::MatrixXd EMFieldComputer::computeSpatialDerivative(
    const Eigen::MatrixXd& field,
    double dx_or_dy,
    int direction,
    bool is_phase_field)
{
    const int Nx = field.rows();
    const int Ny = field.cols();

    Eigen::MatrixXd derivative(Nx, Ny);

    if (direction == 0) {
        // x-derivative: ∂_x f ≈ (f[i+1,j] - f[i-1,j]) / (2dx)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double f_plus = getPeriodicValue(field, i, j, 1, 0);
                double f_minus = getPeriodicValue(field, i, j, -1, 0);
                double diff = f_plus - f_minus;

                // CRITICAL FIX: For phase fields, wrap difference to [-π, π]
                // This handles branch cuts in vortex configurations
                if (is_phase_field) {
                    // Wrap the difference to [-π, π]
                    while (diff > M_PI) diff -= 2.0 * M_PI;
                    while (diff < -M_PI) diff += 2.0 * M_PI;
                }

                derivative(i, j) = diff / (2.0 * dx_or_dy);
            }
        }
    } else {
        // y-derivative: ∂_y f ≈ (f[i,j+1] - f[i,j-1]) / (2dy)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double f_plus = getPeriodicValue(field, i, j, 0, 1);
                double f_minus = getPeriodicValue(field, i, j, 0, -1);
                double diff = f_plus - f_minus;

                // CRITICAL FIX: For phase fields, wrap difference to [-π, π]
                if (is_phase_field) {
                    // Wrap the difference to [-π, π]
                    while (diff > M_PI) diff -= 2.0 * M_PI;
                    while (diff < -M_PI) diff += 2.0 * M_PI;
                }

                derivative(i, j) = diff / (2.0 * dx_or_dy);
            }
        }
    }

    return derivative;
}

Eigen::MatrixXd EMFieldComputer::computePhaseGradientConjugateProduct(
    const Eigen::MatrixXd& R_field,
    const Eigen::MatrixXd& theta_field,
    double dx_or_dy,
    int direction,
    RegularizationType reg_type)
{
    /**
     * Conjugate Product Method for Phase Gradients with Regularization
     *
     * Physics:
     *   Z = R·exp(iθ)  (complex Kuramoto field)
     *   ∇Z = (∇R + i·R·∇θ)·exp(iθ)
     *   Z*·∇Z = R·exp(-iθ)·(∇R + i·R·∇θ)·exp(iθ)
     *         = R·(∇R + i·R·∇θ)
     *   Im(Z*·∇Z) = R²·∇θ
     *
     * Regularization prescriptions:
     *   NONE:      A = Im(Z*·∇Z) / R²  = ∇θ       (unregularized)
     *   R_FACTOR:  A = Im(Z*·∇Z) / R   = R·∇θ     (R regularization)
     *   R2_FACTOR: A = Im(Z*·∇Z)       = R²·∇θ    (R² regularization)
     *
     * Numerical Advantage:
     *   At vortex cores where R → 0:
     *     - Numerator Im(Z*·∇Z) = R²·∇θ → 0
     *     - Denominator depends on prescription
     *     - All prescriptions remain finite via L'Hôpital
     *     - No spurious gradients from 2π branch cuts
     *
     * Implementation:
     *   1. Construct Z = R·cos(θ) + i·R·sin(θ) at center and neighbors
     *   2. Compute ∇Z via centered differences on real and imaginary parts
     *   3. Compute Im(Z*·∇Z) = Z_real·(∇Z)_imag - Z_imag·(∇Z)_real
     *   4. Normalize by prescription-dependent denominator
     */

    const int Nx = R_field.rows();
    const int Ny = R_field.cols();
    const double epsilon = 1e-10;  // Regularization to prevent division by zero

    Eigen::MatrixXd phase_gradient(Nx, Ny);

    if (direction == 0) {
        // x-derivative
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                // Get R and θ at center point
                double R_center = R_field(i, j);
                double theta_center = theta_field(i, j);

                // Construct complex field Z at center
                double Z_real = R_center * std::cos(theta_center);
                double Z_imag = R_center * std::sin(theta_center);

                // Get R and θ at neighbors (with periodic BC)
                double R_plus = getPeriodicValue(R_field, i, j, 1, 0);
                double theta_plus = getPeriodicValue(theta_field, i, j, 1, 0);
                double R_minus = getPeriodicValue(R_field, i, j, -1, 0);
                double theta_minus = getPeriodicValue(theta_field, i, j, -1, 0);

                // Construct complex field Z at neighbors
                double Z_plus_real = R_plus * std::cos(theta_plus);
                double Z_plus_imag = R_plus * std::sin(theta_plus);
                double Z_minus_real = R_minus * std::cos(theta_minus);
                double Z_minus_imag = R_minus * std::sin(theta_minus);

                // Compute gradient of complex field (centered difference)
                double dZ_dx_real = (Z_plus_real - Z_minus_real) / (2.0 * dx_or_dy);
                double dZ_dx_imag = (Z_plus_imag - Z_minus_imag) / (2.0 * dx_or_dy);

                // Compute Im(Z*·∇Z) = Z_real·(∇Z)_imag - Z_imag·(∇Z)_real
                // This gives: Im(Z*·∇Z) = R²·∇θ
                double Im_Z_conj_dZ = Z_real * dZ_dx_imag - Z_imag * dZ_dx_real;

                // Apply regularization prescription
                double R_squared = R_center * R_center;
                double denominator;

                switch (reg_type) {
                    case RegularizationType::NONE:
                        // A = Im(Z*·∇Z) / R² = ∇θ (unregularized)
                        denominator = std::max(R_squared, epsilon);
                        phase_gradient(i, j) = Im_Z_conj_dZ / denominator;
                        break;

                    case RegularizationType::R_FACTOR:
                        // A = Im(Z*·∇Z) / R = R·∇θ (R regularization)
                        denominator = std::max(R_center, epsilon);
                        phase_gradient(i, j) = Im_Z_conj_dZ / denominator;
                        break;

                    case RegularizationType::R2_FACTOR:
                        // A = Im(Z*·∇Z) = R²·∇θ (R² regularization - no normalization)
                        phase_gradient(i, j) = Im_Z_conj_dZ;
                        break;
                }
            }
        }
    } else {
        // y-derivative
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                // Get R and θ at center point
                double R_center = R_field(i, j);
                double theta_center = theta_field(i, j);

                // Construct complex field Z at center
                double Z_real = R_center * std::cos(theta_center);
                double Z_imag = R_center * std::sin(theta_center);

                // Get R and θ at neighbors (with periodic BC)
                double R_plus = getPeriodicValue(R_field, i, j, 0, 1);
                double theta_plus = getPeriodicValue(theta_field, i, j, 0, 1);
                double R_minus = getPeriodicValue(R_field, i, j, 0, -1);
                double theta_minus = getPeriodicValue(theta_field, i, j, 0, -1);

                // Construct complex field Z at neighbors
                double Z_plus_real = R_plus * std::cos(theta_plus);
                double Z_plus_imag = R_plus * std::sin(theta_plus);
                double Z_minus_real = R_minus * std::cos(theta_minus);
                double Z_minus_imag = R_minus * std::sin(theta_minus);

                // Compute gradient of complex field (centered difference)
                double dZ_dy_real = (Z_plus_real - Z_minus_real) / (2.0 * dx_or_dy);
                double dZ_dy_imag = (Z_plus_imag - Z_minus_imag) / (2.0 * dx_or_dy);

                // Compute Im(Z*·∇Z) = Z_real·(∇Z)_imag - Z_imag·(∇Z)_real
                // This gives: Im(Z*·∇Z) = R²·∇θ
                double Im_Z_conj_dZ = Z_real * dZ_dy_imag - Z_imag * dZ_dy_real;

                // Apply regularization prescription
                double R_squared = R_center * R_center;
                double denominator;

                switch (reg_type) {
                    case RegularizationType::NONE:
                        // A = Im(Z*·∇Z) / R² = ∇θ (unregularized)
                        denominator = std::max(R_squared, epsilon);
                        phase_gradient(i, j) = Im_Z_conj_dZ / denominator;
                        break;

                    case RegularizationType::R_FACTOR:
                        // A = Im(Z*·∇Z) / R = R·∇θ (R regularization)
                        denominator = std::max(R_center, epsilon);
                        phase_gradient(i, j) = Im_Z_conj_dZ / denominator;
                        break;

                    case RegularizationType::R2_FACTOR:
                        // A = Im(Z*·∇Z) = R²·∇θ (R² regularization - no normalization)
                        phase_gradient(i, j) = Im_Z_conj_dZ;
                        break;
                }
            }
        }
    }

    return phase_gradient;
}
