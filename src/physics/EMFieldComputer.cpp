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
    double dx, double dy, double dt)
{
    const int Nx = theta_current.rows();
    const int Ny = theta_current.cols();

    EMFields fields(Nx, Ny);

    // Compute temporal derivative: φ = ∂_t θ
    // Use backward difference: ∂_t θ ≈ (θ_current - θ_previous) / dt
    fields.phi = (theta_current - theta_previous) / dt;

    // Compute spatial derivatives: A = ∇θ
    // A_x = ∂_x θ, A_y = ∂_y θ
    fields.A_x = computeSpatialDerivative(theta_current, dx, 0);
    fields.A_y = computeSpatialDerivative(theta_current, dy, 1);

    // Compute field strengths E, B
    computeFieldStrengths(fields, dx, dy);

    // Compute diagnostics
    computeDiagnostics(fields);

    // Compute total field energy
    fields.total_field_energy = computeFieldEnergy(fields, dx, dy);

    return fields;
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
            double E2 = fields.E_x(i, j) * fields.E_x(i, j) +
                       fields.E_y(i, j) * fields.E_y(i, j);
            double B2 = fields.B_z(i, j) * fields.B_z(i, j);

            total_energy += (E2 + B2) * inv_8pi;
        }
    }

    // Multiply by cell volume
    total_energy *= dx * dy;

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
    int direction)
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
                derivative(i, j) = (f_plus - f_minus) / (2.0 * dx_or_dy);
            }
        }
    } else {
        // y-derivative: ∂_y f ≈ (f[i,j+1] - f[i,j-1]) / (2dy)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double f_plus = getPeriodicValue(field, i, j, 0, 1);
                double f_minus = getPeriodicValue(field, i, j, 0, -1);
                derivative(i, j) = (f_plus - f_minus) / (2.0 * dx_or_dy);
            }
        }
    }

    return derivative;
}
