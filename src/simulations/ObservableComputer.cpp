#include "simulations/ObservableComputer.h"
#include "analysis/GeometryAnalyzer.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>

ObservableComputer::Observables ObservableComputer::compute(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double delta,
    double time,
    double E0,
    double norm_tolerance,
    double energy_tolerance) {

    Observables obs;
    obs.time = time;

    // Dirac observables
    obs.norm = computeNorm(dirac);
    obs.norm_error = obs.norm - 1.0;
    obs.energy_kinetic = computeKineticEnergy(dirac);
    obs.energy_potential = computePotentialEnergy(dirac, R_field, delta);
    obs.energy_total = obs.energy_kinetic + obs.energy_potential;

    obs.position_x = computePositionExpectation(dirac, 0);
    obs.position_y = computePositionExpectation(dirac, 1);
    obs.momentum_x = computeMomentumExpectation(dirac, 0);
    obs.momentum_y = computeMomentumExpectation(dirac, 1);

    // Sync field observables
    auto [R_avg, R_max, R_min, R_var] = computeSyncFieldStats(R_field);
    obs.R_avg = R_avg;
    obs.R_max = R_max;
    obs.R_min = R_min;
    obs.R_variance = R_var;

    // Validation
    obs.norm_valid = std::abs(obs.norm_error) < norm_tolerance;

    if (E0 != 0.0) {
        double energy_drift = std::abs(obs.energy_total - E0) / std::abs(E0);
        obs.energy_valid = energy_drift < energy_tolerance;
    } else {
        obs.energy_valid = true; // Can't validate without E0
    }

    return obs;
}

double ObservableComputer::computeNorm(const DiracEvolution& dirac) {
    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    double norm = 0.0;
    for (int i = 0; i < Nx * Ny; ++i) {
        // 4-component spinor: |Ψ|² = Σ_α |ψ_α|²
        for (int alpha = 0; alpha < 4; ++alpha) {
            std::complex<double> psi_val = psi[i * 4 + alpha];
            norm += std::norm(psi_val);  // |ψ|² = ψ*ψ
        }
    }

    // Integrate: ||Ψ||² = ∫|Ψ|² dA = Σ |Ψ|² · dx²
    norm *= dx * dx;

    return norm;
}

double ObservableComputer::computeEnergy(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double delta) {

    double T = computeKineticEnergy(dirac);
    double V = computePotentialEnergy(dirac, R_field, delta);

    return T + V;
}

double ObservableComputer::computeKineticEnergy(const DiracEvolution& dirac) {
    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    double T = 0.0;

    // T = ∫ Ψ†(-iα·∇)Ψ dA
    // Using finite differences for ∇
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Periodic boundaries
            int idx_xp = iy * Nx + ((ix + 1) % Nx);
            int idx_xm = iy * Nx + ((ix - 1 + Nx) % Nx);
            int idx_yp = ((iy + 1) % Ny) * Nx + ix;
            int idx_ym = ((iy - 1 + Ny) % Ny) * Nx + ix;

            for (int alpha = 0; alpha < 4; ++alpha) {
                std::complex<double> psi_center = psi[idx * 4 + alpha];

                // ∂Ψ/∂x ≈ (Ψ(x+dx) - Ψ(x-dx)) / (2dx)
                std::complex<double> psi_xp = psi[idx_xp * 4 + alpha];
                std::complex<double> psi_xm = psi[idx_xm * 4 + alpha];
                std::complex<double> dpsi_dx = (psi_xp - psi_xm) / (2.0 * dx);

                // ∂Ψ/∂y ≈ (Ψ(y+dy) - Ψ(y-dy)) / (2dy)
                std::complex<double> psi_yp = psi[idx_yp * 4 + alpha];
                std::complex<double> psi_ym = psi[idx_ym * 4 + alpha];
                std::complex<double> dpsi_dy = (psi_yp - psi_ym) / (2.0 * dx);

                // -iα·∇Ψ contribution (simplified, using |∇Ψ|²)
                double grad_term = std::norm(dpsi_dx) + std::norm(dpsi_dy);
                T += grad_term;
            }
        }
    }

    T *= dx * dx;  // Spatial integration
    return T;
}

double ObservableComputer::computePotentialEnergy(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double delta) {

    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    double V = 0.0;

    // V = ∫ Ψ†(β·m(x))Ψ dA where m(x) = Δ·R(x)
    for (int i = 0; i < Nx * Ny; ++i) {
        double mass = delta * R_field[i];

        // For 4-component spinor with β matrix (simplified: use |Ψ|² · m)
        for (int alpha = 0; alpha < 4; ++alpha) {
            V += std::norm(psi[i * 4 + alpha]) * mass;
        }
    }

    V *= dx * dx;  // Spatial integration
    return V;
}

std::complex<double> ObservableComputer::computePositionExpectation(
    const DiracEvolution& dirac,
    int component) {

    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    std::complex<double> expectation(0.0, 0.0);

    // <Ψ|x|Ψ> = ∫ Ψ†·x·Ψ dA
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            double position = (component == 0) ? (ix * dx) : (iy * dx);

            for (int alpha = 0; alpha < 4; ++alpha) {
                std::complex<double> psi_val = psi[idx * 4 + alpha];
                expectation += std::conj(psi_val) * position * psi_val;
            }
        }
    }

    expectation *= dx * dx;  // Spatial integration
    return expectation;
}

std::complex<double> ObservableComputer::computeMomentumExpectation(
    const DiracEvolution& dirac,
    int component) {

    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    std::complex<double> expectation(0.0, 0.0);
    std::complex<double> i_unit(0.0, 1.0);

    // <Ψ|p|Ψ> = <Ψ|-i∇|Ψ> = ∫ Ψ†·(-i∇Ψ) dA
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Periodic boundaries
            int idx_p, idx_m;
            if (component == 0) {  // x-direction
                idx_p = iy * Nx + ((ix + 1) % Nx);
                idx_m = iy * Nx + ((ix - 1 + Nx) % Nx);
            } else {  // y-direction
                idx_p = ((iy + 1) % Ny) * Nx + ix;
                idx_m = ((iy - 1 + Ny) % Ny) * Nx + ix;
            }

            for (int alpha = 0; alpha < 4; ++alpha) {
                std::complex<double> psi_center = psi[idx * 4 + alpha];
                std::complex<double> psi_p = psi[idx_p * 4 + alpha];
                std::complex<double> psi_m = psi[idx_m * 4 + alpha];

                // ∇Ψ ≈ (Ψ(+) - Ψ(-)) / (2dx)
                std::complex<double> dpsi = (psi_p - psi_m) / (2.0 * dx);

                // -i∇Ψ
                expectation += std::conj(psi_center) * (-i_unit) * dpsi;
            }
        }
    }

    expectation *= dx * dx;  // Spatial integration
    return expectation;
}

std::tuple<double, double, double, double> ObservableComputer::computeSyncFieldStats(
    const std::vector<double>& R_field) {

    if (R_field.empty()) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    double R_avg = std::accumulate(R_field.begin(), R_field.end(), 0.0) / R_field.size();
    double R_max = *std::max_element(R_field.begin(), R_field.end());
    double R_min = *std::min_element(R_field.begin(), R_field.end());

    // Variance: Var(R) = E[R²] - E[R]²
    double R_sq_sum = 0.0;
    for (double r : R_field) {
        R_sq_sum += r * r;
    }
    double R_variance = (R_sq_sum / R_field.size()) - (R_avg * R_avg);

    return {R_avg, R_max, R_min, R_variance};
}

std::pair<double, double> ObservableComputer::computeRFieldGradient(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double pos_x, double pos_y,
    double dx) {

    // Convert position to grid indices
    int ix = static_cast<int>(std::round(pos_x / dx));
    int iy = static_cast<int>(std::round(pos_y / dx));

    // Periodic boundaries
    ix = (ix % Nx + Nx) % Nx;
    iy = (iy % Ny + Ny) % Ny;

    int idx = iy * Nx + ix;

    // Centered finite differences
    int idx_xp = iy * Nx + ((ix + 1) % Nx);
    int idx_xm = iy * Nx + ((ix - 1 + Nx) % Nx);
    int idx_yp = ((iy + 1) % Ny) * Nx + ix;
    int idx_ym = ((iy - 1 + Ny) % Ny) * Nx + ix;

    double dR_dx = (R_field[idx_xp] - R_field[idx_xm]) / (2.0 * dx);
    double dR_dy = (R_field[idx_yp] - R_field[idx_ym]) / (2.0 * dx);

    return {dR_dx, dR_dy};
}

double ObservableComputer::computeCorrelation(
    const std::vector<double>& series1,
    const std::vector<double>& series2) {

    if (series1.size() != series2.size() || series1.empty()) {
        return 0.0;
    }

    size_t n = series1.size();

    // Compute means
    double mean1 = std::accumulate(series1.begin(), series1.end(), 0.0) / n;
    double mean2 = std::accumulate(series2.begin(), series2.end(), 0.0) / n;

    // Compute covariance and variances
    double cov = 0.0;
    double var1 = 0.0;
    double var2 = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double d1 = series1[i] - mean1;
        double d2 = series2[i] - mean2;
        cov += d1 * d2;
        var1 += d1 * d1;
        var2 += d2 * d2;
    }

    // Pearson correlation coefficient
    double denom = std::sqrt(var1 * var2);
    if (denom < 1e-10) {
        return 0.0;  // Avoid division by zero
    }

    return cov / denom;
}

std::string ObservableComputer::toCSVLine(const Observables& obs) {
    std::ostringstream oss;
    oss << std::setprecision(10);

    oss << obs.time << ","
        << obs.norm << ","
        << obs.norm_error << ","
        << obs.energy_total << ","
        << obs.energy_kinetic << ","
        << obs.energy_potential << ","
        << obs.position_x.real() << ","
        << obs.position_x.imag() << ","
        << obs.position_y.real() << ","
        << obs.position_y.imag() << ","
        << obs.momentum_x.real() << ","
        << obs.momentum_x.imag() << ","
        << obs.momentum_y.real() << ","
        << obs.momentum_y.imag() << ","
        << obs.R_avg << ","
        << obs.R_max << ","
        << obs.R_min << ","
        << obs.R_variance << ","
        << (obs.norm_valid ? 1 : 0) << ","
        << (obs.energy_valid ? 1 : 0);

    return oss.str();
}

std::string ObservableComputer::getCSVHeader() {
    return "time,norm,norm_error,E_total,E_kin,E_pot,"
           "pos_x_re,pos_x_im,pos_y_re,pos_y_im,"
           "mom_x_re,mom_x_im,mom_y_re,mom_y_im,"
           "R_avg,R_max,R_min,R_var,"
           "norm_valid,energy_valid";
}

// ========== Klein-Gordon Observable Implementations ==========

/**
 * IMPORTANT: Klein-Gordon Conserved Quantities
 *
 * The Klein-Gordon equation conserves DIFFERENT quantities than Schrödinger/Dirac:
 *
 * Conserved:
 *   - Klein-Gordon current: j^0 = i(φ*·∂_tφ - φ·∂_tφ*)
 *   - Phase-space norm: ∫[|φ|² + |∂_tφ/m|²]dx
 *   - Energy: E = ∫[|∂_tφ|² + |∇φ|² + m²|φ|²]dx
 *
 * NOT conserved:
 *   - Position-space norm: ||φ||² = ∫|φ|²dx
 *
 * For relativistic wavepackets, ||φ||² can grow significantly (e.g., 1.0 → 1.5)
 * as energy trades between kinetic (∂_tφ) and position (φ) components.
 * This is CORRECT physics, not a numerical error.
 */
ObservableComputer::Observables ObservableComputer::computeKG(
    const KleinGordonEvolution& kg,
    const std::vector<double>& R_field,
    double delta,
    double time,
    double E0,
    double norm_tolerance,
    double energy_tolerance) {

    Observables obs;
    obs.time = time;

    // Klein-Gordon observables
    obs.norm = computeNormKG(kg);
    obs.norm_error = obs.norm - 1.0;
    obs.energy_kinetic = computeKineticEnergyKG(kg);
    obs.energy_potential = computePotentialEnergyKG(kg, R_field, delta);
    obs.energy_total = obs.energy_kinetic + obs.energy_potential;

    obs.position_x = computePositionExpectationKG(kg, 0);
    obs.position_y = computePositionExpectationKG(kg, 1);
    obs.momentum_x = computeMomentumExpectationKG(kg, 0);
    obs.momentum_y = computeMomentumExpectationKG(kg, 1);

    // Sync field observables
    auto [R_avg, R_max, R_min, R_var] = computeSyncFieldStats(R_field);
    obs.R_avg = R_avg;
    obs.R_max = R_max;
    obs.R_min = R_min;
    obs.R_variance = R_var;

    // Validation
    obs.norm_valid = std::abs(obs.norm_error) < norm_tolerance;

    if (E0 != 0.0) {
        double energy_drift = std::abs(obs.energy_total - E0) / std::abs(E0);
        obs.energy_valid = energy_drift < energy_tolerance;
    } else {
        obs.energy_valid = true; // Can't validate without E0
    }

    return obs;
}

double ObservableComputer::computeNormKG(const KleinGordonEvolution& kg) {
    const auto& phi = kg.getField();
    int Nx = kg.getNx();
    int Ny = kg.getNy();
    double dx = kg.getDx();

    double norm = 0.0;
    for (int i = 0; i < Nx * Ny; ++i) {
        // Klein-Gordon: |φ|² = φ*φ
        norm += std::norm(phi[i]);
    }

    // Integrate: ||φ||² = ∫|φ|² dA = Σ |φ|² · dx²
    norm *= dx * dx;

    return norm;
}

double ObservableComputer::computeEnergyKG(
    const KleinGordonEvolution& kg,
    const std::vector<double>& R_field,
    double delta) {

    double T = computeKineticEnergyKG(kg);
    double V = computePotentialEnergyKG(kg, R_field, delta);

    return T + V;
}

double ObservableComputer::computeKineticEnergyKG(const KleinGordonEvolution& kg) {
    const auto& phi_dot = kg.getFieldDot();
    int Nx = kg.getNx();
    int Ny = kg.getNy();
    double dx = kg.getDx();

    double T = 0.0;

    // T = ∫|∂_tφ|² dx
    for (int i = 0; i < Nx * Ny; ++i) {
        T += std::norm(phi_dot[i]);
    }

    T *= dx * dx;  // Spatial integration

    return T;
}

double ObservableComputer::computePotentialEnergyKG(
    const KleinGordonEvolution& kg,
    const std::vector<double>& R_field,
    double delta) {

    const auto& phi = kg.getField();
    int Nx = kg.getNx();
    int Ny = kg.getNy();
    double dx = kg.getDx();

    double V = 0.0;

    // V = ∫[|∇φ|² + m²|φ|²] dx
    // where m(x,y) = Δ·R(x,y)
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Periodic boundaries
            int idx_xp = iy * Nx + ((ix + 1) % Nx);
            int idx_xm = iy * Nx + ((ix - 1 + Nx) % Nx);
            int idx_yp = ((iy + 1) % Ny) * Nx + ix;
            int idx_ym = ((iy - 1 + Ny) % Ny) * Nx + ix;

            std::complex<float> phi_center = phi[idx];

            // Gradient energy: |∇φ|²
            std::complex<float> phi_xp = phi[idx_xp];
            std::complex<float> phi_xm = phi[idx_xm];
            std::complex<float> dphi_dx = (phi_xp - phi_xm) / (2.0f * static_cast<float>(dx));

            std::complex<float> phi_yp = phi[idx_yp];
            std::complex<float> phi_ym = phi[idx_ym];
            std::complex<float> dphi_dy = (phi_yp - phi_ym) / (2.0f * static_cast<float>(dx));

            double grad_energy = std::norm(dphi_dx) + std::norm(dphi_dy);

            // Mass energy: m²|φ|²
            double mass = delta * R_field[idx];
            double mass_energy = mass * mass * std::norm(phi_center);

            V += grad_energy + mass_energy;
        }
    }

    V *= dx * dx;  // Spatial integration
    return V;
}

std::complex<double> ObservableComputer::computePositionExpectationKG(
    const KleinGordonEvolution& kg,
    int component) {

    const auto& phi = kg.getField();
    int Nx = kg.getNx();
    int Ny = kg.getNy();
    double dx = kg.getDx();

    std::complex<double> expectation(0.0, 0.0);
    double norm = 0.0;

    // <x> = ∫ x|φ|² dx / ∫|φ|² dx
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            double position = (component == 0) ? (ix * dx) : (iy * dx);
            double density = std::norm(phi[idx]);

            expectation += position * density;
            norm += density;
        }
    }

    expectation *= dx * dx;  // Spatial integration
    norm *= dx * dx;

    // Normalize by total norm
    if (norm > 1e-10) {
        expectation /= norm;
    }

    return expectation;
}

std::complex<double> ObservableComputer::computeMomentumExpectationKG(
    const KleinGordonEvolution& kg,
    int component) {

    const auto& phi = kg.getField();
    int Nx = kg.getNx();
    int Ny = kg.getNy();
    double dx = kg.getDx();

    std::complex<double> expectation(0.0, 0.0);
    std::complex<double> i_unit(0.0, 1.0);

    // <p> = -i ∫ φ*(∇φ) dx
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Periodic boundaries
            int idx_p, idx_m;
            if (component == 0) {  // x-direction
                idx_p = iy * Nx + ((ix + 1) % Nx);
                idx_m = iy * Nx + ((ix - 1 + Nx) % Nx);
            } else {  // y-direction
                idx_p = ((iy + 1) % Ny) * Nx + ix;
                idx_m = ((iy - 1 + Ny) % Ny) * Nx + ix;
            }

            std::complex<float> phi_center = phi[idx];
            std::complex<float> phi_p = phi[idx_p];
            std::complex<float> phi_m = phi[idx_m];

            // ∇φ ≈ (φ(+) - φ(-)) / (2dx)
            std::complex<double> dphi =
                (std::complex<double>(phi_p) - std::complex<double>(phi_m)) / (2.0 * dx);

            // -i·φ*·∇φ
            expectation += std::conj(std::complex<double>(phi_center)) * (-i_unit) * dphi;
        }
    }

    expectation *= dx * dx;  // Spatial integration
    return expectation;
}

// ========== Phase 3 Vacuum Structure Observable Implementations ==========

std::tuple<double, double, double, double> ObservableComputer::computeEnergyDensityByRegion(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double delta,
    float x_boundary) {

    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    double E_left = 0.0;
    double E_right = 0.0;
    double R_sum_left = 0.0;
    double R_sum_right = 0.0;
    int count_left = 0;
    int count_right = 0;

    // Compute energy and R_avg in each region
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;
            float x = static_cast<float>(ix);

            bool is_left = (x < x_boundary);

            // Periodic boundaries for gradients
            int idx_xp = iy * Nx + ((ix + 1) % Nx);
            int idx_xm = iy * Nx + ((ix - 1 + Nx) % Nx);
            int idx_yp = ((iy + 1) % Ny) * Nx + ix;
            int idx_ym = ((iy - 1 + Ny) % Ny) * Nx + ix;

            // Kinetic energy contribution (gradient terms)
            double T_local = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                std::complex<double> psi_center = psi[idx * 4 + alpha];

                std::complex<double> psi_xp = psi[idx_xp * 4 + alpha];
                std::complex<double> psi_xm = psi[idx_xm * 4 + alpha];
                std::complex<double> dpsi_dx = (psi_xp - psi_xm) / (2.0 * dx);

                std::complex<double> psi_yp = psi[idx_yp * 4 + alpha];
                std::complex<double> psi_ym = psi[idx_ym * 4 + alpha];
                std::complex<double> dpsi_dy = (psi_yp - psi_ym) / (2.0 * dx);

                T_local += std::norm(dpsi_dx) + std::norm(dpsi_dy);
            }
            T_local *= dx * dx;

            // Potential energy contribution
            double V_local = 0.0;
            double mass = delta * R_field[idx];
            for (int alpha = 0; alpha < 4; ++alpha) {
                V_local += std::norm(psi[idx * 4 + alpha]) * mass;
            }
            V_local *= dx * dx;

            double E_local = T_local + V_local;

            // Accumulate by region
            if (is_left) {
                E_left += E_local;
                R_sum_left += R_field[idx];
                count_left++;
            } else {
                E_right += E_local;
                R_sum_right += R_field[idx];
                count_right++;
            }
        }
    }

    // Compute average R in each region
    double R_avg_left = (count_left > 0) ? R_sum_left / count_left : 0.0;
    double R_avg_right = (count_right > 0) ? R_sum_right / count_right : 0.0;

    // Compute energy density: ρ = E / V
    // Volume = count * dx^2 (2D area)
    double V_left = count_left * dx * dx;
    double V_right = count_right * dx * dx;

    double rho_left = (V_left > 1e-10) ? E_left / V_left : 0.0;
    double rho_right = (V_right > 1e-10) ? E_right / V_right : 0.0;

    return {rho_left, rho_right, R_avg_left, R_avg_right};
}

double ObservableComputer::computeDefectForce(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double delta,
    float x1) {

    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    double force = 0.0;

    // Find grid index closest to x1
    int ix1 = static_cast<int>(std::round(x1));
    if (ix1 < 0 || ix1 >= Nx) {
        return 0.0;  // Out of bounds
    }

    // Integrate force along defect line (vary y, fix x = ix1)
    // F = ∫ β⟨Ψ|β|Ψ⟩ · Δ∇R_x dy
    for (int iy = 0; iy < Ny; ++iy) {
        int idx = iy * Nx + ix1;

        // Compute ∇R_x at this point (centered difference)
        int idx_xp = iy * Nx + ((ix1 + 1) % Nx);
        int idx_xm = iy * Nx + ((ix1 - 1 + Nx) % Nx);

        double dR_dx = (R_field[idx_xp] - R_field[idx_xm]) / (2.0 * dx);

        // Compute β expectation value ⟨Ψ|β|Ψ⟩
        // For 4-component Dirac spinor, β matrix has structure:
        // β = diag(1, 1, -1, -1) in standard representation
        // ⟨Ψ|β|Ψ⟩ = |ψ₁|² + |ψ₂|² - |ψ₃|² - |ψ₄|²
        double beta_expectation = 0.0;
        for (int alpha = 0; alpha < 4; ++alpha) {
            double sign = (alpha < 2) ? 1.0 : -1.0;
            beta_expectation += sign * std::norm(psi[idx * 4 + alpha]);
        }

        // Force contribution: F_y = β_expectation · Δ · dR/dx · dy
        force += beta_expectation * delta * dR_dx * dx;
    }

    return force;
}

// ==============================================================================
// Two-Particle Observable Methods (Test 3.4: Antiparticle Separation)
// ==============================================================================

ObservableComputer::TwoParticleObservables ObservableComputer::computeTwoParticle(
    const DiracEvolution& particle,
    const DiracEvolution& antiparticle,
    const std::vector<double>& R_field,
    double delta,
    double time) {

    TwoParticleObservables obs;
    obs.time = time;

    // Compute particle observables
    obs.norm_particle = particle.getNorm();

    float x_p, y_p;
    particle.getCenterOfMass(x_p, y_p);
    obs.position_x_particle = x_p;
    obs.position_y_particle = y_p;

    // Compute antiparticle observables
    obs.norm_antiparticle = antiparticle.getNorm();

    float x_a, y_a;
    antiparticle.getCenterOfMass(x_a, y_a);
    obs.position_x_antiparticle = x_a;
    obs.position_y_antiparticle = y_a;

    // Compute separation distance
    double dx = obs.position_x_particle - obs.position_x_antiparticle;
    double dy = obs.position_y_particle - obs.position_y_antiparticle;
    obs.separation_distance = std::sqrt(dx * dx + dy * dy);

    // Compute forces on each particle: F = -β_sign · Δ · ∇R
    // For particle (β_sign = +1): F_p = -Δ · ∇R (attracted to high-R regions)
    // For antiparticle (β_sign = -1): F_a = +Δ · ∇R (repelled from high-R regions)

    uint32_t Nx = particle.getNx();
    uint32_t Ny = particle.getNy();

    // Particle force (at center of mass)
    int ix_p = static_cast<int>(std::round(x_p));
    int iy_p = static_cast<int>(std::round(y_p));
    if (ix_p >= 0 && ix_p < static_cast<int>(Nx) &&
        iy_p >= 0 && iy_p < static_cast<int>(Ny)) {

        int idx_p = iy_p * Nx + ix_p;
        int idx_xp = iy_p * Nx + ((ix_p + 1) % Nx);
        int idx_xm = iy_p * Nx + ((ix_p - 1 + Nx) % Nx);
        int idx_yp = ((iy_p + 1) % Ny) * Nx + ix_p;
        int idx_ym = ((iy_p - 1 + Ny) % Ny) * Nx + ix_p;

        double dR_dx_p = (R_field[idx_xp] - R_field[idx_xm]) / 2.0;
        double dR_dy_p = (R_field[idx_yp] - R_field[idx_ym]) / 2.0;

        // Particle force: F_p = -β_sign · Δ · ∇R = -Δ · ∇R
        obs.force_particle_x = -delta * dR_dx_p;
        obs.force_particle_y = -delta * dR_dy_p;
    } else {
        obs.force_particle_x = 0.0;
        obs.force_particle_y = 0.0;
    }

    // Antiparticle force (at center of mass)
    int ix_a = static_cast<int>(std::round(x_a));
    int iy_a = static_cast<int>(std::round(y_a));
    if (ix_a >= 0 && ix_a < static_cast<int>(Nx) &&
        iy_a >= 0 && iy_a < static_cast<int>(Ny)) {

        int idx_a = iy_a * Nx + ix_a;
        int idx_xp = iy_a * Nx + ((ix_a + 1) % Nx);
        int idx_xm = iy_a * Nx + ((ix_a - 1 + Nx) % Nx);
        int idx_yp = ((iy_a + 1) % Ny) * Nx + ix_a;
        int idx_ym = ((iy_a - 1 + Ny) % Ny) * Nx + ix_a;

        double dR_dx_a = (R_field[idx_xp] - R_field[idx_xm]) / 2.0;
        double dR_dy_a = (R_field[idx_yp] - R_field[idx_ym]) / 2.0;

        // Antiparticle force: F_a = -β_sign · Δ · ∇R = +Δ · ∇R (opposite sign)
        obs.force_antiparticle_x = delta * dR_dx_a;
        obs.force_antiparticle_y = delta * dR_dy_a;
    } else {
        obs.force_antiparticle_x = 0.0;
        obs.force_antiparticle_y = 0.0;
    }

    // Compute force dot product F_p · F_a
    // Should be < 0 for opposite forces (anti-correlation)
    obs.force_dot_product = obs.force_particle_x * obs.force_antiparticle_x +
                             obs.force_particle_y * obs.force_antiparticle_y;

    // Note: Momentum expectation values not included for brevity
    // Can be added if needed using computeMomentumExpectation()

    return obs;
}

// ==============================================================================
// Phase 4 Time Gradient Observable Implementations
// ==============================================================================

double ObservableComputer::computePhaseAccumulation(const DiracEvolution& dirac) {
    /**
     * Compute global phase φ = arg(⟨Ψ|Ψ⟩) = arg(∫Ψ*Ψ dx)
     *
     * Physics: Tracks overall phase evolution of wavepacket
     * In time-dilation mode (Test 4.1), φ(t) evolves slower where R < 1
     *
     * For 4-component spinor:
     * ⟨Ψ|Ψ⟩ = Σ_x Σ_α ψ*_α(x) ψ_α(x)
     *
     * Returns phase in radians [-π, π]
     */
    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    std::complex<double> global_amplitude(0.0, 0.0);

    // Sum over all grid points and spinor components
    for (int i = 0; i < Nx * Ny; ++i) {
        for (int alpha = 0; alpha < 4; ++alpha) {
            std::complex<double> psi_val = psi[i * 4 + alpha];

            // Accumulate ψ*_α · ψ_α (complex number with phase information)
            global_amplitude += std::conj(psi_val) * psi_val;
        }
    }

    // Extract phase: φ = arg(⟨Ψ|Ψ⟩)
    double phase = std::atan2(global_amplitude.imag(), global_amplitude.real());

    return phase;
}

std::tuple<double, double> ObservableComputer::computeTemporalForce(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    std::complex<double> momentum_prev_x,
    std::complex<double> momentum_prev_y,
    double dt,
    double delta) {
    /**
     * Compute temporal force F_temp = dp/dt - F_spatial
     *
     * Physics decomposition:
     * - F_total = dp/dt (total force from momentum change)
     * - F_spatial = -Δ·∇R (known spatial force from gradient)
     * - F_temporal = F_total - F_spatial (residual, temporal contribution)
     *
     * This tests if time gradients contribute to force beyond spatial effects
     */

    // Compute current momentum
    std::complex<double> momentum_curr_x = computeMomentumExpectation(dirac, 0);
    std::complex<double> momentum_curr_y = computeMomentumExpectation(dirac, 1);

    // Momentum change: dp = p_curr - p_prev
    double dp_x = momentum_curr_x.real() - momentum_prev_x.real();
    double dp_y = momentum_curr_y.real() - momentum_prev_y.real();

    // Total force: F_total = dp/dt
    double F_total_x = dp_x / dt;
    double F_total_y = dp_y / dt;

    // Get particle position
    float pos_x, pos_y;
    dirac.getCenterOfMass(pos_x, pos_y);

    // Compute R-field gradient at particle center
    auto [grad_R_x, grad_R_y] = computeRFieldGradient(
        R_field,
        dirac.getNx(), dirac.getNy(),
        pos_x, pos_y,
        dirac.getDx()
    );

    // Spatial force: F_spatial = -Δ·∇R
    double F_spatial_x = -delta * grad_R_x;
    double F_spatial_y = -delta * grad_R_y;

    // Temporal force: F_temp = F_total - F_spatial (residual)
    double F_temp_x = F_total_x - F_spatial_x;
    double F_temp_y = F_total_y - F_spatial_y;

    return {F_temp_x, F_temp_y};
}

std::tuple<double, double> ObservableComputer::computeGeodesicAcceleration(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double L_domain) {
    /**
     * Compute geodesic acceleration a_geo = -Γ^i_μν v^μ v^ν
     *
     * For metric ds² = -R²dt² + dx² + dy²:
     * Dominant term: a^i ≈ -Γ^i_00 = R ∂R/∂x^i
     *
     * Uses GeometryAnalyzer to compute Christoffel symbols and
     * geodesic acceleration from R-field geometry.
     */

    // Get particle position
    float pos_x_f, pos_y_f;
    dirac.getCenterOfMass(pos_x_f, pos_y_f);
    double pos_x = static_cast<double>(pos_x_f);
    double pos_y = static_cast<double>(pos_y_f);

    // Get particle velocity from momentum
    std::complex<double> px = computeMomentumExpectation(dirac, 0);
    std::complex<double> py = computeMomentumExpectation(dirac, 1);
    double vx = px.real();  // In Planck units, p ≈ v for m ~ 1
    double vy = py.real();

    // Compute geodesic acceleration using GeometryAnalyzer
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    auto accel = GeometryAnalyzer::computeGeodesicAcceleration(
        R_field, Nx, Ny, pos_x, pos_y, vx, vy, L_domain);

    return {accel(0), accel(1)};
}

std::tuple<double, double> ObservableComputer::computeChristoffelAtParticle(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    double L_domain) {
    /**
     * Compute Christoffel symbols Γ^x_00, Γ^y_00 at particle position
     *
     * These are the dominant components for geodesic acceleration:
     * a^x ≈ -Γ^x_00, a^y ≈ -Γ^y_00
     */

    // Get particle position
    float pos_x_f, pos_y_f;
    dirac.getCenterOfMass(pos_x_f, pos_y_f);
    double pos_x = static_cast<double>(pos_x_f);
    double pos_y = static_cast<double>(pos_y_f);

    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    auto gamma = GeometryAnalyzer::computeChristoffelInterpolated(
        R_field, Nx, Ny, pos_x, pos_y, L_domain);

    return {gamma.Gamma_x_tt, gamma.Gamma_y_tt};
}

double ObservableComputer::computeRiemannScalarAtParticle(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field) {
    /**
     * Compute Riemann curvature scalar R at particle position
     *
     * Measures spacetime curvature induced by R-field
     */

    // Get particle position in grid units
    float pos_x_f, pos_y_f;
    dirac.getCenterOfMass(pos_x_f, pos_y_f);

    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    // Convert to grid indices (nearest grid point)
    int ix = static_cast<int>(std::round(pos_x_f));
    int iy = static_cast<int>(std::round(pos_y_f));

    // Clamp to interior for centered stencil
    ix = std::clamp(ix, 1, Nx - 2);
    iy = std::clamp(iy, 1, Ny - 2);

    return GeometryAnalyzer::computeRiemannScalar(
        R_field, Nx, Ny, ix, iy, dx, dx);
}

double ObservableComputer::computeGaussianCurvatureAtParticle(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field) {
    /**
     * Compute Gaussian curvature K at particle position
     *
     * Measures spatial curvature induced by R-field
     */

    // Get particle position in grid units
    float pos_x_f, pos_y_f;
    dirac.getCenterOfMass(pos_x_f, pos_y_f);

    int Nx = dirac.getNx();
    int Ny = dirac.getNy();
    double dx = dirac.getDx();

    // Convert to grid indices (nearest grid point)
    int ix = static_cast<int>(std::round(pos_x_f));
    int iy = static_cast<int>(std::round(pos_y_f));

    // Clamp to interior for centered stencil
    ix = std::clamp(ix, 1, Nx - 2);
    iy = std::clamp(iy, 1, Ny - 2);

    return GeometryAnalyzer::computeGaussianCurvature(
        R_field, Nx, Ny, ix, iy, dx, dx);
}
