#include "simulations/ObservableComputer.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <iomanip>

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
