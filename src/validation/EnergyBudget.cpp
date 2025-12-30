#include "EnergyBudget.h"
#include <cmath>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <algorithm>

EnergyBudget::EnergyBudget(SMFTEngine* engine, const std::string& output_dir)
    : _engine(engine)
    , _output_dir(output_dir)
    , _E0(0.0)
    , _E0_set(false)
{
    // Default grid parameters - will be set when engine is initialized
    _Nx = 128;  // Default size
    _Ny = 128;

    // Default physics parameters - will be set properly when called with TestConfig
    _delta = 1.0f;
    _coupling = 0.1f;
    _K = 1.0f;
    _damping = 0.0f;

    // Default grid spacing - will be set properly from TestConfig
    float L_domain = 100.0f;  // Default domain size in Planck lengths
    _dx = L_domain / _Nx;
    _dy = L_domain / _Ny;
}

void EnergyBudget::setPhysicsParameters(float delta, float coupling, float K, float damping, float L_domain) {
    _delta = delta;
    _coupling = coupling;
    _K = K;
    _damping = damping;

    // Update grid spacing based on domain size
    _dx = L_domain / _Nx;
    _dy = L_domain / _Ny;
}

EnergyBudget::EnergyComponents EnergyBudget::computeComponents(
    DiracEvolution* dirac_evolution,
    double time,
    bool compute_em)
{
    EnergyComponents components = {0};

    // Get current field states from engine
    std::vector<float> theta_field = _engine->getPhaseField();
    std::vector<float> R_field = _engine->getSyncField();

    // 1. Dirac field energy (if Dirac evolution is active)
    if (dirac_evolution != nullptr) {
        // Use ObservableComputer for consistency
        components.T_dirac_kinetic = ObservableComputer::computeKineticEnergy(*dirac_evolution);

        // Convert float R_field to double for ObservableComputer
        std::vector<double> R_field_double(R_field.begin(), R_field.end());
        components.V_dirac_mass = ObservableComputer::computePotentialEnergy(
            *dirac_evolution, R_field_double, _delta);

        // Get Dirac spinor for coupling energy
        auto psi = dirac_evolution->getSpinor();
        components.E_coupling = computeCouplingEnergy(psi, R_field, _coupling, _dx, _dy);
    }

    // 2. Kuramoto field energy
    auto [E_gradient, E_sync] = computeKuramotoEnergy(theta_field, R_field, _K, _dx, _dy);
    components.E_kuramoto_gradient = E_gradient;
    components.E_kuramoto_sync = E_sync;

    // 3. EM field energy (if requested and enabled)
    if (compute_em && dirac_evolution != nullptr) {
        // For now, use the EM field energy from ObservableComputer if available
        // This would need theta_previous for full EM computation
        components.E_em_field = 0.0;  // Placeholder - would need theta history
    }

    // 4. Total energy
    components.E_total = components.T_dirac_kinetic + components.V_dirac_mass +
                         components.E_kuramoto_gradient + components.E_kuramoto_sync +
                         components.E_coupling + components.E_em_field;

    // 5. Dissipation tracking
    components.P_dissipated = computeDissipationRate(theta_field, _damping, _dx, _dy);

    // Track cumulative dissipation
    if (!_time_series.time.empty()) {
        double dt = time - _time_series.time.back();
        double E_prev = _time_series.components.back().E_dissipated_cumulative;
        components.E_dissipated_cumulative = E_prev + components.P_dissipated * dt;

        // Compute energy derivative
        components.dE_dt = (components.E_total - _time_series.components.back().E_total) / dt;

        // Conservation error: |dE/dt + P_dissipated| should be ~0
        components.conservation_error = std::abs(components.dE_dt + components.P_dissipated);
    } else {
        components.E_dissipated_cumulative = 0.0;
        components.dE_dt = 0.0;
        components.conservation_error = 0.0;
    }

    return components;
}

void EnergyBudget::trackTimeSeries(double t, const EnergyComponents& components) {
    _time_series.time.push_back(t);
    _time_series.components.push_back(components);

    // Set initial energy on first call
    if (!_E0_set) {
        _E0 = components.E_total;
        _E0_set = true;
    }

    // Compute derived metrics
    double rel_error = std::abs(components.E_total - _E0) / _E0;
    _time_series.relative_error.push_back(rel_error);

    // Energy drift: E(t) - E₀ + ∫P dt
    // This should be ~0 for perfect conservation with damping
    double drift = components.E_total - _E0 + components.E_dissipated_cumulative;
    _time_series.energy_drift.push_back(drift);
}

bool EnergyBudget::analyzeConservation(const std::string& output_file) {
    std::ofstream report(output_file);
    if (!report.is_open()) {
        std::cerr << "Error: Could not open output file " << output_file << std::endl;
        return false;
    }

    report << "=== ENERGY BUDGET ANALYSIS REPORT ===" << std::endl;
    report << "Simulation parameters:" << std::endl;
    report << "  Grid: " << _Nx << "×" << _Ny << std::endl;
    report << "  Grid spacing: dx=" << _dx << ", dy=" << _dy << std::endl;
    report << "  Delta (mass gap): " << _delta << std::endl;
    report << "  Coupling: " << _coupling << std::endl;
    report << "  Kuramoto K: " << _K << std::endl;
    report << "  Damping γ: " << _damping << std::endl;
    report << std::endl;

    if (_time_series.time.empty()) {
        report << "ERROR: No time series data available" << std::endl;
        return false;
    }

    // Initial and final energies
    const auto& E_initial = _time_series.components.front();
    const auto& E_final = _time_series.components.back();
    double t_final = _time_series.time.back();

    report << "Energy components at t=0:" << std::endl;
    report << std::fixed << std::setprecision(8);
    report << "  T_dirac (kinetic):    " << E_initial.T_dirac_kinetic << std::endl;
    report << "  V_dirac (mass):       " << E_initial.V_dirac_mass << std::endl;
    report << "  E_kuramoto_gradient:  " << E_initial.E_kuramoto_gradient << std::endl;
    report << "  E_kuramoto_sync:      " << E_initial.E_kuramoto_sync << std::endl;
    report << "  E_coupling:           " << E_initial.E_coupling << std::endl;
    report << "  E_total:              " << E_initial.E_total << std::endl;
    report << std::endl;

    report << "Energy components at t=" << t_final << ":" << std::endl;
    report << "  T_dirac (kinetic):    " << E_final.T_dirac_kinetic << std::endl;
    report << "  V_dirac (mass):       " << E_final.V_dirac_mass << std::endl;
    report << "  E_kuramoto_gradient:  " << E_final.E_kuramoto_gradient << std::endl;
    report << "  E_kuramoto_sync:      " << E_final.E_kuramoto_sync << std::endl;
    report << "  E_coupling:           " << E_final.E_coupling << std::endl;
    report << "  E_total:              " << E_final.E_total << std::endl;
    report << std::endl;

    // Energy conservation metrics
    double E_change = E_final.E_total - E_initial.E_total;
    double E_relative_change = E_change / E_initial.E_total;
    double E_dissipated_total = E_final.E_dissipated_cumulative;

    report << "Energy conservation analysis:" << std::endl;
    report << "  Total energy change ΔE:        " << E_change << std::endl;
    report << "  Relative change ΔE/E₀:         " << E_relative_change << std::endl;
    report << "  Energy dissipated (cumulative): " << E_dissipated_total << std::endl;
    report << "  Energy drift (ΔE + ∫P dt):     " << _time_series.energy_drift.back() << std::endl;
    report << std::endl;

    // Analyze conservation error over time
    double max_conservation_error = 0.0;
    double avg_conservation_error = 0.0;
    for (const auto& comp : _time_series.components) {
        max_conservation_error = std::max(max_conservation_error, comp.conservation_error);
        avg_conservation_error += comp.conservation_error;
    }
    avg_conservation_error /= _time_series.components.size();

    report << "Conservation error |dE/dt + P|:" << std::endl;
    report << "  Maximum:    " << max_conservation_error << std::endl;
    report << "  Average:    " << avg_conservation_error << std::endl;
    report << std::endl;

    // Determine conservation status
    bool conserved = false;
    if (_damping < 1e-6) {
        // Undamped case: should conserve to machine precision
        conserved = (std::abs(E_relative_change) < 1e-4);
        report << "Undamped system (γ=" << _damping << ")" << std::endl;
        report << "  Expected: Energy conserved to ~1e-4" << std::endl;
        report << "  Observed: |ΔE/E₀| = " << std::abs(E_relative_change) << std::endl;
        report << "  Status: " << (conserved ? "PASS" : "FAIL") << std::endl;
    } else {
        // Damped case: energy + dissipation should be conserved
        double drift_relative = std::abs(_time_series.energy_drift.back()) / E_initial.E_total;
        conserved = (drift_relative < 1e-3);
        report << "Damped system (γ=" << _damping << ")" << std::endl;
        report << "  Expected: E(t) - E₀ + ∫P dt ≈ 0" << std::endl;
        report << "  Observed drift: " << _time_series.energy_drift.back() << std::endl;
        report << "  Relative drift: " << drift_relative << std::endl;
        report << "  Status: " << (conserved ? "PASS" : "FAIL") << std::endl;

        // Fit exponential decay if significant damping
        if (_damping > 0.01) {
            double gamma_measured = fitExponentialDecay();
            report << std::endl;
            report << "Exponential decay analysis:" << std::endl;
            report << "  Expected decay rate: " << _damping << std::endl;
            report << "  Measured decay rate: " << gamma_measured << std::endl;
            report << "  Relative error: " << std::abs(gamma_measured - _damping) / _damping << std::endl;
        }
    }

    report << std::endl;
    report << "=== OVERALL RESULT: " << (conserved ? "PASS" : "FAIL") << " ===" << std::endl;

    report.close();
    return conserved;
}

void EnergyBudget::writeTimeSeriesCSV(const std::string& csv_file) {
    std::ofstream csv(csv_file);
    if (!csv.is_open()) {
        std::cerr << "Error: Could not open CSV file " << csv_file << std::endl;
        return;
    }

    // Write header
    csv << "time,E_total,T_dirac,V_dirac,E_kuramoto_grad,E_kuramoto_sync,";
    csv << "E_coupling,E_em,P_dissipated,E_dissipated_cum,";
    csv << "dE_dt,conservation_error,relative_error,energy_drift" << std::endl;

    // Write data
    for (size_t i = 0; i < _time_series.time.size(); ++i) {
        const auto& comp = _time_series.components[i];
        csv << std::fixed << std::setprecision(8);
        csv << _time_series.time[i] << ",";
        csv << comp.E_total << ",";
        csv << comp.T_dirac_kinetic << ",";
        csv << comp.V_dirac_mass << ",";
        csv << comp.E_kuramoto_gradient << ",";
        csv << comp.E_kuramoto_sync << ",";
        csv << comp.E_coupling << ",";
        csv << comp.E_em_field << ",";
        csv << comp.P_dissipated << ",";
        csv << comp.E_dissipated_cumulative << ",";
        csv << comp.dE_dt << ",";
        csv << comp.conservation_error << ",";
        csv << _time_series.relative_error[i] << ",";
        csv << _time_series.energy_drift[i] << std::endl;
    }

    csv.close();
}

std::pair<double, double> EnergyBudget::computeKuramotoEnergy(
    const std::vector<float>& theta,
    const std::vector<float>& R_field,
    float K,
    float dx, float dy)
{
    int Nx = std::sqrt(theta.size());
    int Ny = Nx;  // Assuming square grid

    double E_gradient = 0.0;
    double E_sync = 0.0;

    // Compute phase gradient energy: (1/2)∫|∇θ|² dx
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;

            // Compute gradients using centered differences with periodic BC
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            int jp = (j + 1) % Ny;
            int jm = (j - 1 + Ny) % Ny;

            float dtheta_dx = (theta[j * Nx + ip] - theta[j * Nx + im]) / (2.0f * dx);
            float dtheta_dy = (theta[jp * Nx + i] - theta[jm * Nx + i]) / (2.0f * dy);

            E_gradient += 0.5 * (dtheta_dx * dtheta_dx + dtheta_dy * dtheta_dy) * dx * dy;
        }
    }

    // Compute synchronization energy: -K·Σcos(θᵢ - θⱼ)
    // Sum over nearest neighbors
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;
            float theta_i = theta[idx];

            // Right neighbor
            int i_right = (i + 1) % Nx;
            float theta_right = theta[j * Nx + i_right];
            E_sync -= K * std::cos(theta_right - theta_i) * 0.5;  // 0.5 to avoid double counting

            // Top neighbor
            int j_top = (j + 1) % Ny;
            float theta_top = theta[j_top * Nx + i];
            E_sync -= K * std::cos(theta_top - theta_i) * 0.5;
        }
    }

    return {E_gradient, E_sync};
}

double EnergyBudget::computeCouplingEnergy(
    const std::vector<std::complex<double>>& psi,
    const std::vector<float>& R_field,
    float coupling,
    float dx, float dy)
{
    // Assume 4-component spinor
    int Ntotal = psi.size() / 4;

    double E_coupling = 0.0;

    for (int idx = 0; idx < Ntotal; ++idx) {
        // Compute |ψ|² = sum of all 4 spinor components
        double psi_squared = 0.0;
        for (int comp = 0; comp < 4; ++comp) {
            std::complex<double> psi_comp = psi[comp * Ntotal + idx];
            psi_squared += std::norm(psi_comp);
        }

        // Coupling energy: λ·|ψ|²·R
        E_coupling += coupling * psi_squared * R_field[idx] * dx * dy;
    }

    return E_coupling;
}

double EnergyBudget::computeDissipationRate(
    const std::vector<float>& theta,
    float damping,
    float dx, float dy)
{
    if (damping < 1e-10) return 0.0;  // No dissipation if γ ≈ 0

    int Nx = std::sqrt(theta.size());
    int Ny = Nx;  // Assuming square grid

    double P_dissipated = 0.0;

    // P = γ·∫|∇θ|² dx
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            // Compute gradients using centered differences with periodic BC
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;
            int jp = (j + 1) % Ny;
            int jm = (j - 1 + Ny) % Ny;

            float dtheta_dx = (theta[j * Nx + ip] - theta[j * Nx + im]) / (2.0f * dx);
            float dtheta_dy = (theta[jp * Nx + i] - theta[jm * Nx + i]) / (2.0f * dy);

            P_dissipated += damping * (dtheta_dx * dtheta_dx + dtheta_dy * dtheta_dy) * dx * dy;
        }
    }

    return P_dissipated;
}

double EnergyBudget::computePhaseGradientEnergy(
    const std::vector<float>& theta,
    float dx, float dy)
{
    auto [E_gradient, E_sync] = computeKuramotoEnergy(theta, std::vector<float>(), 0.0f, dx, dy);
    return E_gradient;
}

double EnergyBudget::computeSynchronizationEnergy(
    const std::vector<float>& theta,
    float K)
{
    int Nx = std::sqrt(theta.size());
    int Ny = Nx;

    double E_sync = 0.0;

    // Sum over nearest neighbor pairs
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            float theta_i = theta[j * Nx + i];

            // Right neighbor
            int i_right = (i + 1) % Nx;
            float theta_right = theta[j * Nx + i_right];
            E_sync -= K * std::cos(theta_right - theta_i) * 0.5;

            // Top neighbor
            int j_top = (j + 1) % Ny;
            float theta_top = theta[j_top * Nx + i];
            E_sync -= K * std::cos(theta_top - theta_i) * 0.5;
        }
    }

    return E_sync;
}

bool EnergyBudget::validateConservation(double tolerance) {
    if (_time_series.components.empty()) return false;

    // Check maximum conservation error
    double max_error = 0.0;
    for (const auto& comp : _time_series.components) {
        max_error = std::max(max_error, comp.conservation_error);
    }

    return max_error < tolerance;
}

double EnergyBudget::fitExponentialDecay() {
    if (_time_series.time.size() < 10) return 0.0;

    // Simple linear fit to log(E) vs t
    // E(t) = E₀·exp(-γt) => log(E) = log(E₀) - γt

    std::vector<double> log_E;
    for (const auto& comp : _time_series.components) {
        if (comp.E_total > 0) {
            log_E.push_back(std::log(comp.E_total));
        }
    }

    // Linear regression
    double sum_t = 0.0, sum_logE = 0.0, sum_t2 = 0.0, sum_t_logE = 0.0;
    int n = std::min(log_E.size(), _time_series.time.size());

    for (int i = 0; i < n; ++i) {
        double t = _time_series.time[i];
        double y = log_E[i];
        sum_t += t;
        sum_logE += y;
        sum_t2 += t * t;
        sum_t_logE += t * y;
    }

    // Slope = -γ
    double slope = (n * sum_t_logE - sum_t * sum_logE) / (n * sum_t2 - sum_t * sum_t);

    return -slope;  // Return positive decay rate
}

double EnergyBudget::richardsonExtrapolation(
    const std::vector<double>& dt_values,
    const std::string& output_file)
{
    // This would require running multiple simulations with different dt values
    // For now, return a placeholder
    // In a full implementation, this would:
    // 1. Run simulation with dt, dt/2, dt/4
    // 2. Measure energy conservation error for each
    // 3. Extrapolate to dt→0 using Richardson's method

    std::cerr << "Richardson extrapolation not fully implemented" << std::endl;
    return 0.0;
}

double EnergyBudget::compareIntegrators(const std::string& output_file) {
    // This would require implementing Euler stepping as a comparison
    // For now, return a placeholder
    // In a full implementation, this would:
    // 1. Run with Strang splitting (current)
    // 2. Run with first-order Euler
    // 3. Compare energy conservation errors

    std::cerr << "Integrator comparison not fully implemented" << std::endl;
    return 1.0;
}

std::map<std::string, double> EnergyBudget::getConservationMetrics() const {
    std::map<std::string, double> metrics;

    if (!_time_series.components.empty()) {
        const auto& E_initial = _time_series.components.front();
        const auto& E_final = _time_series.components.back();

        metrics["E_initial"] = E_initial.E_total;
        metrics["E_final"] = E_final.E_total;
        metrics["E_change"] = E_final.E_total - E_initial.E_total;
        metrics["E_relative_change"] = (E_final.E_total - E_initial.E_total) / E_initial.E_total;
        metrics["E_dissipated_total"] = E_final.E_dissipated_cumulative;
        metrics["energy_drift"] = _time_series.energy_drift.back();

        // Compute max relative error
        double max_rel_error = 0.0;
        for (double err : _time_series.relative_error) {
            max_rel_error = std::max(max_rel_error, err);
        }
        metrics["max_relative_error"] = max_rel_error;

        // Average conservation error
        double avg_cons_error = 0.0;
        for (const auto& comp : _time_series.components) {
            avg_cons_error += comp.conservation_error;
        }
        avg_cons_error /= _time_series.components.size();
        metrics["avg_conservation_error"] = avg_cons_error;
    }

    return metrics;
}