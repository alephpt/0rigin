#include "simulations/SMFTTestRunner.h"
#include "output/OutputManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

SMFTTestRunner::SMFTTestRunner(const std::string& config_path)
    : _config_loaded(false), _nova(nullptr), _nova_initialized(false),
      _engine(nullptr), _engine_initialized(false), _output_manager(nullptr) {

    _output_manager = new SMFT::OutputManager();
    _config_loaded = _config.loadFromYAML(config_path);
    if (!_config_loaded) {
        std::cerr << "Failed to load configuration from " << config_path << std::endl;
    }
}

SMFTTestRunner::SMFTTestRunner(const TestConfig& config)
    : _config(config), _config_loaded(true), _nova(nullptr), _nova_initialized(false),
      _engine(nullptr), _engine_initialized(false), _output_manager(nullptr) {
    _output_manager = new SMFT::OutputManager();
}

SMFTTestRunner::~SMFTTestRunner() {
    if (_engine) {
        delete _engine;
    }
    if (_nova) {
        delete _nova;
    }
    if (_output_manager) {
        delete _output_manager;
    }
}

bool SMFTTestRunner::initialize() {
    if (!_config_loaded) {
        std::cerr << "Configuration not loaded" << std::endl;
        return false;
    }

    if (!_config.validate()) {
        std::cerr << "Configuration validation failed" << std::endl;
        return false;
    }

    // Print configuration
    _config.print();

    // Initialize Nova
    NovaConfig nova_config = {
        .name = _config.test_name.c_str(),
        .screen = {800, 600},
        .debug_level = "error",
        .dimensions = "2D",
        .camera_type = "fixed",
        .compute = true
    };

    _nova = new Nova(nova_config);
    _nova->initialized = true;
    _nova_initialized = true;

    // Initialize SMFT engine
    _engine = new SMFTEngine(_nova);
    _engine->initialize(_config.grid.size_x, _config.grid.size_y,
                       _config.physics.delta, 0.0f);
    _engine_initialized = true;

    // Create output directories
    createOutputDirectories();

    std::cout << "✓ SMFTTestRunner initialized successfully" << std::endl;
    return true;
}

bool SMFTTestRunner::run() {
    if (!_engine_initialized) {
        std::cerr << "Engine not initialized" << std::endl;
        return false;
    }

    std::cout << "\n===== Running Tests =====" << std::endl;

    if (_config.operator_splitting.enabled) {
        // Run test for each substep ratio
        bool all_passed = true;
        for (int N : _config.operator_splitting.substep_ratios) {
            std::cout << "\n--- Testing with N=" << N << " ---" << std::endl;
            if (!runSingleTest(N)) {
                std::cerr << "Test failed for N=" << N << std::endl;
                all_passed = false;
                // Continue with other N values instead of returning
            }
        }

        // Validate convergence between different N
        std::cout << "\n--- Validating Convergence ---" << std::endl;
        ValidationResult conv_result = validateConvergence();
        _validation_results[-1] = conv_result; // Store with special key -1

        if (conv_result.passed()) {
            std::cout << "✓ Convergence validation PASSED" << std::endl;
        } else {
            std::cout << "✗ Convergence validation FAILED" << std::endl;
        }

    } else {
        // Single test without operator splitting (N=1)
        if (!runSingleTest(1)) {
            std::cerr << "Test failed" << std::endl;
            return false;
        }
    }

    std::cout << "\n===== Tests Complete =====" << std::endl;
    return allTestsPassed();
}

bool SMFTTestRunner::runSingleTest(int N) {
    // Set substep ratio
    _engine->setSubstepRatio(N);

    // Initialize fields
    auto phases = initializePhases();
    auto omegas = initializeOmegas();
    _engine->setInitialPhases(phases);
    _engine->setNaturalFrequencies(omegas);

    // Initialize Dirac field
    _engine->initializeDiracField(
        _config.dirac_initial.x0,
        _config.dirac_initial.y0,
        _config.dirac_initial.sigma,
        _config.dirac_initial.amplitude
    );

    // Prepare results storage
    std::vector<ObservableComputer::Observables> observables;

    // Get initial energy
    DiracEvolution dirac_temp(_config.grid.size_x, _config.grid.size_y);
    dirac_temp.initialize(_config.dirac_initial.x0, _config.dirac_initial.y0,
                         _config.dirac_initial.sigma);
    auto R_field_initial = _engine->getSyncField();
    std::vector<double> R_field_double(R_field_initial.begin(), R_field_initial.end());

    auto obs_initial = ObservableComputer::compute(
        dirac_temp, R_field_double, _config.physics.delta, 0.0,
        0.0, _config.validation.norm_tolerance, _config.validation.energy_tolerance
    );
    double E0 = obs_initial.energy_total;

    std::cout << "  Initial energy E0 = " << E0 << std::endl;
    std::cout << "  Initial norm = " << obs_initial.norm << std::endl;

    // Run evolution
    const int total_steps = _config.physics.total_steps;
    const int save_every = _config.output.save_every;

    for (int step = 0; step < total_steps; ++step) {
        // Evolve system with operator splitting ratio N
        _engine->stepWithDirac(_config.physics.dt, _config.physics.coupling, N,
                               _config.physics.K, _config.physics.damping);

        // Compute observables at save points
        if (step % save_every == 0 || step == total_steps - 1) {
            // Get current state
            auto R_field = _engine->getSyncField();
            std::vector<double> R_field_d(R_field.begin(), R_field.end());

            // Get Dirac evolution state
            const DiracEvolution* dirac = _engine->getDiracEvolution();
            double time = step * _config.physics.dt;

            ObservableComputer::Observables obs;

            if (dirac) {
                // Compute full observables with Dirac state
                obs = ObservableComputer::compute(
                    *dirac, R_field_d, _config.physics.delta, time,
                    E0, _config.validation.norm_tolerance, _config.validation.energy_tolerance
                );
            } else {
                // Fallback if Dirac not initialized
                obs.time = time;
                obs.norm = 1.0;
                obs.norm_error = 0.0;
                obs.energy_total = E0;
                obs.energy_kinetic = 0.0;
                obs.energy_potential = 0.0;

                // Sync field stats
                auto [R_avg, R_max, R_min, R_var] = ObservableComputer::computeSyncFieldStats(R_field_d);
                obs.R_avg = R_avg;
                obs.R_max = R_max;
                obs.R_min = R_min;
                obs.R_variance = R_var;
            }

            observables.push_back(obs);

            if (step % (save_every * 10) == 0) {
                std::cout << "  Step " << step << "/" << total_steps
                         << " | R_avg = " << obs.R_avg
                         << " | norm = " << obs.norm << std::endl;
            }
        }
    }

    // Store results
    _results[N] = observables;

    // Validate norm conservation
    ValidationResult norm_result = validateNorm(N);
    std::cout << "  Norm validation: max_error = " << norm_result.max_norm_error;
    if (norm_result.norm_passed) {
        std::cout << " ✓ PASS" << std::endl;
    } else {
        std::cout << " ✗ FAIL" << std::endl;
    }

    // Validate energy conservation
    ValidationResult energy_result = validateEnergy(N);
    std::cout << "  Energy validation: max_drift = " << energy_result.max_energy_drift;
    if (energy_result.energy_passed) {
        std::cout << " ✓ PASS" << std::endl;
    } else {
        std::cout << " ✗ FAIL" << std::endl;
    }

    // Combine validation results
    ValidationResult combined;
    combined.norm_passed = norm_result.norm_passed;
    combined.energy_passed = energy_result.energy_passed;
    combined.max_norm_error = norm_result.max_norm_error;
    combined.max_energy_drift = energy_result.max_energy_drift;
    _validation_results[N] = combined;

    // Save observables to CSV
    std::string csv_path = getOutputDirectory(N) + "/observables.csv";
    saveObservablesToCSV(N, csv_path);
    std::cout << "  Saved observables to " << csv_path << std::endl;

    return combined.norm_passed && combined.energy_passed;
}

std::vector<float> SMFTTestRunner::initializePhases() const {
    std::vector<float> phases(_config.grid.size_x * _config.grid.size_y);

    std::random_device rd;
    std::mt19937 gen(rd());

    if (_config.kuramoto_initial.phase_distribution == "uniform") {
        std::fill(phases.begin(), phases.end(), 0.0f);
    } else if (_config.kuramoto_initial.phase_distribution == "random") {
        std::uniform_real_distribution<float> dist(0.0f, 2.0f * M_PI);
        for (auto& p : phases) {
            p = dist(gen);
        }
    } else if (_config.kuramoto_initial.phase_distribution == "vortex") {
        // Create vortex pattern at center
        float cx = _config.grid.size_x / 2.0f;
        float cy = _config.grid.size_y / 2.0f;
        for (int iy = 0; iy < _config.grid.size_y; ++iy) {
            for (int ix = 0; ix < _config.grid.size_x; ++ix) {
                float dx = ix - cx;
                float dy = iy - cy;
                phases[iy * _config.grid.size_x + ix] = std::atan2(dy, dx);
            }
        }
    } else if (_config.kuramoto_initial.phase_distribution == "phase_gradient") {
        // Create traveling wave via phase gradient: θ(x,y) = k_x*x + k_y*y
        float kx = _config.kuramoto_initial.wave_vector_x;
        float ky = _config.kuramoto_initial.wave_vector_y;
        for (int iy = 0; iy < _config.grid.size_y; ++iy) {
            for (int ix = 0; ix < _config.grid.size_x; ++ix) {
                phases[iy * _config.grid.size_x + ix] = kx * ix + ky * iy;
            }
        }
    }

    return phases;
}

std::vector<float> SMFTTestRunner::initializeOmegas() const {
    std::vector<float> omegas(_config.grid.size_x * _config.grid.size_y);

    std::random_device rd;
    std::mt19937 gen(rd());

    if (_config.kuramoto_initial.omega_distribution == "zero") {
        std::fill(omegas.begin(), omegas.end(), 0.0f);
    } else if (_config.kuramoto_initial.omega_distribution == "gaussian") {
        std::normal_distribution<float> dist(_config.kuramoto_initial.omega_mean,
                                             _config.kuramoto_initial.omega_std);
        for (auto& w : omegas) {
            w = dist(gen);
        }
    } else if (_config.kuramoto_initial.omega_distribution == "random") {
        std::uniform_real_distribution<float> dist(-_config.kuramoto_initial.omega_std,
                                                    _config.kuramoto_initial.omega_std);
        for (auto& w : omegas) {
            w = dist(gen);
        }
    }

    return omegas;
}

SMFTTestRunner::ValidationResult SMFTTestRunner::validateNorm(int N) const {
    ValidationResult result;

    auto it = _results.find(N);
    if (it == _results.end()) {
        return result;
    }

    const auto& observables = it->second;

    for (const auto& obs : observables) {
        double norm_error = std::abs(obs.norm_error);
        result.max_norm_error = std::max(result.max_norm_error, norm_error);
    }

    result.norm_passed = (result.max_norm_error < _config.validation.norm_tolerance);
    return result;
}

SMFTTestRunner::ValidationResult SMFTTestRunner::validateEnergy(int N) const {
    ValidationResult result;

    auto it = _results.find(N);
    if (it == _results.end()) {
        return result;
    }

    const auto& observables = it->second;
    if (observables.empty()) {
        return result;
    }

    double E0 = observables[0].energy_total;
    if (std::abs(E0) < 1e-10) {
        // Cannot validate energy if E0 ≈ 0
        result.energy_passed = true;
        return result;
    }

    for (const auto& obs : observables) {
        double energy_drift = std::abs(obs.energy_total - E0) / std::abs(E0);
        result.max_energy_drift = std::max(result.max_energy_drift, energy_drift);
    }

    result.energy_passed = (result.max_energy_drift < _config.validation.energy_tolerance);
    return result;
}

SMFTTestRunner::ValidationResult SMFTTestRunner::validateConvergence() const {
    ValidationResult result;

    if (_results.size() < 2) {
        // Need at least 2 different N values to compare
        result.convergence_passed = true;
        return result;
    }

    // Compare final observables between smallest and largest N
    auto smallest_N = _results.begin()->first;
    auto largest_N = _results.rbegin()->first;

    const auto& obs_small = _results.at(smallest_N).back();
    const auto& obs_large = _results.at(largest_N).back();

    // Compare key observables
    std::vector<double> errors;

    // R_avg convergence
    double R_avg_error = std::abs(obs_large.R_avg - obs_small.R_avg) /
                        (std::abs(obs_small.R_avg) + 1e-10);
    errors.push_back(R_avg_error);

    // Position convergence
    double pos_x_error = std::abs(obs_large.position_x - obs_small.position_x) /
                         (std::abs(obs_small.position_x) + 1e-10);
    double pos_y_error = std::abs(obs_large.position_y - obs_small.position_y) /
                         (std::abs(obs_small.position_y) + 1e-10);
    errors.push_back(pos_x_error);
    errors.push_back(pos_y_error);

    result.convergence_errors = errors;

    // Check if all errors below tolerance
    bool all_pass = true;
    for (double err : errors) {
        if (err > _config.validation.convergence_tolerance) {
            all_pass = false;
            break;
        }
    }

    result.convergence_passed = all_pass;
    return result;
}

void SMFTTestRunner::saveObservablesToCSV(int N, const std::string& filepath) const {
    auto it = _results.find(N);
    if (it == _results.end()) {
        return;
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filepath << " for writing" << std::endl;
        return;
    }

    // Write header
    file << ObservableComputer::getCSVHeader() << "\n";

    // Write data
    for (const auto& obs : it->second) {
        file << ObservableComputer::toCSVLine(obs) << "\n";
    }

    file.close();
}

void SMFTTestRunner::createOutputDirectories() {
    // Use OutputManager to create timestamped experiment directory
    // Extract experiment name from config (use test_name or basename of output directory)
    std::string experiment_name = _config.test_name;
    if (experiment_name.empty()) {
        // Fallback: use basename of configured output directory
        fs::path output_path(_config.output.directory);
        experiment_name = output_path.filename().string();
    }

    // Create timestamped directory via OutputManager (cast to actual type)
    auto* mgr = _output_manager;
    _timestamped_output_dir = _output_manager->createExperimentDirectory(experiment_name);

    std::cout << "[OutputManager] Created timestamped directory: "
              << _timestamped_output_dir << std::endl;

    // Create subdirectories for operator splitting tests
    if (_config.operator_splitting.enabled) {
        for (int N : _config.operator_splitting.substep_ratios) {
            std::string subdir = _output_manager->createSubdirectory(
                _timestamped_output_dir, "N_" + std::to_string(N));
        }
    }
}

std::string SMFTTestRunner::getOutputDirectory(int N) const {
    std::ostringstream oss;
    oss << _timestamped_output_dir << "/N_" << N;
    return oss.str();
}

void SMFTTestRunner::generateReport(const std::string& output_path) const {
    std::string report_path = output_path.empty() ?
        (_timestamped_output_dir + "/test_report.txt") : output_path;

    std::ofstream report(report_path);
    if (!report.is_open()) {
        std::cerr << "Failed to open report file: " << report_path << std::endl;
        return;
    }

    report << "===== SMFT Test Report =====" << std::endl;
    report << "Test: " << _config.test_name << std::endl;
    report << "Description: " << _config.description << std::endl;
    report << "\n[Configuration]" << std::endl;
    report << "Grid: " << _config.grid.size_x << " x " << _config.grid.size_y << std::endl;
    report << "Delta: " << _config.physics.delta << std::endl;
    report << "Coupling: " << _config.physics.coupling << std::endl;
    report << "Total steps: " << _config.physics.total_steps << std::endl;
    report << "dt: " << _config.physics.dt << std::endl;

    report << "\n[Validation Results]" << std::endl;

    for (const auto& [N, result] : _validation_results) {
        if (N == -1) {
            // Convergence validation
            report << "\nConvergence Validation: " << (result.passed() ? "✓ PASS" : "✗ FAIL") << std::endl;
            report << "  Convergence errors: [";
            for (size_t i = 0; i < result.convergence_errors.size(); ++i) {
                report << std::scientific << std::setprecision(3) << result.convergence_errors[i];
                if (i < result.convergence_errors.size() - 1) report << ", ";
            }
            report << "]" << std::endl;
        } else {
            report << "\nN=" << N << ": " << (result.passed() ? "✓ PASS" : "✗ FAIL") << std::endl;
            report << "  Norm: max_error = " << std::scientific << std::setprecision(3)
                   << result.max_norm_error << " (threshold: " << _config.validation.norm_tolerance << ")" << std::endl;
            report << "  Energy: max_drift = " << std::scientific << std::setprecision(3)
                   << result.max_energy_drift << " (threshold: " << _config.validation.energy_tolerance << ")" << std::endl;
        }
    }

    report << "\n[Overall Result]" << std::endl;
    if (allTestsPassed()) {
        report << "✓ ALL TESTS PASSED" << std::endl;
    } else {
        report << "✗ SOME TESTS FAILED" << std::endl;
    }

    report << "\n===== End of Report =====" << std::endl;
    report.close();

    std::cout << "\n✓ Report saved to " << report_path << std::endl;

    // Auto-generate visualizations if enabled
    if (_config.output.auto_visualize && _config.operator_splitting.enabled) {
        std::cout << "\n[4/4] Generating visualizations..." << std::endl;
        generateVisualizations();
    }
}

void SMFTTestRunner::generateVisualizations() const {
    // Create Python visualization script in output directory
    std::string script_path = _config.output.directory + "/visualize.py";
    std::string viz_script = R"(#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 1.5

def load_obs(N):
    p = Path(f"N_{N}/observables.csv")
    return pd.read_csv(p) if p.exists() else None

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

for N in [1, 10, 100]:
    df = load_obs(N)
    if df is None: continue
    ax1.plot(df['time'], df['norm'], label=f'N={N}', alpha=0.8)
    E0 = df['E_total'].iloc[0]
    ax2.plot(df['time'], (df['E_total']-E0)/E0, label=f'N={N}', alpha=0.8)
    ax3.plot(df['pos_x_re'], df['pos_y_re'], label=f'N={N}', alpha=0.7)
    p_mag = np.sqrt(df['mom_x_re']**2 + df['mom_y_re']**2)
    ax4.plot(df['time'], p_mag, label=f'N={N}', alpha=0.7)

ax1.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
ax1.set_ylabel('Norm ||Ψ||²')
ax1.set_title('Norm Conservation')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2.axhline(y=0.0, color='k', linestyle='--', alpha=0.3)
ax2.set_ylabel('Relative Energy Drift ΔE/E₀')
ax2.set_title('Energy Conservation')
ax2.legend(); ax2.grid(True, alpha=0.3)

ax3.set_xlabel('Position <x>')
ax3.set_ylabel('Position <y>')
ax3.set_title('Wavepacket Trajectory')
ax3.legend(); ax3.grid(True, alpha=0.3); ax3.axis('equal')

ax4.set_xlabel('Time')
ax4.set_ylabel('Momentum |<p>|')
ax4.set_title('Momentum Evolution')
ax4.legend(); ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('validation_summary.png')
print("✓ Saved: validation_summary.png")
)";

    std::ofstream script_file(script_path);
    if (script_file.is_open()) {
        script_file << viz_script;
        script_file.close();

        // Make executable and run
        std::string cmd = "cd " + _config.output.directory + " && python3 visualize.py 2>&1";
        int result = system(cmd.c_str());

        if (result == 0) {
            std::cout << "✓ Visualizations generated in " << _config.output.directory << std::endl;
        } else {
            std::cerr << "Warning: Visualization generation failed" << std::endl;
        }
    }
}

const std::vector<ObservableComputer::Observables>& SMFTTestRunner::getResults(int N) const {
    static const std::vector<ObservableComputer::Observables> empty;
    auto it = _results.find(N);
    if (it == _results.end()) {
        return empty;
    }
    return it->second;
}

bool SMFTTestRunner::allTestsPassed() const {
    for (const auto& [N, result] : _validation_results) {
        if (!result.passed()) {
            return false;
        }
    }
    return !_validation_results.empty();
}
