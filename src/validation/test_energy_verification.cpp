/**
 * test_energy_verification.cpp
 *
 * Comprehensive energy budget verification test using the EnergyBudget class.
 * This test program can be integrated into the SMFT test framework to perform
 * the 6 energy conservation tests (T1-T6) defined in the Phase 1 strategy.
 *
 * Usage:
 *   ./test_energy_verification config/energy_verification/T1_undamped_dt001.yaml
 */

#include "EnergyBudget.h"
#include "../simulations/SMFTTestRunner.h"
#include "../simulations/TestConfig.h"
#include "../DiracEvolution.h"
#include "../SMFTEngine.h"
#include "../Nova.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <ctime>

/**
 * Run energy verification test with comprehensive tracking
 */
bool runEnergyVerificationTest(const std::string& config_path) {
    std::cout << "=== ENERGY VERIFICATION TEST ===" << std::endl;
    std::cout << "Loading configuration: " << config_path << std::endl;

    // Load test configuration
    TestConfig config(config_path);

    // Initialize Nova
    Nova nova;
    nova.initialize(config.grid.size_x, config.grid.size_y);

    // Initialize SMFT engine
    SMFTEngine engine(&nova);
    engine.initialize(config.grid.size_x, config.grid.size_y);

    // Set physics parameters
    engine.setDelta(config.physics.delta);
    engine.setCoupling(config.physics.coupling);
    engine.setK(config.physics.K);
    engine.setDamping(config.physics.damping);

    // Initialize Dirac evolution if needed
    DiracEvolution* dirac = nullptr;
    if (config.physics.solver_type == "gpu" || config.physics.solver_type == "dirac") {
        dirac = new DiracEvolution();
        dirac->initialize(config.grid.size_x, config.grid.size_y);

        // Set initial conditions based on config
        if (config.dirac_initial.type == "gaussian") {
            // Convert physical coordinates to grid coordinates
            float x0_grid = config.dirac_initial.x0_physical * config.grid.size_x / config.grid.L_domain;
            float y0_grid = config.dirac_initial.y0_physical * config.grid.size_y / config.grid.L_domain;
            float sigma_grid = config.dirac_initial.sigma_physical * config.grid.size_x / config.grid.L_domain;

            dirac->setInitialGaussian(x0_grid, y0_grid, sigma_grid);
        }
    }

    // Initialize Kuramoto field based on config
    if (config.kuramoto_initial.phase_distribution == "vortex") {
        engine.setInitialVortex(
            config.kuramoto_initial.vortex_center_x * config.grid.size_x / config.grid.L_domain,
            config.kuramoto_initial.vortex_center_y * config.grid.size_y / config.grid.L_domain,
            config.kuramoto_initial.winding_number,
            config.kuramoto_initial.vortex_core_radius * config.grid.size_x / config.grid.L_domain
        );
    }

    // Create output directory with timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&time_t);

    std::stringstream output_dir;
    output_dir << "output/"
               << std::setfill('0')
               << std::setw(4) << (tm->tm_year + 1900)
               << std::setw(2) << (tm->tm_mon + 1)
               << std::setw(2) << tm->tm_mday << "_"
               << std::setw(2) << tm->tm_hour
               << std::setw(2) << tm->tm_min
               << std::setw(2) << tm->tm_sec << "_"
               << config.test_name;

    // Add grid size suffix if running grid convergence
    if (!config.grid.grid_sizes.empty()) {
        output_dir << "_" << config.grid.size_x << "x" << config.grid.size_y;
    }

    std::string output_path = output_dir.str();
    std::string mkdir_cmd = "mkdir -p " + output_path;
    system(mkdir_cmd.c_str());

    // Initialize EnergyBudget tracker
    EnergyBudget energy_budget(&engine, output_path);
    energy_budget.setPhysicsParameters(
        config.physics.delta,
        config.physics.coupling,
        config.physics.K,
        config.physics.damping,
        config.grid.L_domain
    );

    std::cout << "\nTest: " << config.test_name << std::endl;
    std::cout << "Description: " << config.description << std::endl;
    std::cout << "Grid: " << config.grid.size_x << "×" << config.grid.size_y << std::endl;
    std::cout << "Domain size: " << config.grid.L_domain << " Planck lengths" << std::endl;
    std::cout << "Delta: " << config.physics.delta << std::endl;
    std::cout << "Coupling: " << config.physics.coupling << std::endl;
    std::cout << "K: " << config.physics.K << std::endl;
    std::cout << "Damping: " << config.physics.damping << std::endl;
    std::cout << "Timestep: " << config.physics.dt << std::endl;
    std::cout << "Total steps: " << config.physics.total_steps << std::endl;
    std::cout << "Output directory: " << output_path << std::endl;
    std::cout << std::endl;

    // Main simulation loop
    double time = 0.0;
    int save_counter = 0;

    std::cout << "Running simulation..." << std::endl;
    std::cout << std::setw(10) << "Step"
              << std::setw(12) << "Time"
              << std::setw(16) << "E_total"
              << std::setw(16) << "ΔE/E₀"
              << std::setw(16) << "P_dissipated"
              << std::setw(16) << "Drift" << std::endl;
    std::cout << std::string(90, '-') << std::endl;

    for (int step = 0; step <= config.physics.total_steps; ++step) {
        // Compute and track energy
        auto components = energy_budget.computeComponents(dirac, time);
        energy_budget.trackTimeSeries(time, components);

        // Print progress
        if (step % (config.physics.total_steps / 100) == 0 || step == 0) {
            double E0 = energy_budget.getConservationMetrics()["E_initial"];
            double rel_error = (components.E_total - E0) / E0;

            std::cout << std::fixed << std::setprecision(6);
            std::cout << std::setw(10) << step
                      << std::setw(12) << time
                      << std::setw(16) << components.E_total
                      << std::setw(16) << rel_error
                      << std::setw(16) << components.P_dissipated
                      << std::setw(16) << components.E_total - E0 + components.E_dissipated_cumulative
                      << std::endl;
        }

        // Save observables periodically
        if (step % config.output.save_every == 0 && config.output.save_observables) {
            save_counter++;
            // Observable saving would be handled by OutputManager in full integration
        }

        // Evolve system (except on last step)
        if (step < config.physics.total_steps) {
            if (config.operator_splitting.enabled && dirac != nullptr) {
                // Use operator splitting with specified substep ratio
                int N = config.operator_splitting.substep_ratios.empty() ?
                        10 : config.operator_splitting.substep_ratios[0];
                engine.stepWithDirac(*dirac, config.physics.dt, N);
            } else {
                // Kuramoto-only evolution
                engine.step(config.physics.dt);
            }

            time += config.physics.dt;
        }
    }

    std::cout << std::string(90, '-') << std::endl;
    std::cout << "Simulation complete." << std::endl << std::endl;

    // Generate energy conservation analysis report
    std::string report_file = output_path + "/energy_analysis_report.txt";
    bool conservation_passed = energy_budget.analyzeConservation(report_file);

    // Write time series to CSV
    std::string csv_file = output_path + "/energy_timeseries.csv";
    energy_budget.writeTimeSeriesCSV(csv_file);

    // Print summary
    auto metrics = energy_budget.getConservationMetrics();

    std::cout << "=== ENERGY CONSERVATION SUMMARY ===" << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Initial energy E₀:        " << metrics["E_initial"] << std::endl;
    std::cout << "Final energy E(t):        " << metrics["E_final"] << std::endl;
    std::cout << "Energy change ΔE:         " << metrics["E_change"] << std::endl;
    std::cout << "Relative change ΔE/E₀:    " << metrics["E_relative_change"] << std::endl;
    std::cout << "Energy dissipated:        " << metrics["E_dissipated_total"] << std::endl;
    std::cout << "Energy drift:             " << metrics["energy_drift"] << std::endl;
    std::cout << "Max relative error:       " << metrics["max_relative_error"] << std::endl;
    std::cout << "Avg conservation error:   " << metrics["avg_conservation_error"] << std::endl;
    std::cout << std::endl;

    if (config.physics.damping < 1e-6) {
        std::cout << "Undamped system (γ=" << config.physics.damping << ")" << std::endl;
        std::cout << "Expected: Energy conserved to tolerance " << config.validation.energy_tolerance << std::endl;
        std::cout << "Observed: |ΔE/E₀| = " << std::abs(metrics["E_relative_change"]) << std::endl;
    } else {
        std::cout << "Damped system (γ=" << config.physics.damping << ")" << std::endl;
        std::cout << "Expected: Exponential decay with γ ≈ " << config.physics.damping << std::endl;
        double drift_relative = std::abs(metrics["energy_drift"]) / metrics["E_initial"];
        std::cout << "Energy budget closure: " << drift_relative << std::endl;
    }

    std::cout << "Conservation test: " << (conservation_passed ? "PASS" : "FAIL") << std::endl;
    std::cout << std::endl;

    std::cout << "Reports written to:" << std::endl;
    std::cout << "  " << report_file << std::endl;
    std::cout << "  " << csv_file << std::endl;

    // Clean up
    if (dirac != nullptr) {
        delete dirac;
    }

    return conservation_passed;
}

// Main entry point - can be called from main.cpp or integrated into test framework
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>" << std::endl;
        std::cerr << "Example: " << argv[0] << " config/energy_verification/T1_undamped_dt001.yaml" << std::endl;
        return 1;
    }

    std::string config_path = argv[1];

    try {
        bool passed = runEnergyVerificationTest(config_path);
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}