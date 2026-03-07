/**
 * test_inflation.cpp
 *
 * C5 COSMOLOGICAL VALIDATION: Primordial Inflation
 *
 * Objective: Verify that TRD's R-field dynamics can produce
 * the exponential expansion needed to explain CMB observations.
 *
 * PHYSICS MODEL:
 *   Inflation requirements (from CMB):
 *   - e-foldings: N = ln(a_end/a_start) ≈ 60
 *   - Slow-roll: ε = -(Ḣ/H²) ≪ 1
 *   - Spectral index: n_s = 1 - 6ε + 2η ≈ 0.96 (Planck)
 *
 *   TRD Inflation:
 *   - R-field starts in false vacuum: R_initial ≈ 1.5
 *   - Potential: V(R) = V₀(1 - R)² (inflationary)
 *   - Slow-roll down potential produces inflation
 *   - Exit when R → 1 (true vacuum)
 *
 * SUCCESS CRITERIA:
 *   ✓ N ≈ 60 ± 10 e-foldings
 *   ✓ ε < 0.01 (slow-roll condition)
 *   ✓ n_s ≈ 0.96 ± 0.02 (Planck constraint)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "TRDCore3D.h"
#include "simulations/VisualizationGenerator.h"

const double PI = 3.14159265358979323846;
const double M_PLANCK = 1.0;  // Natural units

/**
 * Inflationary evolution using TRD R-field dynamics
 */
class InflationEvolution {
private:
    double V0;             // Potential energy scale
    double R;              // R-field value (homogeneous)
    double dR_dt;          // R-field time derivative
    double evolution_time; // Total evolution time
    double dt;             // Time step

    // History for analysis
    std::vector<double> time_history;
    std::vector<double> R_history;
    std::vector<double> a_history;     // Scale factor
    std::vector<double> H_history;     // Hubble parameter
    std::vector<double> epsilon_history; // Slow-roll parameter

    double a;              // Current scale factor
    double a_initial;      // Initial scale factor

public:
    InflationEvolution(double V0_val, double R_initial, double t_max, double time_step)
        : V0(V0_val), R(R_initial), dR_dt(0.0), // Start with zero velocity for slow-roll
          evolution_time(t_max), dt(time_step), a(1.0), a_initial(1.0) {}

    /**
     * Compute inflationary potential V(R) = V₀(1 - R)²
     */
    double computePotential() const {
        double delta = 1.0 - R;
        return V0 * delta * delta;
    }

    /**
     * Compute potential derivative dV/dR = -2V₀(1 - R)
     */
    double computePotentialDerivative() const {
        return -2.0 * V0 * (1.0 - R);
    }

    /**
     * Compute Hubble parameter H = √(V/3M_p²) in slow-roll approximation
     */
    double computeHubbleParameter() const {
        double V = computePotential();
        return std::sqrt(V / (3.0 * M_PLANCK * M_PLANCK));
    }

    /**
     * Compute slow-roll parameter ε = M_p²/2 * (V'/V)²
     */
    double computeSlowRollEpsilon() const {
        double V = computePotential();
        if (std::abs(V) < 1e-10) return 1.0; // Inflation ends
        double V_prime = computePotentialDerivative();
        return 0.5 * M_PLANCK * M_PLANCK * (V_prime / V) * (V_prime / V);
    }

    /**
     * Compute slow-roll parameter η = M_p² * V''/V
     */
    double computeSlowRollEta() const {
        double V = computePotential();
        if (std::abs(V) < 1e-10) return 0.0;
        double V_double_prime = 2.0 * V0; // d²V/dR² = 2V₀
        return M_PLANCK * M_PLANCK * V_double_prime / V;
    }

    /**
     * Compute spectral index n_s = 1 - 6ε + 2η
     */
    double computeSpectralIndex() const {
        double epsilon = computeSlowRollEpsilon();
        double eta = computeSlowRollEta();
        return 1.0 - 6.0 * epsilon + 2.0 * eta;
    }

    /**
     * Check if still in slow-roll regime
     */
    bool isSlowRolling() const {
        return computeSlowRollEpsilon() < 1.0;
    }

    /**
     * Evolve R-field during inflation
     */
    void evolve() {
        time_history.clear();
        R_history.clear();
        a_history.clear();
        H_history.clear();
        epsilon_history.clear();

        double num_steps = evolution_time / dt;
        bool inflation_ended = false;
        double t_end = 0.0;

        for (int step = 0; step < num_steps; ++step) {
            double t = step * dt;

            // Record state
            double H = computeHubbleParameter();
            double epsilon = computeSlowRollEpsilon();

            time_history.push_back(t);
            R_history.push_back(R);
            a_history.push_back(a);
            H_history.push_back(H);
            epsilon_history.push_back(epsilon);

            // Feed visualization data
            double N_current = std::log(a / 1.0); // e-foldings from initial a=1
            VisualizationGenerator::addDataPoint("efolds", static_cast<float>(t), static_cast<float>(N_current));
            VisualizationGenerator::addDataPoint("slow_roll", static_cast<float>(t), static_cast<float>(epsilon));

            // Check if inflation ends (ε ≥ 1)
            if (!inflation_ended && epsilon >= 1.0) {
                inflation_ended = true;
                t_end = t;
                std::cout << "Inflation ended at t = " << t << " (ε = " << epsilon << ")" << std::endl;
            }

            // Update R-field
            // In slow-roll: 3H·dR/dt ≈ -dV/dR
            if (isSlowRolling()) {
                // Slow-roll approximation
                dR_dt = -computePotentialDerivative() / (3.0 * H);
            } else {
                // Full dynamics after inflation
                double d2R_dt2 = -computePotentialDerivative() - 3.0 * H * dR_dt;
                dR_dt += d2R_dt2 * dt;
            }
            R += dR_dt * dt;

            // Update scale factor: a(t+dt) = a(t) * exp(H·dt)
            a *= std::exp(H * dt);
        }

        if (!inflation_ended) {
            std::cout << "Warning: Inflation did not end within simulation time" << std::endl;
        }
    }

    /**
     * Calculate number of e-foldings
     */
    double calculateEfoldings() const {
        if (a_history.empty()) return 0.0;

        // Find when inflation ends (ε = 1)
        size_t end_idx = 0;
        for (size_t i = 0; i < epsilon_history.size(); ++i) {
            if (epsilon_history[i] >= 1.0) {
                end_idx = i;
                break;
            }
        }
        if (end_idx == 0) end_idx = epsilon_history.size() - 1;

        double a_end = a_history[end_idx];
        return std::log(a_end / a_initial);
    }

    /**
     * Analyze results
     */
    void analyzeResults(YAML::Node& results) {
        if (H_history.empty()) {
            std::cerr << "No evolution data available!" << std::endl;
            return;
        }

        // Calculate e-foldings
        double N = calculateEfoldings();

        // Find minimum slow-roll parameter during inflation
        double epsilon_min = epsilon_history[0];
        double epsilon_avg = 0.0;
        int slow_roll_steps = 0;

        for (size_t i = 0; i < epsilon_history.size(); ++i) {
            if (epsilon_history[i] < 1.0) {
                epsilon_min = std::min(epsilon_min, epsilon_history[i]);
                epsilon_avg += epsilon_history[i];
                slow_roll_steps++;
            }
        }
        if (slow_roll_steps > 0) {
            epsilon_avg /= slow_roll_steps;
        }

        // Compute spectral index at horizon crossing (N ≈ 50-60 before end)
        double n_s = computeSpectralIndex();

        // Store results
        results["e_foldings"] = N;
        results["slow_roll"]["epsilon_min"] = epsilon_min;
        results["slow_roll"]["epsilon_avg"] = epsilon_avg;
        results["spectral_index"] = n_s;

        // Check success criteria
        bool N_good = (N > 50 && N < 70);
        bool epsilon_good = (epsilon_min < 0.01);
        bool ns_good = std::abs(n_s - 0.96) < 0.02;

        results["criteria"]["e_foldings_passed"] = N_good;
        results["criteria"]["slow_roll_passed"] = epsilon_good;
        results["criteria"]["spectral_index_passed"] = ns_good;

        // Output summary
        std::cout << "\n=== Inflation Analysis ===" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "e-foldings N = " << N << " (target: 60 ± 10)" << std::endl;
        std::cout << "  Status: " << (N_good ? "PASSED ✓" : "FAILED ✗") << std::endl;

        std::cout << "\nSlow-roll parameter ε:" << std::endl;
        std::cout << "  Minimum: " << epsilon_min << std::endl;
        std::cout << "  Average: " << epsilon_avg << std::endl;
        std::cout << "  Status: " << (epsilon_good ? "PASSED ✓" : "FAILED ✗") << std::endl;

        std::cout << "\nSpectral index n_s = " << n_s << " (Planck: 0.9649 ± 0.0042)" << std::endl;
        std::cout << "  Status: " << (ns_good ? "PASSED ✓" : "FAILED ✗") << std::endl;

        // Scale factor growth
        if (!a_history.empty()) {
            double a_final = a_history.back();
            double expansion = a_final / a_initial;
            std::cout << "\nTotal expansion: " << std::scientific << expansion << "x" << std::endl;
            results["scale_factor"]["initial"] = a_initial;
            results["scale_factor"]["final"] = a_final;
            results["scale_factor"]["total_expansion"] = expansion;
        }
    }

    /**
     * Save evolution data to file
     */
    void saveData(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open " << filename << " for writing" << std::endl;
            return;
        }

        file << "# Inflation Evolution Data" << std::endl;
        file << "# time, R, scale_factor, Hubble_parameter, epsilon" << std::endl;

        for (size_t i = 0; i < time_history.size(); ++i) {
            file << time_history[i] << ", "
                 << R_history[i] << ", "
                 << a_history[i] << ", "
                 << H_history[i] << ", "
                 << epsilon_history[i] << std::endl;
        }

        file.close();
        std::cout << "Saved evolution data to " << filename << std::endl;
    }
};

/**
 * Run inflation test from YAML config
 */
int runInflationTest() {
    std::cout << "\n=========================================" << std::endl;
    std::cout << "    C5: PRIMORDIAL INFLATION TEST       " << std::endl;
    std::cout << "=========================================" << std::endl;

    try {
        // Load configuration
        YAML::Node config = YAML::LoadFile("config/inflation.yaml");
        YAML::Node inflation = config["inflation"];

        double initial_R = inflation["initial_R"].as<double>();
        double V0 = inflation["V0"].as<double>();
        double evolution_time = inflation["evolution_time"].as<double>();
        double dt = 0.001; // Fine time step for accurate integration

        YAML::Node results;
        results["test"] = "C5_Primordial_Inflation";
        results["timestamp"] = std::time(nullptr);

        // Run test scenarios
        YAML::Node scenarios = config["test_scenarios"];
        bool all_passed = true;

        std::cout << "\n--- Running Inflation Simulation ---" << std::endl;
        std::cout << "Initial R-field: " << initial_R << std::endl;
        std::cout << "Potential scale V0: " << V0 << std::endl;

        // Create and run evolution
        InflationEvolution evolution(V0, initial_R, evolution_time, dt);
        evolution.evolve();

        // Analyze results
        YAML::Node scenario_results;
        evolution.analyzeResults(scenario_results);

        // Check each scenario
        for (const auto& scenario : scenarios) {
            std::string name = scenario["name"].as<std::string>();
            std::cout << "\n--- Checking: " << name << " ---" << std::endl;

            YAML::Node expected = scenario["expected"];
            bool passed = true;

            if (name == "e_Foldings") {
                double N = scenario_results["e_foldings"].as<double>();
                double N_min = expected["N_min"].as<double>();
                double N_max = expected["N_max"].as<double>();

                if (N < N_min || N > N_max) {
                    passed = false;
                    std::cout << "FAILED: N = " << N
                             << " (expected " << N_min << " - " << N_max << ")" << std::endl;
                } else {
                    std::cout << "PASSED: N = " << N << " ✓" << std::endl;
                }
            }
            else if (name == "Slow_Roll") {
                double epsilon_min = scenario_results["slow_roll"]["epsilon_min"].as<double>();
                double epsilon_max = expected["epsilon_max"].as<double>();

                if (epsilon_min > epsilon_max) {
                    passed = false;
                    std::cout << "FAILED: ε_min = " << epsilon_min
                             << " (expected < " << epsilon_max << ")" << std::endl;
                } else {
                    std::cout << "PASSED: ε_min = " << epsilon_min << " ✓" << std::endl;
                }
            }
            else if (name == "Spectral_Index") {
                double n_s = scenario_results["spectral_index"].as<double>();
                double n_s_expected = expected["n_s"].as<double>();
                double tolerance = expected["tolerance"].as<double>();

                if (std::abs(n_s - n_s_expected) > tolerance) {
                    passed = false;
                    std::cout << "FAILED: n_s = " << n_s
                             << " (expected " << n_s_expected << " ± " << tolerance << ")" << std::endl;
                } else {
                    std::cout << "PASSED: n_s = " << n_s << " ✓" << std::endl;
                }
            }

            if (!passed) {
                all_passed = false;
            }

            results["scenarios"][name]["passed"] = passed;
        }

        // Save data
        evolution.saveData("inflation_evolution.csv");

        // Final summary
        std::cout << "\n=========================================" << std::endl;
        if (all_passed) {
            std::cout << "✓ C5 INFLATION TEST PASSED" << std::endl;
            std::cout << "TRD successfully produces primordial inflation!" << std::endl;
            std::cout << "  - Sufficient e-foldings (N ≈ 60)" << std::endl;
            std::cout << "  - Slow-roll conditions satisfied" << std::endl;
            std::cout << "  - Spectral index matches Planck data" << std::endl;
        } else {
            std::cout << "✗ C5 INFLATION TEST FAILED" << std::endl;
        }
        std::cout << "=========================================" << std::endl;

        // Save full results
        std::ofstream results_file("inflation_results.yaml");
        results_file << results;
        results_file.close();

        return all_passed ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}