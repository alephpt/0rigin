/**
 * test_finite_temperature.cpp
 *
 * F3: Finite Temperature Effects - Thermal Phase Transitions
 *
 * Tests whether TRD phase dynamics reproduce thermal phase transitions
 * through stochastic evolution with thermal noise.
 *
 * Physics:
 *   Langevin dynamics: dθ/dt = f(θ) + sqrt(2kT/γ)·η(t)
 *     f(θ): Deterministic Kuramoto coupling
 *     η(t): Gaussian white noise <η(t)η(t')> = δ(t-t')
 *     T: Temperature in TRD units
 *     γ: Damping coefficient
 *
 *   Order parameter: R = |<e^(iθ)>|
 *     R ≈ 1: Synchronized (low temperature)
 *     R ≈ 0: Desynchronized (high temperature)
 *
 *   Phase transition at critical temperature T_c:
 *     T < T_c: Synchronized phase (ordered)
 *     T > T_c: Desynchronized phase (disordered)
 *
 *   Fluctuation-dissipation theorem: <noise²> = 2γkT
 *
 * TRD Implementation:
 *   Sweep temperature from T=0.1 to T=5.0 TRD units
 *   Measure equilibrium R(T) for each temperature
 *   Identify critical temperature where R drops sharply
 *   Compare to Kuramoto model: T_c ∝ K (coupling strength)
 *
 * Golden Key: 1 TRD unit = 246 GeV
 *   Critical temperature: T_c ~ K/2 for uniform frequency distribution
 *   Quality gate: T_c within 20% of expected value, R transition >50%
 */

#include "TRDCore3D.h"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>
#include <random>
#include <algorithm>

// Physical constants
const float PI = 3.14159265358979323846f;
const float BOLTZMANN = 1.0f;  // k_B = 1 in natural units

// Global config path (for main.cpp integration)
extern std::string g_test_config_path;

/**
 * Thermal Bath Configuration
 *
 * Adds Gaussian white noise to simulate thermal fluctuations
 */
struct ThermalBathConfig {
    float temperature = 1.0f;      // Current temperature (TRD units)
    float damping = 1.0f;          // Damping coefficient γ
    float noise_variance = 0.0f;   // <noise²> = 2γkT (computed)
    bool fluctuation_dissipation = true;  // Enforce F-D theorem
};

/**
 * Stochastic Evolution Engine
 *
 * Extends TRDCore3D with thermal noise for Langevin dynamics
 */
class StochasticEvolution {
public:
    StochasticEvolution(TRDCore3D& core, uint32_t seed = 42)
        : _core(core), _rng(seed), _normal_dist(0.0f, 1.0f) {
        _N_total = core.getTotalPoints();
        _noise_field.resize(_N_total, 0.0f);
    }

    /**
     * Evolve with thermal noise using Euler-Maruyama method
     *
     * dθ/dt = f(θ) + sqrt(2γkT)·η(t)
     *
     * Discretized:
     *   θ(t+dt) = θ(t) + f(θ)·dt + sqrt(2γkT·dt)·N(0,1)
     *
     * Note: We use the symplectic base integrator for f(θ) and add
     * noise correction to preserve fluctuation-dissipation theorem
     */
    void evolveWithNoise(float dt, const ThermalBathConfig& bath) {
        // Compute noise variance from fluctuation-dissipation theorem
        float noise_std = sqrt(2.0f * bath.damping * BOLTZMANN *
                               bath.temperature * dt);

        // Generate noise field
        for (uint32_t idx = 0; idx < _N_total; ++idx) {
            _noise_field[idx] = _normal_dist(_rng) * noise_std;
        }

        // Deterministic evolution (symplectic)
        _core.evolveSymplecticCPU(dt);

        // Add noise perturbation
        auto& theta = _core.getTheta();
        for (uint32_t idx = 0; idx < _N_total; ++idx) {
            theta[idx] += _noise_field[idx];

            // Wrap phase to [-π, π]
            while (theta[idx] > PI) theta[idx] -= 2.0f * PI;
            while (theta[idx] < -PI) theta[idx] += 2.0f * PI;
        }
    }

    /**
     * Equilibrate system at given temperature
     *
     * Run until R field stabilizes (equilibrium)
     *
     * Note: Initial state is set by computePhaseDiagram based on temperature:
     *   T < 0.5: Ordered initial state (all phases aligned)
     *   T ≥ 0.5: Random initial state
     */
    float equilibrate(float temperature, float dt, uint32_t equilibration_steps,
                     const ThermalBathConfig& bath_template) {
        ThermalBathConfig bath = bath_template;
        bath.temperature = temperature;

        // Track R convergence
        std::vector<float> R_history;
        R_history.reserve(100);

        for (uint32_t step = 0; step < equilibration_steps; ++step) {
            evolveWithNoise(dt, bath);

            // Compute R every 10 steps
            if (step % 10 == 0) {
                _core.computeRField();
                float R_avg = _core.getAverageR();
                R_history.push_back(R_avg);
            }
        }

        // Return average R over last 50% of equilibration
        size_t start_idx = R_history.size() / 2;
        float R_equilibrium = 0.0f;
        for (size_t i = start_idx; i < R_history.size(); ++i) {
            R_equilibrium += R_history[i];
        }
        R_equilibrium /= (R_history.size() - start_idx);

        return R_equilibrium;
    }

private:
    TRDCore3D& _core;
    uint32_t _N_total;

    // Random number generation
    std::mt19937 _rng;
    std::normal_distribution<float> _normal_dist;

    // Noise field
    std::vector<float> _noise_field;
};

/**
 * Phase Diagram Computer
 *
 * Sweeps temperature and measures order parameter R(T)
 */
struct PhaseDiagramPoint {
    float temperature;
    float order_parameter;  // R(T)
};

class PhaseDiagramCalculator {
public:
    static std::vector<PhaseDiagramPoint> computePhaseDiagram(
        const std::vector<float>& temperatures,
        TRDCore3D& core,
        const ThermalBathConfig& bath_config,
        float dt,
        uint32_t equilibration_steps,
        uint32_t seed) {

        std::vector<PhaseDiagramPoint> phase_diagram;
        phase_diagram.reserve(temperatures.size());

        std::cout << "\nComputing phase diagram (temperature sweep):" << std::endl;

        for (size_t i = 0; i < temperatures.size(); ++i) {
            float T = temperatures[i];

            // Temperature-dependent initialization
            if (T < 0.5f) {
                // Low temperature: Initialize from ordered state (all phases aligned)
                auto& theta = core.getTheta();
                uint32_t N_total = core.getTotalPoints();
                for (uint32_t idx = 0; idx < N_total; ++idx) {
                    theta[idx] = 0.0f;  // Synchronized phases
                }
                std::cout << "  T = " << std::fixed << std::setprecision(3)
                         << T << " (ordered init) ... " << std::flush;
            } else {
                // High temperature: Random initial condition
                core.initializeRandom(seed + i);
                std::cout << "  T = " << std::fixed << std::setprecision(3)
                         << T << " (random init) ... " << std::flush;
            }

            // Create stochastic evolution engine
            StochasticEvolution stochastic(core, seed + i + 1000);

            // Equilibrate at this temperature
            float R_eq = stochastic.equilibrate(T, dt, equilibration_steps,
                                                bath_config);

            std::cout << "R = " << std::setprecision(4) << R_eq << std::endl;

            phase_diagram.push_back({T, R_eq});
        }

        return phase_diagram;
    }

    /**
     * Find critical temperature where R drops below threshold
     *
     * Uses linear interpolation between points
     */
    static float findCriticalTemperature(
        const std::vector<PhaseDiagramPoint>& phase_diagram,
        float R_threshold = 0.5f) {

        // Find crossing point where R < threshold
        for (size_t i = 1; i < phase_diagram.size(); ++i) {
            if (phase_diagram[i-1].order_parameter >= R_threshold &&
                phase_diagram[i].order_parameter < R_threshold) {

                // Linear interpolation
                float T1 = phase_diagram[i-1].temperature;
                float T2 = phase_diagram[i].temperature;
                float R1 = phase_diagram[i-1].order_parameter;
                float R2 = phase_diagram[i].order_parameter;

                float T_c = T1 + (R_threshold - R1) * (T2 - T1) / (R2 - R1);
                return T_c;
            }
        }

        // No transition found
        return -1.0f;
    }
};

/**
 * Main test function
 */
int runFiniteTemperatureTest() {
    std::string config_path = g_test_config_path.empty() ?
        "config/finite_temperature.yaml" : g_test_config_path;

    std::cout << "\n═══════════════════════════════════════════════" << std::endl;
    std::cout << "F3: FINITE TEMPERATURE EFFECTS" << std::endl;
    std::cout << "═══════════════════════════════════════════════" << std::endl;
    std::cout << "Config: " << config_path << std::endl;

    // Load configuration
    YAML::Node config = YAML::LoadFile(config_path);

    // Grid parameters
    uint32_t grid_size = config["physics"]["grid_size"].as<uint32_t>();
    float dx = config["physics"]["spatial_step"].as<float>();
    float dt = config["physics"]["time_step"].as<float>();
    float coupling = config["physics"]["coupling_strength"].as<float>();

    // Thermal bath parameters
    auto thermal_config = config["physics"]["thermal_bath"];
    std::vector<float> temperatures;
    for (auto T_node : thermal_config["temperatures"]) {
        temperatures.push_back(T_node.as<float>());
    }
    float damping = thermal_config["damping"].as<float>();

    // Evolution parameters
    uint32_t equilibration_steps =
        config["simulation"]["equilibration_steps"].as<uint32_t>();
    uint32_t seed = config["simulation"]["random_seed"].as<uint32_t>();

    // Quality gates
    auto phase_transition = config["physics"]["phase_transitions"][0];
    float T_c_expected =
        phase_transition["critical_temperature_expected"].as<float>();
    float tolerance = 0.20f;  // 20% tolerance

    std::cout << "\nPhysics Parameters:" << std::endl;
    std::cout << "  Grid: " << grid_size << "³" << std::endl;
    std::cout << "  Coupling K = " << coupling << std::endl;
    std::cout << "  Temperature range: [" << temperatures.front()
             << ", " << temperatures.back() << "]" << std::endl;
    std::cout << "  Damping γ = " << damping << std::endl;
    std::cout << "  Expected T_c ≈ " << T_c_expected << std::endl;

    // Initialize TRD core
    TRDCore3D::Config trd_config;
    trd_config.Nx = grid_size;
    trd_config.Ny = grid_size;
    trd_config.Nz = grid_size;
    trd_config.dx = dx;
    trd_config.dt = dt;
    trd_config.coupling_strength = coupling;
    trd_config.mode = TRDCore3D::IntegrationMode::SYMPLECTIC;

    TRDCore3D core;
    core.initialize(trd_config);

    // Thermal bath configuration
    ThermalBathConfig bath;
    bath.damping = damping;
    bath.fluctuation_dissipation = true;

    // Compute phase diagram
    auto phase_diagram = PhaseDiagramCalculator::computePhaseDiagram(
        temperatures, core, bath, dt, equilibration_steps, seed);

    // Find critical temperature
    float T_c_measured = PhaseDiagramCalculator::findCriticalTemperature(
        phase_diagram, 0.5f);

    // Calculate transition sharpness
    float R_high_T = phase_diagram.back().order_parameter;  // Highest T
    float R_low_T = phase_diagram.front().order_parameter;   // Lowest T
    float transition_sharpness = (R_low_T - R_high_T) / R_low_T;

    // Write phase diagram to CSV
    std::string output_dir = config["simulation"]["output_directory"].as<std::string>();
    std::string csv_path = output_dir + "/phase_diagram.csv";
    std::ofstream csv_file(csv_path);
    csv_file << "# Phase diagram: Temperature vs Order Parameter R\n";
    csv_file << "# Golden Key: 1 TRD unit = 246 GeV\n";
    csv_file << "Temperature,OrderParameter\n";
    for (const auto& point : phase_diagram) {
        csv_file << point.temperature << "," << point.order_parameter << "\n";
    }
    csv_file.close();

    // Results
    std::cout << "\n═══════════════════════════════════════════════" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "═══════════════════════════════════════════════" << std::endl;

    std::cout << "\nPhase Diagram:" << std::endl;
    std::cout << "  T_low = " << temperatures.front()
             << " → R = " << std::fixed << std::setprecision(4)
             << R_low_T << " (synchronized)" << std::endl;
    std::cout << "  T_high = " << temperatures.back()
             << " → R = " << R_high_T << " (desynchronized)" << std::endl;

    std::cout << "\nCritical Temperature:" << std::endl;
    if (T_c_measured > 0) {
        std::cout << "  T_c (measured) = " << std::setprecision(3)
                 << T_c_measured << std::endl;
        std::cout << "  T_c (expected) = " << T_c_expected << std::endl;
        float error = fabs(T_c_measured - T_c_expected) / T_c_expected;
        std::cout << "  Relative error = " << std::setprecision(2)
                 << (error * 100.0f) << "%" << std::endl;
    } else {
        std::cout << "  WARNING: No clear transition found!" << std::endl;
    }

    std::cout << "\nTransition Sharpness:" << std::endl;
    std::cout << "  ΔR/R_low = " << std::setprecision(2)
             << (transition_sharpness * 100.0f) << "%" << std::endl;

    std::cout << "\nFluctuation-Dissipation Theorem:" << std::endl;
    std::cout << "  <noise²> = 2γkT" << std::endl;
    std::cout << "  Example: T = " << temperatures[temperatures.size()/2]
             << " → σ² = " << std::setprecision(4)
             << (2.0f * damping * BOLTZMANN * temperatures[temperatures.size()/2])
             << std::endl;

    // Quality gates
    std::cout << "\n═══════════════════════════════════════════════" << std::endl;
    std::cout << "QUALITY GATES" << std::endl;
    std::cout << "═══════════════════════════════════════════════" << std::endl;

    bool all_passed = true;

    // Gate 1: Critical temperature accuracy
    bool gate1 = false;
    if (T_c_measured > 0) {
        float error = fabs(T_c_measured - T_c_expected) / T_c_expected;
        gate1 = (error < tolerance);
        std::cout << "\n[Gate 1] Critical Temperature Accuracy" << std::endl;
        std::cout << "  Requirement: |T_c - T_c_expected|/T_c_expected < "
                 << (tolerance * 100) << "%" << std::endl;
        std::cout << "  Measured: " << std::setprecision(2) << (error * 100)
                 << "%" << std::endl;
        std::cout << "  Status: " << (gate1 ? "✅ PASS" : "❌ FAIL") << std::endl;
    } else {
        std::cout << "\n[Gate 1] Critical Temperature Accuracy" << std::endl;
        std::cout << "  Status: ❌ FAIL (no transition detected)" << std::endl;
    }
    all_passed &= gate1;

    // Gate 2: Transition sharpness
    bool gate2 = (transition_sharpness > 0.5f);
    std::cout << "\n[Gate 2] Phase Transition Sharpness" << std::endl;
    std::cout << "  Requirement: ΔR/R_low > 50%" << std::endl;
    std::cout << "  Measured: " << std::setprecision(2)
             << (transition_sharpness * 100.0f) << "%" << std::endl;
    std::cout << "  Status: " << (gate2 ? "✅ PASS" : "❌ FAIL") << std::endl;
    all_passed &= gate2;

    // Gate 3: Ordered state at low temperature
    bool gate3 = (R_low_T > 0.8f);
    std::cout << "\n[Gate 3] Ordered State (Low T)" << std::endl;
    std::cout << "  Requirement: R(T_low) > 0.8" << std::endl;
    std::cout << "  Measured: " << std::setprecision(4) << R_low_T << std::endl;
    std::cout << "  Status: " << (gate3 ? "✅ PASS" : "❌ FAIL") << std::endl;
    all_passed &= gate3;

    // Gate 4: Disordered state at high temperature
    bool gate4 = (R_high_T < 0.3f);
    std::cout << "\n[Gate 4] Disordered State (High T)" << std::endl;
    std::cout << "  Requirement: R(T_high) < 0.3" << std::endl;
    std::cout << "  Measured: " << std::setprecision(4) << R_high_T << std::endl;
    std::cout << "  Status: " << (gate4 ? "✅ PASS" : "❌ FAIL") << std::endl;
    all_passed &= gate4;

    // Final verdict
    std::cout << "\n═══════════════════════════════════════════════" << std::endl;
    if (all_passed) {
        std::cout << "✅ ALL QUALITY GATES PASSED" << std::endl;
        std::cout << "\nTRD successfully reproduces thermal phase transitions!" << std::endl;
        std::cout << "Phase diagram saved to: " << csv_path << std::endl;
    } else {
        std::cout << "❌ SOME QUALITY GATES FAILED" << std::endl;
        std::cout << "\nReview phase diagram for analysis: " << csv_path << std::endl;
    }
    std::cout << "═══════════════════════════════════════════════\n" << std::endl;

    return all_passed ? 0 : 1;
}
