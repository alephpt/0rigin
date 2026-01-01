#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

/**
 * TestConfig - YAML-driven test configuration for unified SMFT testing
 *
 * Supports loading test parameters, initial conditions, validation tolerances,
 * and output specifications from YAML files.
 *
 * Example YAML:
 * ---
 * test_name: "timesync_validation"
 * description: "Validate operator splitting with different N ratios"
 *
 * grid:
 *   size_x: 64
 *   size_y: 64
 *
 * physics:
 *   delta: 2.5
 *   coupling: 0.1
 *   dt: 0.01
 *   total_steps: 100
 *
 * initial_conditions:
 *   dirac:
 *     type: "gaussian"
 *     x0: 32.0
 *     y0: 32.0
 *     sigma: 3.0
 *   kuramoto:
 *     phase_distribution: "uniform"
 *     omega_distribution: "zero"
 *
 * operator_splitting:
 *   enabled: true
 *   substep_ratios: [1, 10, 100]
 *
 * validation:
 *   norm_tolerance: 1.0e-4
 *   energy_tolerance: 1.0e-2
 *   convergence_tolerance: 0.05
 *
 * output:
 *   directory: "output/timesync_validation"
 *   save_every: 10
 *   formats: ["csv", "binary"]
 */
class TestConfig {
public:
    // Grid configuration
    struct GridConfig {
        int size_x = 64;
        int size_y = 64;
    };

    // Physics parameters
    struct PhysicsConfig {
        float delta = 2.5f;           // Mass gap parameter
        float coupling = 0.1f;         // Kuramoto-Dirac coupling
        float dt = 0.01f;              // Timestep
        int total_steps = 100;         // Number of steps to run
        float K = 1.0f;                // Kuramoto coupling strength
        float damping = 0.1f;          // Phase damping
    };

    // Dirac initial condition
    struct DiracInitialCondition {
        std::string type = "gaussian"; // Type: "gaussian", "plane_wave", etc.
        float x0 = 32.0f;              // Center x
        float y0 = 32.0f;              // Center y
        float sigma = 3.0f;            // Width parameter
        float amplitude = 1.0f;        // Initial amplitude
    };

    // Kuramoto initial condition
    struct KuramotoInitialCondition {
        std::string phase_distribution = "uniform";  // "uniform", "random", "vortex", "phase_gradient"
        std::string omega_distribution = "zero";     // "zero", "gaussian", "random"
        float omega_mean = 0.0f;
        float omega_std = 0.1f;
        // Wave vector for phase_gradient initialization
        float wave_vector_x = 0.0f;
        float wave_vector_y = 0.0f;
    };

    // Operator splitting configuration
    struct OperatorSplittingConfig {
        bool enabled = false;
        std::vector<int> substep_ratios = {1, 10, 100};
    };

    // Validation tolerances
    struct ValidationConfig {
        float norm_tolerance = 1.0e-4f;      // ||Ψ||² - 1 tolerance
        float energy_tolerance = 1.0e-2f;    // |ΔE/E₀| tolerance
        float convergence_tolerance = 0.05f; // Convergence check (5%)
    };

    // Output configuration
    struct OutputConfig {
        std::string directory = "output/test";
        int save_every = 10;                      // Save observables every N steps
        std::vector<std::string> formats = {"csv"}; // Output formats
        bool auto_visualize = false;               // Auto-generate plots after test
    };

    // Full test configuration
    std::string test_name;
    std::string description;

    GridConfig grid;
    PhysicsConfig physics;
    DiracInitialCondition dirac_initial;
    KuramotoInitialCondition kuramoto_initial;
    OperatorSplittingConfig operator_splitting;
    ValidationConfig validation;
    OutputConfig output;

    /**
     * Load configuration from YAML file
     * @param yaml_path Path to YAML configuration file
     * @return true if successful, false otherwise
     */
    bool loadFromYAML(const std::string& yaml_path);

    /**
     * Save configuration to YAML file
     * @param yaml_path Path to output YAML file
     * @return true if successful, false otherwise
     */
    bool saveToYAML(const std::string& yaml_path) const;

    /**
     * Create a default configuration with standard parameters
     * @param test_name Name for the test
     * @return Default test configuration
     */
    static TestConfig createDefault(const std::string& test_name = "default_test");

    /**
     * Validate configuration parameters
     * @return true if valid, false otherwise
     */
    bool validate() const;

    /**
     * Print configuration to stdout
     */
    void print() const;

private:
    // Helper methods for YAML parsing
    void parseGrid(const YAML::Node& node);
    void parsePhysics(const YAML::Node& node);
    void parseDiracInitial(const YAML::Node& node);
    void parseKuramotoInitial(const YAML::Node& node);
    void parseOperatorSplitting(const YAML::Node& node);
    void parseValidation(const YAML::Node& node);
    void parseOutput(const YAML::Node& node);
};

#endif // TEST_CONFIG_H
