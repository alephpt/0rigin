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
        // Physical length scales (Planck units: ℏ = c = G = 1, Δ = m_Planck = 1)
        float L_domain = 100.0f;  // Domain size in Planck lengths (default: 100 ℓ_P)

        // Grid convergence testing: array of grid sizes to test
        std::vector<int> grid_sizes;  // e.g., [64, 128, 256] - tests all if non-empty
    };

    // Physics parameters
    struct PhysicsConfig {
        std::string solver_type = "gpu";  // Solver: "gpu" (default), "dirac" (CPU), or "klein_gordon"
        float delta = 2.5f;           // Mass gap parameter
        float coupling = 0.1f;         // Kuramoto-Dirac coupling
        float dt = 0.01f;              // Timestep
        int total_steps = 100;         // Number of steps to run
        float K = 1.0f;                // Kuramoto coupling strength
        float damping = 0.1f;          // Phase damping
        float noise_strength = 0.0f;   // Langevin noise amplitude (σ, 0 = deterministic)

        // Phase transition scan parameters (Test 3.3)
        std::vector<float> noise_scan;  // Array of noise values to scan (e.g., [0.0, 0.05, 0.1, ...])

        // Electromagnetic coupling parameters (Phase 5)
        bool em_coupling_enabled = false;     // Enable EM coupling from Kuramoto phase
        float em_coupling_strength = 1.0f;    // Effective charge q
        std::string em_coupling_type = "perturbative";  // "perturbative" or "full_gauge"
        std::string em_regularization = "none";  // "none", "R", "R2" - regularization for A = ∇θ prescription
    };

    // Dirac initial condition
    struct DiracInitialCondition {
        std::string type = "gaussian"; // Type: "gaussian", "boosted_gaussian", "plane_wave", "linear_defects", "domain_split", "two_particle"
        float x0 = 32.0f;              // Center x (DEPRECATED: use x0_physical)
        float y0 = 32.0f;              // Center y (DEPRECATED: use y0_physical)
        float sigma = 3.0f;            // Width parameter (DEPRECATED: use sigma_physical)
        float amplitude = 1.0f;        // Initial amplitude

        // Grid-independent physical parameters (Planck lengths)
        float x0_physical = -1.0f;     // Center x in physical units (-1 = use x0 grid units for backward compatibility)
        float y0_physical = -1.0f;     // Center y in physical units (-1 = use y0 grid units for backward compatibility)
        float sigma_physical = -1.0f;  // Width in physical units (-1 = use sigma grid units for backward compatibility)

        // Relativistic boost parameters (Scenario 2.3)
        std::vector<float> boost_velocities;  // Array of velocities to test (c = 1 in Planck units)
        float boost_vx = 0.0f;                // Single velocity boost in x-direction
        float boost_vy = 0.0f;                // Single velocity boost in y-direction

        // Phase 3 IC parameters (Test 3.1: Casimir Force)
        float defect_separation = 10.0f;      // Distance between linear defects [ℓ_P]
        float defect_width = 1.0f;            // Transition width of defect [ℓ_P]
        int winding_number = 1;               // Winding orientation for defects

        // Phase 3 IC parameters (Test 3.2: Vacuum Energy)
        float split_x = 50.0f;                // x-coordinate of domain boundary [ℓ_P]
        float R_left = 1.0f;                  // Order parameter left side
        float R_right = 0.0f;                 // Order parameter right side
        float transition_width = 1.0f;        // Width of domain boundary transition [ℓ_P]

        // Phase 3 IC parameters (Test 3.4: Two-particle mode)
        bool two_particle_mode = false;       // Enable dual particle evolution
        struct ParticleConfig {
            float x0 = 40.0f;
            float y0 = 50.0f;
            float sigma = 3.0f;
            float boost_vx = 0.0f;
            float boost_vy = 0.0f;
        };
        ParticleConfig particle_1;            // First particle config
        ParticleConfig particle_2;            // Second particle config
    };

    // Kuramoto initial condition
    struct KuramotoInitialCondition {
        std::string phase_distribution = "uniform";  // "uniform", "random", "vortex", "phase_gradient", "multi_vortex"
        std::string omega_distribution = "zero";     // "zero", "gaussian", "random"
        float omega_mean = 0.0f;
        float omega_std = 0.1f;
        // Wave vector for phase_gradient initialization
        float wave_vector_x = 0.0f;
        float wave_vector_y = 0.0f;

        // Single vortex configuration (grid-independent physical parameters)
        int winding_number = 1;              // Topological charge
        float vortex_core_radius = 3.0f;     // Core radius in Planck lengths (grid-independent)
        float vortex_center_x = 50.0f;       // Center x in physical units (default: domain center)
        float vortex_center_y = 50.0f;       // Center y in physical units (default: domain center)

        // Vortex pair configuration (Sprint 2: Multi-Defect Interactions)
        float vortex_separation = 10.0f;     // Separation between vortex pair in Planck lengths

        // Multi-vortex configuration (Phase 5/6: EM tests)
        struct VortexConfig {
            float center_x = 50.0f;          // Center x in Planck lengths
            float center_y = 50.0f;          // Center y in Planck lengths
            int winding = 1;                 // Topological charge (+1, -1, etc.)
            float core_radius = 3.0f;        // Core radius in Planck lengths
        };
        std::vector<VortexConfig> vortices;  // Array of vortices for multi_vortex type
    };

    // Test particle configuration (Lorentz force validation)
    struct TestParticleConfig {
        bool enabled = false;
        double x0 = 60.0;          // Initial position x [Planck lengths]
        double y0 = 50.0;          // Initial position y [Planck lengths]
        double vx0 = 0.0;          // Initial velocity vx [c = 1]
        double vy0 = 0.1;          // Initial velocity vy [c = 1]
        double charge = -1.0;      // Particle charge [natural units]
        double mass = 1.0;         // Particle mass [natural units]
        int record_every = 10;     // Record trajectory every N steps
        bool use_uniform_B = false; // Use uniform B-field instead of from phase
        double uniform_B_z = 0.1;  // Uniform B-field strength
    };

    // Operator splitting configuration
    struct OperatorSplittingConfig {
        bool enabled = false;
        std::vector<int> substep_ratios = {1, 10, 100};
    };

    // Validation tolerances
    struct ValidationConfig {
        // Global validation (always enforced)
        float norm_tolerance = 5.0e-3f;      // ||Ψ||² - 1 tolerance (0.5%)
        float energy_tolerance = 1.0e-2f;    // |ΔE/E₀| tolerance (1%)
        float convergence_tolerance = 0.05f; // Convergence check (5%)
        bool enforce_R_bounds = true;        // Check 0 ≤ R ≤ 1
        bool enforce_causality = true;       // Check v ≤ c
        bool check_numerical_stability = true; // Check for NaN/Inf

        // Scenario-specific validation
        std::string scenario = "none";       // "none", "defect_localization", "traveling_wave", "relativistic_mass"

        // Scenario: traveling_wave and relativistic_mass
        bool require_vortex = false;         // Require W = ±1
        float winding_tolerance = 0.2f;      // |W - 1| < tolerance
        bool require_core = false;           // Require R_min < 0.5
        float core_R_threshold = 0.5f;       // R_min threshold
        bool require_boost = false;          // Require initial momentum
        float initial_momentum_tolerance = 0.05f; // 5%
        bool validate_gamma_factor = false;  // Validate γ_measured vs theory
        float gamma_tolerance = 0.05f;       // 5%

        // EM validation (when em_coupling enabled)
        bool validate_maxwell_equations = false;  // Verify ∇·E, ∇×B, etc.
        float maxwell_tolerance = 1.0e-6f;        // Maxwell equation residual tolerance
        bool validate_flux_quantization = false;  // Verify ∮A·dl = 2πW
        float flux_quantization_tolerance = 0.1f; // 10% tolerance for flux quantization
        int expected_winding_number = 1;          // Expected vortex winding W
        bool validate_gauge_invariance = false;   // Test θ → θ + α invariance
        float gauge_shift_angle = 0.785398f;      // Gauge shift angle (π/4 rad)
        float gauge_invariance_tolerance = 1.0e-10f; // Field strength difference tolerance

        // Validation timing
        bool validate_initial_state = true;
        bool validate_during_evolution = false;
        int validation_interval = 100;       // Steps between runtime checks
        bool validate_final_state = true;
        bool fail_on_critical = true;        // Abort simulation on critical failure
        bool verbose = true;                 // Print detailed validation output
    };

    // Output configuration
    struct OutputConfig {
        std::string directory = "output/test";
        int save_every = 10;                      // Save observables every N steps
        std::vector<std::string> formats = {"csv"}; // Output formats
        bool auto_visualize = false;               // Auto-generate plots after test
        bool save_spatial_snapshots = false;       // Save full spatial fields (theta, R) at snapshots
        std::vector<int> snapshot_steps = {};      // Timesteps at which to save spatial snapshots (empty = auto)

        // Plotting configuration (NEW)
        bool enable_plots = false;                 // Enable C++ plotting via matplotlib-cpp
        int plot_dpi = 300;                        // DPI for saved figures
        std::vector<std::string> plot_formats = {"png"}; // Plot output formats (png, pdf, svg)
        bool plot_6panel = true;                   // Generate 6-panel observable plot
        bool plot_conservation = true;             // Generate conservation law plots
        bool plot_spatial_fields = true;           // Generate spatial field plots
        bool plot_gamma_validation = false;        // Generate gamma factor validation plot (Scenario 2.3)
    };

    // Analysis configuration
    struct AnalysisConfig {
        std::string mode = "simulation"; // "simulation", "dispersion", "stability"
        bool track_trajectory = false;
        bool measure_velocity = false;
        bool compute_effective_mass = false;
        bool test_gamma_factor = false;
        bool compute_energy = false;
        bool lorentz_invariants = false;

        // Boosted frame analysis (Phase 2.5B)
        bool perform_lorentz_transform = false;       // Transform to boosted frame
        bool measure_R_in_boosted_frame = false;      // Measure R'(x',t') in particle rest frame
        bool test_lorentz_covariance = false;         // Test if m' = m in all frames
        bool compute_frame_invariants = false;        // Verify E'² - p'² = E² - p²

        // Dispersion analysis parameters
        int dispersion_k_steps = 50;  // Number of k points to scan
        float dispersion_max_k = 3.14159f; // Max k (pi/dx)
    };

    // Full test configuration
    std::string test_name;
    std::string description;

    GridConfig grid;
    PhysicsConfig physics;
    DiracInitialCondition dirac_initial;
    KuramotoInitialCondition kuramoto_initial;
    TestParticleConfig test_particle;
    OperatorSplittingConfig operator_splitting;
    ValidationConfig validation;
    OutputConfig output;
    AnalysisConfig analysis;

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
    void parseAnalysis(const YAML::Node& node);
    void parseEMCoupling(const YAML::Node& node);  // Parse em_coupling section (Phase 5/6)
    void parseTestParticle(const YAML::Node& node);  // Parse test_particle section (Lorentz force)
};

#endif // TEST_CONFIG_H
