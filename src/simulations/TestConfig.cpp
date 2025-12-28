#include "TestConfig.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

bool TestConfig::loadFromYAML(const std::string& yaml_path) {
    try {
        YAML::Node config = YAML::LoadFile(yaml_path);

        // Load basic info
        if (config["test_name"]) {
            test_name = config["test_name"].as<std::string>();
        }
        if (config["description"]) {
            description = config["description"].as<std::string>();
        }

        // Load sections
        if (config["grid"]) parseGrid(config["grid"]);
        if (config["physics"]) parsePhysics(config["physics"]);
        if (config["initial_conditions"]) {
            if (config["initial_conditions"]["dirac"]) {
                parseDiracInitial(config["initial_conditions"]["dirac"]);
            }
            if (config["initial_conditions"]["kuramoto"]) {
                parseKuramotoInitial(config["initial_conditions"]["kuramoto"]);
            }
        }
        if (config["operator_splitting"]) parseOperatorSplitting(config["operator_splitting"]);
        if (config["validation"]) parseValidation(config["validation"]);
        if (config["output"]) parseOutput(config["output"]);
        if (config["analysis"]) parseAnalysis(config["analysis"]);

        // Parse EM coupling section (Phase 5/6) - can override physics section
        if (config["em_coupling"]) parseEMCoupling(config["em_coupling"]);

        return validate();

    } catch (const YAML::Exception& e) {
        std::cerr << "YAML parsing error: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << std::endl;
        return false;
    }
}

bool TestConfig::saveToYAML(const std::string& yaml_path) const {
    try {
        YAML::Emitter out;
        out << YAML::BeginMap;

        // Basic info
        out << YAML::Key << "test_name" << YAML::Value << test_name;
        out << YAML::Key << "description" << YAML::Value << description;

        // Grid
        out << YAML::Key << "grid" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "size_x" << YAML::Value << grid.size_x;
        out << YAML::Key << "size_y" << YAML::Value << grid.size_y;
        out << YAML::EndMap;

        // Physics
        out << YAML::Key << "physics" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "delta" << YAML::Value << physics.delta;
        out << YAML::Key << "coupling" << YAML::Value << physics.coupling;
        out << YAML::Key << "dt" << YAML::Value << physics.dt;
        out << YAML::Key << "total_steps" << YAML::Value << physics.total_steps;
        out << YAML::Key << "K" << YAML::Value << physics.K;
        out << YAML::Key << "damping" << YAML::Value << physics.damping;
        out << YAML::EndMap;

        // Initial conditions
        out << YAML::Key << "initial_conditions" << YAML::Value << YAML::BeginMap;

        out << YAML::Key << "dirac" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "type" << YAML::Value << dirac_initial.type;
        out << YAML::Key << "x0" << YAML::Value << dirac_initial.x0;
        out << YAML::Key << "y0" << YAML::Value << dirac_initial.y0;
        out << YAML::Key << "sigma" << YAML::Value << dirac_initial.sigma;
        out << YAML::Key << "amplitude" << YAML::Value << dirac_initial.amplitude;
        out << YAML::EndMap;

        out << YAML::Key << "kuramoto" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "phase_distribution" << YAML::Value << kuramoto_initial.phase_distribution;
        out << YAML::Key << "omega_distribution" << YAML::Value << kuramoto_initial.omega_distribution;
        out << YAML::Key << "omega_mean" << YAML::Value << kuramoto_initial.omega_mean;
        out << YAML::Key << "omega_std" << YAML::Value << kuramoto_initial.omega_std;
        out << YAML::EndMap;

        out << YAML::EndMap; // initial_conditions

        // Operator splitting
        out << YAML::Key << "operator_splitting" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "enabled" << YAML::Value << operator_splitting.enabled;
        out << YAML::Key << "substep_ratios" << YAML::Value << YAML::Flow << operator_splitting.substep_ratios;
        out << YAML::EndMap;

        // Validation
        out << YAML::Key << "validation" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "norm_tolerance" << YAML::Value << validation.norm_tolerance;
        out << YAML::Key << "energy_tolerance" << YAML::Value << validation.energy_tolerance;
        out << YAML::Key << "convergence_tolerance" << YAML::Value << validation.convergence_tolerance;
        out << YAML::EndMap;

        // Output
        out << YAML::Key << "output" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "directory" << YAML::Value << output.directory;
        out << YAML::Key << "save_every" << YAML::Value << output.save_every;
        out << YAML::Key << "formats" << YAML::Value << YAML::Flow << output.formats;
        out << YAML::EndMap;

        // Analysis
        out << YAML::Key << "analysis" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "mode" << YAML::Value << analysis.mode;
        out << YAML::Key << "track_trajectory" << YAML::Value << analysis.track_trajectory;
        out << YAML::Key << "measure_velocity" << YAML::Value << analysis.measure_velocity;
        out << YAML::Key << "compute_effective_mass" << YAML::Value << analysis.compute_effective_mass;
        out << YAML::Key << "test_gamma_factor" << YAML::Value << analysis.test_gamma_factor;
        out << YAML::Key << "compute_energy" << YAML::Value << analysis.compute_energy;
        out << YAML::Key << "lorentz_invariants" << YAML::Value << analysis.lorentz_invariants;
        out << YAML::Key << "dispersion_k_steps" << YAML::Value << analysis.dispersion_k_steps;
        out << YAML::Key << "dispersion_max_k" << YAML::Value << analysis.dispersion_max_k;
        out << YAML::EndMap;

        out << YAML::EndMap;

        // Write to file
        std::ofstream fout(yaml_path);
        if (!fout.is_open()) {
            std::cerr << "Failed to open file for writing: " << yaml_path << std::endl;
            return false;
        }
        fout << out.c_str();
        fout.close();

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error saving config: " << e.what() << std::endl;
        return false;
    }
}

TestConfig TestConfig::createDefault(const std::string& test_name) {
    TestConfig config;
    config.test_name = test_name;
    config.description = "Default SMFT test configuration";

    // Defaults already set in struct definitions
    return config;
}

bool TestConfig::validate() const {
    // Check grid
    if (grid.size_x <= 0 || grid.size_y <= 0) {
        std::cerr << "Invalid grid size" << std::endl;
        return false;
    }

    // Check physics parameters
    if (physics.dt <= 0.0f) {
        std::cerr << "Invalid timestep dt" << std::endl;
        return false;
    }
    if (physics.total_steps <= 0) {
        std::cerr << "Invalid total_steps" << std::endl;
        return false;
    }

    // Check Dirac initial condition
    // Skip grid-based validation if using physical coordinates
    if (dirac_initial.x0_physical < 0) {
        // Using deprecated grid coordinates - validate against grid size
        if (dirac_initial.x0 < 0 || dirac_initial.x0 >= grid.size_x) {
            std::cerr << "Dirac x0 out of bounds" << std::endl;
            return false;
        }
        if (dirac_initial.y0 < 0 || dirac_initial.y0 >= grid.size_y) {
            std::cerr << "Dirac y0 out of bounds" << std::endl;
            return false;
        }
        if (dirac_initial.sigma <= 0.0f) {
            std::cerr << "Invalid Dirac sigma" << std::endl;
            return false;
        }
    } else {
        // Using physical coordinates - validate against domain size
        if (dirac_initial.x0_physical < 0 || dirac_initial.x0_physical > grid.L_domain) {
            std::cerr << "Dirac x0_physical out of domain bounds" << std::endl;
            return false;
        }
        if (dirac_initial.y0_physical < 0 || dirac_initial.y0_physical > grid.L_domain) {
            std::cerr << "Dirac y0_physical out of domain bounds" << std::endl;
            return false;
        }
        if (dirac_initial.sigma_physical <= 0.0f) {
            std::cerr << "Invalid Dirac sigma_physical" << std::endl;
            return false;
        }
    }

    // Check operator splitting
    if (operator_splitting.enabled && operator_splitting.substep_ratios.empty()) {
        std::cerr << "Operator splitting enabled but no substep ratios specified" << std::endl;
        return false;
    }

    // Check validation tolerances
    if (validation.norm_tolerance <= 0.0f || validation.energy_tolerance <= 0.0f) {
        std::cerr << "Invalid validation tolerances" << std::endl;
        return false;
    }

    return true;
}

void TestConfig::print() const {
    std::cout << "\n===== Test Configuration: " << test_name << " =====" << std::endl;
    std::cout << "Description: " << description << std::endl;

    std::cout << "\n[Grid]" << std::endl;
    std::cout << "  size: " << grid.size_x << " x " << grid.size_y << std::endl;

    std::cout << "\n[Physics]" << std::endl;
    std::cout << "  solver_type: " << physics.solver_type << std::endl;
    std::cout << "  delta: " << physics.delta << std::endl;
    std::cout << "  coupling: " << physics.coupling << std::endl;
    std::cout << "  dt: " << physics.dt << std::endl;
    std::cout << "  total_steps: " << physics.total_steps << std::endl;
    std::cout << "  K (Kuramoto): " << physics.K << std::endl;
    std::cout << "  damping: " << physics.damping << std::endl;
    std::cout << "  noise_strength: " << physics.noise_strength << std::endl;
    if (!physics.noise_scan.empty()) {
        std::cout << "  noise_scan: [";
        for (size_t i = 0; i < physics.noise_scan.size(); ++i) {
            std::cout << physics.noise_scan[i];
            if (i < physics.noise_scan.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    // EM coupling parameters (Phase 5/6)
    std::cout << "  em_coupling_enabled: " << (physics.em_coupling_enabled ? "true" : "false") << std::endl;
    if (physics.em_coupling_enabled) {
        std::cout << "  em_coupling_strength: " << physics.em_coupling_strength << std::endl;
        std::cout << "  em_coupling_type: " << physics.em_coupling_type << std::endl;
    }

    std::cout << "\n[Initial Conditions - Dirac]" << std::endl;
    std::cout << "  type: " << dirac_initial.type << std::endl;
    std::cout << "  center: (" << dirac_initial.x0 << ", " << dirac_initial.y0 << ")" << std::endl;
    std::cout << "  sigma: " << dirac_initial.sigma << std::endl;
    std::cout << "  amplitude: " << dirac_initial.amplitude << std::endl;

    std::cout << "\n[Initial Conditions - Kuramoto]" << std::endl;
    std::cout << "  phase_distribution: " << kuramoto_initial.phase_distribution << std::endl;
    std::cout << "  omega_distribution: " << kuramoto_initial.omega_distribution << std::endl;
    std::cout << "  omega_mean: " << kuramoto_initial.omega_mean << std::endl;
    std::cout << "  omega_std: " << kuramoto_initial.omega_std << std::endl;

    std::cout << "\n[Operator Splitting]" << std::endl;
    std::cout << "  enabled: " << (operator_splitting.enabled ? "yes" : "no") << std::endl;
    if (operator_splitting.enabled) {
        std::cout << "  substep_ratios: [";
        for (size_t i = 0; i < operator_splitting.substep_ratios.size(); ++i) {
            std::cout << operator_splitting.substep_ratios[i];
            if (i < operator_splitting.substep_ratios.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    std::cout << "\n[Validation]" << std::endl;
    std::cout << "  norm_tolerance: " << validation.norm_tolerance << std::endl;
    std::cout << "  energy_tolerance: " << validation.energy_tolerance << std::endl;
    std::cout << "  convergence_tolerance: " << validation.convergence_tolerance << std::endl;

    std::cout << "\n[Output]" << std::endl;
    std::cout << "  directory: " << output.directory << std::endl;
    std::cout << "  save_every: " << output.save_every << " steps" << std::endl;
    std::cout << "  formats: [";
    for (size_t i = 0; i < output.formats.size(); ++i) {
        std::cout << output.formats[i];
        if (i < output.formats.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

// Private parsing methods

void TestConfig::parseGrid(const YAML::Node& node) {
    if (node["size_x"]) grid.size_x = node["size_x"].as<int>();
    if (node["size_y"]) grid.size_y = node["size_y"].as<int>();
    if (node["L_domain"]) grid.L_domain = node["L_domain"].as<float>();

    // Parse grid_sizes array for convergence testing
    if (node["grid_sizes"]) {
        grid.grid_sizes.clear();
        for (const auto& size : node["grid_sizes"]) {
            grid.grid_sizes.push_back(size.as<int>());
        }
    }
}

void TestConfig::parsePhysics(const YAML::Node& node) {
    if (node["solver_type"]) physics.solver_type = node["solver_type"].as<std::string>();
    if (node["delta"]) physics.delta = node["delta"].as<float>();
    if (node["coupling"]) physics.coupling = node["coupling"].as<float>();
    if (node["dt"]) physics.dt = node["dt"].as<float>();
    if (node["total_steps"]) physics.total_steps = node["total_steps"].as<int>();
    if (node["K"]) physics.K = node["K"].as<float>();
    if (node["damping"]) physics.damping = node["damping"].as<float>();
    if (node["noise_strength"]) physics.noise_strength = node["noise_strength"].as<float>();

    // Parse noise_scan array for phase transition testing (Test 3.3)
    if (node["noise_scan"]) {
        physics.noise_scan.clear();
        for (const auto& sigma : node["noise_scan"]) {
            physics.noise_scan.push_back(sigma.as<float>());
        }
    }

    // Parse electromagnetic coupling parameters (Phase 5)
    if (node["em_coupling_enabled"]) physics.em_coupling_enabled = node["em_coupling_enabled"].as<bool>();
    if (node["em_coupling_strength"]) physics.em_coupling_strength = node["em_coupling_strength"].as<float>();
    if (node["em_coupling_type"]) physics.em_coupling_type = node["em_coupling_type"].as<std::string>();
}

void TestConfig::parseDiracInitial(const YAML::Node& node) {
    if (node["type"]) dirac_initial.type = node["type"].as<std::string>();
    if (node["x0"]) dirac_initial.x0 = node["x0"].as<float>();
    if (node["y0"]) dirac_initial.y0 = node["y0"].as<float>();
    if (node["sigma"]) dirac_initial.sigma = node["sigma"].as<float>();
    if (node["amplitude"]) dirac_initial.amplitude = node["amplitude"].as<float>();

    // Grid-independent physical parameters
    if (node["x0_physical"]) dirac_initial.x0_physical = node["x0_physical"].as<float>();
    if (node["y0_physical"]) dirac_initial.y0_physical = node["y0_physical"].as<float>();
    if (node["sigma_physical"]) dirac_initial.sigma_physical = node["sigma_physical"].as<float>();

    // Relativistic boost parameters (Scenario 2.3)
    if (node["boost_velocities"]) {
        dirac_initial.boost_velocities.clear();
        for (const auto& v : node["boost_velocities"]) {
            dirac_initial.boost_velocities.push_back(v.as<float>());
        }
    }
    if (node["boost_vx"]) dirac_initial.boost_vx = node["boost_vx"].as<float>();
    if (node["boost_vy"]) dirac_initial.boost_vy = node["boost_vy"].as<float>();

    // Phase 3 IC parameters (Test 3.1: Casimir Force)
    if (node["defect_separation"]) dirac_initial.defect_separation = node["defect_separation"].as<float>();
    if (node["defect_width"]) dirac_initial.defect_width = node["defect_width"].as<float>();
    if (node["winding_number"]) dirac_initial.winding_number = node["winding_number"].as<int>();

    // Phase 3 IC parameters (Test 3.2: Vacuum Energy)
    if (node["split_x"]) dirac_initial.split_x = node["split_x"].as<float>();
    if (node["R_left"]) dirac_initial.R_left = node["R_left"].as<float>();
    if (node["R_right"]) dirac_initial.R_right = node["R_right"].as<float>();
    if (node["transition_width"]) dirac_initial.transition_width = node["transition_width"].as<float>();

    // Phase 3 IC parameters (Test 3.4: Two-particle mode)
    if (node["two_particle_mode"]) dirac_initial.two_particle_mode = node["two_particle_mode"].as<bool>();
    if (node["particle_1"]) {
        const auto& p1 = node["particle_1"];
        if (p1["x0"]) dirac_initial.particle_1.x0 = p1["x0"].as<float>();
        if (p1["y0"]) dirac_initial.particle_1.y0 = p1["y0"].as<float>();
        if (p1["sigma"]) dirac_initial.particle_1.sigma = p1["sigma"].as<float>();
        if (p1["boost_vx"]) dirac_initial.particle_1.boost_vx = p1["boost_vx"].as<float>();
        if (p1["boost_vy"]) dirac_initial.particle_1.boost_vy = p1["boost_vy"].as<float>();
    }
    if (node["particle_2"]) {
        const auto& p2 = node["particle_2"];
        if (p2["x0"]) dirac_initial.particle_2.x0 = p2["x0"].as<float>();
        if (p2["y0"]) dirac_initial.particle_2.y0 = p2["y0"].as<float>();
        if (p2["sigma"]) dirac_initial.particle_2.sigma = p2["sigma"].as<float>();
        if (p2["boost_vx"]) dirac_initial.particle_2.boost_vx = p2["boost_vx"].as<float>();
        if (p2["boost_vy"]) dirac_initial.particle_2.boost_vy = p2["boost_vy"].as<float>();
    }
}

void TestConfig::parseKuramotoInitial(const YAML::Node& node) {
    if (node["phase_distribution"]) {
        kuramoto_initial.phase_distribution = node["phase_distribution"].as<std::string>();
    }
    if (node["omega_distribution"]) {
        kuramoto_initial.omega_distribution = node["omega_distribution"].as<std::string>();
    }
    if (node["omega_mean"]) kuramoto_initial.omega_mean = node["omega_mean"].as<float>();
    if (node["omega_std"]) kuramoto_initial.omega_std = node["omega_std"].as<float>();
    if (node["wave_vector_x"]) kuramoto_initial.wave_vector_x = node["wave_vector_x"].as<float>();
    if (node["wave_vector_y"]) kuramoto_initial.wave_vector_y = node["wave_vector_y"].as<float>();

    // Vortex configuration (grid-independent)
    if (node["winding_number"]) kuramoto_initial.winding_number = node["winding_number"].as<int>();
    if (node["vortex_core_radius"]) kuramoto_initial.vortex_core_radius = node["vortex_core_radius"].as<float>();
    if (node["vortex_center_x"]) kuramoto_initial.vortex_center_x = node["vortex_center_x"].as<float>();
    if (node["vortex_center_y"]) kuramoto_initial.vortex_center_y = node["vortex_center_y"].as<float>();
}

void TestConfig::parseOperatorSplitting(const YAML::Node& node) {
    if (node["enabled"]) operator_splitting.enabled = node["enabled"].as<bool>();
    if (node["substep_ratios"]) {
        operator_splitting.substep_ratios = node["substep_ratios"].as<std::vector<int>>();
    }
}

void TestConfig::parseValidation(const YAML::Node& node) {
    // Global validation
    if (node["norm_tolerance"]) validation.norm_tolerance = node["norm_tolerance"].as<float>();
    if (node["energy_tolerance"]) validation.energy_tolerance = node["energy_tolerance"].as<float>();
    if (node["convergence_tolerance"]) validation.convergence_tolerance = node["convergence_tolerance"].as<float>();
    if (node["enforce_R_bounds"]) validation.enforce_R_bounds = node["enforce_R_bounds"].as<bool>();
    if (node["enforce_causality"]) validation.enforce_causality = node["enforce_causality"].as<bool>();
    if (node["check_numerical_stability"]) validation.check_numerical_stability = node["check_numerical_stability"].as<bool>();

    // Scenario-specific validation
    if (node["scenario"]) validation.scenario = node["scenario"].as<std::string>();
    if (node["require_vortex"]) validation.require_vortex = node["require_vortex"].as<bool>();
    if (node["winding_tolerance"]) validation.winding_tolerance = node["winding_tolerance"].as<float>();
    if (node["require_core"]) validation.require_core = node["require_core"].as<bool>();
    if (node["core_R_threshold"]) validation.core_R_threshold = node["core_R_threshold"].as<float>();
    if (node["require_boost"]) validation.require_boost = node["require_boost"].as<bool>();
    if (node["initial_momentum_tolerance"]) validation.initial_momentum_tolerance = node["initial_momentum_tolerance"].as<float>();
    if (node["validate_gamma_factor"]) validation.validate_gamma_factor = node["validate_gamma_factor"].as<bool>();
    if (node["gamma_tolerance"]) validation.gamma_tolerance = node["gamma_tolerance"].as<float>();

    // Validation timing
    if (node["validate_initial_state"]) validation.validate_initial_state = node["validate_initial_state"].as<bool>();
    if (node["validate_during_evolution"]) validation.validate_during_evolution = node["validate_during_evolution"].as<bool>();
    if (node["validation_interval"]) validation.validation_interval = node["validation_interval"].as<int>();
    if (node["validate_final_state"]) validation.validate_final_state = node["validate_final_state"].as<bool>();
    if (node["fail_on_critical"]) validation.fail_on_critical = node["fail_on_critical"].as<bool>();
    if (node["verbose"]) validation.verbose = node["verbose"].as<bool>();
}

void TestConfig::parseOutput(const YAML::Node& node) {
    if (node["directory"]) output.directory = node["directory"].as<std::string>();
    if (node["save_every"]) output.save_every = node["save_every"].as<int>();
    if (node["formats"]) output.formats = node["formats"].as<std::vector<std::string>>();
    if (node["auto_visualize"]) output.auto_visualize = node["auto_visualize"].as<bool>();
    if (node["save_spatial_snapshots"]) output.save_spatial_snapshots = node["save_spatial_snapshots"].as<bool>();
    if (node["snapshot_steps"]) output.snapshot_steps = node["snapshot_steps"].as<std::vector<int>>();

    // Plotting configuration
    if (node["enable_plots"]) output.enable_plots = node["enable_plots"].as<bool>();
    if (node["plot_dpi"]) output.plot_dpi = node["plot_dpi"].as<int>();
    if (node["plot_formats"]) output.plot_formats = node["plot_formats"].as<std::vector<std::string>>();
    if (node["plot_6panel"]) output.plot_6panel = node["plot_6panel"].as<bool>();
    if (node["plot_conservation"]) output.plot_conservation = node["plot_conservation"].as<bool>();
    if (node["plot_spatial_fields"]) output.plot_spatial_fields = node["plot_spatial_fields"].as<bool>();
    if (node["plot_gamma_validation"]) output.plot_gamma_validation = node["plot_gamma_validation"].as<bool>();
}

void TestConfig::parseAnalysis(const YAML::Node& node) {
    if (node["mode"]) analysis.mode = node["mode"].as<std::string>();
    if (node["track_trajectory"]) analysis.track_trajectory = node["track_trajectory"].as<bool>();
    if (node["measure_velocity"]) analysis.measure_velocity = node["measure_velocity"].as<bool>();
    if (node["compute_effective_mass"]) analysis.compute_effective_mass = node["compute_effective_mass"].as<bool>();
    if (node["test_gamma_factor"]) analysis.test_gamma_factor = node["test_gamma_factor"].as<bool>();
    if (node["compute_energy"]) analysis.compute_energy = node["compute_energy"].as<bool>();
    if (node["lorentz_invariants"]) analysis.lorentz_invariants = node["lorentz_invariants"].as<bool>();

    // Boosted frame analysis (Phase 2.5B)
    if (node["perform_lorentz_transform"]) analysis.perform_lorentz_transform = node["perform_lorentz_transform"].as<bool>();
    if (node["measure_R_in_boosted_frame"]) analysis.measure_R_in_boosted_frame = node["measure_R_in_boosted_frame"].as<bool>();
    if (node["test_lorentz_covariance"]) analysis.test_lorentz_covariance = node["test_lorentz_covariance"].as<bool>();
    if (node["compute_frame_invariants"]) analysis.compute_frame_invariants = node["compute_frame_invariants"].as<bool>();

    // Dispersion analysis
    if (node["dispersion_k_steps"]) analysis.dispersion_k_steps = node["dispersion_k_steps"].as<int>();
    if (node["dispersion_max_k"]) analysis.dispersion_max_k = node["dispersion_max_k"].as<float>();
}

void TestConfig::parseEMCoupling(const YAML::Node& node) {
    // Parse EM coupling configuration from em_coupling section (Phase 5/6)
    // This allows configs to use em_coupling: section instead of physics.em_coupling_*

    if (node["enabled"]) {
        physics.em_coupling_enabled = node["enabled"].as<bool>();
    }

    if (node["coupling_strength"]) {
        physics.em_coupling_strength = node["coupling_strength"].as<float>();
    }

    if (node["field_type"]) {
        physics.em_coupling_type = node["field_type"].as<std::string>();
    }

    // Optional: Parse field-specific parameters for future use
    // These are not yet stored in PhysicsConfig but can be logged/used in test runner
    if (node["gauge_transformation_enabled"]) {
        // For Test 3: Gauge invariance
        // Store or log this parameter - implementation pending
    }

    if (node["gauge_function"]) {
        // For Test 3: χ(x,y) functional form
    }

    if (node["k_x"] || node["k_y"]) {
        // For Test 3: Gauge gradient wave vector
    }

    if (node["B_z"]) {
        // For Test 1/4: External magnetic field
    }
}
