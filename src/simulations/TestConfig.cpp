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
    std::cout << "  delta: " << physics.delta << std::endl;
    std::cout << "  coupling: " << physics.coupling << std::endl;
    std::cout << "  dt: " << physics.dt << std::endl;
    std::cout << "  total_steps: " << physics.total_steps << std::endl;
    std::cout << "  K (Kuramoto): " << physics.K << std::endl;
    std::cout << "  damping: " << physics.damping << std::endl;

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
    if (node["delta"]) physics.delta = node["delta"].as<float>();
    if (node["coupling"]) physics.coupling = node["coupling"].as<float>();
    if (node["dt"]) physics.dt = node["dt"].as<float>();
    if (node["total_steps"]) physics.total_steps = node["total_steps"].as<int>();
    if (node["K"]) physics.K = node["K"].as<float>();
    if (node["damping"]) physics.damping = node["damping"].as<float>();
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
    if (node["norm_tolerance"]) validation.norm_tolerance = node["norm_tolerance"].as<float>();
    if (node["energy_tolerance"]) validation.energy_tolerance = node["energy_tolerance"].as<float>();
    if (node["convergence_tolerance"]) {
        validation.convergence_tolerance = node["convergence_tolerance"].as<float>();
    }
}

void TestConfig::parseOutput(const YAML::Node& node) {
    if (node["directory"]) output.directory = node["directory"].as<std::string>();
    if (node["save_every"]) output.save_every = node["save_every"].as<int>();
    if (node["formats"]) output.formats = node["formats"].as<std::vector<std::string>>();
    if (node["auto_visualize"]) output.auto_visualize = node["auto_visualize"].as<bool>();
    if (node["save_spatial_snapshots"]) output.save_spatial_snapshots = node["save_spatial_snapshots"].as<bool>();
    if (node["snapshot_steps"]) output.snapshot_steps = node["snapshot_steps"].as<std::vector<int>>();
}
