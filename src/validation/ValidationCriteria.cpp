#include "validation/ValidationCriteria.h"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace Validation {

std::string ValidationReport::toString() const {
    std::ostringstream oss;

    oss << "================================================================================\n";
    oss << "VALIDATION REPORT\n";
    oss << "================================================================================\n\n";

    // Global results
    if (!global_results.empty()) {
        oss << "GLOBAL REQUIREMENTS (Apply to ALL simulations):\n";
        oss << "--------------------------------------------------------------------------------\n";
        for (const auto& result : global_results) {
            oss << (result.passed ? "✓" : "✗") << " " << result.name << "\n";
            oss << "  " << result.message << "\n";
            if (!result.passed && result.is_critical) {
                oss << "  ⚠ CRITICAL FAILURE - Simulation INVALID\n";
            }
            oss << "\n";
        }
        oss << "Global status: " << (global_pass ? "✅ PASS" : "❌ FAIL") << "\n\n";
    }

    // Scenario-specific results
    if (!scenario_results.empty()) {
        oss << "SCENARIO-SPECIFIC REQUIREMENTS:\n";
        oss << "--------------------------------------------------------------------------------\n";
        for (const auto& result : scenario_results) {
            oss << (result.passed ? "✓" : "✗") << " " << result.name << "\n";
            oss << "  " << result.message << "\n";
            if (!result.passed && result.is_critical) {
                oss << "  ⚠ CRITICAL FAILURE - Test-specific requirement violated\n";
            }
            oss << "\n";
        }
        oss << "Scenario status: " << (scenario_pass ? "✅ PASS" : "❌ FAIL") << "\n\n";
    }

    // Overall summary
    oss << "================================================================================\n";
    oss << "OVERALL STATUS: " << (overall_pass ? "✅ PASS" : "❌ FAIL") << "\n";
    oss << summary << "\n";
    oss << "================================================================================\n";

    return oss.str();
}

void ValidationReport::saveToFile(const std::string& filepath) const {
    std::ofstream file(filepath);
    if (file.is_open()) {
        file << toString();
        file.close();
    }
}

ValidationConfig ValidationConfig::fromYAML(const std::map<std::string, std::string>& yaml_data) {
    ValidationConfig config;

    // Parse global criteria
    if (yaml_data.count("norm_tolerance")) {
        config.global.norm_tolerance = std::stod(yaml_data.at("norm_tolerance"));
    }
    if (yaml_data.count("energy_tolerance")) {
        config.global.energy_tolerance = std::stod(yaml_data.at("energy_tolerance"));
    }
    if (yaml_data.count("R_bounds_check")) {
        config.global.enforce_R_bounds = (yaml_data.at("R_bounds_check") == "true");
    }
    if (yaml_data.count("causality_check")) {
        config.global.enforce_causality = (yaml_data.at("causality_check") == "true");
    }

    // Parse scenario type
    if (yaml_data.count("scenario")) {
        const std::string& scenario_str = yaml_data.at("scenario");
        if (scenario_str == "defect_localization") {
            config.scenario = ScenarioType::DEFECT_LOCALIZATION;
        } else if (scenario_str == "traveling_wave") {
            config.scenario = ScenarioType::TRAVELING_WAVE;
        } else if (scenario_str == "relativistic_mass") {
            config.scenario = ScenarioType::RELATIVISTIC_MASS;
        } else if (scenario_str == "breakdown_investigation") {
            config.scenario = ScenarioType::BREAKDOWN_INVESTIGATION;
        }
    }

    // Parse scenario-specific criteria
    if (yaml_data.count("gamma_tolerance")) {
        config.traveling_wave.gamma_tolerance = std::stod(yaml_data.at("gamma_tolerance"));
    }
    if (yaml_data.count("require_vortex")) {
        config.traveling_wave.require_vortex = (yaml_data.at("require_vortex") == "true");
    }
    if (yaml_data.count("require_core")) {
        config.traveling_wave.require_core = (yaml_data.at("require_core") == "true");
    }
    if (yaml_data.count("require_boost")) {
        config.traveling_wave.require_initial_boost = (yaml_data.at("require_boost") == "true");
    }

    // Parse validation timing
    if (yaml_data.count("validation_interval")) {
        config.validation_interval = std::stoi(yaml_data.at("validation_interval"));
    }

    return config;
}

} // namespace Validation
