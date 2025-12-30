#ifndef VALIDATION_CRITERIA_H
#define VALIDATION_CRITERIA_H

#include <string>
#include <vector>
#include <map>

/**
 * ValidationCriteria - Defines global and scenario-specific validation requirements
 *
 * Implements the validation framework from VALIDATION_FRAMEWORK.md:
 * - Global requirements (apply to ALL simulations)
 * - Scenario-specific requirements (vary by test type)
 */

namespace Validation {

// ============================================================================
// GLOBAL REQUIREMENTS (Always enforced for ALL simulations)
// ============================================================================

struct GlobalCriteria {
    // 1. Probability Conservation: d/dt ∫|ψ|² dx dy = 0
    double norm_tolerance = 0.005;  // |Δ||ψ||²| < 0.5%

    // 2. Energy Conservation: dE_total/dt = 0
    double energy_tolerance = 0.01;  // |ΔE/E₀| < 1%

    // 3. Order Parameter Bounds: 0 ≤ R(x,y,t) ≤ 1
    bool enforce_R_bounds = true;

    // 4. Gauge Invariance: ⟨ψ|O|ψ⟩ independent of θ → θ + α(t)
    bool check_gauge_invariance = false;  // Advanced check (optional)

    // 5. Causality: |dx/dt| ≤ c (c=1 in Planck units)
    bool enforce_causality = true;
    double c_light = 1.0;

    // 6. Numerical Stability: All fields finite and non-divergent
    bool check_numerical_stability = true;
};

// ============================================================================
// SCENARIO-SPECIFIC REQUIREMENTS
// ============================================================================

enum class ScenarioType {
    NONE,
    DEFECT_LOCALIZATION,       // Scenario 2.1
    TRAVELING_WAVE,            // Scenario 2.2
    RELATIVISTIC_MASS,         // Scenario 2.3
    BREAKDOWN_INVESTIGATION,   // Scenario 2.4
    VORTEX_PAIR_SEPARATION,    // Test 2.7 (Sprint 2)
    ANNIHILATION_DYNAMICS      // Test 2.8 (Sprint 2)
};

struct DefectLocalizationCriteria {
    // Required:
    bool require_vortex = true;         // W = ±1
    double winding_tolerance = 0.2;     // |W - 1| < 0.2

    bool require_core = true;           // R_min < 0.5
    double core_R_threshold = 0.5;

    bool check_core_stability = true;   // |dx_core/dt| < 0.01
    double core_velocity_max = 0.01;

    // NOT required:
    bool require_particle = false;      // Vacuum only
    bool require_boost = false;
};

struct TravelingWaveCriteria {
    // Required:
    bool require_vortex = true;         // W = ±1
    double winding_tolerance = 0.2;

    bool require_core = true;           // R_min < 0.5
    double core_R_threshold = 0.5;

    bool require_gaussian = true;       // ψ(x,y,t=0) boosted Gaussian
    double gaussian_width_tolerance = 0.2;

    bool require_initial_boost = true;  // ⟨p⟩(t=0) = γmv ± 5%
    double initial_momentum_tolerance = 0.05;  // 5%

    bool check_particle_tracking = true; // Particle tracks vortex motion

    bool validate_gamma_factor = true;   // γ_measured within 5% of theory
    double gamma_tolerance = 0.05;       // 5%

    // NOT required:
    bool require_multiple_vortices = false;
    bool require_grid_convergence = false;
};

struct RelativisticMassCriteria {
    // Includes all TravelingWave requirements:
    TravelingWaveCriteria base_criteria;

    // Additional requirements:
    bool require_grid_convergence = true;
    double grid_convergence_tolerance = 0.02;  // <2% error between grids

    bool require_N_convergence = true;
    double N_convergence_tolerance = 0.05;     // N=10 optimal

    // NOT required:
    bool require_long_time_evolution = false;  // Testing instantaneous mass
};

// ============================================================================
// Sprint 2: Multi-Vortex Interaction Criteria
// ============================================================================

struct VortexPairSeparationCriteria {
    // Required:
    bool require_W_total_zero = true;          // W_total = 0 (vortex-antivortex pair)
    double W_total_tolerance = 0.1;

    bool check_two_cores = true;               // Detect two distinct R-field cores
    double separation_expected = 10.0;         // Expected separation [ℓ_P]
    double separation_tolerance = 5.0;         // ±5 ℓ_P

    bool validate_force_law = false;           // F ∝ 1/d² (requires time series)
    bool validate_binding_energy = false;      // E_int ∝ -log(d/a) (requires time series)

    // NOT required:
    bool require_particle = false;             // Vacuum only
    bool require_long_evolution = false;       // Instantaneous separation scan
};

struct AnnihilationDynamicsCriteria {
    // Required:
    bool require_W_total_zero = true;          // W_total = 0 throughout
    double W_total_tolerance = 0.1;

    bool check_R_field_smoothing = true;       // R_min → 1.0 after annihilation
    double R_final_threshold = 0.9;            // R_min > 0.9 at t_final
    double R_smoothing_tolerance = 0.1;

    bool validate_separation_decrease = false; // d(t) decreases (requires tracking)
    bool validate_em_energy_release = false;   // ΔE_EM during annihilation (requires EM field)

    // NOT required:
    bool require_particle = false;             // Vacuum only
};

// ============================================================================
// VALIDATION RESULTS
// ============================================================================

struct CriterionResult {
    std::string name;
    bool passed;
    double measured_value;
    double expected_value;
    double tolerance;
    std::string message;
    bool is_critical;  // If true, failure invalidates entire simulation
};

struct ValidationReport {
    // Global checks
    std::vector<CriterionResult> global_results;
    bool global_pass = true;

    // Scenario-specific checks
    std::vector<CriterionResult> scenario_results;
    bool scenario_pass = true;

    // Overall
    bool overall_pass = true;
    std::string summary;

    // Serialize to string
    std::string toString() const;

    // Save to file
    void saveToFile(const std::string& filepath) const;
};

// ============================================================================
// VALIDATION CONFIGURATION
// ============================================================================

struct ValidationConfig {
    // Global criteria (always enforced)
    GlobalCriteria global;

    // Scenario type (determines which specific checks to run)
    ScenarioType scenario = ScenarioType::NONE;

    // Scenario-specific criteria
    DefectLocalizationCriteria defect_localization;
    TravelingWaveCriteria traveling_wave;
    RelativisticMassCriteria relativistic_mass;

    // Sprint 2: Multi-vortex criteria
    VortexPairSeparationCriteria vortex_pair_separation;
    AnnihilationDynamicsCriteria annihilation_dynamics;

    // When to check
    bool validate_initial_state = true;
    bool validate_during_evolution = true;  // Check every N steps
    int validation_interval = 100;          // Steps between checks
    bool validate_final_state = true;

    // Output
    bool verbose = true;
    bool fail_on_critical = true;  // Stop simulation on critical failure
    std::string report_directory = "validation_reports";

    // Load from YAML
    static ValidationConfig fromYAML(const std::map<std::string, std::string>& yaml_data);
};

} // namespace Validation

#endif // VALIDATION_CRITERIA_H
