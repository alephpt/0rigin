#include "validation/GlobalValidator.h"
#include <cmath>
#include <sstream>
#include <algorithm>

namespace Validation {

// ============================================================================
// PUBLIC METHODS
// ============================================================================

ValidationReport GlobalValidator::validateInitialState(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const GlobalCriteria& criteria) {

    ValidationReport report;
    report.global_pass = true;

    // Check 1: Probability normalization
    double norm = ObservableComputer::computeNorm(dirac);
    CriterionResult norm_check = checkProbabilityConservation(norm, 1.0, criteria.norm_tolerance);
    report.global_results.push_back(norm_check);
    if (!norm_check.passed && norm_check.is_critical) {
        report.global_pass = false;
    }

    // Check 3: Order parameter bounds
    if (criteria.enforce_R_bounds) {
        CriterionResult bounds_check = checkOrderParameterBounds(R_field);
        report.global_results.push_back(bounds_check);
        if (!bounds_check.passed && bounds_check.is_critical) {
            report.global_pass = false;
        }
    }

    // Check 6: Numerical stability
    if (criteria.check_numerical_stability) {
        // Create dummy observables for initial check
        ObservableComputer::Observables obs;
        obs.norm = norm;
        obs.energy_total = 0.0;
        obs.position_x = std::complex<double>(0.0, 0.0);
        obs.position_y = std::complex<double>(0.0, 0.0);
        obs.momentum_x = std::complex<double>(0.0, 0.0);
        obs.momentum_y = std::complex<double>(0.0, 0.0);
        obs.R_avg = 0.0;

        CriterionResult stability_check = checkNumericalStability(obs, R_field);
        report.global_results.push_back(stability_check);
        if (!stability_check.passed && stability_check.is_critical) {
            report.global_pass = false;
        }
    }

    report.overall_pass = report.global_pass;

    std::ostringstream summary;
    summary << "Initial state validation: "
            << (report.overall_pass ? "PASS" : "FAIL")
            << " (" << report.global_results.size() << " checks)";
    report.summary = summary.str();

    return report;
}

ValidationReport GlobalValidator::validateRuntimeState(
    const ObservableComputer::Observables& current_obs,
    const ObservableComputer::Observables& initial_obs,
    const std::vector<double>& R_field,
    const GlobalCriteria& criteria,
    double time) {

    ValidationReport report;
    report.global_pass = true;

    // Check 1: Probability conservation
    CriterionResult norm_check = checkProbabilityConservation(
        current_obs.norm, initial_obs.norm, criteria.norm_tolerance);
    report.global_results.push_back(norm_check);
    if (!norm_check.passed && norm_check.is_critical) {
        report.global_pass = false;
    }

    // Check 2: Energy conservation
    CriterionResult energy_check = checkEnergyConservation(
        current_obs.energy_total, initial_obs.energy_total, criteria.energy_tolerance);
    report.global_results.push_back(energy_check);
    if (!energy_check.passed && energy_check.is_critical) {
        report.global_pass = false;
    }

    // Check 3: Order parameter bounds
    if (criteria.enforce_R_bounds) {
        CriterionResult bounds_check = checkOrderParameterBounds(R_field);
        report.global_results.push_back(bounds_check);
        if (!bounds_check.passed && bounds_check.is_critical) {
            report.global_pass = false;
        }
    }

    // Check 5: Causality
    if (criteria.enforce_causality) {
        CriterionResult causality_check = checkCausality(current_obs, criteria.c_light);
        report.global_results.push_back(causality_check);
        if (!causality_check.passed && causality_check.is_critical) {
            report.global_pass = false;
        }
    }

    // Check 6: Numerical stability
    if (criteria.check_numerical_stability) {
        CriterionResult stability_check = checkNumericalStability(current_obs, R_field);
        report.global_results.push_back(stability_check);
        if (!stability_check.passed && stability_check.is_critical) {
            report.global_pass = false;
        }
    }

    report.overall_pass = report.global_pass;

    std::ostringstream summary;
    summary << "Runtime validation at t=" << time << ": "
            << (report.overall_pass ? "PASS" : "FAIL");
    report.summary = summary.str();

    return report;
}

ValidationReport GlobalValidator::validateFinalState(
    const ObservableComputer::Observables& final_obs,
    const ObservableComputer::Observables& initial_obs,
    const std::vector<double>& R_field,
    const GlobalCriteria& criteria) {

    // Use runtime validation for final state
    return validateRuntimeState(final_obs, initial_obs, R_field, criteria, final_obs.time);
}

// ============================================================================
// PRIVATE METHODS - Individual Checks
// ============================================================================

CriterionResult GlobalValidator::checkProbabilityConservation(
    double norm,
    double norm_initial,
    double tolerance) {

    CriterionResult result;
    result.name = "Probability Conservation";
    result.measured_value = norm;
    result.expected_value = norm_initial;
    result.tolerance = tolerance;
    result.is_critical = true;

    double drift = std::abs(norm - norm_initial);
    result.passed = (drift < tolerance);

    std::ostringstream msg;
    msg << "||ψ||² = " << norm
        << " (drift: " << (drift * 100.0) << "% from initial " << norm_initial << ")";
    if (result.passed) {
        msg << " ✓ PASS";
    } else {
        msg << " ✗ FAIL (exceeds tolerance " << (tolerance * 100.0) << "%)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult GlobalValidator::checkEnergyConservation(
    double E_current,
    double E_initial,
    double tolerance) {

    CriterionResult result;
    result.name = "Energy Conservation";
    result.measured_value = E_current;
    result.expected_value = E_initial;
    result.tolerance = tolerance;
    result.is_critical = true;

    double drift = 0.0;
    if (std::abs(E_initial) > 1e-10) {
        drift = std::abs(E_current - E_initial) / std::abs(E_initial);
    }
    result.passed = (drift < tolerance);

    std::ostringstream msg;
    msg << "E_total = " << E_current
        << " (drift: " << (drift * 100.0) << "% from initial " << E_initial << ")";
    if (result.passed) {
        msg << " ✓ PASS";
    } else {
        msg << " ✗ FAIL (exceeds tolerance " << (tolerance * 100.0) << "%)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult GlobalValidator::checkOrderParameterBounds(
    const std::vector<double>& R_field) {

    CriterionResult result;
    result.name = "Order Parameter Bounds";
    result.is_critical = true;

    // Find min and max
    double R_min = *std::min_element(R_field.begin(), R_field.end());
    double R_max = *std::max_element(R_field.begin(), R_field.end());

    result.measured_value = R_min;  // Store min for diagnostic
    result.expected_value = 0.0;
    result.tolerance = 1.0;

    // Check if all values in [0, 1]
    bool all_in_bounds = (R_min >= 0.0) && (R_max <= 1.0);
    result.passed = all_in_bounds;

    std::ostringstream msg;
    msg << "R ∈ [" << R_min << ", " << R_max << "]";
    if (result.passed) {
        msg << " ✓ PASS (within [0,1])";
    } else {
        msg << " ✗ FAIL (out of bounds [0,1])";
    }
    result.message = msg.str();

    return result;
}

CriterionResult GlobalValidator::checkGaugeInvariance(
    const ObservableComputer::Observables& obs) {

    CriterionResult result;
    result.name = "Gauge Invariance";
    result.is_critical = false;  // Advanced check, not critical
    result.passed = true;  // TODO: Implement actual gauge invariance test
    result.message = "Gauge invariance check not implemented (advanced feature)";

    return result;
}

CriterionResult GlobalValidator::checkCausality(
    const ObservableComputer::Observables& obs,
    double c_light) {

    CriterionResult result;
    result.name = "Causality";
    result.is_critical = true;

    // Compute velocity: v = |p| / E
    double p_mag = std::sqrt(
        obs.momentum_x.real() * obs.momentum_x.real() +
        obs.momentum_y.real() * obs.momentum_y.real()
    );

    double velocity = 0.0;
    if (std::abs(obs.energy_total) > 1e-10) {
        velocity = p_mag / std::abs(obs.energy_total);
    }

    result.measured_value = velocity;
    result.expected_value = c_light;
    result.tolerance = c_light * 0.01;  // Allow 1% margin

    result.passed = (velocity <= c_light * 1.01);  // Small margin for numerical errors

    std::ostringstream msg;
    msg << "v = |p|/E = " << velocity << " (c = " << c_light << ")";
    if (result.passed) {
        msg << " ✓ PASS (v ≤ c)";
    } else {
        msg << " ✗ FAIL (superluminal!)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult GlobalValidator::checkNumericalStability(
    const ObservableComputer::Observables& obs,
    const std::vector<double>& R_field) {

    CriterionResult result;
    result.name = "Numerical Stability";
    result.is_critical = true;

    // Check for NaN or Inf in observables
    bool stable = std::isfinite(obs.norm) &&
                  std::isfinite(obs.energy_total) &&
                  std::isfinite(obs.position_x.real()) &&
                  std::isfinite(obs.position_y.real()) &&
                  std::isfinite(obs.momentum_x.real()) &&
                  std::isfinite(obs.momentum_y.real()) &&
                  std::isfinite(obs.R_avg);

    // Check R field
    if (stable) {
        for (double r : R_field) {
            if (!std::isfinite(r)) {
                stable = false;
                break;
            }
        }
    }

    result.passed = stable;
    result.message = stable ? "All fields finite ✓ PASS" : "NaN or Inf detected ✗ FAIL";

    return result;
}

} // namespace Validation
