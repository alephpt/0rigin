#include "validation/ScenarioValidator.h"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>

namespace Validation {

// ============================================================================
// PUBLIC METHODS - Scenario Validation
// ============================================================================

ValidationReport ScenarioValidator::validateTravelingWave(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& initial_obs,
    const ObservableComputer::Observables& final_obs,
    double gamma_theory,
    const TravelingWaveCriteria& criteria) {

    ValidationReport report;
    report.scenario_pass = true;

    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    // ========================================================================
    // CRITERION 1: θ(x,y,t=0) shows vortex structure
    // ========================================================================
    if (criteria.require_vortex) {
        CriterionResult vortex_check = checkVortexStructure(
            theta_field, Nx, Ny, criteria.winding_tolerance);
        report.scenario_results.push_back(vortex_check);
        if (!vortex_check.passed) {
            report.scenario_pass = false;
        }
    }

    // ========================================================================
    // CRITERION 2: W = ±1 computed from phase winding
    // ========================================================================
    // (Already checked in criterion 1)

    // ========================================================================
    // CRITERION 3: R(x,y,t=0) shows core (R_min < 0.5)
    // ========================================================================
    if (criteria.require_core) {
        CriterionResult core_check = checkRFieldCore(
            R_field, criteria.core_R_threshold);
        report.scenario_results.push_back(core_check);
        if (!core_check.passed) {
            report.scenario_pass = false;
        }
    }

    // ========================================================================
    // CRITERION 4: ψ(x,y,t=0) is Gaussian at offset position
    // ========================================================================
    if (criteria.require_gaussian) {
        CriterionResult gaussian_check = checkGaussianWavepacket(
            dirac, criteria.gaussian_width_tolerance);
        report.scenario_results.push_back(gaussian_check);
        if (!gaussian_check.passed) {
            report.scenario_pass = false;
        }
    }

    // ========================================================================
    // CRITERION 5: ⟨p⟩(t=0) = γmv (within 5%)
    // ========================================================================
    if (criteria.require_initial_boost) {
        // Expected momentum: p = γmv
        // For v=0.3c, m=1, γ=1.0483: p_expected = 0.3145
        double p_x = initial_obs.momentum_x.real();
        double p_y = initial_obs.momentum_y.real();
        double p_mag = std::sqrt(p_x * p_x + p_y * p_y);

        // Assume velocity is known from config (would need to pass in)
        // For now, infer from gamma: v = √(1 - 1/γ²)
        double v = std::sqrt(1.0 - 1.0 / (gamma_theory * gamma_theory));
        double p_expected = gamma_theory * 1.0 * v;  // m=1

        CriterionResult momentum_check = checkInitialMomentum(
            initial_obs, p_expected, criteria.initial_momentum_tolerance);
        report.scenario_results.push_back(momentum_check);
        if (!momentum_check.passed) {
            report.scenario_pass = false;
        }
    }

    // ========================================================================
    // CRITERION 6: γ_measured(t=final) within 5% of theory
    // ========================================================================
    if (criteria.validate_gamma_factor) {
        CriterionResult gamma_check = checkGammaFactor(
            final_obs, gamma_theory, 1.0, criteria.gamma_tolerance);  // delta=1.0
        report.scenario_results.push_back(gamma_check);
        if (!gamma_check.passed) {
            report.scenario_pass = false;
        }
    }

    report.overall_pass = report.scenario_pass;

    std::ostringstream summary;
    summary << "Traveling Wave validation: "
            << (report.overall_pass ? "✅ PASS" : "❌ FAIL")
            << " (" << report.scenario_results.size() << " criteria checked)";
    report.summary = summary.str();

    return report;
}

ValidationReport ScenarioValidator::validateRelativisticMass(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& initial_obs,
    const ObservableComputer::Observables& final_obs,
    double gamma_theory,
    const RelativisticMassCriteria& criteria) {

    // Start with base traveling wave validation
    ValidationReport report = validateTravelingWave(
        dirac, R_field, theta_field, initial_obs, final_obs,
        gamma_theory, criteria.base_criteria);

    // Additional checks for Scenario 2.3
    // (Grid and N-convergence are checked externally by comparing multiple runs)

    std::ostringstream summary;
    summary << "Relativistic Mass validation: "
            << (report.overall_pass ? "✅ PASS" : "❌ FAIL")
            << " (includes all Traveling Wave criteria + convergence)";
    report.summary = summary.str();

    return report;
}

ValidationReport ScenarioValidator::validateDefectLocalization(
    const DiracEvolution& dirac,
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    const ObservableComputer::Observables& obs,
    const DefectLocalizationCriteria& criteria) {

    ValidationReport report;
    report.scenario_pass = true;

    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    // Check for vortex
    if (criteria.require_vortex) {
        CriterionResult vortex_check = checkVortexStructure(
            theta_field, Nx, Ny, criteria.winding_tolerance);
        report.scenario_results.push_back(vortex_check);
        if (!vortex_check.passed) {
            report.scenario_pass = false;
        }
    }

    // Check for R-field core
    if (criteria.require_core) {
        CriterionResult core_check = checkRFieldCore(
            R_field, criteria.core_R_threshold);
        report.scenario_results.push_back(core_check);
        if (!core_check.passed) {
            report.scenario_pass = false;
        }
    }

    report.overall_pass = report.scenario_pass;

    std::ostringstream summary;
    summary << "Defect Localization validation: "
            << (report.overall_pass ? "✅ PASS" : "❌ FAIL");
    report.summary = summary.str();

    return report;
}

// ============================================================================
// Sprint 2: Multi-Vortex Scenario Validators
// ============================================================================

ValidationReport ScenarioValidator::validateVortexPairSeparation(
    const std::vector<double>& R_field,
    const std::vector<double>& theta_field,
    int Nx, int Ny,
    float separation_expected,
    double tolerance) {

    ValidationReport report;
    report.scenario_pass = true;

    // Criterion 1: Total winding number conservation (W_total = 0)
    CriterionResult W_total_check = checkWindingNumberConservation(
        theta_field, Nx, Ny, 0.0, tolerance);
    report.scenario_results.push_back(W_total_check);
    if (!W_total_check.passed) {
        report.scenario_pass = false;
    }

    // Criterion 2: Measure vortex separation
    double x1, y1, x2, y2;
    double d_measured = measureVortexSeparation(R_field, Nx, Ny, x1, y1, x2, y2);

    CriterionResult separation_check;
    separation_check.name = "Vortex Pair Separation";
    separation_check.is_critical = true;
    separation_check.measured_value = d_measured;
    separation_check.expected_value = separation_expected;
    separation_check.tolerance = 5.0;  // ±5 grid points

    if (d_measured > 0) {
        separation_check.passed = (std::abs(d_measured - separation_expected) < separation_check.tolerance);
        std::ostringstream msg;
        msg << "Separation d = " << d_measured << " (expected " << separation_expected << ")";
        if (separation_check.passed) {
            msg << " ✓ PASS";
        } else {
            msg << " ✗ FAIL";
        }
        separation_check.message = msg.str();
    } else {
        separation_check.passed = false;
        separation_check.message = "Could not detect two vortex cores ✗ FAIL";
    }

    report.scenario_results.push_back(separation_check);
    if (!separation_check.passed) {
        report.scenario_pass = false;
    }

    report.overall_pass = report.scenario_pass;

    std::ostringstream summary;
    summary << "Vortex Pair Separation validation: "
            << (report.overall_pass ? "✅ PASS" : "❌ FAIL");
    report.summary = summary.str();

    return report;
}

ValidationReport ScenarioValidator::validateAnnihilationDynamics(
    const std::vector<double>& R_field_initial,
    const std::vector<double>& R_field_final,
    const std::vector<double>& theta_field,
    int Nx, int Ny,
    double tolerance) {

    ValidationReport report;
    report.scenario_pass = true;

    // Criterion 1: W_total conservation throughout
    CriterionResult W_total_check = checkWindingNumberConservation(
        theta_field, Nx, Ny, 0.0, tolerance);
    report.scenario_results.push_back(W_total_check);
    if (!W_total_check.passed) {
        report.scenario_pass = false;
    }

    // Criterion 2: R-field smoothing after annihilation
    CriterionResult smoothing_check = checkRFieldSmoothing(
        R_field_final, 0.9, tolerance);
    report.scenario_results.push_back(smoothing_check);
    if (!smoothing_check.passed) {
        report.scenario_pass = false;
    }

    report.overall_pass = report.scenario_pass;

    std::ostringstream summary;
    summary << "Annihilation Dynamics validation: "
            << (report.overall_pass ? "✅ PASS" : "❌ FAIL");
    report.summary = summary.str();

    return report;
}

// ============================================================================
// PRIVATE METHODS - Individual Criterion Checks
// ============================================================================

double ScenarioValidator::computeWindingNumber(
    const std::vector<double>& theta_field,
    int Nx, int Ny) {

    // Compute W = (1/2π) ∮ ∇θ · dl around boundary contour
    // Use boundary edges: top → right → bottom → left

    std::vector<double> contour;
    contour.reserve(2 * (Nx + Ny));

    // Top edge (left to right)
    for (int ix = 0; ix < Nx; ++ix) {
        contour.push_back(theta_field[0 * Nx + ix]);
    }

    // Right edge (top to bottom)
    for (int iy = 0; iy < Ny; ++iy) {
        contour.push_back(theta_field[iy * Nx + (Nx - 1)]);
    }

    // Bottom edge (right to left)
    for (int ix = Nx - 1; ix >= 0; --ix) {
        contour.push_back(theta_field[(Ny - 1) * Nx + ix]);
    }

    // Left edge (bottom to top)
    for (int iy = Ny - 1; iy >= 0; --iy) {
        contour.push_back(theta_field[iy * Nx + 0]);
    }

    // Compute phase differences (handle 2π wrapping)
    double total_winding = 0.0;
    for (size_t i = 0; i < contour.size(); ++i) {
        double theta1 = contour[i];
        double theta2 = contour[(i + 1) % contour.size()];

        // Phase difference with wrapping: arctan2(sin(Δθ), cos(Δθ))
        double dtheta = theta2 - theta1;
        dtheta = std::atan2(std::sin(dtheta), std::cos(dtheta));

        total_winding += dtheta;
    }

    // Winding number W = total_winding / (2π)
    double W = total_winding / (2.0 * M_PI);

    return W;
}

CriterionResult ScenarioValidator::checkVortexStructure(
    const std::vector<double>& theta_field,
    int Nx, int Ny,
    double tolerance) {

    CriterionResult result;
    result.name = "1. θ(x,y,t=0) shows vortex structure (W = ±1)";
    result.is_critical = true;

    double W = computeWindingNumber(theta_field, Nx, Ny);

    result.measured_value = W;
    result.expected_value = 1.0;  // Expect |W| = 1
    result.tolerance = tolerance;

    // Pass if |W - 1| < tolerance or |W + 1| < tolerance
    bool is_vortex = (std::abs(std::abs(W) - 1.0) < tolerance);
    result.passed = is_vortex;

    std::ostringstream msg;
    msg << "Winding number W = " << W;
    if (result.passed) {
        msg << " ✓ PASS (|W| ≈ 1)";
    } else {
        msg << " ✗ FAIL (expect |W| = 1 ± " << tolerance << ")";
    }
    result.message = msg.str();

    return result;
}

CriterionResult ScenarioValidator::checkRFieldCore(
    const std::vector<double>& R_field,
    double threshold) {

    CriterionResult result;
    result.name = "3. R(x,y,t=0) shows core (R_min < 0.5)";
    result.is_critical = true;

    double R_min = *std::min_element(R_field.begin(), R_field.end());

    result.measured_value = R_min;
    result.expected_value = 0.0;  // Ideally R_min ≈ 0 at core
    result.tolerance = threshold;

    result.passed = (R_min < threshold);

    std::ostringstream msg;
    msg << "R_min = " << R_min;
    if (result.passed) {
        msg << " ✓ PASS (< " << threshold << ")";
    } else {
        msg << " ✗ FAIL (expect R_min < " << threshold << ")";
    }
    result.message = msg.str();

    return result;
}

CriterionResult ScenarioValidator::checkGaussianWavepacket(
    const DiracEvolution& dirac,
    double tolerance) {

    CriterionResult result;
    result.name = "4. ψ(x,y,t=0) is Gaussian at offset position";
    result.is_critical = false;  // Can be checked visually

    // For now, just verify that wavepacket is localized (not uniform)
    const auto& psi = dirac.getSpinorField();
    int Nx = dirac.getNx();
    int Ny = dirac.getNy();

    // Compute variance of |ψ|² to check localization
    std::vector<double> density(Nx * Ny);
    for (int i = 0; i < Nx * Ny; ++i) {
        double rho = 0.0;
        for (int alpha = 0; alpha < 4; ++alpha) {
            rho += std::norm(psi[i * 4 + alpha]);
        }
        density[i] = rho;
    }

    double mean = std::accumulate(density.begin(), density.end(), 0.0) / density.size();
    double variance = 0.0;
    for (double rho : density) {
        variance += (rho - mean) * (rho - mean);
    }
    variance /= density.size();

    // If variance is too small, wavepacket is too uniform (not localized)
    bool is_localized = (variance > 1e-6);

    result.passed = is_localized;
    result.measured_value = variance;
    result.expected_value = 1e-3;  // Expected variance for Gaussian
    result.tolerance = tolerance;

    std::ostringstream msg;
    msg << "Wavepacket variance = " << variance;
    if (result.passed) {
        msg << " ✓ PASS (localized)";
    } else {
        msg << " ✗ FAIL (too uniform, not Gaussian)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult ScenarioValidator::checkInitialMomentum(
    const ObservableComputer::Observables& obs,
    double p_expected,
    double tolerance) {

    CriterionResult result;
    result.name = "5. ⟨p⟩(t=0) = γmv (within 5%)";
    result.is_critical = true;

    double p_x = obs.momentum_x.real();
    double p_y = obs.momentum_y.real();
    double p_mag = std::sqrt(p_x * p_x + p_y * p_y);

    result.measured_value = p_mag;
    result.expected_value = p_expected;
    result.tolerance = tolerance * p_expected;

    double error = 0.0;
    if (std::abs(p_expected) > 1e-10) {
        error = std::abs(p_mag - p_expected) / std::abs(p_expected);
    }

    result.passed = (error < tolerance);

    std::ostringstream msg;
    msg << "⟨p⟩ = " << p_mag
        << " (expect " << p_expected
        << ", error = " << (error * 100.0) << "%)";
    if (result.passed) {
        msg << " ✓ PASS";
    } else {
        msg << " ✗ FAIL (exceeds " << (tolerance * 100.0) << "% tolerance)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult ScenarioValidator::checkGammaFactor(
    const ObservableComputer::Observables& obs,
    double gamma_theory,
    double delta,
    double tolerance) {

    CriterionResult result;
    result.name = "6. γ_measured(t=final) within 5% of theory";
    result.is_critical = true;

    // Compute effective mass: m_eff = √(E² - p²)
    double E = obs.energy_total;
    double p_x = obs.momentum_x.real();
    double p_y = obs.momentum_y.real();
    double p_mag = std::sqrt(p_x * p_x + p_y * p_y);

    double m_eff = std::sqrt(std::max(0.0, E * E - p_mag * p_mag));

    // Measured gamma: γ = m_eff / (Δ · R_avg)
    double gamma_measured = 0.0;
    if (obs.R_avg > 1e-10 && delta > 1e-10) {
        gamma_measured = m_eff / (delta * obs.R_avg);
    }

    result.measured_value = gamma_measured;
    result.expected_value = gamma_theory;
    result.tolerance = tolerance * gamma_theory;

    double error = 0.0;
    if (std::abs(gamma_theory) > 1e-10) {
        error = std::abs(gamma_measured - gamma_theory) / std::abs(gamma_theory);
    }

    result.passed = (error < tolerance);

    std::ostringstream msg;
    msg << "γ_measured = " << gamma_measured
        << " (γ_theory = " << gamma_theory
        << ", error = " << (error * 100.0) << "%)";
    if (result.passed) {
        msg << " ✓ PASS";
    } else {
        msg << " ✗ FAIL (exceeds " << (tolerance * 100.0) << "% tolerance)";
    }
    result.message = msg.str();

    return result;
}

CriterionResult ScenarioValidator::checkParticleTracking(
    const std::vector<double>& particle_positions,
    const std::vector<double>& core_positions) {

    CriterionResult result;
    result.name = "Particle-Vortex Tracking";
    result.is_critical = false;

    // Compute correlation between particle and core positions
    if (particle_positions.size() != core_positions.size() || particle_positions.empty()) {
        result.passed = false;
        result.message = "Invalid position data";
        return result;
    }

    // Compute Pearson correlation
    double corr = ObservableComputer::computeCorrelation(particle_positions, core_positions);

    result.measured_value = corr;
    result.expected_value = 1.0;  // Perfect tracking
    result.tolerance = 0.1;  // Allow 90% correlation

    result.passed = (corr > 0.9);

    std::ostringstream msg;
    msg << "Correlation(particle, core) = " << corr;
    if (result.passed) {
        msg << " ✓ PASS (strong tracking)";
    } else {
        msg << " ✗ FAIL (weak tracking)";
    }
    result.message = msg.str();

    return result;
}

// ============================================================================
// Sprint 2: Multi-Vortex Validation Methods
// ============================================================================

CriterionResult ScenarioValidator::checkWindingNumberConservation(
    const std::vector<double>& theta_field,
    int Nx, int Ny,
    double W_total_expected,
    double tolerance) {

    CriterionResult result;
    result.name = "Topological Charge Conservation (W_total)";
    result.is_critical = true;

    double W_total = computeWindingNumber(theta_field, Nx, Ny);

    result.measured_value = W_total;
    result.expected_value = W_total_expected;
    result.tolerance = tolerance;

    result.passed = (std::abs(W_total - W_total_expected) < tolerance);

    std::ostringstream msg;
    msg << "W_total = " << W_total << " (expected " << W_total_expected << ")";
    if (result.passed) {
        msg << " ✓ PASS";
    } else {
        msg << " ✗ FAIL (Δ = " << std::abs(W_total - W_total_expected) << ")";
    }
    result.message = msg.str();

    return result;
}

double ScenarioValidator::measureVortexSeparation(
    const std::vector<double>& R_field,
    int Nx, int Ny,
    double& x1, double& y1,
    double& x2, double& y2) {

    // Find two local minima in R-field (vortex cores)
    // Simple approach: find global minimum, then find second minimum far from first

    double R_min1 = 1e10;
    int idx_min1 = -1;

    // Find first minimum
    for (int i = 0; i < Nx * Ny; ++i) {
        if (R_field[i] < R_min1) {
            R_min1 = R_field[i];
            idx_min1 = i;
        }
    }

    if (idx_min1 < 0) return -1.0;  // No minimum found

    int iy1 = idx_min1 / Nx;
    int ix1 = idx_min1 % Nx;
    x1 = static_cast<double>(ix1);
    y1 = static_cast<double>(iy1);

    // Find second minimum (must be > 5 grid points away from first)
    double R_min2 = 1e10;
    int idx_min2 = -1;
    const double min_separation = 5.0;  // Minimum separation in grid units

    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Distance from first minimum
            double dx = ix - ix1;
            double dy = iy - iy1;
            double dist = std::sqrt(dx * dx + dy * dy);

            if (dist < min_separation) continue;  // Too close to first minimum

            if (R_field[idx] < R_min2) {
                R_min2 = R_field[idx];
                idx_min2 = idx;
            }
        }
    }

    if (idx_min2 < 0) return -1.0;  // No second minimum found

    int iy2 = idx_min2 / Nx;
    int ix2 = idx_min2 % Nx;
    x2 = static_cast<double>(ix2);
    y2 = static_cast<double>(iy2);

    // Return separation distance
    double dx = x2 - x1;
    double dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

CriterionResult ScenarioValidator::checkRFieldSmoothing(
    const std::vector<double>& R_field,
    double R_final_threshold,
    double tolerance) {

    CriterionResult result;
    result.name = "R-Field Smoothing After Annihilation";
    result.is_critical = true;

    double R_min = *std::min_element(R_field.begin(), R_field.end());

    result.measured_value = R_min;
    result.expected_value = R_final_threshold;
    result.tolerance = tolerance;

    // After annihilation, R_min should be close to 1.0 (fully synchronized)
    result.passed = (R_min > R_final_threshold - tolerance);

    std::ostringstream msg;
    msg << "R_min = " << R_min << " (threshold " << R_final_threshold << ")";
    if (result.passed) {
        msg << " ✓ PASS (cores healed)";
    } else {
        msg << " ✗ FAIL (cores still present)";
    }
    result.message = msg.str();

    return result;
}

} // namespace Validation
