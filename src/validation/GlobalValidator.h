#ifndef GLOBAL_VALIDATOR_H
#define GLOBAL_VALIDATOR_H

#include "ValidationCriteria.h"
#include "simulations/ObservableComputer.h"
#include "DiracEvolution.h"
#include <vector>

namespace Validation {

/**
 * GlobalValidator - Enforces the 6 global requirements for ALL simulations
 *
 * These checks MUST pass for any SMFT simulation to be considered valid.
 * Violations indicate fundamental physics or numerical errors.
 */
class GlobalValidator {
public:
    /**
     * Validate global requirements at initialization
     *
     * Checks:
     * 1. Probability normalization: ||ψ||² ≈ 1
     * 2. Order parameter bounds: 0 ≤ R ≤ 1
     * 3. Causality (if boosted): ||v|| < c
     * 4. Numerical stability: All fields finite
     */
    static ValidationReport validateInitialState(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        const std::vector<double>& theta_field,
        const GlobalCriteria& criteria);

    /**
     * Validate global requirements during evolution
     *
     * Called every N steps to monitor:
     * 1. Probability conservation
     * 2. Energy conservation
     * 3. R bounds
     * 4. Numerical stability
     */
    static ValidationReport validateRuntimeState(
        const ObservableComputer::Observables& current_obs,
        const ObservableComputer::Observables& initial_obs,
        const std::vector<double>& R_field,
        const GlobalCriteria& criteria,
        double time);

    /**
     * Validate global requirements at end of simulation
     *
     * Final check of conservation laws and bounds
     */
    static ValidationReport validateFinalState(
        const ObservableComputer::Observables& final_obs,
        const ObservableComputer::Observables& initial_obs,
        const std::vector<double>& R_field,
        const GlobalCriteria& criteria);

private:
    // Individual global checks

    /**
     * Check 1: Probability Conservation
     * ||ψ||² should remain ≈ 1 throughout simulation
     */
    static CriterionResult checkProbabilityConservation(
        double norm,
        double norm_initial,
        double tolerance);

    /**
     * Check 2: Energy Conservation
     * E_total should remain constant (isolated system)
     */
    static CriterionResult checkEnergyConservation(
        double E_current,
        double E_initial,
        double tolerance);

    /**
     * Check 3: Order Parameter Bounds
     * 0 ≤ R(x,y,t) ≤ 1 for all points
     */
    static CriterionResult checkOrderParameterBounds(
        const std::vector<double>& R_field);

    /**
     * Check 4: Gauge Invariance (advanced, optional)
     * Observables should be independent of global phase shift
     */
    static CriterionResult checkGaugeInvariance(
        const ObservableComputer::Observables& obs);

    /**
     * Check 5: Causality
     * Particle velocity ||v|| = ||⟨p⟩|| / E should be < c
     */
    static CriterionResult checkCausality(
        const ObservableComputer::Observables& obs,
        double c_light);

    /**
     * Check 6: Numerical Stability
     * All fields should be finite (no NaN, no Inf)
     */
    static CriterionResult checkNumericalStability(
        const ObservableComputer::Observables& obs,
        const std::vector<double>& R_field);
};

} // namespace Validation

#endif // GLOBAL_VALIDATOR_H
