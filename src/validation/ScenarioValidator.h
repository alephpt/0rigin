#ifndef SCENARIO_VALIDATOR_H
#define SCENARIO_VALIDATOR_H

#include "ValidationCriteria.h"
#include "simulations/ObservableComputer.h"
#include "DiracEvolution.h"
#include <vector>

namespace Validation {

/**
 * ScenarioValidator - Enforces scenario-specific validation requirements
 *
 * Different simulation scenarios have different physics requirements.
 * This validator applies the appropriate checks based on scenario type.
 */
class ScenarioValidator {
public:
    /**
     * Validate Scenario 2.1: Defect Localization
     *
     * Required:
     * - Vortex present: W = ±1
     * - R-field core: R_min < 0.5
     * - Core position stable: |dx_core/dt| < 0.01
     *
     * NOT required:
     * - Particle present (vacuum only)
     * - Boosted conditions
     */
    static ValidationReport validateDefectLocalization(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        const std::vector<double>& theta_field,
        const ObservableComputer::Observables& obs,
        const DefectLocalizationCriteria& criteria);

    /**
     * Validate Scenario 2.2: Traveling Wave Surfing
     *
     * Required (at t=0):
     * - θ(x,y,t=0) vortex: W = ±1
     * - R(x,y,t=0) core: R_min < 0.5
     * - ψ(x,y,t=0) boosted Gaussian at offset
     * - ⟨p⟩(t=0) = γmv ± 5%
     *
     * Required (at t=final):
     * - Particle tracks vortex motion
     * - γ_measured within 5% of theory
     */
    static ValidationReport validateTravelingWave(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        const std::vector<double>& theta_field,
        const ObservableComputer::Observables& initial_obs,
        const ObservableComputer::Observables& final_obs,
        double gamma_theory,
        const TravelingWaveCriteria& criteria);

    /**
     * Validate Scenario 2.3: Relativistic Mass Validation
     *
     * Includes all Scenario 2.2 requirements, plus:
     * - Grid convergence demonstrated
     * - Operator splitting convergence (N=1, 10, 100)
     */
    static ValidationReport validateRelativisticMass(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        const std::vector<double>& theta_field,
        const ObservableComputer::Observables& initial_obs,
        const ObservableComputer::Observables& final_obs,
        double gamma_theory,
        const RelativisticMassCriteria& criteria);

    /**
     * Validate Test 2.7: Vortex Pair Separation Scan (Sprint 2)
     *
     * Required:
     * - W_total = 0 (topological charge conservation)
     * - Two distinct cores in R-field
     * - Binding energy E_int(d) ∝ -log(d/a)
     * - Interaction force F(d) ∝ 1/d²
     */
    static ValidationReport validateVortexPairSeparation(
        const std::vector<double>& R_field,
        const std::vector<double>& theta_field,
        int Nx, int Ny,
        float separation_expected,
        double tolerance);

    /**
     * Validate Test 2.8: Annihilation Dynamics (Sprint 2)
     *
     * Required:
     * - W_total = 0 throughout evolution
     * - Vortex separation d(t) decreases monotonically
     * - R-field smoothing: R_min → 1.0 after annihilation
     * - Energy conservation with EM energy release
     */
    static ValidationReport validateAnnihilationDynamics(
        const std::vector<double>& R_field_initial,
        const std::vector<double>& R_field_final,
        const std::vector<double>& theta_field,
        int Nx, int Ny,
        double tolerance);

private:
    // Individual scenario-specific checks

    /**
     * Compute topological winding number W = (1/2π) ∮ ∇θ · dl
     * Returns W ∈ ℤ (should be ±1 for single vortex)
     */
    static double computeWindingNumber(
        const std::vector<double>& theta_field,
        int Nx, int Ny);

    /**
     * Check if vortex is present: |W - 1| < tolerance
     */
    static CriterionResult checkVortexStructure(
        const std::vector<double>& theta_field,
        int Nx, int Ny,
        double tolerance);

    /**
     * Check if R-field has vortex core: R_min < threshold
     */
    static CriterionResult checkRFieldCore(
        const std::vector<double>& R_field,
        double threshold);

    /**
     * Check if wavepacket is Gaussian-shaped at offset position
     */
    static CriterionResult checkGaussianWavepacket(
        const DiracEvolution& dirac,
        double tolerance);

    /**
     * Check initial momentum: ⟨p⟩(t=0) = γmv ± tolerance
     */
    static CriterionResult checkInitialMomentum(
        const ObservableComputer::Observables& obs,
        double p_expected,
        double tolerance);

    /**
     * Check final gamma factor: γ_measured = γ_theory ± tolerance
     *
     * Measured from late-time data:
     * m_eff = √(E² - p²)
     * γ_measured = m_eff / (Δ · R_avg)
     */
    static CriterionResult checkGammaFactor(
        const ObservableComputer::Observables& obs,
        double gamma_theory,
        double delta,
        double tolerance);

    /**
     * Check if particle tracks vortex motion
     * Correlate ⟨x⟩(t) with x_core(t)
     */
    static CriterionResult checkParticleTracking(
        const std::vector<double>& particle_positions,
        const std::vector<double>& core_positions);

    // Sprint 2: Multi-vortex validation helpers

    /**
     * Check total winding number conservation for vortex pair
     * W_total = W₁ + W₂ = 0 for vortex-antivortex pair
     */
    static CriterionResult checkWindingNumberConservation(
        const std::vector<double>& theta_field,
        int Nx, int Ny,
        double W_total_expected,
        double tolerance);

    /**
     * Detect vortex cores in R-field and measure separation
     * Returns separation distance d in grid units
     */
    static double measureVortexSeparation(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double& x1, double& y1,
        double& x2, double& y2);

    /**
     * Check if R-field has smoothed out after annihilation
     * R_min should approach R_max ≈ 1.0
     */
    static CriterionResult checkRFieldSmoothing(
        const std::vector<double>& R_field,
        double R_final_threshold,
        double tolerance);
};

} // namespace Validation

#endif // SCENARIO_VALIDATOR_H
