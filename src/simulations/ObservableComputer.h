#ifndef OBSERVABLE_COMPUTER_H
#define OBSERVABLE_COMPUTER_H

#include "DiracEvolution.h"
#include "KleinGordonEvolution.h"
#include <complex>
#include <vector>
#include <map>
#include <string>

/**
 * ObservableComputer - Centralized computation of all physical observables
 *
 * Provides consistent, validated implementations of:
 * - Norm conservation: ||Ψ||²
 * - Energy (kinetic + potential)
 * - Position/momentum expectation values
 * - Sync field statistics (R_avg, R_max, etc.)
 *
 * Used by unified test framework to ensure quantitative validation.
 */
class ObservableComputer {
public:
    struct Observables {
        double time;

        // Dirac observables
        double norm;              // ||Ψ||² - should be ≈ 1.0
        double norm_error;        // ||Ψ||² - 1.0
        double energy_total;
        double energy_kinetic;
        double energy_potential;

        // Position expectation <Ψ|x|Ψ>, <Ψ|y|Ψ>
        std::complex<double> position_x;
        std::complex<double> position_y;

        // Momentum expectation <Ψ|p_x|Ψ>, <Ψ|p_y|Ψ>
        std::complex<double> momentum_x;
        std::complex<double> momentum_y;

        // Acceleration (computed from velocity time derivative)
        std::complex<double> acceleration_x;
        std::complex<double> acceleration_y;

        // Geodesic acceleration (from spacetime geometry)
        double geodesic_acceleration_x;
        double geodesic_acceleration_y;

        // Christoffel symbols at particle position
        double Gamma_x_tt;
        double Gamma_y_tt;

        // Curvature diagnostics
        double riemann_scalar;
        double gaussian_curvature;

        // Kuramoto sync field observables
        double R_avg;
        double R_max;
        double R_min;
        double R_variance;

        // Validation flags
        bool norm_valid;          // |norm_error| < tolerance
        bool energy_valid;        // |ΔE/E₀| < tolerance
    };

    /**
     * Two-particle observables (Test 3.4: Antiparticle Separation)
     */
    struct TwoParticleObservables {
        double time;

        // Particle observables
        double norm_particle;
        double position_x_particle;
        double position_y_particle;
        double momentum_x_particle;
        double momentum_y_particle;

        // Antiparticle observables
        double norm_antiparticle;
        double position_x_antiparticle;
        double position_y_antiparticle;
        double momentum_x_antiparticle;
        double momentum_y_antiparticle;

        // Separation dynamics
        double separation_distance;  // |r_particle - r_antiparticle|
        double force_particle_x;     // Force on particle
        double force_particle_y;
        double force_antiparticle_x; // Force on antiparticle
        double force_antiparticle_y;
        double force_dot_product;    // F_p · F_a (should be < 0 for opposite forces)
    };

    /**
     * Compute all observables for current state
     *
     * @param dirac Dirac evolution state
     * @param R_field Kuramoto sync field (size Nx*Ny)
     * @param delta Vacuum potential Δ
     * @param time Current simulation time
     * @param E0 Initial energy (for relative energy conservation check)
     * @param em_field_energy EM field energy contribution (default 0.0)
     * @param norm_tolerance Tolerance for ||Ψ||² ≈ 1 (default 1e-4)
     * @param energy_tolerance Tolerance for ΔE/E₀ (default 1e-2)
     * @return Observables struct with all computed values
     */
    static Observables compute(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta,
        double time,
        double E0 = 0.0,
        double em_field_energy = 0.0,
        double norm_tolerance = 1e-4,
        double energy_tolerance = 1e-2);

    /**
     * Compute norm ||Ψ||² = ∫|Ψ|² dA
     * For 4-component spinor: ||Ψ||² = Σ_i (|ψ₁|² + |ψ₂|² + |ψ₃|² + |ψ₄|²)
     */
    static double computeNorm(const DiracEvolution& dirac);

    /**
     * Compute total energy E = T + V
     * T = kinetic energy (Dirac operator terms)
     * V = potential energy (mass coupling)
     */
    static double computeEnergy(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta);

    /**
     * Compute kinetic energy from Dirac operator
     * T = ∫ Ψ†(-iα·∇)Ψ dA
     */
    static double computeKineticEnergy(const DiracEvolution& dirac);

    /**
     * Compute potential energy from mass field
     * V = ∫ Ψ†(β·m(x))Ψ dA where m(x) = Δ·R(x)
     */
    static double computePotentialEnergy(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta);

    /**
     * Compute position expectation value <Ψ|x|Ψ> or <Ψ|y|Ψ>
     * component: 0 for x, 1 for y
     */
    static std::complex<double> computePositionExpectation(
        const DiracEvolution& dirac,
        int component);

    /**
     * Compute momentum expectation value <Ψ|p_x|Ψ> or <Ψ|p_y|Ψ>
     * Uses p = -i∇ operator
     * component: 0 for p_x, 1 for p_y
     */
    static std::complex<double> computeMomentumExpectation(
        const DiracEvolution& dirac,
        int component);

    // ========== Two-Particle Observable Methods (Test 3.4) ==========

    /**
     * Compute two-particle observables for particle-antiparticle separation test
     *
     * @param particle Particle DiracEvolution (β = +1)
     * @param antiparticle Antiparticle DiracEvolution (β = -1)
     * @param R_field Kuramoto sync field (size Nx*Ny)
     * @param delta Vacuum potential Δ
     * @param time Current simulation time
     * @return TwoParticleObservables struct with all computed values
     */
    static TwoParticleObservables computeTwoParticle(
        const DiracEvolution& particle,
        const DiracEvolution& antiparticle,
        const std::vector<double>& R_field,
        double delta,
        double time);

    // ========== Klein-Gordon Observable Methods ==========

    /**
     * Compute all observables for Klein-Gordon scalar field
     *
     * @param kg Klein-Gordon evolution state
     * @param R_field Kuramoto sync field (size Nx*Ny)
     * @param delta Vacuum potential Δ
     * @param time Current simulation time
     * @param E0 Initial energy (for relative energy conservation check)
     * @param em_field_energy EM field energy contribution (default 0.0)
     * @param norm_tolerance Tolerance for ||φ||² ≈ 1 (default 1e-4)
     * @param energy_tolerance Tolerance for ΔE/E₀ (default 1e-2)
     * @return Observables struct with all computed values
     */
    static Observables computeKG(
        const KleinGordonEvolution& kg,
        const std::vector<double>& R_field,
        double delta,
        double time,
        double E0 = 0.0,
        double em_field_energy = 0.0,
        double norm_tolerance = 1e-4,
        double energy_tolerance = 1e-2);

    /**
     * Compute norm ||φ||² = ∫|φ|² dA for Klein-Gordon scalar field
     */
    static double computeNormKG(const KleinGordonEvolution& kg);

    /**
     * Compute total energy for Klein-Gordon field
     * E = ∫[|∂_tφ|² + |∇φ|² + m²|φ|²] dx
     */
    static double computeEnergyKG(
        const KleinGordonEvolution& kg,
        const std::vector<double>& R_field,
        double delta);

    /**
     * Compute kinetic energy for Klein-Gordon field
     * T = ∫|∂_tφ|² dx
     */
    static double computeKineticEnergyKG(const KleinGordonEvolution& kg);

    /**
     * Compute potential energy for Klein-Gordon field
     * V = ∫[|∇φ|² + m²|φ|²] dx
     */
    static double computePotentialEnergyKG(
        const KleinGordonEvolution& kg,
        const std::vector<double>& R_field,
        double delta);

    /**
     * Compute position expectation value <x> or <y> for Klein-Gordon field
     * <x> = ∫ x|φ|² dx / ∫|φ|² dx
     * component: 0 for x, 1 for y
     */
    static std::complex<double> computePositionExpectationKG(
        const KleinGordonEvolution& kg,
        int component);

    /**
     * Compute momentum expectation value <p_x> or <p_y> for Klein-Gordon field
     * <p> = -i ∫ φ*(∇φ) dx
     * component: 0 for p_x, 1 for p_y
     */
    static std::complex<double> computeMomentumExpectationKG(
        const KleinGordonEvolution& kg,
        int component);

    /**
     * Compute sync field statistics
     * Returns {R_avg, R_max, R_min, R_variance}
     */
    static std::tuple<double, double, double, double> computeSyncFieldStats(
        const std::vector<double>& R_field);

    /**
     * Compute R-field gradient at particle center of mass
     * Returns (∇R_x, ∇R_y) using centered finite differences
     */
    static std::pair<double, double> computeRFieldGradient(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double pos_x, double pos_y,
        double dx = 1.0);

    /**
     * Compute Pearson correlation coefficient between two time series
     * Returns ρ ∈ [-1, 1]
     */
    static double computeCorrelation(
        const std::vector<double>& series1,
        const std::vector<double>& series2);

    // ========== Phase 4 Time Gradient Observables ==========

    /**
     * Compute global phase accumulation φ = arg(⟨Ψ|Ψ⟩)
     * Used for Test 4.1: Time Dilation Measurement
     *
     * Physics: Tracks overall phase evolution of wavepacket
     * In time-dilation mode, φ(t) evolves slower where R < 1
     *
     * @param dirac Dirac evolution state
     * @return Global phase in radians [-π, π]
     */
    static double computePhaseAccumulation(const DiracEvolution& dirac);

    /**
     * Compute temporal force from momentum change and R-field gradient
     * Used for Test 4.2: Temporal Force Measurement
     *
     * Physics: F_temp = dp/dt - F_spatial = dp/dt + Δ·∇R
     * Measures residual force after accounting for spatial gradient
     *
     * @param dirac Dirac evolution state
     * @param R_field Current synchronization field
     * @param momentum_prev Previous momentum expectation (p_x, p_y)
     * @param dt Timestep
     * @param delta Vacuum potential Δ
     * @return (F_temp_x, F_temp_y) tuple
     */
    static std::tuple<double, double> computeTemporalForce(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        std::complex<double> momentum_prev_x,
        std::complex<double> momentum_prev_y,
        double dt,
        double delta);

    /**
     * Compute geodesic acceleration from R-field geometry
     * Used for Test 4.3: Geodesic Deviation Analysis
     *
     * Physics: a_geo = -Γ^i_μν v^μ v^ν (geodesic equation)
     * For metric ds² = -R²dt² + dx² + dy²:
     *   a^x ≈ -Γ^x_00 = R ∂R/∂x (dominant term)
     *   a^y ≈ -Γ^y_00 = R ∂R/∂y
     *
     * @param dirac Dirac evolution state
     * @param R_field Current synchronization field
     * @param L_domain Domain size (Planck units)
     * @return (a_geo_x, a_geo_y) tuple
     */
    static std::tuple<double, double> computeGeodesicAcceleration(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double L_domain);

    /**
     * Compute Christoffel symbols at particle position
     * Used for Test 4.3: Geodesic Deviation Analysis
     *
     * @param dirac Dirac evolution state
     * @param R_field Current synchronization field
     * @param L_domain Domain size (Planck units)
     * @return (Γ^x_00, Γ^y_00) tuple
     */
    static std::tuple<double, double> computeChristoffelAtParticle(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double L_domain);

    /**
     * Compute Riemann curvature scalar at particle position
     * Used for Test 4.3: Geodesic Deviation Analysis
     *
     * @param dirac Dirac evolution state
     * @param R_field Current synchronization field
     * @return Riemann scalar R
     */
    static double computeRiemannScalarAtParticle(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field);

    /**
     * Compute Gaussian curvature at particle position
     * Used for Test 4.3: Geodesic Deviation Analysis
     *
     * @param dirac Dirac evolution state
     * @param R_field Current synchronization field
     * @return Gaussian curvature K
     */
    static double computeGaussianCurvatureAtParticle(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field);

    // ========== Phase 5 EM Field Observables ==========

    /**
     * Electromagnetic field observables structure
     *
     * Tracks EM field properties extracted from Kuramoto phase gradients
     * for validation of electromagnetic coupling hypothesis (Phase 5).
     */
    struct EMObservables {
        double field_energy;        // ∫ (E² + B²)/(8π) d²x
        double max_E_magnitude;     // max |E|
        double max_B_magnitude;     // max |B|
        double total_flux;          // ∫ |S| d²x (Poynting)
        double avg_lorentz_force;   // ⟨|F_L|⟩

        // Validation metrics
        double charge_density_rms;  // RMS(ρ)
        double current_density_rms; // RMS(j)
        double maxwell_violation;   // |∇·E - 4πρ| + |∇×B - 4πj/c|
    };

    /**
     * Compute EM field observables from Kuramoto phase field
     *
     * Extracts gauge potential A_μ from phase gradients, computes field strengths,
     * and validates against Maxwell equations and particle dynamics.
     *
     * Physics:
     *   - A_μ = ∂_μ θ (gauge potential from phase)
     *   - E = -∇φ - ∂_t A (electric field)
     *   - B = ∇×A (magnetic field)
     *   - F_L = ρE + J×B (Lorentz force)
     *
     * Used for Phase 5 (Scenario 2.6B): Electromagnetic Coupling Validation
     *
     * @param theta_current Current Kuramoto phase field θ(t) [Nx*Ny]
     * @param theta_previous Previous phase field θ(t-dt) [Nx*Ny]
     * @param psi Dirac spinor field (4-component) [4*Nx*Ny]
     * @param Nx Grid width
     * @param Ny Grid height
     * @param dx Grid spacing x
     * @param dy Grid spacing y
     * @param dt Timestep
     * @return EMObservables struct with field diagnostics and validation metrics
     */
    static EMObservables computeEMObservables(
        const std::vector<float>& theta_current,
        const std::vector<float>& theta_previous,
        const std::vector<std::complex<double>>& psi,
        int Nx, int Ny,
        double dx, double dy, double dt);

    // ========== Phase 3 Vacuum Structure Observables ==========

    /**
     * Compute energy density by region (Test 3.2: Vacuum Energy)
     *
     * Divides domain into left/right regions and computes:
     * - ρ_left = E_left / V_left
     * - ρ_right = E_right / V_right
     * - ⟨R⟩_left, ⟨R⟩_right
     *
     * @param dirac Dirac evolution state
     * @param R_field Order parameter field
     * @param delta Vacuum potential Δ
     * @param x_boundary Region boundary (grid units)
     * @return {ρ_left, ρ_right, R_avg_left, R_avg_right}
     */
    static std::tuple<double, double, double, double> computeEnergyDensityByRegion(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta,
        float x_boundary);

    /**
     * Compute force between parallel defect sheets (Test 3.1: Casimir Force)
     *
     * Measures force per unit length between linear defects via:
     * F = ∫ β⟨Ψ|β|Ψ⟩ · Δ∇R dy
     *
     * Integration along defect line 1 (x = x1)
     *
     * @param dirac Dirac evolution state
     * @param R_field Order parameter field
     * @param delta Vacuum potential Δ
     * @param x1 Position of first defect (grid units)
     * @return Force magnitude (positive = repulsive, negative = attractive)
     */
    static double computeDefectForce(
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta,
        float x1);

    /**
     * Write observables to CSV line
     * Format: time,norm,norm_error,E_total,E_kin,E_pot,<x_re>,<x_im>,...
     */
    static std::string toCSVLine(const Observables& obs);

    /**
     * Get CSV header
     */
    static std::string getCSVHeader();
};

#endif // OBSERVABLE_COMPUTER_H
