/**
 * EMFieldComputer.h
 *
 * Electromagnetic field extraction from Kuramoto phase field
 *
 * Physics Foundation:
 *   Hypothesis: Kuramoto phase θ(x,y,t) encodes electromagnetic potential
 *
 *   GAUGE POTENTIAL (A_μ = ∂_μ θ):
 *     φ = ∂_t θ     (scalar potential / temporal phase gradient)
 *     A = ∇θ        (vector potential / spatial phase gradients)
 *
 *   FIELD STRENGTHS (F_μν from A_μ):
 *     E = -∇φ - ∂_t A  (electric field from scalar and vector potential time derivative)
 *     B = ∇×A           (magnetic field from vector potential curl)
 *
 *   2D SIMPLIFICATION:
 *     B = B_z ẑ = (∂_x A_y - ∂_y A_x) ẑ
 *     ∇×E = (0, 0, ∂_x E_y - ∂_y E_x)  (curl reduces to scalar in 2D)
 *
 *   KEY PHYSICS INSIGHT:
 *     For smooth phase θ with ∂_t(∇θ) ≈ ∇(∂_t θ):
 *       E ≈ -∇(∂_t θ) - ∂_t(∇θ) ≈ 0  (fields vanish in smooth regions)
 *
 *     Non-zero EM fields ONLY at topological defects:
 *       - Vortex cores (where θ has 2π winding)
 *       - Phase singularities (unbounded gradients)
 *       - Temporal discontinuities
 *
 * Numerical Methods:
 *   - Spatial derivatives: Centered finite differences (2nd order O(h²))
 *     ∂_x f ≈ (f[i+1,j] - f[i-1,j]) / (2dx)
 *   - Temporal derivatives: Backward difference (1st order O(dt))
 *     ∂_t f ≈ (f(t) - f(t-dt)) / dt
 *   - Boundary conditions: Periodic (wraparound indexing)
 *
 * Purpose: Phase 5 (Scenario 2.6B) - Electromagnetic Coupling
 *   - Extract gauge potential A_μ from phase gradients
 *   - Compute field strengths E, B for minimal coupling
 *   - Measure field energy density ∫(E² + B²) dV
 *   - Validate Maxwell equations on discrete grid
 *   - Test if EM coupling explains particle dynamics
 *   - Extract effective fine structure constant α_eff = q²/(ħc)
 *
 * Usage Pattern:
 *   1. computeFromPhase() → computes potentials and fields in one call
 *   2. Optional: computeFieldStrengths() for separate phase/field updates
 *   3. computeFieldEnergy() for energy analysis
 *   4. computeDiagnostics() for max/avg field statistics
 *   5. computePoyntingVector() for energy flux analysis
 *   6. computeChargeCurrent() for source term computation
 *
 * Example:
 *   ```cpp
 *   // Initialize phase fields
 *   Eigen::MatrixXd theta_t0 = ...; // phase at t
 *   Eigen::MatrixXd theta_t1 = ...; // phase at t+dt
 *   double dx = 0.1, dy = 0.1, dt = 0.01;
 *
 *   // Extract EM fields
 *   auto fields = EMFieldComputer::computeFromPhase(theta_t1, theta_t0, dx, dy, dt);
 *
 *   // Analyze field energy
 *   double U_field = EMFieldComputer::computeFieldEnergy(fields, dx, dy);
 *   std::cout << "Total field energy: " << U_field << std::endl;
 *
 *   // Check Poynting vector (energy flux)
 *   auto poynting = EMFieldComputer::computePoyntingVector(fields, dx, dy);
 *
 *   // Validate Maxwell equations (with charge/current sources)
 *   Eigen::MatrixXd rho = ...; // charge density
 *   Eigen::MatrixXd J_x = ..., J_y = ...;
 *   bool is_valid = EMFieldComputer::validateMaxwellEquations(
 *       fields, rho, J_x, J_y, dx, dy);
 *   ```
 */

#pragma once

#include <eigen3/Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

class EMFieldComputer {
public:
    // ========================================================================
    // DATA STRUCTURES
    // ========================================================================

    /**
     * Electromagnetic field data structure
     *
     * Stores computed gauge potential and field strengths on discrete grid.
     * All fields are [Nx × Ny] matrices corresponding to grid points.
     *
     * Units: Natural (Planck) units where ħ = c = G = Δx = 1
     *        Fields normalized to order parameter scale
     */
    struct EMFields {
        // GAUGE POTENTIAL: A_μ = ∂_μ θ
        // ────────────────────────────────
        Eigen::MatrixXd phi;        // Scalar potential φ = ∂_t θ [Nx × Ny]
                                    // Time evolution rate of phase
                                    // Units: [radians/time]

        Eigen::MatrixXd A_x;        // Vector potential x-component A_x = ∂_x θ [Nx × Ny]
                                    // Spatial gradient in x direction
                                    // Units: [radians/length]

        Eigen::MatrixXd A_y;        // Vector potential y-component A_y = ∂_y θ [Nx × Ny]
                                    // Spatial gradient in y direction
                                    // Units: [radians/length]

        // FIELD STRENGTHS: F_μν computed from A_μ
        // ─────────────────────────────────────────
        Eigen::MatrixXd E_x;        // Electric field x-component E_x [Nx × Ny]
                                    // Computed: E_x = -∂_x φ - ∂_t A_x
                                    // Units: [field_strength]

        Eigen::MatrixXd E_y;        // Electric field y-component E_y [Nx × Ny]
                                    // Computed: E_y = -∂_y φ - ∂_t A_y
                                    // Units: [field_strength]

        Eigen::MatrixXd B_z;        // Magnetic field z-component (2D only) [Nx × Ny]
                                    // Computed: B_z = ∂_x A_y - ∂_y A_x
                                    // Units: [field_strength]

        // DIAGNOSTICS: Summary statistics
        // ────────────────────────────────
        double total_field_energy;  // Total electromagnetic energy ∫(E² + B²)/(8π) dA
                                    // Units: [energy]

        double max_E;               // Peak electric field magnitude max(|E|)
                                    // Units: [field_strength]

        double max_B;               // Peak magnetic field magnitude max(|B_z|)
                                    // Units: [field_strength]

        double avg_E;               // Average electric field magnitude <|E|>
                                    // Units: [field_strength]

        double avg_B;               // Average magnetic field magnitude <|B_z|>
                                    // Units: [field_strength]

        // TEMPORAL TRACKING (for time derivatives)
        // ─────────────────────────────────────────
        bool has_previous;          // Whether fields have been computed before
                                    // (allows accurate time derivatives)

        // Constructor
        EMFields(int Nx, int Ny)
            : phi(Nx, Ny), A_x(Nx, Ny), A_y(Nx, Ny),
              E_x(Nx, Ny), E_y(Nx, Ny), B_z(Nx, Ny),
              total_field_energy(0.0), max_E(0.0), max_B(0.0),
              avg_E(0.0), avg_B(0.0), has_previous(false)
        {
            phi.setZero();
            A_x.setZero();
            A_y.setZero();
            E_x.setZero();
            E_y.setZero();
            B_z.setZero();
        }
    };

    /**
     * Poynting vector components: S = (1/(4π)) E × B
     *
     * Represents energy flux (power per unit area) in electromagnetic field.
     * In 2D with B = B_z ẑ:
     *   S_x = (1/(4π)) E_y B_z
     *   S_y = (1/(4π)) (-E_x B_z)
     *   |S| = (1/(4π)) |E| |B| (for orthogonal E, B)
     */
    struct PoyntingVector {
        Eigen::MatrixXd S_x;        // x-component of Poynting vector [Nx × Ny]
        Eigen::MatrixXd S_y;        // y-component of Poynting vector [Nx × Ny]
        double total_flux;          // Total energy flux ∫|S| dA
        double max_flux;            // Peak energy flux magnitude
        double avg_flux;            // Average energy flux magnitude

        PoyntingVector(int Nx, int Ny)
            : S_x(Nx, Ny), S_y(Nx, Ny),
              total_flux(0.0), max_flux(0.0), avg_flux(0.0)
        {
            S_x.setZero();
            S_y.setZero();
        }
    };

    /**
     * Maxwell equation validation summary
     *
     * Contains residuals for each Maxwell equation showing how well
     * the computed fields satisfy fundamental electromagnetic laws.
     */
    struct MaxwellValidation {
        // RMS residuals (should be O(10^-2) for well-resolved fields)
        double gauss_law_rms;       // |∇·E - 4πρ|_rms
        double ampere_law_rms;      // |∇×B - (4π/c)J - (1/c)∂_t E|_rms
        double faraday_law_rms;     // |∇×E + (1/c)∂_t B|_rms
        double no_monopole_rms;     // |∇·B|_rms (should be ~ 0)
        double coulomb_gauge_rms;   // |∇·A|_rms (constraint in Coulomb gauge)

        // Peak residuals
        double gauss_law_max;
        double ampere_law_max;
        double faraday_law_max;

        // Overall quality metrics
        bool all_equations_satisfied; // true if all RMS residuals < threshold
        double overall_error;         // Weighted sum of residuals

        MaxwellValidation()
            : gauss_law_rms(0.0), ampere_law_rms(0.0), faraday_law_rms(0.0),
              no_monopole_rms(0.0), coulomb_gauge_rms(0.0),
              gauss_law_max(0.0), ampere_law_max(0.0), faraday_law_max(0.0),
              all_equations_satisfied(false), overall_error(0.0)
        {}
    };

    // ========================================================================
    // PUBLIC METHODS: CORE COMPUTATION
    // ========================================================================

    /**
     * Extract electromagnetic potential and fields from Kuramoto phase
     *
     * Main entry point: computes gauge potential A_μ from phase gradients,
     * then derives field strengths E, B.
     *
     * Physics:
     *   1. A_μ = ∂_μ θ  (gauge potential from phase derivatives)
     *   2. E = -∇φ - ∂_t A  (electric field)
     *   3. B = ∇×A  (magnetic field)
     *   4. U = ∫(E²+B²)/(8π) dV  (field energy)
     *
     * Implementation:
     *   - Spatial derivatives: Centered differences with periodic BC
     *   - Temporal derivative: Backward difference
     *   - Automatic diagnostics computation
     *
     * @param[in] theta_current    Current phase field θ(t) [Nx × Ny]
     *                              Range: [0, 2π) (radians)
     *                              Represents local synchronization phase
     *
     * @param[in] theta_previous   Previous phase field θ(t - dt) [Nx × Ny]
     *                              Required for accurate temporal derivatives
     *                              If unavailable, use backward-only estimate
     *
     * @param[in] dx               Grid spacing in x-direction [Planck lengths]
     *                              Typical: 0.1 to 1.0 ℓ_P
     *
     * @param[in] dy               Grid spacing in y-direction [Planck lengths]
     *
     * @param[in] dt               Timestep [Planck time]
     *                              Used for temporal derivative ∂_t θ ≈ Δθ/Δt
     *
     * @return EMFields            Struct with potentials A_μ, fields E,B,
     *                              and diagnostics (energy, max values, etc.)
     *
     * Complexity: O(Nx × Ny) - single pass over grid
     * Accuracy: O(h²) spatial, O(dt) temporal
     *
     * Notes:
     *   - Assumes periodic boundaries (phase wraps at domain edges)
     *   - First timestep may have large temporal derivative error
     *   - For vortex detection, look for peaks in E and B
     *   - For smooth background: E ≈ B ≈ 0 (should be ~10^-3)
     */
    static EMFields computeFromPhase(
        const Eigen::MatrixXd& theta_current,
        const Eigen::MatrixXd& theta_previous,
        double dx, double dy, double dt);

    // ========================================================================
    // PUBLIC METHODS: FIELD STRENGTH COMPUTATION
    // ========================================================================

    /**
     * Compute electric and magnetic field strengths from gauge potential
     *
     * Physics:
     *   E = -∇φ - ∂_t A
     *     - ∇φ term: static potential gradient (e.g., from temperature)
     *     - ∂_t A term: induction from changing vector potential
     *   B = ∇×A = ∂_x A_y - ∂_y A_x  (2D curl)
     *
     * Implementation:
     *   - Uses already-computed potentials in EMFields struct
     *   - Computes field strengths in-place on input struct
     *   - Assumes periodic boundaries for spatial derivatives
     *
     * Mathematical Detail:
     *   For smooth phase (∂_t(∇θ) ≈ ∇(∂_t θ)):
     *     E = -∇(∂_t θ) - ∂_t(∇θ)
     *       ≈ -∇(∂_t θ) - ∂_t(∇θ)  (exact for smooth fields)
     *       ≈ 0  (if time and space derivatives commute)
     *
     *   Large E, B indicate:
     *     - Vortex cores (where θ is singular)
     *     - Rapid phase changes (temporal shock)
     *     - Fine structures (unresolved gradients)
     *
     * @param[in,out] fields       EMFields struct with A_μ already computed
     *                              Updated in-place: E_x, E_y, B_z filled
     *
     * @param[in] dx               Grid spacing in x-direction
     *
     * @param[in] dy               Grid spacing in y-direction
     *
     * Complexity: O(Nx × Ny) - 4 spatial derivatives computed
     */
    static void computeFieldStrengths(EMFields& fields, double dx, double dy);

    // ========================================================================
    // PUBLIC METHODS: ENERGY & FLUX
    // ========================================================================

    /**
     * Compute total electromagnetic field energy
     *
     * Physics (Gaussian units, natural units ħ=c=1):
     *   Energy density:   u(r) = (E² + B²) / (8π)
     *   Total energy:     U = ∫ u(r) dV = ∫ (E² + B²)/(8π) dx dy
     *
     * Interpretation:
     *   - E² term: electrostatic energy from charge/potential
     *   - B² term: magnetic energy from currents/loops
     *   - Factor 1/(8π): Gaussian unit conversion
     *
     * Numerical Integration:
     *   U ≈ Σ_i,j u(i,j) · dx·dy  (trapezoid rule, O(h²) accurate)
     *
     * @param[in] fields           EMFields struct with E, B computed
     *
     * @param[in] dx               Grid spacing x (for volume element)
     *
     * @param[in] dy               Grid spacing y (for volume element)
     *
     * @return Total energy U [energy units]
     *
     * Complexity: O(Nx × Ny) - single pass over all grid points
     *
     * Typical Values:
     *   - Vacuum background: U ~ 10^-6
     *   - Single vortex: U ~ 0.01 to 0.1
     *   - Vortex pair: U ~ 0.02 to 0.2
     */
    static double computeFieldEnergy(const EMFields& fields, double dx, double dy);

    /**
     * Compute Poynting vector (energy flux)
     *
     * Physics:
     *   Poynting vector: S = (1/(4π)) E × B  [energy/(area·time)]
     *   In 2D with B = B_z ẑ:
     *     S_x = (1/(4π)) E_y · B_z
     *     S_y = (1/(4π)) (-E_x · B_z)
     *     |S| = (1/(4π)) |E| |B| sin(angle between E and B)
     *
     * Interpretation:
     *   - |S| = power per unit area flowing through surface
     *   - Direction: perpendicular to both E and B
     *   - ∫∮ S·dA = power radiated from volume
     *   - ∂u/∂t + ∇·S = -J·E  (Poynting theorem, energy conservation)
     *
     * @param[in] fields           EMFields struct with E, B computed
     *
     * @param[in] dx               Grid spacing x
     *
     * @param[in] dy               Grid spacing y
     *
     * @return PoyntingVector struct with S_x, S_y and statistics
     *
     * Complexity: O(Nx × Ny)
     *
     * Notes:
     *   - S = 0 in smooth background regions (E ≈ 0 or B ≈ 0)
     *   - Large |S| at vortex boundaries (energy flowing in/out)
     *   - Can detect vortex motion by tracking flux patterns
     */
    static PoyntingVector computePoyntingVector(const EMFields& fields,
                                                double dx, double dy);

    // ========================================================================
    // PUBLIC METHODS: FORCE & SOURCE TERMS
    // ========================================================================

    /**
     * Compute Lorentz force density at grid point
     *
     * Physics:
     *   F = ρE + J×B
     *
     * where:
     *   ρ = charge density (e.g., from spinor norm: Ψ†Ψ)
     *   J = current density (e.g., from spinor: Ψ†α Ψ)
     *   E = electric field
     *   B = magnetic field
     *
     * In 2D with B = B_z ẑ:
     *   J×B = (J_x, J_y, 0) × (0, 0, B_z)
     *       = (J_y·B_z, -J_x·B_z, 0)
     *
     * Used For:
     *   - Minimal coupling: ∂_t p = F_Lorentz + other forces
     *   - Validation: Does EM coupling explain observed dynamics?
     *   - Discovery: Extract effective coupling strength α_eff
     *
     * @param[in] fields           EMFields struct with E, B
     *
     * @param[in] charge_density   Charge density ρ at point (scalar)
     *                              Typical: 0.1 to 1.0
     *
     * @param[in] current_density  Current density J = (J_x, J_y) at point
     *                              Typical: 0.01 to 0.1 (usually small)
     *
     * @param[in] ix, iy           Grid indices for field lookup
     *
     * @return Lorentz force F = (F_x, F_y) at point [force units]
     *
     * Complexity: O(1) - single point calculation
     */
    static Eigen::Vector2d computeLorentzForce(
        const EMFields& fields,
        double charge_density,
        const Eigen::Vector2d& current_density,
        int ix, int iy);

    /**
     * Compute charge density from matter field
     *
     * Physics:
     *   For spinor field Ψ(r) = (ψ_1, ψ_2, ψ_3, ψ_4):
     *   ρ(r) = Ψ†(r)Ψ(r) = Σ_{α=1}^4 |ψ_α(r)|²
     *
     * @param[in] spinor_field     4-component spinor [4 × N_points]
     *                              Interleaved: [ψ_1, ψ_2, ψ_3, ψ_4, ...]
     *
     * @param[in] Nx, Ny           Grid dimensions
     *
     * @return Charge density ρ [Nx × Ny]
     *
     * Complexity: O(Nx × Ny)
     */
    static Eigen::MatrixXd computeChargeDensity(
        const std::vector<std::complex<double>>& spinor_field,
        int Nx, int Ny);

    /**
     * Compute current density from matter field
     *
     * Physics (Dirac equation):
     *   J^μ = Ψ†γ^μ Ψ  (probability current / vector current)
     *
     * For 2D spatial part (α_x, α_y Pauli matrices):
     *   J_x(r) = Re[Ψ†(r) α_x Ψ(r)]
     *   J_y(r) = Re[Ψ†(r) α_y Ψ(r)]
     *
     * where α_x = [[0,1],[1,0]], α_y = [[0,-i],[i,0]] (Pauli matrices)
     *
     * @param[in] spinor_field     4-component spinor
     *
     * @param[in] Nx, Ny           Grid dimensions
     *
     * @param[out] J_x             Current density x-component [Nx × Ny]
     *
     * @param[out] J_y             Current density y-component [Nx × Ny]
     *
     * Complexity: O(Nx × Ny)
     */
    static void computeCurrentDensity(
        const std::vector<std::complex<double>>& spinor_field,
        int Nx, int Ny,
        Eigen::MatrixXd& J_x,
        Eigen::MatrixXd& J_y);

    // ========================================================================
    // PUBLIC METHODS: VALIDATION & DIAGNOSTICS
    // ========================================================================

    /**
     * Compute field diagnostics (max, average, statistics)
     *
     * Updates the diagnostic fields in EMFields struct:
     *   - max_E, max_B: peak field strengths
     *   - avg_E, avg_B: average field magnitudes
     *   - Optionally: correlation with known defects
     *
     * Used To:
     *   - Monitor field evolution
     *   - Detect vortex creation/annihilation
     *   - Quantify noise in potential extraction
     *
     * @param[in,out] fields       EMFields struct to update diagnostics
     *
     * Complexity: O(Nx × Ny) - single pass over all points
     */
    static void computeDiagnostics(EMFields& fields);

    /**
     * Validate Maxwell equations on discrete grid
     *
     * Physics:
     *   ∇·E = 4πρ           (Gauss's law)
     *   ∇·B = 0             (No monopoles - automatic in 2D)
     *   ∇×E = -(1/c)∂_t B   (Faraday's law)
     *   ∇×B = (4π/c)J + (1/c)∂_t E   (Ampere-Maxwell law)
     *   ∇·A = 0             (Coulomb gauge constraint)
     *
     * Returns residuals showing how well computed fields satisfy laws.
     * Small residuals (RMS < 0.01) indicate well-resolved physics.
     *
     * @param[in] fields           Current EM fields
     *
     * @param[in] fields_prev      Previous EM fields (for time derivatives)
     *
     * @param[in] charge_density   Charge density ρ [Nx × Ny]
     *                              From spinor: ρ = Ψ†Ψ
     *
     * @param[in] current_x, current_y  Current density J [Nx × Ny]
     *                              From spinor: J = Ψ†α Ψ
     *
     * @param[in] dx, dy           Grid spacing
     *
     * @param[in] dt               Timestep (for time derivatives)
     *
     * @return MaxwellValidation struct with residuals and quality metrics
     *
     * Complexity: O(Nx × Ny) - multiple passes (5 differential operators)
     *
     * Interpretation:
     *   - gauss_law_rms ~ 10^-3: Good agreement with source term
     *   - ampere_law_rms ~ 10^-2: Acceptable Ampere's law
     *   - faraday_law_rms ~ 10^-2: Acceptable Faraday's law
     *   - couloumb_gauge_rms ~ 0: Should be near zero (constraint)
     *
     * If residuals are large:
     *     - Check grid resolution (dx too large?)
     *     - Check timestep (dt too large?)
     *     - Look for unresolved vortex cores
     *     - May indicate numerical instability
     */
    static MaxwellValidation validateMaxwellEquations(
        const EMFields& fields,
        const EMFields& fields_prev,
        const Eigen::MatrixXd& charge_density,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y,
        double dx, double dy, double dt);

    /**
     * Compute energy conservation check
     *
     * Physics (Poynting's theorem):
     *   ∂u/∂t + ∇·S = -J·E
     *
     * where:
     *   u = (E² + B²)/(8π) is field energy density
     *   S = (1/(4π)) E×B is Poynting vector (energy flux)
     *   J·E is power dissipation/work by charges
     *
     * LHS: Rate of field energy change + outward flux
     * RHS: Work done by electromagnetic force on charges
     *
     * @param[in] fields_curr      EM fields at time t
     *
     * @param[in] fields_prev      EM fields at time t-dt
     *
     * @param[in] current_x, current_y  Current density J
     *
     * @param[in] E_prev_x, E_prev_y    Previous electric field
     *
     * @param[in] dt, dx, dy       Time and space steps
     *
     * @return RMS residual of Poynting's theorem
     *         Should be small if energy conservation holds
     *
     * Complexity: O(Nx × Ny)
     */
    static double validateEnergyConservation(
        const EMFields& fields_curr,
        const EMFields& fields_prev,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y,
        const Eigen::MatrixXd& E_prev_x,
        const Eigen::MatrixXd& E_prev_y,
        double dt, double dx, double dy);

    // ========================================================================
    // PUBLIC METHODS: FINE STRUCTURE CONSTANT EXTRACTION
    // ========================================================================

    /**
     * Extract effective fine structure constant from field coupling
     *
     * Physics:
     *   Fine structure constant: α = e²/(ħc) in SI units
     *   In natural units (ħ = c = 1): α = e²
     *   Measured value: α ≈ 1/137
     *
     * Hypothesis:
     *   EM coupling strength in SMFT encodes α_eff = effective fine structure
     *   Can be measured from ratio of EM force to kinetic energy
     *   Or from coupling strength in Dirac equation: ∂_t Ψ = (-iγ^μ A_μ) Ψ + ...
     *
     * Method:
     *   α_eff ≈ |F_EM|² / (kinetic energy)  (order-of-magnitude estimate)
     *   or: α_eff ≈ |E·B| / (particle energy)²
     *
     * @param[in] fields           EM fields with E, B computed
     *
     * @param[in] kinetic_energy   Average kinetic energy of matter field
     *                              From Dirac: KE ~ p²/m for free particles
     *
     * @return Effective fine structure constant α_eff
     *
     * Notes:
     *   - Purely experimental/exploratory at this stage
     *   - Will be refined based on simulation results
     *   - Look for consistency: α_eff should be ~constant across runs
     *   - If α_eff ≈ 1/137: potential deep physics discovered!
     */
    static double computeEffectiveAlpha(
        const EMFields& fields,
        double kinetic_energy);

    // ========================================================================
    // PRIVATE METHODS: NUMERICAL UTILITIES
    // ========================================================================

private:

    /**
     * Apply periodic boundary conditions for finite difference
     *
     * Wraps indices to [0, N) range, enabling center-of-grid field access
     * at domain boundaries.
     *
     * @param[in] field            2D field [Nx × Ny]
     *
     * @param[in] i, j             Base indices
     *
     * @param[in] offset_i, offset_j  Index offsets (+1, -1, etc.)
     *
     * @return field(i_wrapped, j_wrapped) where indices wrap periodically
     *
     * Complexity: O(1)
     */
    static double getPeriodicValue(
        const Eigen::MatrixXd& field,
        int i, int j,
        int offset_i, int offset_j);

    /**
     * Compute spatial derivative with centered finite differences
     *
     * Physics:
     *   ∂_x f ≈ (f[i+1,j] - f[i-1,j]) / (2dx)   [O(h²) accurate]
     *   ∂_y f ≈ (f[i,j+1] - f[i,j-1]) / (2dy)
     *
     * Uses periodic boundary conditions for wraparound at edges.
     *
     * Accuracy:
     *   Truncation error: O(h²) where h = dx or dy
     *   For smooth fields: error ~ h²·(d³f/dx³)
     *   For vortex cores: error can be larger (unresolved singularities)
     *
     * @param[in] field            2D field to differentiate
     *
     * @param[in] dx_or_dy         Grid spacing (dx for x-deriv, dy for y-deriv)
     *
     * @param[in] direction        0 for x-derivative, 1 for y-derivative
     *
     * @return Derivative field ∂f/∂(x or y) [Nx × Ny]
     *
     * Complexity: O(Nx × Ny)
     */
    static Eigen::MatrixXd computeSpatialDerivative(
        const Eigen::MatrixXd& field,
        double dx_or_dy,
        int direction);

    /**
     * Compute divergence of vector field with periodic boundaries
     *
     * Physics:
     *   ∇·F = ∂_x F_x + ∂_y F_y
     *
     * Uses centered finite differences on both components.
     *
     * @param[in] F_x, F_y        Vector field components [Nx × Ny]
     *
     * @param[in] dx, dy          Grid spacing
     *
     * @return Divergence field ∇·F [Nx × Ny]
     *
     * Complexity: O(Nx × Ny)
     */
    static Eigen::MatrixXd computeDivergence(
        const Eigen::MatrixXd& F_x,
        const Eigen::MatrixXd& F_y,
        double dx, double dy);

    /**
     * Compute curl z-component (2D only)
     *
     * Physics:
     *   (∇×F)_z = ∂_x F_y - ∂_y F_x
     *
     * In 2D, only z-component survives; x,y components are zero.
     *
     * @param[in] F_x, F_y        Vector field components [Nx × Ny]
     *
     * @param[in] dx, dy          Grid spacing
     *
     * @return Curl z-component (∇×F)_z [Nx × Ny]
     *
     * Complexity: O(Nx × Ny)
     */
    static Eigen::MatrixXd computeCurl(
        const Eigen::MatrixXd& F_x,
        const Eigen::MatrixXd& F_y,
        double dx, double dy);

    /**
     * Compute RMS (root-mean-square) of field
     *
     * RMS(f) = √(<f²>) = √(Σ_i |f_i|² / N)
     *
     * Used for error metrics (residuals of Maxwell equations).
     *
     * @param[in] field            Field values [Nx × Ny]
     *
     * @return RMS = sqrt(<field²>)
     *
     * Complexity: O(Nx × Ny)
     */
    static double computeRMS(const Eigen::MatrixXd& field);

    /**
     * Compute Pearson correlation coefficient
     *
     * Physics:
     *   ρ(A,B) = <(A - <A>)(B - <B>)> / (σ_A σ_B)
     *   Range: [-1, +1]
     *
     * Used for comparing predicted vs observed forces.
     *
     * @param[in] field_A, field_B  Fields to correlate [Nx × Ny]
     *
     * @return Correlation ρ ∈ [-1, 1]
     *
     * Complexity: O(Nx × Ny)
     */
    static double computeCorrelation(
        const Eigen::MatrixXd& field_A,
        const Eigen::MatrixXd& field_B);

}; // class EMFieldComputer
