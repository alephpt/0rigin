/**
 * EMFieldComputer.h
 *
 * LOW-LEVEL ELECTROMAGNETIC FIELD EXTRACTION FROM KURAMOTO PHASE
 *
 * ARCHITECTURE NOTE:
 *   This class provides ONLY low-level field computation from phase gradients.
 *   For higher-level validation, observables, and Maxwell equation checking,
 *   use EMObservables.h which builds on top of this class.
 *
 *   Layered Design:
 *     EMFieldComputer (low-level)  → Field extraction from phase
 *     EMObservables (high-level)   → Validation & analysis using EMFieldComputer
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
 *   - Compute Poynting vector (energy flux)
 *   - Compute Lorentz force at grid points
 *
 * Usage Pattern:
 *   1. computeFromPhase() → computes potentials and fields in one call
 *   2. Optional: computeFieldStrengths() for separate phase/field updates
 *   3. computeFieldEnergy() for energy analysis
 *   4. computeDiagnostics() for max/avg field statistics
 *   5. computePoyntingVector() for energy flux analysis
 *   6. computeLorentzForce() for force at grid point
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
 *   // For validation & observables, use EMObservables.h
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

    /**
     * Convenience wrapper: Extract EM fields from std::vector<float> phase data
     *
     * Same physics as Eigen version, but accepts CPU vectors (row-major layout).
     * Useful for integration with SMFTEngine which uses std::vector<float> internally.
     *
     * @param[in] theta_current    Current phase field [Nx*Ny flat vector]
     * @param[in] theta_previous   Previous phase field [Nx*Ny flat vector]
     * @param[in] Nx               Grid width
     * @param[in] Ny               Grid height
     * @param[in] dx               Grid spacing x
     * @param[in] dy               Grid spacing y
     * @param[in] dt               Timestep
     *
     * @return EMFields            Potentials and fields (stored as Eigen matrices)
     */
    static EMFields computeFromPhase(
        const std::vector<float>& theta_current,
        const std::vector<float>& theta_previous,
        int Nx, int Ny,
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
     * Note: For charge/current density computation from spinor fields,
     *       use EMObservables::computeChargeDensity() and
     *       EMObservables::computeCurrentDensity()
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

    // ========================================================================
    // PUBLIC METHODS: DIAGNOSTICS
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
     * Note: For Maxwell equation validation, energy conservation checks,
     *       and fine structure constant extraction, use EMObservables class.
     *
     * @param[in,out] fields       EMFields struct to update diagnostics
     *
     * Complexity: O(Nx × Ny) - single pass over all points
     */
    static void computeDiagnostics(EMFields& fields);

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

}; // class EMFieldComputer
