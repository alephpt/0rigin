#ifndef OBSERVABLE_COMPUTER_H
#define OBSERVABLE_COMPUTER_H

#include "DiracEvolution.h"
#include <complex>
#include <vector>
#include <map>
#include <string>

// Forward declaration
class TRDEngine;

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
        double time = 0.0;

        // Dirac observables
        double norm = 0.0;              // ||Ψ||² - should be ≈ 1.0
        double norm_error = 0.0;        // ||Ψ||² - 1.0
        double energy_total = 0.0;
        double energy_kinetic = 0.0;
        double energy_potential = 0.0;

        // Position expectation <Ψ|x|Ψ>, <Ψ|y|Ψ>
        std::complex<double> position_x{0.0, 0.0};
        std::complex<double> position_y{0.0, 0.0};

        // Momentum expectation <Ψ|p_x|Ψ>, <Ψ|p_y|Ψ>
        std::complex<double> momentum_x{0.0, 0.0};
        std::complex<double> momentum_y{0.0, 0.0};

        // Kuramoto sync field observables
        double R_avg = -999.0;  // SENTINEL: If this stays -999, assignment failed!
        double R_max = 0.0;
        double R_min = 0.0;
        double R_variance = 0.0;

        // EM field observables (Stückelberg gauge-restored)
        double EM_B_max = 0.0;          // Maximum magnetic field strength |B_z|
        double EM_B_rms = 0.0;          // RMS magnetic field
        double EM_energy = 0.0;         // Total EM field energy

        // Validation flags
        bool norm_valid = false;          // |norm_error| < tolerance
        bool energy_valid = false;        // |ΔE/E₀| < tolerance
    };

    // HACK: Global workaround for pointer corruption bug
    static thread_local Observables* g_result_hack;

    /**
     * Compute all observables for current state
     *
     * CRITICAL: Uses output parameter to avoid compiler bug where return-by-value
     * corrupts struct fields (R_avg was becoming 0 in caller despite correct value in function).
     *
     * @param result [OUT] Observables struct to populate
     * @param dirac Dirac evolution state
     * @param R_field Kuramoto sync field (size Nx*Ny)
     * @param delta Vacuum potential Δ
     * @param time Current simulation time
     * @param E0 Initial energy (for relative energy conservation check)
     * @param norm_tolerance Tolerance for ||Ψ||² ≈ 1 (default 1e-4)
     * @param energy_tolerance Tolerance for ΔE/E₀ (default 1e-2)
     * @param engine Optional TRDEngine pointer for EM observables (default nullptr)
     */
    static void compute(
        Observables* result,
        const DiracEvolution& dirac,
        const std::vector<double>& R_field,
        double delta,
        double time,
        double E0 = 0.0,
        double norm_tolerance = 1e-4,
        double energy_tolerance = 1e-2,
        const TRDEngine* engine = nullptr);

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

    /**
     * Compute sync field statistics
     * Returns {R_avg, R_max, R_min, R_variance}
     */
    static std::tuple<double, double, double, double> computeSyncFieldStats(
        const std::vector<double>& R_field);

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
