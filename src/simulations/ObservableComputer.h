#ifndef OBSERVABLE_COMPUTER_H
#define OBSERVABLE_COMPUTER_H

#include "DiracEvolution.h"
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
     * Compute all observables for current state
     *
     * @param dirac Dirac evolution state
     * @param R_field Kuramoto sync field (size Nx*Ny)
     * @param delta Vacuum potential Δ
     * @param time Current simulation time
     * @param E0 Initial energy (for relative energy conservation check)
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
