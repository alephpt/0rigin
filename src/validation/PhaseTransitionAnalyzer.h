/**
 * PhaseTransitionAnalyzer.h
 *
 * Analysis tools for Test 3.3: Noise-Driven Phase Transition
 *
 * Physics: Order parameter ⟨R⟩ undergoes phase transition at critical noise σ_c.
 * Expect ⟨R⟩ ∝ (σ_c - σ)^β with β ≈ 0.125 (2D Ising universality class).
 *
 * Methods:
 * - findCriticalPoint: Locate σ_c from noise scan data
 * - fitCriticalExponent: Extract β from power-law fit near σ_c
 * - computeOrderParameter: Spatial average ⟨R⟩ with optional variance
 */

#pragma once

#include <vector>
#include <string>
#include <cstdint>

class PhaseTransitionAnalyzer {
public:
    /**
     * Data point for noise scan
     */
    struct DataPoint {
        float sigma;      // Noise amplitude
        float R_mean;     // Spatial average ⟨R⟩
        float R_variance; // Variance σ²(R)
        float R_min;      // Min R value
        float R_max;      // Max R value
    };

    /**
     * Critical point fit result
     */
    struct CriticalPointResult {
        float sigma_c;           // Critical noise (estimated)
        float sigma_c_error;     // Uncertainty in σ_c
        float beta;              // Critical exponent
        float beta_error;        // Uncertainty in β
        float R2;                // Goodness of fit (R²)
        bool success;            // Fit converged
        std::string message;     // Status/error message
    };

    /**
     * Compute spatial average order parameter ⟨R⟩ from R field
     *
     * @param R_field Synchronization field R(x,y)
     * @param Nx Grid width
     * @param Ny Grid height
     * @return Spatial average ⟨R⟩ = (1/N)·Σ R_i
     */
    static float computeOrderParameter(const std::vector<float>& R_field,
                                        uint32_t Nx, uint32_t Ny);

    /**
     * Compute order parameter statistics (mean, variance, min, max)
     *
     * @param R_field Synchronization field R(x,y)
     * @param Nx Grid width
     * @param Ny Grid height
     * @param mean Output: spatial average ⟨R⟩
     * @param variance Output: variance σ²(R)
     * @param min_val Output: min R value
     * @param max_val Output: max R value
     */
    static void computeOrderParameterStats(const std::vector<float>& R_field,
                                             uint32_t Nx, uint32_t Ny,
                                             float& mean, float& variance,
                                             float& min_val, float& max_val);

    /**
     * Find critical point σ_c from noise scan data
     *
     * Method: Locate maximum derivative |d⟨R⟩/dσ|
     * This marks the steepest drop in order parameter
     *
     * @param data Noise scan data points (sorted by increasing σ)
     * @return Critical point fit result
     */
    static CriticalPointResult findCriticalPoint(const std::vector<DataPoint>& data);

    /**
     * Fit critical exponent β from power-law near σ_c
     *
     * Model: ⟨R⟩ ∝ (σ_c - σ)^β for σ < σ_c
     *
     * Uses log-log linear regression:
     * log(⟨R⟩) = β·log(σ_c - σ) + const
     *
     * @param data Noise scan data points
     * @param sigma_c Critical noise (from findCriticalPoint)
     * @param fit_range_fraction Fraction of data below σ_c to use (default: 0.5)
     * @return Critical exponent β with error estimate
     */
    static CriticalPointResult fitCriticalExponent(const std::vector<DataPoint>& data,
                                                     float sigma_c,
                                                     float fit_range_fraction = 0.5f);

    /**
     * Detect phase transition type (first-order vs second-order)
     *
     * First-order: Discontinuous jump in ⟨R⟩
     * Second-order: Continuous but singular d⟨R⟩/dσ
     *
     * @param data Noise scan data points
     * @return "first_order", "second_order", or "none"
     */
    static std::string detectTransitionType(const std::vector<DataPoint>& data);

    /**
     * Generate validation report for phase transition test
     *
     * @param data Noise scan data points
     * @param output_path Path to write report file
     * @return true if report written successfully
     */
    static bool writeReport(const std::vector<DataPoint>& data,
                             const std::string& output_path);

private:
    /**
     * Compute numerical derivative d⟨R⟩/dσ using centered differences
     *
     * @param data Noise scan data (sorted by σ)
     * @return Derivative at each point
     */
    static std::vector<float> computeDerivative(const std::vector<DataPoint>& data);

    /**
     * Linear regression: y = a·x + b
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @param a Output: slope
     * @param b Output: intercept
     * @param R2 Output: coefficient of determination
     * @return true if regression successful
     */
    static bool linearRegression(const std::vector<float>& x,
                                   const std::vector<float>& y,
                                   float& a, float& b, float& R2);
};
