/**
 * FiniteSizeScaling.h
 *
 * Finite-size scaling analysis tools for phase transition universality class determination
 *
 * Physics: Near a critical point, observables follow scaling forms:
 * - ⟨R⟩(σ, L) = L^(-β/ν) f_R((σ - σ_c)L^(1/ν))
 * - χ(σ, L) = L^(γ/ν) f_χ((σ - σ_c)L^(1/ν))
 * - U_L(σ, L) = g_U((σ - σ_c)L^(1/ν))
 *
 * Where f_R, f_χ, g_U are universal scaling functions.
 * Data from different L should collapse when plotted with proper scaling.
 */

#pragma once

#include <vector>
#include <map>
#include <string>
#include <cstdint>

class FiniteSizeScaling {
public:
    /**
     * Data point for FSS analysis
     */
    struct FSSDataPoint {
        int L;                // System size
        float sigma;          // Noise amplitude
        float R_mean;         // Spatial average ⟨R⟩
        float R_squared;      // ⟨R²⟩
        float R_fourth;       // ⟨R⁴⟩
        float R_variance;     // Var(R) = ⟨R²⟩ - ⟨R⟩²
        int num_samples;      // Number of samples averaged
    };

    /**
     * FSS analysis result
     */
    struct FSSResult {
        // Critical point
        float sigma_c;             // True thermodynamic critical point
        float sigma_c_error;       // Uncertainty in σ_c

        // Critical exponents
        float beta;                // Order parameter exponent ⟨R⟩ ~ (σ_c - σ)^β
        float beta_error;          // Uncertainty in β
        float nu;                  // Correlation length exponent ξ ~ |σ - σ_c|^(-ν)
        float nu_error;            // Uncertainty in ν
        float gamma;               // Susceptibility exponent χ ~ |σ - σ_c|^(-γ)
        float gamma_error;         // Uncertainty in γ
        float eta;                 // Anomalous dimension (from hyperscaling)
        float eta_error;           // Uncertainty in η

        // Quality metrics
        float R_squared;           // Quality of data collapse (R²)
        float chi_squared;         // χ² goodness of fit
        int degrees_of_freedom;    // Statistical DOF

        // Universality class identification
        std::string universality_class;  // "2D_Ising", "2D_XY", "Novel", etc.
        float confidence_level;          // Statistical confidence in classification

        // Status
        bool success;              // Analysis converged successfully
        std::string message;       // Status/error message
    };

    /**
     * Binder cumulant result
     */
    struct BinderResult {
        float U_L;                 // Binder cumulant value
        float U_L_error;           // Statistical error
        float sigma_c_crossing;    // Critical point from crossing
        float crossing_error;      // Uncertainty in crossing point
    };

    /**
     * Constructor
     */
    FiniteSizeScaling();

    /**
     * Add data point to FSS analysis
     *
     * @param L System size
     * @param sigma Noise amplitude
     * @param R_field Synchronization field R(x,y)
     */
    void addDataPoint(int L, float sigma, const std::vector<float>& R_field);

    /**
     * Add pre-computed data point
     *
     * @param point FSS data point with computed moments
     */
    void addDataPoint(const FSSDataPoint& point);

    /**
     * Compute Binder cumulant U_L = 1 - ⟨R⁴⟩/(3⟨R²⟩²)
     *
     * At criticality, U_L becomes L-independent (universal value)
     * Crossing of U_L curves for different L gives σ_c
     *
     * @param R_field Synchronization field R(x,y)
     * @return Binder cumulant value
     */
    static float computeBinderCumulant(const std::vector<float>& R_field);

    /**
     * Compute susceptibility χ = N·(⟨R²⟩ - ⟨R⟩²)
     *
     * @param R_field Synchronization field R(x,y)
     * @param L System size (assuming square grid L×L)
     * @return Magnetic susceptibility
     */
    static float computeSusceptibility(const std::vector<float>& R_field, int L);

    /**
     * Find critical point from Binder cumulant crossing
     *
     * Different system sizes have U_L curves that cross at σ_c
     *
     * @return Critical point and error estimate
     */
    BinderResult findCriticalPointFromBinder() const;

    /**
     * Perform full FSS analysis with critical exponents
     *
     * Fits data to scaling forms and extracts β, ν, γ, η
     * Uses optimization to find best data collapse
     *
     * @param sigma_c_guess Initial guess for critical point (optional)
     * @return FSS analysis results with exponents
     */
    FSSResult fitCriticalExponents(float sigma_c_guess = -1.0f) const;

    /**
     * Perform data collapse analysis
     *
     * Transforms data to scaling variables:
     * X = (σ - σ_c)·L^(1/ν)
     * Y = ⟨R⟩·L^(β/ν)
     *
     * Good collapse indicates correct exponents
     *
     * @param sigma_c Critical point
     * @param beta_over_nu Ratio β/ν
     * @param one_over_nu Exponent 1/ν
     * @return Quality of collapse (R² value)
     */
    float evaluateDataCollapse(float sigma_c, float beta_over_nu, float one_over_nu) const;

    /**
     * Identify universality class from exponents
     *
     * Compares measured exponents to known 2D universality classes
     *
     * @param result FSS analysis results
     * @return Updated result with universality classification
     */
    FSSResult identifyUniversalityClass(const FSSResult& result) const;

    /**
     * Write comprehensive FSS report
     *
     * Includes:
     * - Critical exponents with error bars
     * - Data collapse plots
     * - Universality class identification
     * - Quality metrics
     *
     * @param output_dir Output directory for report
     * @return true if report written successfully
     */
    bool writeDataCollapseReport(const std::string& output_dir) const;

    /**
     * Write data for external plotting (Python/gnuplot)
     *
     * @param output_dir Output directory
     * @return true if data written successfully
     */
    bool writeDataForPlotting(const std::string& output_dir) const;

    /**
     * Clear all stored data
     */
    void clear();

    /**
     * Get number of data points for given L
     */
    size_t getDataCount(int L) const;

    /**
     * Get total number of data points
     */
    size_t getTotalDataCount() const;

    /**
     * Get list of system sizes
     */
    std::vector<int> getSystemSizes() const;

private:
    // Data storage: map from L to vector of data points
    std::map<int, std::vector<FSSDataPoint>> data_by_L_;

    /**
     * Compute moments of R field
     *
     * @param R_field Synchronization field
     * @param mean Output: ⟨R⟩
     * @param second Output: ⟨R²⟩
     * @param fourth Output: ⟨R⁴⟩
     */
    static void computeMoments(const std::vector<float>& R_field,
                                float& mean, float& second, float& fourth);

    /**
     * Optimize scaling parameters for best data collapse
     *
     * Uses nonlinear optimization to find σ_c, β/ν, 1/ν
     *
     * @param initial_sigma_c Initial guess for σ_c
     * @return Optimized FSS result
     */
    FSSResult optimizeDataCollapse(float initial_sigma_c) const;

    /**
     * Bootstrap error estimation
     *
     * Resample data to estimate uncertainties in exponents
     *
     * @param result FSS result to add errors to
     * @param num_bootstrap Number of bootstrap samples
     * @return Result with error estimates
     */
    FSSResult bootstrapErrors(const FSSResult& result, int num_bootstrap = 100) const;

    /**
     * Check hyperscaling relations
     *
     * In 2D: 2β + γ = 2 (with α = 0)
     *        γ = ν(2 - η)
     *
     * @param result FSS results
     * @return true if hyperscaling satisfied within errors
     */
    bool checkHyperscaling(const FSSResult& result) const;

    /**
     * Linear regression for log-log fits
     */
    static bool linearRegression(const std::vector<float>& x,
                                  const std::vector<float>& y,
                                  float& slope, float& intercept,
                                  float& R2, float& slope_error);

    /**
     * Find intersection of two curves (for Binder crossing)
     */
    static bool findIntersection(const std::vector<float>& x1,
                                  const std::vector<float>& y1,
                                  const std::vector<float>& x2,
                                  const std::vector<float>& y2,
                                  float& x_cross, float& y_cross);
};