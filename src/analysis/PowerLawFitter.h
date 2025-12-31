#ifndef POWER_LAW_FITTER_H
#define POWER_LAW_FITTER_H

#include <vector>
#include <tuple>

/**
 * PowerLawFitter - Fit power law relationships F(x) = A·x^α
 *
 * Purpose: Extract scaling exponents from force-distance data
 * for Phase 3 Test 3.1 (Casimir-like Force).
 *
 * Uses log-log linear regression: log(F) = log(A) + α·log(x)
 *
 * Expected: F(d) ∝ d^(-2) for Casimir-like force
 */
class PowerLawFitter {
public:
    struct FitResult {
        double exponent;        // Power law exponent α
        double prefactor;       // Prefactor A
        double r_squared;       // Goodness of fit (R²)
        double exponent_error;  // Uncertainty in α
        bool success;           // Fit converged successfully
    };

    /**
     * Fit power law F(x) = A·x^α using log-log linear regression
     *
     * @param x Independent variable (e.g., separation distance d)
     * @param y Dependent variable (e.g., force F)
     * @return FitResult with exponent, prefactor, and quality metrics
     */
    static FitResult fit(const std::vector<double>& x, const std::vector<double>& y);

    /**
     * Perform linear regression on (x, y) data
     *
     * Fits: y = m·x + b
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @return {slope m, intercept b, R²}
     */
    static std::tuple<double, double, double> linearRegression(
        const std::vector<double>& x,
        const std::vector<double>& y);

    /**
     * Compute coefficient of determination R²
     *
     * Measures goodness of fit:
     * - R² = 1: perfect fit
     * - R² = 0: no correlation
     *
     * @param y_data Observed data
     * @param y_fit Fitted model predictions
     * @return R² value in [0, 1]
     */
    static double computeRSquared(
        const std::vector<double>& y_data,
        const std::vector<double>& y_fit);

    /**
     * Filter out invalid data points (NaN, Inf, zero, negative)
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @return {filtered_x, filtered_y}
     */
    static std::pair<std::vector<double>, std::vector<double>> filterValidData(
        const std::vector<double>& x,
        const std::vector<double>& y);
};

#endif // POWER_LAW_FITTER_H
