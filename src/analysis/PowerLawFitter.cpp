#include "analysis/PowerLawFitter.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

PowerLawFitter::FitResult PowerLawFitter::fit(
    const std::vector<double>& x,
    const std::vector<double>& y) {

    FitResult result;
    result.success = false;

    // Filter out invalid data
    auto [x_valid, y_valid] = filterValidData(x, y);

    if (x_valid.size() < 3) {
        // Need at least 3 points for meaningful fit
        return result;
    }

    // Transform to log-log space: log(y) = log(A) + α·log(x)
    std::vector<double> log_x, log_y;
    log_x.reserve(x_valid.size());
    log_y.reserve(y_valid.size());

    for (size_t i = 0; i < x_valid.size(); ++i) {
        if (x_valid[i] > 0 && y_valid[i] > 0) {
            log_x.push_back(std::log(x_valid[i]));
            log_y.push_back(std::log(y_valid[i]));
        }
    }

    if (log_x.size() < 3) {
        return result;
    }

    // Perform linear regression in log-log space
    auto [slope, intercept, r_squared] = linearRegression(log_x, log_y);

    // Extract power law parameters
    result.exponent = slope;           // α = slope in log-log plot
    result.prefactor = std::exp(intercept);  // A = exp(intercept)
    result.r_squared = r_squared;

    // Estimate uncertainty in exponent using residual standard error
    // σ_α ≈ s / sqrt(Σ(x - x̄)²) where s is residual standard error
    double mean_log_x = std::accumulate(log_x.begin(), log_x.end(), 0.0) / log_x.size();
    double sum_sq_dev = 0.0;
    double sum_residuals_sq = 0.0;

    for (size_t i = 0; i < log_x.size(); ++i) {
        double x_dev = log_x[i] - mean_log_x;
        sum_sq_dev += x_dev * x_dev;

        double y_pred = slope * log_x[i] + intercept;
        double residual = log_y[i] - y_pred;
        sum_residuals_sq += residual * residual;
    }

    int df = static_cast<int>(log_x.size()) - 2;  // Degrees of freedom
    if (df > 0 && sum_sq_dev > 1e-10) {
        double residual_std = std::sqrt(sum_residuals_sq / df);
        result.exponent_error = residual_std / std::sqrt(sum_sq_dev);
    } else {
        result.exponent_error = 0.0;
    }

    result.success = true;
    return result;
}

std::tuple<double, double, double> PowerLawFitter::linearRegression(
    const std::vector<double>& x,
    const std::vector<double>& y) {

    if (x.size() != y.size() || x.empty()) {
        return {0.0, 0.0, 0.0};
    }

    size_t n = x.size();

    // Compute means
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;

    // Compute slope and intercept using least squares
    double numerator = 0.0;
    double denominator = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        numerator += dx * dy;
        denominator += dx * dx;
    }

    if (std::abs(denominator) < 1e-10) {
        // Degenerate case: all x values identical
        return {0.0, mean_y, 0.0};
    }

    double slope = numerator / denominator;
    double intercept = mean_y - slope * mean_x;

    // Compute fitted values and R²
    std::vector<double> y_fit;
    y_fit.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        y_fit.push_back(slope * x[i] + intercept);
    }

    double r_squared = computeRSquared(y, y_fit);

    return {slope, intercept, r_squared};
}

double PowerLawFitter::computeRSquared(
    const std::vector<double>& y_data,
    const std::vector<double>& y_fit) {

    if (y_data.size() != y_fit.size() || y_data.empty()) {
        return 0.0;
    }

    size_t n = y_data.size();

    // Compute mean of observed data
    double mean_y = std::accumulate(y_data.begin(), y_data.end(), 0.0) / n;

    // Compute total sum of squares (SS_tot) and residual sum of squares (SS_res)
    double ss_tot = 0.0;
    double ss_res = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double residual = y_data[i] - y_fit[i];
        double deviation = y_data[i] - mean_y;

        ss_res += residual * residual;
        ss_tot += deviation * deviation;
    }

    if (ss_tot < 1e-10) {
        // Degenerate case: all y values identical
        return 1.0;
    }

    // R² = 1 - (SS_res / SS_tot)
    double r_squared = 1.0 - (ss_res / ss_tot);

    // Clamp to [0, 1] (negative R² can occur for bad fits)
    return std::max(0.0, std::min(1.0, r_squared));
}

std::pair<std::vector<double>, std::vector<double>> PowerLawFitter::filterValidData(
    const std::vector<double>& x,
    const std::vector<double>& y) {

    std::vector<double> x_valid, y_valid;

    if (x.size() != y.size()) {
        return {x_valid, y_valid};
    }

    for (size_t i = 0; i < x.size(); ++i) {
        bool x_valid_val = std::isfinite(x[i]) && x[i] > 0;
        bool y_valid_val = std::isfinite(y[i]) && y[i] != 0.0;

        if (x_valid_val && y_valid_val) {
            x_valid.push_back(x[i]);
            y_valid.push_back(std::abs(y[i]));  // Use magnitude for force
        }
    }

    return {x_valid, y_valid};
}
