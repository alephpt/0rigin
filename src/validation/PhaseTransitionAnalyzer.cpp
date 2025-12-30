/**
 * PhaseTransitionAnalyzer.cpp
 *
 * Implementation of phase transition analysis tools
 */

#include "PhaseTransitionAnalyzer.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

float PhaseTransitionAnalyzer::computeOrderParameter(const std::vector<float>& R_field,
                                                       uint32_t Nx, uint32_t Ny) {
    if (R_field.empty()) return 0.0f;

    float sum = std::accumulate(R_field.begin(), R_field.end(), 0.0f);
    return sum / static_cast<float>(Nx * Ny);
}

void PhaseTransitionAnalyzer::computeOrderParameterStats(const std::vector<float>& R_field,
                                                           uint32_t Nx, uint32_t Ny,
                                                           float& mean, float& variance,
                                                           float& min_val, float& max_val) {
    if (R_field.empty()) {
        mean = variance = min_val = max_val = 0.0f;
        return;
    }

    // Mean
    mean = computeOrderParameter(R_field, Nx, Ny);

    // Min/Max
    auto minmax = std::minmax_element(R_field.begin(), R_field.end());
    min_val = *minmax.first;
    max_val = *minmax.second;

    // Variance: Var(R) = ⟨R²⟩ - ⟨R⟩²
    float R_squared_sum = 0.0f;
    for (float R : R_field) {
        R_squared_sum += R * R;
    }
    float R_squared_mean = R_squared_sum / static_cast<float>(Nx * Ny);
    variance = R_squared_mean - mean * mean;
}

PhaseTransitionAnalyzer::CriticalPointResult
PhaseTransitionAnalyzer::findCriticalPoint(const std::vector<DataPoint>& data) {
    CriticalPointResult result;
    result.success = false;

    if (data.size() < 3) {
        result.message = "Insufficient data points (need at least 3)";
        return result;
    }

    // Compute derivative d⟨R⟩/dσ
    std::vector<float> dR_dsigma = computeDerivative(data);

    // Find maximum |dR/dσ| (steepest drop)
    float max_derivative = 0.0f;
    size_t critical_idx = 0;

    for (size_t i = 0; i < dR_dsigma.size(); ++i) {
        float abs_deriv = std::abs(dR_dsigma[i]);
        if (abs_deriv > max_derivative) {
            max_derivative = abs_deriv;
            critical_idx = i;
        }
    }

    // Critical point estimate
    result.sigma_c = data[critical_idx].sigma;
    result.sigma_c_error = 0.0f;  // Simple estimate: no error for now

    // Estimate error from neighboring points
    if (critical_idx > 0 && critical_idx < data.size() - 1) {
        float dsigma = (data[critical_idx + 1].sigma - data[critical_idx - 1].sigma) / 2.0f;
        result.sigma_c_error = dsigma;
    }

    result.success = true;
    result.message = "Critical point found at maximum |dR/dσ|";

    return result;
}

PhaseTransitionAnalyzer::CriticalPointResult
PhaseTransitionAnalyzer::fitCriticalExponent(const std::vector<DataPoint>& data,
                                               float sigma_c,
                                               float fit_range_fraction) {
    CriticalPointResult result;
    result.success = false;
    result.sigma_c = sigma_c;

    // Filter data: only σ < σ_c (ordered phase)
    std::vector<float> log_delta_sigma;
    std::vector<float> log_R;

    for (const auto& point : data) {
        if (point.sigma < sigma_c && point.R_mean > 0.0f) {
            float delta = sigma_c - point.sigma;
            if (delta > 1e-6f) {  // Avoid log(0)
                log_delta_sigma.push_back(std::log(delta));
                log_R.push_back(std::log(point.R_mean));
            }
        }
    }

    if (log_delta_sigma.size() < 3) {
        result.message = "Insufficient data below σ_c for power-law fit";
        return result;
    }

    // Use only the closest fit_range_fraction of points to σ_c
    size_t num_fit_points = static_cast<size_t>(log_delta_sigma.size() * fit_range_fraction);
    num_fit_points = std::max(num_fit_points, size_t(3));  // At least 3 points
    num_fit_points = std::min(num_fit_points, log_delta_sigma.size());

    // Take the largest log(δσ) values (closest to σ_c)
    std::vector<float> x_fit(log_delta_sigma.end() - num_fit_points, log_delta_sigma.end());
    std::vector<float> y_fit(log_R.end() - num_fit_points, log_R.end());

    // Linear regression: log(R) = β·log(σ_c - σ) + const
    float beta, intercept, R2;
    if (!linearRegression(x_fit, y_fit, beta, intercept, R2)) {
        result.message = "Linear regression failed";
        return result;
    }

    result.beta = beta;
    result.R2 = R2;
    result.success = true;

    // Estimate error from R² (rough heuristic)
    result.beta_error = beta * std::sqrt(1.0f - R2) / std::sqrt(static_cast<float>(num_fit_points));

    result.message = "Power-law fit successful";

    return result;
}

std::string PhaseTransitionAnalyzer::detectTransitionType(const std::vector<DataPoint>& data) {
    if (data.size() < 3) return "none";

    // Compute derivative
    std::vector<float> dR_dsigma = computeDerivative(data);

    // Check for discontinuity (first-order)
    // Look for large jump in ⟨R⟩ between consecutive points
    float max_jump = 0.0f;
    for (size_t i = 1; i < data.size(); ++i) {
        float jump = std::abs(data[i].R_mean - data[i-1].R_mean);
        max_jump = std::max(max_jump, jump);
    }

    // Heuristic: first-order if max jump > 0.3
    // (this should be tuned based on system size and noise level)
    if (max_jump > 0.3f) {
        return "first_order";
    }

    // Check for continuous but diverging derivative (second-order)
    float max_deriv = 0.0f;
    for (float deriv : dR_dsigma) {
        max_deriv = std::max(max_deriv, std::abs(deriv));
    }

    // Heuristic: second-order if |dR/dσ|_max > 5
    if (max_deriv > 5.0f) {
        return "second_order";
    }

    return "none";
}

bool PhaseTransitionAnalyzer::writeReport(const std::vector<DataPoint>& data,
                                           const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "[PhaseTransitionAnalyzer] Failed to open " << output_path << std::endl;
        return false;
    }

    out << "========================================" << std::endl;
    out << "Phase Transition Analysis Report" << std::endl;
    out << "Test 3.3: Noise-Driven Phase Transition" << std::endl;
    out << "========================================\n" << std::endl;

    // Find critical point
    auto crit_result = findCriticalPoint(data);

    out << "[Critical Point Detection]" << std::endl;
    if (crit_result.success) {
        out << "  σ_c = " << std::fixed << std::setprecision(4)
            << crit_result.sigma_c << " ± " << crit_result.sigma_c_error << std::endl;
        out << "  Status: " << crit_result.message << std::endl;
    } else {
        out << "  FAILED: " << crit_result.message << std::endl;
    }
    out << std::endl;

    // Fit critical exponent
    if (crit_result.success) {
        auto fit_result = fitCriticalExponent(data, crit_result.sigma_c);
        out << "[Critical Exponent]" << std::endl;
        if (fit_result.success) {
            out << "  β = " << std::fixed << std::setprecision(4)
                << fit_result.beta << " ± " << fit_result.beta_error << std::endl;
            out << "  R² = " << fit_result.R2 << std::endl;
            out << "  Expected (2D Ising): β ≈ 0.125" << std::endl;

            float deviation = std::abs(fit_result.beta - 0.125f);
            if (deviation < 0.05f) {
                out << "  PASS: Consistent with 2D Ising universality (|Δβ| < 0.05)" << std::endl;
            } else {
                out << "  WARNING: Deviation from 2D Ising: Δβ = " << deviation << std::endl;
            }
        } else {
            out << "  FAILED: " << fit_result.message << std::endl;
        }
        out << std::endl;
    }

    // Detect transition type
    std::string trans_type = detectTransitionType(data);
    out << "[Transition Type]" << std::endl;
    out << "  Detected: " << trans_type << std::endl;
    out << "  Expected: second_order (continuous phase transition)" << std::endl;
    out << std::endl;

    // Data table
    out << "[Noise Scan Data]" << std::endl;
    out << std::setw(12) << "σ" << " | "
        << std::setw(12) << "⟨R⟩" << " | "
        << std::setw(12) << "σ²(R)" << " | "
        << std::setw(12) << "R_min" << " | "
        << std::setw(12) << "R_max" << std::endl;
    out << std::string(68, '-') << std::endl;

    for (const auto& point : data) {
        out << std::fixed << std::setprecision(6)
            << std::setw(12) << point.sigma << " | "
            << std::setw(12) << point.R_mean << " | "
            << std::setw(12) << point.R_variance << " | "
            << std::setw(12) << point.R_min << " | "
            << std::setw(12) << point.R_max << std::endl;
    }
    out << std::endl;

    out << "========================================" << std::endl;
    out << "Report generation complete" << std::endl;
    out << "========================================" << std::endl;

    out.close();
    return true;
}

// Private helper methods

std::vector<float> PhaseTransitionAnalyzer::computeDerivative(const std::vector<DataPoint>& data) {
    std::vector<float> derivative;
    if (data.size() < 2) return derivative;

    derivative.resize(data.size());

    // Forward difference for first point
    derivative[0] = (data[1].R_mean - data[0].R_mean) /
                     (data[1].sigma - data[0].sigma);

    // Centered differences for interior points
    for (size_t i = 1; i < data.size() - 1; ++i) {
        derivative[i] = (data[i+1].R_mean - data[i-1].R_mean) /
                        (data[i+1].sigma - data[i-1].sigma);
    }

    // Backward difference for last point
    size_t n = data.size() - 1;
    derivative[n] = (data[n].R_mean - data[n-1].R_mean) /
                     (data[n].sigma - data[n-1].sigma);

    return derivative;
}

bool PhaseTransitionAnalyzer::linearRegression(const std::vector<float>& x,
                                                const std::vector<float>& y,
                                                float& a, float& b, float& R2) {
    if (x.size() != y.size() || x.size() < 2) {
        return false;
    }

    size_t n = x.size();

    // Compute means
    float x_mean = std::accumulate(x.begin(), x.end(), 0.0f) / n;
    float y_mean = std::accumulate(y.begin(), y.end(), 0.0f) / n;

    // Compute covariance and variance
    float cov_xy = 0.0f;
    float var_x = 0.0f;
    float var_y = 0.0f;

    for (size_t i = 0; i < n; ++i) {
        float dx = x[i] - x_mean;
        float dy = y[i] - y_mean;
        cov_xy += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    if (var_x < 1e-10f) {
        return false;  // Degenerate case
    }

    // Slope and intercept
    a = cov_xy / var_x;
    b = y_mean - a * x_mean;

    // R² (coefficient of determination)
    float ss_res = 0.0f;  // Residual sum of squares
    for (size_t i = 0; i < n; ++i) {
        float y_pred = a * x[i] + b;
        float residual = y[i] - y_pred;
        ss_res += residual * residual;
    }

    float ss_tot = var_y;  // Total sum of squares (already computed)
    R2 = 1.0f - (ss_res / (ss_tot + 1e-10f));  // Avoid division by zero

    return true;
}
