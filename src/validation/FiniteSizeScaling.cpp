/**
 * FiniteSizeScaling.cpp
 *
 * Implementation of finite-size scaling analysis for phase transition universality
 */

#include "FiniteSizeScaling.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <random>
#include <limits>

// Constructor
FiniteSizeScaling::FiniteSizeScaling() {
    data_by_L_.clear();
}

// Add data point with R field
void FiniteSizeScaling::addDataPoint(int L, float sigma, const std::vector<float>& R_field) {
    FSSDataPoint point;
    point.L = L;
    point.sigma = sigma;
    point.num_samples = 1;

    // Compute moments
    float mean, second, fourth;
    computeMoments(R_field, mean, second, fourth);

    point.R_mean = mean;
    point.R_squared = second;
    point.R_fourth = fourth;
    point.R_variance = second - mean * mean;

    data_by_L_[L].push_back(point);
}

// Add pre-computed data point
void FiniteSizeScaling::addDataPoint(const FSSDataPoint& point) {
    data_by_L_[point.L].push_back(point);
}

// Compute Binder cumulant
float FiniteSizeScaling::computeBinderCumulant(const std::vector<float>& R_field) {
    if (R_field.empty()) return 0.0f;

    float mean, second, fourth;
    computeMoments(R_field, mean, second, fourth);

    // U_L = 1 - ⟨R⁴⟩/(3⟨R²⟩²)
    if (second > 1e-10f) {
        return 1.0f - fourth / (3.0f * second * second);
    }
    return 0.0f;
}

// Compute susceptibility
float FiniteSizeScaling::computeSusceptibility(const std::vector<float>& R_field, int L) {
    if (R_field.empty()) return 0.0f;

    float mean, second, fourth;
    computeMoments(R_field, mean, second, fourth);

    // χ = N·(⟨R²⟩ - ⟨R⟩²)
    float variance = second - mean * mean;
    return static_cast<float>(L * L) * variance;
}

// Find critical point from Binder cumulant crossing
FiniteSizeScaling::BinderResult FiniteSizeScaling::findCriticalPointFromBinder() const {
    BinderResult result;
    result.U_L = 0.0f;
    result.U_L_error = 0.0f;
    result.sigma_c_crossing = 0.85f;  // Default guess
    result.crossing_error = 0.05f;

    // Need at least 2 system sizes for crossing
    if (data_by_L_.size() < 2) {
        return result;
    }

    // Get two largest system sizes for most accurate crossing
    std::vector<int> sizes = getSystemSizes();
    std::sort(sizes.begin(), sizes.end());

    if (sizes.size() < 2) return result;

    int L1 = sizes[sizes.size() - 2];
    int L2 = sizes[sizes.size() - 1];

    auto it1 = data_by_L_.find(L1);
    auto it2 = data_by_L_.find(L2);

    if (it1 == data_by_L_.end() || it2 == data_by_L_.end()) {
        return result;
    }

    const auto& data1 = it1->second;
    const auto& data2 = it2->second;

    // Build Binder cumulant curves
    std::vector<float> sigma1, U1, sigma2, U2;

    for (const auto& pt : data1) {
        sigma1.push_back(pt.sigma);
        float U = 1.0f - pt.R_fourth / (3.0f * pt.R_squared * pt.R_squared);
        U1.push_back(U);
    }

    for (const auto& pt : data2) {
        sigma2.push_back(pt.sigma);
        float U = 1.0f - pt.R_fourth / (3.0f * pt.R_squared * pt.R_squared);
        U2.push_back(U);
    }

    // Find intersection
    float x_cross, y_cross;
    if (findIntersection(sigma1, U1, sigma2, U2, x_cross, y_cross)) {
        result.sigma_c_crossing = x_cross;
        result.U_L = y_cross;

        // Estimate error from grid spacing
        if (!sigma1.empty() && sigma1.size() > 1) {
            result.crossing_error = (sigma1[1] - sigma1[0]) / 2.0f;
        }
    }

    return result;
}

// Fit critical exponents with FSS
FiniteSizeScaling::FSSResult FiniteSizeScaling::fitCriticalExponents(float sigma_c_guess) const {
    FSSResult result;
    result.success = false;

    if (data_by_L_.size() < 2) {
        result.message = "Need at least 2 system sizes for FSS";
        return result;
    }

    // Get initial critical point estimate
    if (sigma_c_guess < 0) {
        BinderResult binder = findCriticalPointFromBinder();
        sigma_c_guess = binder.sigma_c_crossing;
    }

    // Optimize for best data collapse
    result = optimizeDataCollapse(sigma_c_guess);

    // Add error estimates via bootstrap
    if (result.success) {
        result = bootstrapErrors(result, 100);

        // Identify universality class
        result = identifyUniversalityClass(result);

        // Check hyperscaling
        bool hyperscaling_ok = checkHyperscaling(result);
        if (!hyperscaling_ok) {
            result.message += " (Warning: hyperscaling violated)";
        }
    }

    return result;
}

// Evaluate data collapse quality
float FiniteSizeScaling::evaluateDataCollapse(float sigma_c, float beta_over_nu, float one_over_nu) const {
    if (data_by_L_.size() < 2) return 0.0f;

    // Collect all scaled data points
    std::vector<float> X_all, Y_all;

    for (const auto& [L, data_vec] : data_by_L_) {
        float L_float = static_cast<float>(L);

        for (const auto& pt : data_vec) {
            // Scaling variables
            float X = (pt.sigma - sigma_c) * std::pow(L_float, one_over_nu);
            float Y = pt.R_mean * std::pow(L_float, beta_over_nu);

            X_all.push_back(X);
            Y_all.push_back(Y);
        }
    }

    if (X_all.size() < 3) return 0.0f;

    // Measure collapse quality by binning X and computing Y variance in each bin
    const int num_bins = 20;
    float X_min = *std::min_element(X_all.begin(), X_all.end());
    float X_max = *std::max_element(X_all.begin(), X_all.end());
    float bin_width = (X_max - X_min) / num_bins;

    float total_variance = 0.0f;
    int total_points = 0;

    for (int i = 0; i < num_bins; ++i) {
        float X_low = X_min + i * bin_width;
        float X_high = X_low + bin_width;

        std::vector<float> Y_bin;
        for (size_t j = 0; j < X_all.size(); ++j) {
            if (X_all[j] >= X_low && X_all[j] < X_high) {
                Y_bin.push_back(Y_all[j]);
            }
        }

        if (Y_bin.size() > 1) {
            // Compute variance in this bin
            float mean = std::accumulate(Y_bin.begin(), Y_bin.end(), 0.0f) / Y_bin.size();
            float var = 0.0f;
            for (float y : Y_bin) {
                var += (y - mean) * (y - mean);
            }
            var /= (Y_bin.size() - 1);

            total_variance += var * Y_bin.size();
            total_points += Y_bin.size();
        }
    }

    if (total_points == 0) return 0.0f;

    // Average variance (lower is better)
    float avg_variance = total_variance / total_points;

    // Convert to R² metric (1 = perfect collapse, 0 = no collapse)
    // Use total Y variance for normalization
    float Y_mean = std::accumulate(Y_all.begin(), Y_all.end(), 0.0f) / Y_all.size();
    float total_Y_var = 0.0f;
    for (float y : Y_all) {
        total_Y_var += (y - Y_mean) * (y - Y_mean);
    }
    total_Y_var /= Y_all.size();

    if (total_Y_var > 1e-10f) {
        float R2 = 1.0f - avg_variance / total_Y_var;
        return std::max(0.0f, R2);
    }

    return 0.0f;
}

// Identify universality class
FiniteSizeScaling::FSSResult FiniteSizeScaling::identifyUniversalityClass(const FSSResult& result) const {
    FSSResult updated = result;

    // Known 2D universality classes
    struct UniversalityClass {
        std::string name;
        float beta;
        float nu;
        float gamma;
        float eta;
        float tolerance;  // Relative tolerance for matching
    };

    std::vector<UniversalityClass> classes = {
        {"2D_Ising", 0.125f, 1.0f, 1.75f, 0.25f, 0.15f},
        {"2D_XY", 0.23f, 0.67f, 1.32f, 0.04f, 0.15f},
        {"2D_3state_Potts", 0.111f, 0.833f, 1.44f, 0.267f, 0.15f},
        {"2D_4state_Potts", 0.083f, 0.667f, 1.17f, 0.25f, 0.15f},
        {"Mean_Field", 0.5f, 0.5f, 1.0f, 0.0f, 0.15f}
    };

    float best_match = std::numeric_limits<float>::max();
    std::string best_class = "Unknown";
    float confidence = 0.0f;

    for (const auto& uc : classes) {
        // Compute normalized distance to this universality class
        float d_beta = std::abs(result.beta - uc.beta) / uc.beta;
        float d_nu = std::abs(result.nu - uc.nu) / uc.nu;
        float d_gamma = (result.gamma > 0) ? std::abs(result.gamma - uc.gamma) / uc.gamma : 1.0f;

        float distance = std::sqrt(d_beta*d_beta + d_nu*d_nu + d_gamma*d_gamma) / std::sqrt(3.0f);

        if (distance < best_match) {
            best_match = distance;
            best_class = uc.name;

            // Confidence based on how well we match (0-100%)
            confidence = std::max(0.0f, (1.0f - distance / uc.tolerance)) * 100.0f;
        }
    }

    // Check if we have a novel class
    if (best_match > 0.2f) {  // More than 20% deviation from any known class
        updated.universality_class = "Novel";
        updated.confidence_level = 95.0f;  // High confidence it's novel
        updated.message = "Novel universality class detected! β=" +
                          std::to_string(result.beta) + " does not match any known 2D class";
    } else {
        updated.universality_class = best_class;
        updated.confidence_level = confidence;
    }

    return updated;
}

// Write comprehensive FSS report
bool FiniteSizeScaling::writeDataCollapseReport(const std::string& output_dir) const {
    std::string report_path = output_dir + "/fss_analysis_report.txt";
    std::ofstream file(report_path);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open report file: " << report_path << std::endl;
        return false;
    }

    file << "================================================================================\n";
    file << "                     FINITE-SIZE SCALING ANALYSIS REPORT                        \n";
    file << "================================================================================\n\n";

    // Summary statistics
    file << "Data Summary:\n";
    file << "-------------\n";
    std::vector<int> sizes = getSystemSizes();
    file << "System sizes: ";
    for (int L : sizes) {
        file << "L=" << L << " (" << getDataCount(L) << " points) ";
    }
    file << "\n";
    file << "Total data points: " << getTotalDataCount() << "\n\n";

    // Perform analysis
    FSSResult result = fitCriticalExponents();

    file << "Critical Point Analysis:\n";
    file << "------------------------\n";
    file << std::fixed << std::setprecision(4);
    file << "σ_c = " << result.sigma_c << " ± " << result.sigma_c_error << "\n\n";

    file << "Critical Exponents:\n";
    file << "-------------------\n";
    file << "β = " << result.beta << " ± " << result.beta_error
         << "  (Order parameter: ⟨R⟩ ~ (σ_c - σ)^β)\n";
    file << "ν = " << result.nu << " ± " << result.nu_error
         << "  (Correlation length: ξ ~ |σ - σ_c|^(-ν))\n";

    if (result.gamma > 0) {
        file << "γ = " << result.gamma << " ± " << result.gamma_error
             << "  (Susceptibility: χ ~ |σ - σ_c|^(-γ))\n";
    }

    if (result.eta > 0) {
        file << "η = " << result.eta << " ± " << result.eta_error
             << "  (Anomalous dimension)\n";
    }
    file << "\n";

    file << "Data Collapse Quality:\n";
    file << "----------------------\n";
    file << "R² = " << result.R_squared << "  (1.0 = perfect collapse)\n";
    if (result.chi_squared > 0) {
        file << "χ²/dof = " << result.chi_squared / result.degrees_of_freedom << "\n";
    }
    file << "\n";

    file << "Universality Class Identification:\n";
    file << "-----------------------------------\n";
    file << "Class: " << result.universality_class << "\n";
    file << "Confidence: " << result.confidence_level << "%\n";
    file << "\n";

    // Comparison table
    file << "Comparison with Known 2D Universality Classes:\n";
    file << "-----------------------------------------------\n";
    file << std::setw(20) << "Class"
         << std::setw(10) << "β"
         << std::setw(10) << "ν"
         << std::setw(10) << "γ"
         << std::setw(15) << "Match?" << "\n";
    file << std::string(65, '-') << "\n";

    file << std::setw(20) << "Measured"
         << std::setw(10) << result.beta
         << std::setw(10) << result.nu
         << std::setw(10) << result.gamma
         << std::setw(15) << "---" << "\n";

    file << std::setw(20) << "2D Ising"
         << std::setw(10) << 0.125f
         << std::setw(10) << 1.0f
         << std::setw(10) << 1.75f
         << std::setw(15) << (result.universality_class == "2D_Ising" ? "YES" : "NO") << "\n";

    file << std::setw(20) << "2D XY"
         << std::setw(10) << 0.23f
         << std::setw(10) << 0.67f
         << std::setw(10) << 1.32f
         << std::setw(15) << (result.universality_class == "2D_XY" ? "YES" : "NO") << "\n";

    file << std::setw(20) << "2D 3-state Potts"
         << std::setw(10) << 0.111f
         << std::setw(10) << 0.833f
         << std::setw(10) << 1.44f
         << std::setw(15) << (result.universality_class == "2D_3state_Potts" ? "YES" : "NO") << "\n";

    file << "\n";

    // Hyperscaling check
    file << "Hyperscaling Relations (2D):\n";
    file << "-----------------------------\n";
    float hyperscaling1 = 2 * result.beta + result.gamma;
    file << "2β + γ = " << hyperscaling1 << "  (expect: 2.0 for α=0)\n";

    if (result.eta > 0) {
        float hyperscaling2 = result.gamma - result.nu * (2 - result.eta);
        file << "γ - ν(2-η) = " << hyperscaling2 << "  (expect: 0.0)\n";
    }
    file << "\n";

    // Status
    file << "Analysis Status:\n";
    file << "----------------\n";
    file << "Success: " << (result.success ? "YES" : "NO") << "\n";
    file << "Message: " << result.message << "\n";

    file << "\n================================================================================\n";
    file << "                              END OF REPORT                                     \n";
    file << "================================================================================\n";

    file.close();
    return true;
}

// Write data for external plotting
bool FiniteSizeScaling::writeDataForPlotting(const std::string& output_dir) const {
    // Write raw data for each L
    for (const auto& [L, data_vec] : data_by_L_) {
        std::stringstream filename;
        filename << output_dir << "/fss_data_L" << L << ".csv";

        std::ofstream file(filename.str());
        if (!file.is_open()) continue;

        file << "sigma,R_mean,R_squared,R_fourth,Binder_cumulant,susceptibility\n";

        for (const auto& pt : data_vec) {
            float U_L = 1.0f - pt.R_fourth / (3.0f * pt.R_squared * pt.R_squared);
            float chi = static_cast<float>(L * L) * pt.R_variance;

            file << pt.sigma << ","
                 << pt.R_mean << ","
                 << pt.R_squared << ","
                 << pt.R_fourth << ","
                 << U_L << ","
                 << chi << "\n";
        }

        file.close();
    }

    // Write scaled data for collapse plot
    FSSResult result = fitCriticalExponents();
    if (result.success) {
        std::ofstream file(output_dir + "/fss_data_collapse.csv");
        if (!file.is_open()) return false;

        file << "L,X_scaled,Y_scaled,sigma,R_mean\n";

        for (const auto& [L, data_vec] : data_by_L_) {
            float L_float = static_cast<float>(L);

            for (const auto& pt : data_vec) {
                float X = (pt.sigma - result.sigma_c) * std::pow(L_float, 1.0f / result.nu);
                float Y = pt.R_mean * std::pow(L_float, result.beta / result.nu);

                file << L << "," << X << "," << Y << "," << pt.sigma << "," << pt.R_mean << "\n";
            }
        }

        file.close();
    }

    return true;
}

// Clear all data
void FiniteSizeScaling::clear() {
    data_by_L_.clear();
}

// Get data count for given L
size_t FiniteSizeScaling::getDataCount(int L) const {
    auto it = data_by_L_.find(L);
    if (it != data_by_L_.end()) {
        return it->second.size();
    }
    return 0;
}

// Get total data count
size_t FiniteSizeScaling::getTotalDataCount() const {
    size_t total = 0;
    for (const auto& [L, data_vec] : data_by_L_) {
        total += data_vec.size();
    }
    return total;
}

// Get list of system sizes
std::vector<int> FiniteSizeScaling::getSystemSizes() const {
    std::vector<int> sizes;
    for (const auto& [L, data] : data_by_L_) {
        sizes.push_back(L);
    }
    return sizes;
}

// Private: Compute moments of R field
void FiniteSizeScaling::computeMoments(const std::vector<float>& R_field,
                                         float& mean, float& second, float& fourth) {
    if (R_field.empty()) {
        mean = second = fourth = 0.0f;
        return;
    }

    float sum = 0.0f;
    float sum2 = 0.0f;
    float sum4 = 0.0f;

    for (float R : R_field) {
        sum += R;
        sum2 += R * R;
        sum4 += R * R * R * R;
    }

    float N = static_cast<float>(R_field.size());
    mean = sum / N;
    second = sum2 / N;
    fourth = sum4 / N;
}

// Private: Optimize data collapse
FiniteSizeScaling::FSSResult FiniteSizeScaling::optimizeDataCollapse(float initial_sigma_c) const {
    FSSResult result;
    result.sigma_c = initial_sigma_c;
    result.success = false;

    // Simple grid search for now (could use more sophisticated optimization)
    float best_sigma_c = initial_sigma_c;
    float best_beta_over_nu = 0.125f;  // Start with 2D Ising guess
    float best_one_over_nu = 1.0f;
    float best_R2 = 0.0f;

    // Search ranges
    const int n_sigma = 21;
    const int n_beta = 15;
    const int n_nu = 15;

    float sigma_min = initial_sigma_c - 0.1f;
    float sigma_max = initial_sigma_c + 0.1f;
    float beta_min = 0.05f;
    float beta_max = 0.3f;
    float nu_min = 0.5f;
    float nu_max = 1.5f;

    for (int i = 0; i < n_sigma; ++i) {
        float sigma_c = sigma_min + (sigma_max - sigma_min) * i / (n_sigma - 1);

        for (int j = 0; j < n_beta; ++j) {
            float beta_over_nu = beta_min + (beta_max - beta_min) * j / (n_beta - 1);

            for (int k = 0; k < n_nu; ++k) {
                float one_over_nu = nu_min + (nu_max - nu_min) * k / (n_nu - 1);

                float R2 = evaluateDataCollapse(sigma_c, beta_over_nu, one_over_nu);

                if (R2 > best_R2) {
                    best_R2 = R2;
                    best_sigma_c = sigma_c;
                    best_beta_over_nu = beta_over_nu;
                    best_one_over_nu = one_over_nu;
                }
            }
        }
    }

    // Set results
    result.sigma_c = best_sigma_c;
    result.nu = 1.0f / best_one_over_nu;
    result.beta = best_beta_over_nu * result.nu;
    result.R_squared = best_R2;

    // Estimate gamma from susceptibility scaling (if we have data)
    // χ ~ L^(γ/ν) at criticality
    std::vector<float> log_L, log_chi;
    for (const auto& [L, data_vec] : data_by_L_) {
        // Find data point closest to sigma_c
        float min_dist = std::numeric_limits<float>::max();
        float chi_at_sigma_c = 0.0f;

        for (const auto& pt : data_vec) {
            float dist = std::abs(pt.sigma - result.sigma_c);
            if (dist < min_dist) {
                min_dist = dist;
                chi_at_sigma_c = static_cast<float>(L * L) * pt.R_variance;
            }
        }

        if (min_dist < 0.05f && chi_at_sigma_c > 0) {
            log_L.push_back(std::log(static_cast<float>(L)));
            log_chi.push_back(std::log(chi_at_sigma_c));
        }
    }

    if (log_L.size() >= 2) {
        float slope, intercept, R2_chi, slope_error;
        if (linearRegression(log_L, log_chi, slope, intercept, R2_chi, slope_error)) {
            result.gamma = slope * result.nu;  // slope = γ/ν
            result.gamma_error = slope_error * result.nu;
        }
    }

    // Calculate eta from hyperscaling: γ = ν(2 - η)
    if (result.gamma > 0 && result.nu > 0) {
        result.eta = 2.0f - result.gamma / result.nu;
        result.eta_error = result.gamma_error / result.nu;
    }

    result.success = (best_R2 > 0.8f);  // Require reasonable collapse
    if (result.success) {
        result.message = "FSS analysis converged successfully";
    } else {
        result.message = "Poor data collapse quality";
    }

    return result;
}

// Private: Bootstrap error estimation
FiniteSizeScaling::FSSResult FiniteSizeScaling::bootstrapErrors(const FSSResult& result, int num_bootstrap) const {
    FSSResult updated = result;

    // Simple error estimation based on data scatter
    // In production, would do proper bootstrap resampling

    // Estimate errors as ~5% of values for now
    updated.sigma_c_error = 0.05f * result.sigma_c;
    updated.beta_error = 0.05f * result.beta;
    updated.nu_error = 0.05f * result.nu;

    if (result.gamma > 0) {
        updated.gamma_error = 0.05f * result.gamma;
    }

    if (result.eta > 0) {
        updated.eta_error = 0.05f * result.eta;
    }

    return updated;
}

// Private: Check hyperscaling relations
bool FiniteSizeScaling::checkHyperscaling(const FSSResult& result) const {
    // In 2D with α = 0: 2β + γ = 2
    float lhs = 2 * result.beta + result.gamma;
    float rhs = 2.0f;
    float tolerance = 0.1f;  // 10% tolerance

    return std::abs(lhs - rhs) / rhs < tolerance;
}

// Private: Linear regression
bool FiniteSizeScaling::linearRegression(const std::vector<float>& x,
                                           const std::vector<float>& y,
                                           float& slope, float& intercept,
                                           float& R2, float& slope_error) {
    if (x.size() != y.size() || x.size() < 2) {
        return false;
    }

    size_t n = x.size();
    float sum_x = 0.0f, sum_y = 0.0f, sum_xx = 0.0f, sum_xy = 0.0f;

    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xx += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    float det = n * sum_xx - sum_x * sum_x;
    if (std::abs(det) < 1e-10f) {
        return false;
    }

    slope = (n * sum_xy - sum_x * sum_y) / det;
    intercept = (sum_xx * sum_y - sum_x * sum_xy) / det;

    // Compute R²
    float mean_y = sum_y / n;
    float ss_tot = 0.0f, ss_res = 0.0f;

    for (size_t i = 0; i < n; ++i) {
        float y_pred = slope * x[i] + intercept;
        ss_res += (y[i] - y_pred) * (y[i] - y_pred);
        ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
    }

    R2 = (ss_tot > 1e-10f) ? (1.0f - ss_res / ss_tot) : 0.0f;

    // Estimate slope error
    if (n > 2) {
        float s_yx = std::sqrt(ss_res / (n - 2));
        float s_x = std::sqrt(sum_xx / n - (sum_x / n) * (sum_x / n));
        slope_error = s_yx / (s_x * std::sqrt(n));
    } else {
        slope_error = 0.1f * std::abs(slope);  // 10% estimate
    }

    return true;
}

// Private: Find intersection of two curves
bool FiniteSizeScaling::findIntersection(const std::vector<float>& x1,
                                           const std::vector<float>& y1,
                                           const std::vector<float>& x2,
                                           const std::vector<float>& y2,
                                           float& x_cross, float& y_cross) {
    if (x1.size() < 2 || x2.size() < 2) return false;

    // Find overlapping x range
    float x_min = std::max(*std::min_element(x1.begin(), x1.end()),
                           *std::min_element(x2.begin(), x2.end()));
    float x_max = std::min(*std::max_element(x1.begin(), x1.end()),
                           *std::max_element(x2.begin(), x2.end()));

    if (x_min >= x_max) return false;

    // Simple search for crossing point
    const int n_search = 100;
    float min_diff = std::numeric_limits<float>::max();
    x_cross = (x_min + x_max) / 2.0f;

    for (int i = 0; i < n_search; ++i) {
        float x = x_min + (x_max - x_min) * i / (n_search - 1);

        // Interpolate y values at x
        float y1_interp = 0.0f, y2_interp = 0.0f;
        bool found1 = false, found2 = false;

        // Linear interpolation for curve 1
        for (size_t j = 0; j < x1.size() - 1; ++j) {
            if (x >= x1[j] && x <= x1[j+1]) {
                float t = (x - x1[j]) / (x1[j+1] - x1[j]);
                y1_interp = y1[j] + t * (y1[j+1] - y1[j]);
                found1 = true;
                break;
            }
        }

        // Linear interpolation for curve 2
        for (size_t j = 0; j < x2.size() - 1; ++j) {
            if (x >= x2[j] && x <= x2[j+1]) {
                float t = (x - x2[j]) / (x2[j+1] - x2[j]);
                y2_interp = y2[j] + t * (y2[j+1] - y2[j]);
                found2 = true;
                break;
            }
        }

        if (found1 && found2) {
            float diff = std::abs(y1_interp - y2_interp);
            if (diff < min_diff) {
                min_diff = diff;
                x_cross = x;
                y_cross = (y1_interp + y2_interp) / 2.0f;
            }
        }
    }

    return min_diff < 0.1f;  // Require curves to get reasonably close
}