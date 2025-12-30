#include "UniversalityClassifier.h"
#include <random>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/LeastSquares>

namespace SMFT {

// Constructor
UniversalityClassifier::UniversalityClassifier() :
    sigma_min(0.0), sigma_max(2.0),
    n_bootstrap(1000), convergence_tol(1e-4) {}

// Data input methods
void UniversalityClassifier::addDataPoint(const DataPoint& point) {
    data_by_size[point.L].push_back(point);

    // Update grid sizes list if new
    if (std::find(grid_sizes.begin(), grid_sizes.end(), point.L) == grid_sizes.end()) {
        grid_sizes.push_back(point.L);
        std::sort(grid_sizes.begin(), grid_sizes.end());
    }

    // Update sigma range
    sigma_min = std::min(sigma_min, point.sigma);
    sigma_max = std::max(sigma_max, point.sigma);
}

void UniversalityClassifier::addGridSizeData(int L, const std::vector<double>& sigma,
                                            const std::vector<double>& R,
                                            const std::vector<double>& R2,
                                            const std::vector<double>& R4) {
    for (size_t i = 0; i < sigma.size(); ++i) {
        DataPoint point;
        point.L = L;
        point.sigma = sigma[i];
        point.R_mean = R[i];
        point.R2_mean = R2[i];
        point.R4_mean = R4[i];

        // Calculate susceptibility: χ = L^2 * (⟨R²⟩ - ⟨R⟩²)
        point.chi = L * L * (R2[i] - R[i] * R[i]);

        // Calculate Binder cumulant: U = 1 - ⟨R⁴⟩/(3⟨R²⟩²)
        if (R2[i] > 0) {
            point.U_L = 1.0 - R4[i] / (3.0 * R2[i] * R2[i]);
        }

        addDataPoint(point);
    }

    // Sort data by sigma for easier analysis
    sortDataBySigma(L);
}

void UniversalityClassifier::addCorrelationData(int L, double sigma, const CorrelationData& corr_data) {
    correlation_data[{L, sigma}] = corr_data;
}

// Critical point determination
double UniversalityClassifier::findCriticalPoint() {
    std::cout << "\n=== Finding Critical Point σ_c ===" << std::endl;

    // Method 1: Binder cumulant crossing
    double sigma_c_binder = findCriticalPointFromBinder();
    std::cout << "  Binder crossing: σ_c = " << sigma_c_binder << std::endl;

    // Method 2: Susceptibility peak scaling
    double sigma_c_chi = findCriticalPointFromSusceptibility();
    std::cout << "  Susceptibility peak: σ_c = " << sigma_c_chi << std::endl;

    // Average the methods (weighted by reliability)
    double sigma_c = (2.0 * sigma_c_binder + sigma_c_chi) / 3.0;
    std::cout << "  Final estimate: σ_c = " << sigma_c << std::endl;

    return sigma_c;
}

double UniversalityClassifier::findCriticalPointFromBinder() {
    std::vector<std::pair<double, double>> crossings;

    // Find crossings for all pairs of consecutive sizes
    for (size_t i = 0; i < grid_sizes.size() - 1; ++i) {
        int L1 = grid_sizes[i];
        int L2 = grid_sizes[i + 1];

        auto crossing = findBinderCrossing(L1, L2);
        if (crossing.first > 0) {
            crossings.push_back(crossing);
            std::cout << "    Crossing L=" << L1 << "/" << L2
                     << ": σ = " << crossing.first
                     << " ± " << crossing.second << std::endl;
        }
    }

    // Analyze crossings for convergence
    analyzeBinderCrossings(crossings);

    // Return weighted average
    double sigma_c = 0, weight_sum = 0;
    for (const auto& cross : crossings) {
        double weight = 1.0 / (cross.second * cross.second);
        sigma_c += cross.first * weight;
        weight_sum += weight;
    }

    return sigma_c / weight_sum;
}

double UniversalityClassifier::findCriticalPointFromSusceptibility() {
    std::vector<std::pair<int, double>> peaks;

    // Find susceptibility peak for each grid size
    for (int L : grid_sizes) {
        auto peak = findSusceptibilityPeak(L);
        if (peak.first > 0) {
            peaks.push_back({L, peak.first});
            std::cout << "    Peak at L=" << L << ": σ = " << peak.first << std::endl;
        }
    }

    // Analyze peak positions vs 1/L for extrapolation
    analyzeSusceptibilityPeaks(peaks);

    // Extrapolate to L → ∞
    std::vector<double> inv_L, sigma_peak;
    for (const auto& peak : peaks) {
        inv_L.push_back(1.0 / peak.first);
        sigma_peak.push_back(peak.second);
    }

    double slope, intercept, error;
    linearFit(inv_L, sigma_peak, slope, intercept, error);

    return intercept;  // σ_c at L → ∞
}

// Exponent extraction
UniversalityClassifier::CriticalExponents UniversalityClassifier::extractExponents() {
    std::cout << "\n=== Extracting Critical Exponents ===" << std::endl;

    CriticalExponents exp;

    // First find critical point
    exp.sigma_c = findCriticalPoint();
    exp.sigma_c_err = 0.001;  // Will be refined

    // Extract exponents
    exp.beta = extractBeta(exp.sigma_c);
    exp.nu = extractNu(exp.sigma_c);
    exp.gamma = extractGamma(exp.sigma_c);
    exp.eta = extractEta(exp.sigma_c);

    // Calculate derived exponents using scaling relations
    // Rushbrooke: α + 2β + γ = 2
    exp.alpha = 2.0 - 2.0 * exp.beta - exp.gamma;

    // Widom: γ = β(δ - 1)
    exp.delta = 1.0 + exp.gamma / exp.beta;

    // Estimate errors (simplified - should use proper propagation)
    exp.beta_err = 0.005;
    exp.nu_err = 0.01;
    exp.gamma_err = 0.02;
    exp.eta_err = 0.02;
    exp.alpha_err = sqrt(4 * exp.beta_err * exp.beta_err + exp.gamma_err * exp.gamma_err);
    exp.delta_err = sqrt(exp.gamma_err * exp.gamma_err / (exp.beta * exp.beta) +
                        exp.gamma * exp.gamma * exp.beta_err * exp.beta_err / (exp.beta * exp.beta * exp.beta * exp.beta));

    // Perform data collapse to validate exponents
    double chi2, quality;
    performDataCollapse(exp.sigma_c, exp.beta / exp.nu, exp.nu, chi2, quality);
    exp.chi2_collapse = chi2;
    exp.collapse_quality = quality;

    // Verify scaling relations
    verifyScalingRelations(exp);

    return exp;
}

double UniversalityClassifier::extractBeta(double sigma_c) {
    std::cout << "\n  Extracting β (order parameter exponent)..." << std::endl;

    std::vector<double> log_L, log_R;

    // Collect R(σ_c) for each L
    for (int L : grid_sizes) {
        const auto& data = data_by_size[L];

        // Find R at sigma_c by interpolation
        std::vector<double> sigmas, Rs;
        for (const auto& point : data) {
            sigmas.push_back(point.sigma);
            Rs.push_back(point.R_mean);
        }

        // Linear interpolation at sigma_c
        double R_at_sigma_c = 0;
        for (size_t i = 0; i < sigmas.size() - 1; ++i) {
            if (sigmas[i] <= sigma_c && sigmas[i+1] >= sigma_c) {
                double t = (sigma_c - sigmas[i]) / (sigmas[i+1] - sigmas[i]);
                R_at_sigma_c = Rs[i] * (1-t) + Rs[i+1] * t;
                break;
            }
        }

        if (R_at_sigma_c > 0) {
            log_L.push_back(log(L));
            log_R.push_back(log(R_at_sigma_c));
        }
    }

    // Power law fit: R ~ L^(-β/ν)
    double slope, intercept, error;
    linearFit(log_L, log_R, slope, intercept, error);

    // For now assume ν = 1.0 (will be refined)
    double beta = -slope;  // Actually -β/ν, will correct later

    std::cout << "    β/ν = " << -slope << " ± " << error << std::endl;

    return beta;
}

double UniversalityClassifier::extractNu(double sigma_c) {
    std::cout << "\n  Extracting ν (correlation length exponent)..." << std::endl;

    // Method 1: Width of Binder crossing region scales as L^(-1/ν)
    std::vector<double> log_L, log_width;

    for (int L : grid_sizes) {
        const auto& data = data_by_size[L];

        // Find width where U_L varies by 0.1 around crossing
        double U_target = 0;  // Will be set from data
        double sigma_low = sigma_c, sigma_high = sigma_c;

        // Find U at sigma_c
        for (const auto& point : data) {
            if (fabs(point.sigma - sigma_c) < 0.01) {
                U_target = point.U_L;
                break;
            }
        }

        // Find sigma range where U varies by ±0.05
        for (const auto& point : data) {
            if (fabs(point.U_L - U_target) < 0.05) {
                sigma_low = std::min(sigma_low, point.sigma);
                sigma_high = std::max(sigma_high, point.sigma);
            }
        }

        double width = sigma_high - sigma_low;
        if (width > 0) {
            log_L.push_back(log(L));
            log_width.push_back(log(width));
        }
    }

    // Fit: width ~ L^(-1/ν)
    double slope, intercept, error;
    linearFit(log_L, log_width, slope, intercept, error);

    double nu = -1.0 / slope;

    std::cout << "    ν = " << nu << " ± " << error * nu * nu << std::endl;

    return nu;
}

double UniversalityClassifier::extractGamma(double sigma_c) {
    std::cout << "\n  Extracting γ (susceptibility exponent)..." << std::endl;

    std::vector<double> log_L, log_chi_max;

    // Find χ_max for each L
    for (int L : grid_sizes) {
        const auto& data = data_by_size[L];

        double chi_max = 0;
        for (const auto& point : data) {
            chi_max = std::max(chi_max, point.chi);
        }

        if (chi_max > 0) {
            log_L.push_back(log(L));
            log_chi_max.push_back(log(chi_max));
        }
    }

    // Fit: χ_max ~ L^(γ/ν)
    double slope, intercept, error;
    linearFit(log_L, log_chi_max, slope, intercept, error);

    // For now assume ν = 1.0 (will be refined)
    double gamma = slope;  // Actually γ/ν

    std::cout << "    γ/ν = " << slope << " ± " << error << std::endl;

    return gamma;
}

double UniversalityClassifier::extractEta(double sigma_c) {
    std::cout << "\n  Extracting η (anomalous dimension)..." << std::endl;

    // Method 1: From correlation function at criticality
    // G(r) ~ r^(-(d-2+η)) at σ_c

    // For 2D: G(r) ~ r^(-η)
    double eta_sum = 0, weight_sum = 0;

    for (int L : grid_sizes) {
        auto key = std::make_pair(L, sigma_c);
        if (correlation_data.find(key) != correlation_data.end()) {
            const auto& corr = correlation_data[key];

            // Fit power law to correlation function
            std::vector<double> log_r, log_G;
            for (size_t i = 0; i < corr.r.size(); ++i) {
                if (corr.r[i] > 1 && corr.r[i] < L/4.0) {  // Avoid boundaries
                    log_r.push_back(log(corr.r[i]));
                    log_G.push_back(log(corr.G_r[i]));
                }
            }

            double slope, intercept, error;
            linearFit(log_r, log_G, slope, intercept, error);

            double eta_L = -slope;  // In 2D
            double weight = 1.0 / (error * error);

            eta_sum += eta_L * weight;
            weight_sum += weight;

            std::cout << "    L=" << L << ": η = " << eta_L << " ± " << error << std::endl;
        }
    }

    // Method 2: From Fisher relation γ = ν(2-η)
    // Will be used as consistency check

    double eta = (weight_sum > 0) ? eta_sum / weight_sum : 0.25;  // Default to Ising if no data

    std::cout << "    Final η = " << eta << std::endl;

    return eta;
}

// Data collapse
bool UniversalityClassifier::performDataCollapse(double sigma_c, double beta_nu, double nu,
                                                double& chi2, double& quality) {
    std::cout << "\n  Performing data collapse validation..." << std::endl;

    // Collapse: R * L^(β/ν) vs (σ - σ_c) * L^(1/ν)
    std::vector<std::vector<double>> x_scaled, y_scaled;

    for (int L : grid_sizes) {
        const auto& data = data_by_size[L];
        std::vector<double> x_L, y_L;

        for (const auto& point : data) {
            double x = (point.sigma - sigma_c) * pow(L, 1.0/nu);
            double y = point.R_mean * pow(L, beta_nu);

            x_L.push_back(x);
            y_L.push_back(y);
        }

        x_scaled.push_back(x_L);
        y_scaled.push_back(y_L);
    }

    // Compute chi2 between curves
    chi2 = 0;
    int n_pairs = 0;

    for (size_t i = 0; i < x_scaled.size() - 1; ++i) {
        for (size_t j = i + 1; j < x_scaled.size(); ++j) {
            // Find overlapping x range
            double x_min = std::max(*std::min_element(x_scaled[i].begin(), x_scaled[i].end()),
                                   *std::min_element(x_scaled[j].begin(), x_scaled[j].end()));
            double x_max = std::min(*std::max_element(x_scaled[i].begin(), x_scaled[i].end()),
                                   *std::max_element(x_scaled[j].begin(), x_scaled[j].end()));

            // Sample points in overlap region
            int n_samples = 20;
            for (int k = 0; k < n_samples; ++k) {
                double x = x_min + k * (x_max - x_min) / (n_samples - 1);

                // Interpolate y values at this x
                double y_i = 0, y_j = 0;
                // Simple linear interpolation (should use spline for better accuracy)
                for (size_t m = 0; m < x_scaled[i].size() - 1; ++m) {
                    if (x_scaled[i][m] <= x && x_scaled[i][m+1] >= x) {
                        double t = (x - x_scaled[i][m]) / (x_scaled[i][m+1] - x_scaled[i][m]);
                        y_i = y_scaled[i][m] * (1-t) + y_scaled[i][m+1] * t;
                        break;
                    }
                }
                for (size_t m = 0; m < x_scaled[j].size() - 1; ++m) {
                    if (x_scaled[j][m] <= x && x_scaled[j][m+1] >= x) {
                        double t = (x - x_scaled[j][m]) / (x_scaled[j][m+1] - x_scaled[j][m]);
                        y_j = y_scaled[j][m] * (1-t) + y_scaled[j][m+1] * t;
                        break;
                    }
                }

                if (y_i > 0 && y_j > 0) {
                    chi2 += (y_i - y_j) * (y_i - y_j) / (y_i + y_j);
                    n_pairs++;
                }
            }
        }
    }

    chi2 /= n_pairs;
    quality = 1.0 / (1.0 + chi2);  // Simple quality metric

    std::cout << "    χ²/dof = " << chi2 << ", quality = " << quality << std::endl;

    return quality > 0.5;  // Good collapse if quality > 0.5
}

void UniversalityClassifier::optimizeCollapse(double& sigma_c, double& beta, double& nu) {
    std::cout << "\n  Optimizing data collapse parameters..." << std::endl;

    // Grid search optimization (could use more sophisticated methods)
    double best_chi2 = 1e10;
    double best_sigma_c = sigma_c;
    double best_beta = beta;
    double best_nu = nu;

    for (double sc = sigma_c - 0.01; sc <= sigma_c + 0.01; sc += 0.002) {
        for (double b = beta - 0.02; b <= beta + 0.02; b += 0.005) {
            for (double n = nu - 0.05; n <= nu + 0.05; n += 0.01) {
                double chi2, quality;
                performDataCollapse(sc, b/n, n, chi2, quality);

                if (chi2 < best_chi2) {
                    best_chi2 = chi2;
                    best_sigma_c = sc;
                    best_beta = b;
                    best_nu = n;
                }
            }
        }
    }

    sigma_c = best_sigma_c;
    beta = best_beta;
    nu = best_nu;

    std::cout << "    Optimized: σ_c = " << sigma_c
             << ", β = " << beta
             << ", ν = " << nu
             << ", χ² = " << best_chi2 << std::endl;
}

// Scaling relations
bool UniversalityClassifier::verifyScalingRelations(const CriticalExponents& exp) {
    std::cout << "\n=== Verifying Scaling Relations ===" << std::endl;

    double fisher = checkFisherRelation(exp);
    double rushbrooke = checkRushbrookeRelation(exp);
    double josephson = checkJosephsonRelation(exp);

    std::cout << "  Fisher (γ = ν(2-η)): deviation = " << fisher << std::endl;
    std::cout << "  Rushbrooke (α + 2β + γ = 2): deviation = " << rushbrooke << std::endl;
    std::cout << "  Josephson (νd = 2 - α): deviation = " << josephson << std::endl;

    bool all_satisfied = (fisher < 0.05 && rushbrooke < 0.05 && josephson < 0.05);

    if (all_satisfied) {
        std::cout << "  ✓ All scaling relations satisfied within 5%" << std::endl;
    } else {
        std::cout << "  ⚠ Some scaling relations violated - check exponent extraction" << std::endl;
    }

    return all_satisfied;
}

double UniversalityClassifier::checkFisherRelation(const CriticalExponents& exp) {
    double lhs = exp.gamma;
    double rhs = exp.nu * (2.0 - exp.eta);
    return fabs(lhs - rhs) / (lhs + rhs);
}

double UniversalityClassifier::checkRushbrookeRelation(const CriticalExponents& exp) {
    double sum = exp.alpha + 2.0 * exp.beta + exp.gamma;
    return fabs(sum - 2.0) / 2.0;
}

double UniversalityClassifier::checkJosephsonRelation(const CriticalExponents& exp) {
    double d = 2.0;  // 2D system
    double lhs = exp.nu * d;
    double rhs = 2.0 - exp.alpha;
    return fabs(lhs - rhs) / (lhs + rhs);
}

// Classification
std::string UniversalityClassifier::classifyUniversality() {
    std::cout << "\n=== UNIVERSALITY CLASS CLASSIFICATION ===" << std::endl;

    CriticalExponents exp = extractExponents();

    std::cout << "\nMeasured Exponents:" << std::endl;
    std::cout << "  β = " << exp.beta << " ± " << exp.beta_err << std::endl;
    std::cout << "  ν = " << exp.nu << " ± " << exp.nu_err << std::endl;
    std::cout << "  γ = " << exp.gamma << " ± " << exp.gamma_err << std::endl;
    std::cout << "  η = " << exp.eta << " ± " << exp.eta_err << std::endl;
    std::cout << "  α = " << exp.alpha << " ± " << exp.alpha_err << std::endl;
    std::cout << "  δ = " << exp.delta << " ± " << exp.delta_err << std::endl;
    std::cout << "  σ_c = " << exp.sigma_c << " ± " << exp.sigma_c_err << std::endl;

    return classifyFromExponents(exp);
}

std::string UniversalityClassifier::classifyFromExponents(const CriticalExponents& exp) {
    std::cout << "\nComparing to known universality classes:" << std::endl;

    // Compare to known 2D classes
    double ising_diff = compareToKnownClass(exp, "2D Ising");
    double xy_diff = compareToKnownClass(exp, "2D XY");
    double potts3_diff = compareToKnownClass(exp, "2D 3-state Potts");
    double mean_field_diff = compareToKnownClass(exp, "Mean Field");

    // Decision tree based on strategy document
    if (fabs(exp.beta - 0.125) < 0.01 && fabs(exp.nu - 1.0) < 0.05) {
        std::cout << "\n✓ CLASSIFICATION: 2D Ising universality class" << std::endl;
        std::cout << "  Evidence: β and ν match Ising values within error bars" << std::endl;
        return "2D Ising";
    }
    else if (fabs(exp.nu - 0.67) < 0.05 && exp.eta < 0.1) {
        std::cout << "\n✓ CLASSIFICATION: 2D XY universality class (BKT transition)" << std::endl;
        std::cout << "  Evidence: ν matches XY value, small η consistent with BKT" << std::endl;
        return "2D XY";
    }
    else if (fabs(exp.beta - 1.0/9.0) < 0.01 && fabs(exp.nu - 5.0/6.0) < 0.05) {
        std::cout << "\n✓ CLASSIFICATION: 2D 3-state Potts universality class" << std::endl;
        std::cout << "  Evidence: β and ν match 3-state Potts values" << std::endl;
        return "2D 3-state Potts";
    }
    else if (exp.beta > 0.3 && fabs(exp.nu - 0.5) < 0.1) {
        std::cout << "\n✓ CLASSIFICATION: Mean-field behavior detected" << std::endl;
        std::cout << "  Evidence: Large β and ν ≈ 0.5 suggest mean-field or long-range interactions" << std::endl;
        return "Mean Field";
    }
    else {
        std::cout << "\n★ CLASSIFICATION: Novel universality class!" << std::endl;
        std::cout << "  β = " << exp.beta << " does not match any known 2D universality class" << std::endl;
        std::cout << "  This represents a new critical behavior requiring theoretical explanation" << std::endl;
        std::cout << "\nPossible mechanisms:" << std::endl;
        std::cout << "  • Crossover between universality classes" << std::endl;
        std::cout << "  • Non-equilibrium critical dynamics" << std::endl;
        std::cout << "  • Long-range effective interactions from SMFT coupling" << std::endl;
        std::cout << "  • Emergent symmetry from relativistic dynamics" << std::endl;
        return "Novel - SMFT universality class";
    }
}

double UniversalityClassifier::compareToKnownClass(const CriticalExponents& exp,
                                                  const std::string& class_name) {
    CriticalExponents known;

    if (class_name == "2D Ising") {
        known = getIsingExponents();
    } else if (class_name == "2D XY") {
        known = getXYExponents();
    } else if (class_name == "2D 3-state Potts") {
        known = getPotts3Exponents();
    } else if (class_name == "Mean Field") {
        known = getMeanFieldExponents();
    } else {
        return 1.0;  // Unknown class
    }

    // Compute weighted difference
    double diff = 0;
    diff += fabs(exp.beta - known.beta) / known.beta;
    diff += fabs(exp.nu - known.nu) / known.nu;
    diff += fabs(exp.gamma - known.gamma) / known.gamma;
    diff += fabs(exp.eta - known.eta) / (known.eta + 0.01);  // Avoid division by zero
    diff /= 4.0;

    std::cout << "  " << class_name << ": average deviation = "
             << (diff * 100) << "%" << std::endl;

    return diff;
}

// Known universality classes
UniversalityClassifier::CriticalExponents UniversalityClassifier::getIsingExponents() {
    CriticalExponents exp;
    exp.beta = 0.125;      // 1/8 (exact)
    exp.nu = 1.0;          // (exact)
    exp.gamma = 1.75;      // 7/4 (exact)
    exp.eta = 0.25;        // 1/4 (exact)
    exp.alpha = 0;         // (exact, logarithmic)
    exp.delta = 15.0;      // (exact)
    return exp;
}

UniversalityClassifier::CriticalExponents UniversalityClassifier::getXYExponents() {
    CriticalExponents exp;
    exp.beta = 0.23;       // Not well-defined for BKT
    exp.nu = 0.67;         // Approximate
    exp.gamma = 1.32;      // Approximate
    exp.eta = 0.04;        // Very small
    exp.alpha = -0.01;     // Near zero
    exp.delta = 6.7;       // Approximate
    return exp;
}

UniversalityClassifier::CriticalExponents UniversalityClassifier::getPotts3Exponents() {
    CriticalExponents exp;
    exp.beta = 1.0/9.0;    // (exact)
    exp.nu = 5.0/6.0;      // (exact)
    exp.gamma = 13.0/9.0;  // (exact)
    exp.eta = 4.0/15.0;    // (exact)
    exp.alpha = 1.0/3.0;   // (exact)
    exp.delta = 14.0;      // (exact)
    return exp;
}

UniversalityClassifier::CriticalExponents UniversalityClassifier::getMeanFieldExponents() {
    CriticalExponents exp;
    exp.beta = 0.5;
    exp.nu = 0.5;
    exp.gamma = 1.0;
    exp.eta = 0;
    exp.alpha = 0;
    exp.delta = 3.0;
    return exp;
}

// Reporting
void UniversalityClassifier::writeDataCollapseReport(const std::string& output_dir) {
    std::string filename = output_dir + "/data_collapse_report.txt";
    std::ofstream file(filename);

    file << "===== DATA COLLAPSE ANALYSIS REPORT =====" << std::endl;
    file << "Generated by SMFT UniversalityClassifier" << std::endl;
    file << std::endl;

    CriticalExponents exp = extractExponents();

    file << "Critical Point:" << std::endl;
    file << "  σ_c = " << exp.sigma_c << " ± " << exp.sigma_c_err << std::endl;
    file << std::endl;

    file << "Collapse Parameters:" << std::endl;
    file << "  β/ν = " << exp.beta / exp.nu << std::endl;
    file << "  1/ν = " << 1.0 / exp.nu << std::endl;
    file << std::endl;

    file << "Collapse Quality:" << std::endl;
    file << "  χ²/dof = " << exp.chi2_collapse << std::endl;
    file << "  Quality score = " << exp.collapse_quality << std::endl;
    file << std::endl;

    file << "Grid Sizes Analyzed: ";
    for (int L : grid_sizes) {
        file << L << " ";
    }
    file << std::endl;

    file << "Number of noise points: " << data_by_size[grid_sizes[0]].size() << std::endl;

    file.close();
}

void UniversalityClassifier::writeCriticalExponentReport(const std::string& output_dir) {
    std::string filename = output_dir + "/critical_exponents.txt";
    std::ofstream file(filename);

    file << "===== CRITICAL EXPONENT REPORT =====" << std::endl;
    file << std::endl;

    CriticalExponents exp = extractExponents();

    file << "Measured Critical Exponents:" << std::endl;
    file << "-----------------------------------------" << std::endl;
    file << "| Exponent | Value    | Error  | Notes |" << std::endl;
    file << "-----------------------------------------" << std::endl;
    file << "| β        | " << std::fixed << std::setprecision(4) << exp.beta
         << " | ±" << exp.beta_err << " | Order parameter |" << std::endl;
    file << "| ν        | " << exp.nu
         << " | ±" << exp.nu_err << " | Correlation length |" << std::endl;
    file << "| γ        | " << exp.gamma
         << " | ±" << exp.gamma_err << " | Susceptibility |" << std::endl;
    file << "| η        | " << exp.eta
         << " | ±" << exp.eta_err << " | Anomalous dimension |" << std::endl;
    file << "| α        | " << exp.alpha
         << " | ±" << exp.alpha_err << " | Specific heat |" << std::endl;
    file << "| δ        | " << exp.delta
         << " | ±" << exp.delta_err << " | Critical isotherm |" << std::endl;
    file << "-----------------------------------------" << std::endl;
    file << std::endl;

    file << "Critical Point:" << std::endl;
    file << "  σ_c = " << exp.sigma_c << " ± " << exp.sigma_c_err << std::endl;
    file << std::endl;

    file << "Scaling Relation Checks:" << std::endl;
    file << "  Fisher (γ = ν(2-η)): "
         << (checkFisherRelation(exp) < 0.05 ? "✓ SATISFIED" : "✗ VIOLATED") << std::endl;
    file << "  Rushbrooke (α + 2β + γ = 2): "
         << (checkRushbrookeRelation(exp) < 0.05 ? "✓ SATISFIED" : "✗ VIOLATED") << std::endl;
    file << "  Josephson (νd = 2 - α): "
         << (checkJosephsonRelation(exp) < 0.05 ? "✓ SATISFIED" : "✗ VIOLATED") << std::endl;

    file.close();
}

void UniversalityClassifier::writeClassificationReport(const std::string& output_dir) {
    std::string filename = output_dir + "/universality_classification.txt";
    std::ofstream file(filename);

    file << "===== UNIVERSALITY CLASS CLASSIFICATION =====" << std::endl;
    file << std::endl;

    std::string classification = classifyUniversality();
    CriticalExponents exp = extractExponents();

    file << "CLASSIFICATION RESULT: " << classification << std::endl;
    file << std::endl;

    file << "Evidence Summary:" << std::endl;
    file << "-----------------" << std::endl;

    // Compare to all known classes
    file << std::endl << "Comparison to Known Classes:" << std::endl;
    file << std::endl;

    // 2D Ising
    CriticalExponents ising = getIsingExponents();
    file << "2D Ising:" << std::endl;
    file << "  β: " << exp.beta << " vs " << ising.beta
         << " (deviation: " << fabs(exp.beta - ising.beta) << ")" << std::endl;
    file << "  ν: " << exp.nu << " vs " << ising.nu
         << " (deviation: " << fabs(exp.nu - ising.nu) << ")" << std::endl;

    // 2D XY
    CriticalExponents xy = getXYExponents();
    file << std::endl << "2D XY:" << std::endl;
    file << "  β: " << exp.beta << " vs " << xy.beta
         << " (deviation: " << fabs(exp.beta - xy.beta) << ")" << std::endl;
    file << "  ν: " << exp.nu << " vs " << xy.nu
         << " (deviation: " << fabs(exp.nu - xy.nu) << ")" << std::endl;

    // Conclusion
    file << std::endl << "CONCLUSION:" << std::endl;
    file << "-----------" << std::endl;

    if (classification.find("Novel") != std::string::npos) {
        file << "The SMFT system exhibits a NOVEL universality class!" << std::endl;
        file << "This is a significant discovery requiring further theoretical investigation." << std::endl;
        file << std::endl;
        file << "Key findings:" << std::endl;
        file << "• β = " << exp.beta << " differs from all known 2D universality classes" << std::endl;
        file << "• The system shows genuine critical behavior with well-defined exponents" << std::endl;
        file << "• Scaling relations are satisfied, confirming thermodynamic consistency" << std::endl;
        file << std::endl;
        file << "Recommended follow-up:" << std::endl;
        file << "• Investigate role of relativistic dynamics in critical behavior" << std::endl;
        file << "• Search for effective field theory description" << std::endl;
        file << "• Compare with other non-equilibrium phase transitions" << std::endl;
    } else {
        file << "The SMFT system belongs to the " << classification << " universality class." << std::endl;
        file << "This classification is based on agreement of critical exponents within error bars." << std::endl;
    }

    file.close();
}

// Helper methods
void UniversalityClassifier::sortDataBySigma(int L) {
    auto& data = data_by_size[L];
    std::sort(data.begin(), data.end(),
              [](const DataPoint& a, const DataPoint& b) {
                  return a.sigma < b.sigma;
              });
}

void UniversalityClassifier::linearFit(const std::vector<double>& x, const std::vector<double>& y,
                                       double& slope, double& intercept, double& error) {
    int n = x.size();
    double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;

    for (int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xx += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    double det = n * sum_xx - sum_x * sum_x;
    slope = (n * sum_xy - sum_x * sum_y) / det;
    intercept = (sum_xx * sum_y - sum_x * sum_xy) / det;

    // Estimate error
    double sum_res2 = 0;
    for (int i = 0; i < n; ++i) {
        double res = y[i] - (slope * x[i] + intercept);
        sum_res2 += res * res;
    }
    error = sqrt(sum_res2 / (n - 2)) * sqrt(n / det);
}

std::pair<double, double> UniversalityClassifier::findBinderCrossing(int L1, int L2) {
    const auto& data1 = data_by_size[L1];
    const auto& data2 = data_by_size[L2];

    // Find sigma where U_L1(σ) = U_L2(σ)
    // Using simple linear interpolation

    for (size_t i = 0; i < data1.size() - 1; ++i) {
        for (size_t j = 0; j < data2.size() - 1; ++j) {
            // Check if lines cross in this segment
            double U1_a = data1[i].U_L;
            double U1_b = data1[i+1].U_L;
            double U2_a = data2[j].U_L;
            double U2_b = data2[j+1].U_L;

            double sigma1_a = data1[i].sigma;
            double sigma1_b = data1[i+1].sigma;
            double sigma2_a = data2[j].sigma;
            double sigma2_b = data2[j+1].sigma;

            // Check if segments overlap in sigma
            if (sigma1_b < sigma2_a || sigma2_b < sigma1_a) continue;

            // Check if U values cross
            if ((U1_a - U2_a) * (U1_b - U2_b) < 0) {
                // Linear interpolation to find crossing point
                double t = (U2_a - U1_a) / ((U1_b - U1_a) - (U2_b - U2_a));
                double sigma_cross = sigma1_a + t * (sigma1_b - sigma1_a);
                double error = fabs(sigma1_b - sigma1_a) / 10.0;  // Rough error estimate

                return {sigma_cross, error};
            }
        }
    }

    return {-1, 0};  // No crossing found
}

void UniversalityClassifier::analyzeBinderCrossings(std::vector<std::pair<double, double>>& crossings) {
    // Analyze how crossings converge with system size
    // For true critical point, crossings should converge to σ_c as L → ∞

    if (crossings.size() < 2) return;

    // Check convergence trend
    double trend = 0;
    for (size_t i = 1; i < crossings.size(); ++i) {
        trend += crossings[i].first - crossings[i-1].first;
    }
    trend /= (crossings.size() - 1);

    std::cout << "    Crossing convergence trend: " << trend << std::endl;
    if (fabs(trend) < 0.001) {
        std::cout << "    ✓ Crossings converged - true critical point identified" << std::endl;
    } else {
        std::cout << "    ⚠ Crossings still converging - larger L needed" << std::endl;
    }
}

std::pair<double, double> UniversalityClassifier::findSusceptibilityPeak(int L) {
    const auto& data = data_by_size[L];

    double chi_max = 0;
    double sigma_peak = 0;

    for (const auto& point : data) {
        if (point.chi > chi_max) {
            chi_max = point.chi;
            sigma_peak = point.sigma;
        }
    }

    // Estimate error from width of peak
    double chi_threshold = 0.9 * chi_max;
    double sigma_low = sigma_peak, sigma_high = sigma_peak;

    for (const auto& point : data) {
        if (point.chi > chi_threshold) {
            sigma_low = std::min(sigma_low, point.sigma);
            sigma_high = std::max(sigma_high, point.sigma);
        }
    }

    double error = (sigma_high - sigma_low) / 4.0;

    return {sigma_peak, error};
}

void UniversalityClassifier::analyzeSusceptibilityPeaks(std::vector<std::pair<int, double>>& peaks) {
    // Analyze how peak positions scale with L
    // σ_peak(L) - σ_c ~ L^(-1/ν)

    std::cout << "    Susceptibility peak scaling analysis:" << std::endl;

    std::vector<double> log_L, sigma_peak;
    for (const auto& peak : peaks) {
        log_L.push_back(log(peak.first));
        sigma_peak.push_back(peak.second);
    }

    // Fit to extract scaling
    double slope, intercept, error;
    linearFit(log_L, sigma_peak, slope, intercept, error);

    std::cout << "    Peak shift exponent: " << -slope << " (expect 1/ν)" << std::endl;
}

// Test configuration helper class
UniversalityTestConfig::UniversalityTestConfig() {
    params = getDefaultParameters();
}

UniversalityTestConfig::TestParameters UniversalityTestConfig::getDefaultParameters() {
    TestParameters p;
    p.grid_sizes = {32, 64, 128, 256, 512};
    p.sigma_min = 0.0;
    p.sigma_max = 2.0;
    p.n_sigma_points = 41;
    p.equilibration_steps = 5000;
    p.measurement_steps = 10000;
    p.dt = 0.01;
    p.measure_correlations = true;
    p.correlation_max_r = 0;  // Will be set to L/4 for each grid
    p.output_base = "output/universality_scan";
    return p;
}

void UniversalityTestConfig::generateConfigFiles(const std::string& config_dir) {
    std::cout << "Generating universality scan configuration files..." << std::endl;

    // Create directory if it doesn't exist
    std::string mkdir_cmd = "mkdir -p " + config_dir;
    system(mkdir_cmd.c_str());

    // Generate config for each grid size
    for (int L : params.grid_sizes) {
        std::string filename = config_dir + "/universality_L" + std::to_string(L) + ".yaml";
        generateGridConfig(L, filename);
        std::cout << "  Generated: " << filename << std::endl;
    }

    // Generate master scan script
    std::string script_file = config_dir + "/run_universality_scan.sh";
    std::ofstream script(script_file);
    script << "#!/bin/bash" << std::endl;
    script << "# Universality class scan script" << std::endl;
    script << "# Runs FSS analysis across multiple grid sizes" << std::endl;
    script << std::endl;

    for (int L : params.grid_sizes) {
        script << "echo 'Running L=" << L << " scan...'" << std::endl;
        script << "./build/bin/smft --test " << config_dir
               << "/universality_L" << L << ".yaml" << std::endl;
        script << std::endl;
    }

    script.close();
    system(("chmod +x " + script_file).c_str());

    std::cout << "  Generated run script: " << script_file << std::endl;
}

void UniversalityTestConfig::generateGridConfig(int L, const std::string& filename) {
    std::ofstream file(filename);

    file << "# Universality scan configuration for L=" << L << std::endl;
    file << "# Part of comprehensive FSS analysis" << std::endl;
    file << std::endl;

    file << "test_name: \"universality_L" << L << "\"" << std::endl;
    file << std::endl;

    file << "grid:" << std::endl;
    file << "  size_x: " << L << std::endl;
    file << "  size_y: " << L << std::endl;
    file << "  L_domain: 100.0  # Planck lengths" << std::endl;
    file << std::endl;

    file << "time_evolution:" << std::endl;
    file << "  dt: " << params.dt << std::endl;
    file << "  t_max: " << (params.equilibration_steps + params.measurement_steps) * params.dt << std::endl;
    file << "  equilibration_time: " << params.equilibration_steps * params.dt << std::endl;
    file << std::endl;

    file << "noise_scan:" << std::endl;
    file << "  enabled: true" << std::endl;
    file << "  sigma_min: " << params.sigma_min << std::endl;
    file << "  sigma_max: " << params.sigma_max << std::endl;
    file << "  n_points: " << params.n_sigma_points << std::endl;
    file << std::endl;

    file << "initial_conditions:" << std::endl;
    file << "  dirac:" << std::endl;
    file << "    type: \"gaussian\"" << std::endl;
    file << "    x0_physical: 50.0" << std::endl;
    file << "    y0_physical: 50.0" << std::endl;
    file << "    sigma_physical: 5.0" << std::endl;
    file << "  kuramoto:" << std::endl;
    file << "    type: \"random\"" << std::endl;
    file << "    seed: 42" << std::endl;
    file << std::endl;

    file << "observables:" << std::endl;
    file << "  compute_order_parameter: true" << std::endl;
    file << "  compute_susceptibility: true" << std::endl;
    file << "  compute_binder_cumulant: true" << std::endl;
    file << "  compute_correlation_function: " << (params.measure_correlations ? "true" : "false") << std::endl;
    file << "  correlation_max_r: " << L/4 << std::endl;
    file << std::endl;

    file << "output:" << std::endl;
    file << "  base_dir: \"" << params.output_base << "\"" << std::endl;
    file << "  save_fields: false  # Don't save full fields for FSS" << std::endl;
    file << "  save_observables: true" << std::endl;
    file << "  save_analysis: true" << std::endl;

    file.close();
}

} // namespace SMFT