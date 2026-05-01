/**
 * test_njl_gap_curve.cpp
 *
 * Tests the NJL-gap-equation isomorphism claim:
 *   In NJL, the dynamical fermion mass M_dyn satisfies a self-consistent gap
 *   equation; near the critical coupling G_c it grows as
 *      M_dyn ∝ √(G − G_c)        (mean-field critical exponent 1/2).
 *
 *   The TRD analogue with R = ⟨q̄q⟩ / (some scale) and K = G/N_bond:
 *      m_eff(K) = Δ · R(K)   for fixed Δ.
 *   So the same fit M_eff ∝ √(K − K_c) above K_c is the TRD prediction.
 *
 * Method: read R(K) from the Kuramoto phase-diagram CSV (test_kuramoto_phase_diagram
 * must run first), compute m_eff(K) = Δ · R_mean(K), and fit
 *      m_eff² = a·(K − K_c) + b
 * over K > K_c (mean-field analytic).  Slope a is the prefactor; b should be
 * ≈ 0.  The test reports the fit and a quality measure χ² / dof.
 *
 * Why this is separate from the sweep test: the sweep is slow (14 K values × 3
 * seeds × 3000 steps × 64³ ≈ 20 minutes), and the gap analysis is post-processing
 * — no need to rerun.  If the CSV is missing, this test fails with a clear
 * message asking the user to run the Kuramoto sweep first.
 *
 * Output: output/njl_gap_curve/m_vs_K.csv and a residuals YAML.
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct SweepRow {
    float K, R_mean, R_std;
    int seeds;
};

bool read_kuramoto_csv(const std::string& path, std::vector<SweepRow>& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        SweepRow row;
        char comma;
        ss >> row.K >> comma >> row.R_mean >> comma >> row.R_std >> comma >> row.seeds;
        if (ss) out.push_back(row);
    }
    return !out.empty();
}

}  // namespace

int runNJLGapCurveTest() {
    std::cout << "\n=== NJL Gap Curve M(K) ===\n";
    std::cout << "Tests structural isomorphism: m_eff(K) = Δ·R(K)\n";
    std::cout << "should fit M ∝ √(K − K_c) above K_c (NJL mean-field exponent).\n\n";

    std::vector<SweepRow> data;
    const std::string csv_in = "output/kuramoto_phase_diagram/R_vs_K.csv";
    if (!read_kuramoto_csv(csv_in, data)) {
        std::cerr << "ERROR: " << csv_in << " not found.\n";
        std::cerr << "Run config/kuramoto_phase_diagram.yaml first.\n";
        return 2;
    }

    // Analytic mean-field K_c for Gaussian g(ω, σ=0.1) on lattice z=6:
    const float sigma_omega = 0.1f;
    const float z = 6.0f;
    const float K_c_analytic = sigma_omega * std::sqrt(8.0f / static_cast<float>(M_PI)) / z;

    // Choose Δ such that m_eff = Δ·R has natural unit scale; here Δ = 1.
    const float Delta = 1.0f;

    // Fit m_eff² = a·(K - K_c_fit) for K > K_c_fit. We treat K_c_fit as a
    // free parameter and find the value that minimizes residuals.
    auto fit_for_K_c = [&](float K_c_fit, double& slope, double& intercept,
                           double& chi2, int& dof) {
        // linear regression of y = m² vs x = K - K_c_fit, restricted to K > K_c_fit
        double sx=0, sy=0, sxx=0, sxy=0; int n=0;
        for (auto& r : data) {
            if (r.K <= K_c_fit + 1e-6f) continue;
            const double x = r.K - K_c_fit;
            const double m = static_cast<double>(Delta) * static_cast<double>(r.R_mean);
            const double y = m * m;
            sx += x; sy += y; sxx += x*x; sxy += x*y; ++n;
        }
        if (n < 2) { slope=0; intercept=0; chi2=0; dof=0; return; }
        const double denom = n*sxx - sx*sx;
        slope = (n*sxy - sx*sy) / denom;
        intercept = (sy - slope*sx) / n;
        chi2 = 0.0;
        for (auto& r : data) {
            if (r.K <= K_c_fit + 1e-6f) continue;
            const double x = r.K - K_c_fit;
            const double m = static_cast<double>(Delta) * static_cast<double>(r.R_mean);
            const double y = m * m;
            const double pred = slope * x + intercept;
            chi2 += (y - pred) * (y - pred);
        }
        dof = n - 2;
    };

    // Scan K_c_fit and pick the minimum-chi² value.
    float best_K_c = K_c_analytic;
    double best_slope=0, best_intercept=0, best_chi2=1e30;
    int best_dof = 0;
    const int N_scan = 200;
    const float K_lo = K_c_analytic * 0.1f;
    const float K_hi = K_c_analytic * 5.0f;
    for (int s = 0; s < N_scan; ++s) {
        const float K_c_try = K_lo + (K_hi - K_lo) * static_cast<float>(s) / static_cast<float>(N_scan - 1);
        double slope, intercept, chi2; int dof;
        fit_for_K_c(K_c_try, slope, intercept, chi2, dof);
        if (dof < 2 || slope <= 0.0) continue;
        if (chi2 < best_chi2) {
            best_chi2 = chi2; best_slope = slope; best_intercept = intercept;
            best_dof = dof; best_K_c = K_c_try;
        }
    }

    std::cout << "Analytic mean-field K_c (Gaussian, σ=0.1, z=6): " << std::fixed
              << std::setprecision(5) << K_c_analytic << "\n";
    std::cout << "Best-fit K_c (NJL form):                       " << best_K_c << "\n";
    std::cout << "Slope a in m² = a·(K − K_c):                   " << std::scientific
              << std::setprecision(5) << best_slope << "\n";
    std::cout << "Intercept b (should be ≈ 0):                    " << best_intercept << "\n";
    std::cout << "χ² / dof:                                       " << best_chi2
              << " / " << best_dof << "\n";

    // Write outputs
    std::filesystem::create_directories("output/njl_gap_curve");
    std::ofstream csv("output/njl_gap_curve/m_vs_K.csv");
    csv << "K,m_eff,m_eff_squared,prediction_m_squared\n";
    for (auto& r : data) {
        const double m = Delta * r.R_mean;
        double pred = 0.0;
        if (r.K > best_K_c) pred = best_slope * (r.K - best_K_c) + best_intercept;
        csv << std::fixed << std::setprecision(6) << r.K << ","
            << std::scientific << std::setprecision(8) << m << ","
            << m * m << "," << pred << "\n";
    }

    std::ofstream y("output/njl_gap_curve/summary.yaml");
    y << "# NJL gap-curve fit summary\n";
    y << "input_csv: " << csv_in << "\n";
    y << "Delta: " << Delta << "\n";
    y << "K_c_analytic_mean_field: " << std::fixed << std::setprecision(6) << K_c_analytic << "\n";
    y << "K_c_best_fit: " << best_K_c << "\n";
    y << "slope_a: " << std::scientific << std::setprecision(8) << best_slope << "\n";
    y << "intercept_b: " << best_intercept << "\n";
    y << "chi2: " << best_chi2 << "\n";
    y << "dof: " << best_dof << "\n";
    y << "interpretation: |\n";
    y << "  Fit assumes M_dyn^2 = a·(K - K_c) above K_c (NJL mean-field).\n";
    y << "  Slope a > 0 with intercept ≈ 0 supports the isomorphism.\n";
    y << "  Best-fit K_c close to analytic mean-field K_c is corroborating.\n";

    std::cout << "\nWrote output/njl_gap_curve/m_vs_K.csv\n";
    std::cout << "Wrote output/njl_gap_curve/summary.yaml\n";
    return 0;
}
