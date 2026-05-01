/**
 * test_kuramoto_phase_diagram.cpp
 *
 * R(K) sweep on the 3D lattice Kuramoto model with our actual disorder
 * distribution g(ω) = N(0, σ=0.1).  Locates the synchronization transition
 * numerically and compares to the mean-field critical coupling for a
 * Gaussian g.
 *
 * Why this test exists: docs/paper/TRD_Paper.md §2.1 asserts that R is the
 * order parameter of the vacuum and behaves like the Higgs modulus, but the
 * paper has never actually plotted R(K) and never measured K_c.  The first
 * external critic gave a Lorentzian-K_c formula that does not apply to our
 * Gaussian distribution.  The right answer is the numerical phase diagram.
 *
 * Mean-field K_c for Gaussian g(ω, σ) on a lattice with coordination z=6:
 *   K_c = (2 / (π z g(0))) = (2 / (π z)) · σ · √(2π) = σ · √(8π) / z
 * For σ=0.1, z=6:  K_c ≈ 0.0836.
 *
 * That is the *infinite-system mean-field* prediction.  On a finite 64³
 * lattice with random ω one expects a finite-size-rounded transition.
 * The test reports both the analytic K_c and the measured K at which
 * R first crosses 0.5.
 *
 * Output: output/kuramoto_phase_diagram/R_vs_K.csv
 */

#include "TRDCore3D.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

struct SweepPoint {
    float K;
    float R_mean;
    float R_std;
    int seeds_used;
};

float run_to_steady_state(float K, uint32_t seed,
                          uint32_t Nx, uint32_t Ny, uint32_t Nz,
                          float dt, int n_burnin, int n_average) {
    TRDCore3D core;
    TRDCore3D::Config cfg;
    cfg.Nx = Nx; cfg.Ny = Ny; cfg.Nz = Nz;
    cfg.dx = 1.0f;
    cfg.dt = dt;
    cfg.coupling_strength = K;
    cfg.mode = TRDCore3D::IntegrationMode::SYMPLECTIC;
    cfg.local_R_radius = 0;  // global R
    core.initialize(cfg);
    core.initializeRandom(seed);

    // Burn-in
    for (int n = 0; n < n_burnin; ++n) {
        core.evolveSymplecticCPU(dt);
    }
    // Average R over the last n_average steps
    double R_sum = 0.0;
    for (int n = 0; n < n_average; ++n) {
        core.evolveSymplecticCPU(dt);
        core.computeRField();
        R_sum += core.getAverageR();
    }
    return static_cast<float>(R_sum / std::max(1, n_average));
}

}  // namespace

int runKuramotoPhaseDiagramTest() {
    std::cout << "\n=== Kuramoto R(K) Phase Diagram ===\n";
    std::cout << "Sweeping coupling K to locate the synchronization transition\n";
    std::cout << "with our actual g(ω) = N(0, σ=0.1) and lattice z=6.\n\n";

    // Sweep configuration
    const std::vector<float> K_values = {
        0.005f, 0.01f, 0.02f, 0.04f, 0.06f, 0.08f, 0.10f,
        0.15f, 0.20f, 0.30f, 0.50f, 1.00f, 2.00f, 4.00f
    };
    const std::vector<uint32_t> seeds = {42, 137, 1729};
    const uint32_t Nx = 64, Ny = 64, Nz = 64;
    const float dt = 0.01f;
    const int n_burnin = 2000;
    const int n_average = 1000;

    // Analytic mean-field K_c for Gaussian g(ω) with σ:
    //   K_c = 2 / (π · z · g(0))     where g(0) = 1/(σ√(2π))
    //   K_c = σ · √(8/π) / z
    const float sigma_omega = 0.1f;  // matches TRDCore3D::initializeRandom default
    const float z = 6.0f;
    const float K_c_analytic = sigma_omega * std::sqrt(8.0f / static_cast<float>(M_PI)) / z;
    std::cout << "Analytic mean-field K_c (Gaussian, σ=" << sigma_omega
              << ", z=" << z << "): " << std::fixed << std::setprecision(5)
              << K_c_analytic << "\n\n";

    std::cout << std::left << std::setw(10) << "K"
              << std::setw(14) << "R_mean"
              << std::setw(14) << "R_std"
              << std::setw(8)  << "seeds"
              << "\n";
    std::cout << std::string(46, '-') << "\n";

    std::vector<SweepPoint> sweep;
    for (float K : K_values) {
        std::vector<float> R_per_seed;
        for (uint32_t seed : seeds) {
            float R = run_to_steady_state(K, seed, Nx, Ny, Nz, dt, n_burnin, n_average);
            R_per_seed.push_back(R);
        }
        double mean = 0.0;
        for (float r : R_per_seed) mean += r;
        mean /= R_per_seed.size();
        double var = 0.0;
        for (float r : R_per_seed) var += (r - mean) * (r - mean);
        var /= std::max<size_t>(1, R_per_seed.size() - 1);
        SweepPoint p;
        p.K = K;
        p.R_mean = static_cast<float>(mean);
        p.R_std = static_cast<float>(std::sqrt(var));
        p.seeds_used = static_cast<int>(R_per_seed.size());
        sweep.push_back(p);

        std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(4) << K
                  << std::setw(14) << std::setprecision(6) << p.R_mean
                  << std::setw(14) << std::setprecision(6) << p.R_std
                  << std::setw(8)  << p.seeds_used
                  << "\n";
    }

    // Locate the K at which R first crosses 0.5 (a conventional threshold)
    float K_c_measured = -1.0f;
    for (size_t i = 1; i < sweep.size(); ++i) {
        if (sweep[i-1].R_mean < 0.5f && sweep[i].R_mean >= 0.5f) {
            float a = sweep[i-1].R_mean, b = sweep[i].R_mean;
            float Ka = sweep[i-1].K,    Kb = sweep[i].K;
            K_c_measured = Ka + (0.5f - a) * (Kb - Ka) / (b - a);
            break;
        }
    }

    std::cout << "\nMeasured K_c (R first crosses 0.5): ";
    if (K_c_measured > 0.0f) {
        std::cout << std::fixed << std::setprecision(5) << K_c_measured
                  << "  (analytic mean-field: " << K_c_analytic << ")\n";
        std::cout << "Ratio measured/analytic: " << (K_c_measured / K_c_analytic) << "\n";
    } else {
        std::cout << "NOT REACHED in swept range\n";
    }

    // Write CSV
    std::filesystem::create_directories("output/kuramoto_phase_diagram");
    std::ofstream csv("output/kuramoto_phase_diagram/R_vs_K.csv");
    csv << "K,R_mean,R_std,seeds_used\n";
    for (auto& p : sweep) {
        csv << std::fixed << std::setprecision(6)
            << p.K << "," << p.R_mean << "," << p.R_std << "," << p.seeds_used << "\n";
    }
    csv.close();

    // Write summary YAML
    std::ofstream y("output/kuramoto_phase_diagram/summary.yaml");
    y << "# Kuramoto phase-diagram summary\n";
    y << "lattice: { Nx: " << Nx << ", Ny: " << Ny << ", Nz: " << Nz << " }\n";
    y << "disorder: { distribution: gaussian, sigma: " << sigma_omega << " }\n";
    y << "coordination_number: " << z << "\n";
    y << "K_c_analytic_mean_field: " << K_c_analytic << "\n";
    y << "K_c_measured_R_crosses_half: ";
    if (K_c_measured > 0.0f) y << K_c_measured << "\n";
    else y << "null  # transition not reached in swept range\n";
    y << "n_burnin_steps: " << n_burnin << "\n";
    y << "n_average_steps: " << n_average << "\n";
    y << "n_seeds: " << seeds.size() << "\n";
    y.close();

    std::cout << "\nWrote output/kuramoto_phase_diagram/R_vs_K.csv\n";
    std::cout << "Wrote output/kuramoto_phase_diagram/summary.yaml\n";
    return 0;
}
