// test/test_proca_parameter_sweep.cpp
// Systematic Proca parameter sweep: find ANY combination that produces B > 0.01
#include "physics/ProcaEM.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

struct TestResult {
    float alpha;
    float g;
    int steps;
    float B_max;
    float B_rms;
    bool pass;
};

int main() {
    const int nx = 64;
    const int ny = 64;
    const float dx = 1.0f;
    const float dt = 0.01f;

    // Create phase field with vortex
    std::vector<float> theta(nx * ny);
    std::vector<float> R(nx * ny, 1.0f); // Fully synchronized

    float vortex_x = nx / 2.0f;
    float vortex_y = ny / 2.0f;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            float x = i - vortex_x;
            float y = j - vortex_y;
            theta[j * nx + i] = std::atan2(y, x); // Winding W=1
        }
    }

    // Test parameters
    std::vector<float> alphas = {0.1f, 1.0f, 10.0f, 100.0f, 1000.0f};
    std::vector<float> gs = {0.1f, 0.01f, 0.001f};
    std::vector<int> step_counts = {1000, 5000, 10000};

    std::vector<TestResult> results;

    std::cout << "=== PROCA PARAMETER SWEEP ===" << std::endl;
    std::cout << "Grid: " << nx << "×" << ny << std::endl;
    std::cout << "Vortex at: (" << vortex_x << ", " << vortex_y << ")" << std::endl;
    std::cout << "Search space: " << alphas.size() << " × " << gs.size() << " × " << step_counts.size()
              << " = " << (alphas.size() * gs.size() * step_counts.size()) << " combinations" << std::endl;
    std::cout << std::string(80, '=') << std::endl << std::endl;

    bool found_working = false;

    for (float alpha : alphas) {
        for (float g : gs) {
            for (int steps : step_counts) {
                // Create Proca with these parameters
                physics::ProcaEM proca(nx, ny, dx, g, alpha);

                // Evolve
                for (int step = 0; step < steps; ++step) {
                    proca.computePotentials(theta.data(), R.data(), nx, ny, dx, dt);
                }

                // Measure B-field
                proca.computeFieldStrengths();

                float B_max = 0.0f;
                float B_rms = 0.0f;
                int count = 0;

                int r_min = 5, r_max = 15;
                for (int j = ny/2 - r_max; j <= ny/2 + r_max; ++j) {
                    for (int i = nx/2 - r_max; i <= nx/2 + r_max; ++i) {
                        if (i >= 0 && i < nx && j >= 0 && j < ny) {
                            float dx_pos = i - nx/2.0f;
                            float dy_pos = j - ny/2.0f;
                            float r = std::sqrt(dx_pos*dx_pos + dy_pos*dy_pos);

                            if (r >= r_min && r <= r_max) {
                                auto F = proca.getFieldAt(i, j);
                                B_max = std::max(B_max, std::abs(F.Bz));
                                B_rms += F.Bz * F.Bz;
                                ++count;
                            }
                        }
                    }
                }

                if (count > 0) {
                    B_rms = std::sqrt(B_rms / count);
                }

                bool pass = (B_max > 0.01f) || (B_rms > 0.001f);
                results.push_back({alpha, g, steps, B_max, B_rms, pass});

                std::cout << "α=" << std::setw(7) << alpha
                          << " g=" << std::setw(6) << g
                          << " steps=" << std::setw(5) << steps
                          << " → B_max=" << std::scientific << std::setprecision(3) << B_max
                          << " B_rms=" << std::scientific << std::setprecision(3) << B_rms
                          << (pass ? " ✅ PASS" : " ❌")
                          << std::endl;

                if (pass && !found_working) {
                    std::cout << "    *** FOUND WORKING PARAMETERS! ***" << std::endl;
                    found_working = true;
                }
            }
        }
    }

    // Summary
    std::cout << std::endl << std::string(80, '=') << std::endl;
    std::cout << "=== SUMMARY ===" << std::endl;
    std::cout << std::string(80, '=') << std::endl << std::endl;

    int pass_count = 0;
    for (const auto& r : results) {
        if (r.pass) {
            std::cout << "✅ PASS: α=" << std::setw(7) << r.alpha
                     << " g=" << std::setw(6) << r.g
                     << " steps=" << std::setw(5) << r.steps
                     << " B_max=" << std::scientific << std::setprecision(3) << r.B_max
                     << " B_rms=" << std::scientific << std::setprecision(3) << r.B_rms
                     << std::endl;
            pass_count++;
        }
    }

    // Find best result
    TestResult best = results[0];
    for (const auto& r : results) {
        if (r.B_max > best.B_max) {
            best = r;
        }
    }

    std::cout << std::endl;
    std::cout << "Total passing: " << pass_count << " / " << results.size() << std::endl;
    std::cout << "Best B_max: " << std::scientific << std::setprecision(3) << best.B_max
              << " (α=" << best.alpha << ", g=" << best.g << ", steps=" << best.steps << ")" << std::endl;

    // Recommendation
    std::cout << std::endl << std::string(80, '=') << std::endl;
    std::cout << "=== GO/NO-GO RECOMMENDATION ===" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    if (pass_count > 0) {
        std::cout << "✅ SUCCESS: Found " << pass_count << " working parameter combinations" << std::endl;
        std::cout << "   Proca mechanism CAN generate magnetic fields with proper tuning" << std::endl;
        std::cout << "   Recommendation: Proceed to Week 3 with optimized parameters" << std::endl;
        std::cout << "   Suggested: α=" << best.alpha << ", g=" << best.g << ", steps=" << best.steps << std::endl;
        return 0;
    } else {
        std::cout << "❌ FAILURE: No parameter combination achieved B > 0.01" << std::endl;
        std::cout << "   Best achieved: B_max = " << std::scientific << std::setprecision(3) << best.B_max << std::endl;
        std::cout << "   Proca mechanism does not work with tested parameters" << std::endl;
        std::cout << "   Recommendation: Pivot to Stückelberg or Option A" << std::endl;
        return 1;
    }
}
