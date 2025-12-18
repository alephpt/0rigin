/**
 * @file test_dirac_coupling.cpp
 * @brief Dirac Coupling Experiment: 5 Success Criteria Validation
 *
 * From Dirac-Anomaly.md Section VI:
 * 1. Localization at defects (O > 0.7)
 * 2. Stabilization of deep defects (ΔS > 10%)
 * 3. Discrete energy levels (2-5 peaks in histogram)
 * 4. Particle count ∈ [10, 200]
 * 5. Long-time stability (τ > 500 steps)
 *
 * Outputs: build/output/dirac_coupling/
 */

#include "../lib/Nova/Nova.h"
#include "../src/MSFTEngine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <algorithm>

void create_output_dirs() {
    mkdir("build/output", 0755);
    mkdir("build/output/dirac_coupling", 0755);
}

void save_field_2d(const std::vector<float>& field, int Nx, int Ny,
                   const std::string& filename) {
    std::ofstream f(filename);
    if (!f) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            f << field[y * Nx + x];
            if (x < Nx - 1) f << " ";
        }
        f << "\n";
    }
    f.close();
}

void save_timeseries(const std::vector<float>& data, const std::string& filename) {
    std::ofstream f(filename);
    if (!f) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        f << i << " " << data[i] << "\n";
    }
    f.close();
}

// Criterion 1: Localization at defects
float compute_overlap(const std::vector<float>& R_field,
                      const std::vector<float>& psi_density) {
    float numerator = 0.0f;
    float denominator = 0.0f;

    for (size_t i = 0; i < R_field.size(); i++) {
        float defect_indicator = 1.0f - R_field[i];  // High where R is low
        numerator += defect_indicator * psi_density[i];
        denominator += psi_density[i];
    }

    return (denominator > 1e-10) ? (numerator / denominator) : 0.0f;
}

// Criterion 2: Defect counting
struct DefectStats {
    int total_defects;
    int deep_defects;  // R_core < 0.3
};

DefectStats count_defects(const std::vector<float>& R_field, float threshold = 0.5) {
    DefectStats stats = {0, 0};

    for (float r : R_field) {
        if (r < threshold) {
            stats.total_defects++;
            if (r < 0.3) {
                stats.deep_defects++;
            }
        }
    }

    return stats;
}

// Criterion 3: Energy histogram
std::vector<float> compute_local_energies(const std::vector<float>& psi_density,
                                          const std::vector<float>& R_field,
                                          int Nx, int Ny, float Delta) {
    std::vector<float> energies;

    // Find peaks in psi_density (localized particles)
    float threshold = 0.1f * (*std::max_element(psi_density.begin(), psi_density.end()));

    for (int y = 1; y < Ny - 1; y++) {
        for (int x = 1; x < Nx - 1; x++) {
            int idx = y * Nx + x;
            float rho = psi_density[idx];

            // Check if this is a local maximum (peak)
            bool is_peak = true;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int neighbor_idx = (y + dy) * Nx + (x + dx);
                    if (psi_density[neighbor_idx] > rho) {
                        is_peak = false;
                        break;
                    }
                }
                if (!is_peak) break;
            }

            // If peak above threshold, compute local energy
            if (is_peak && rho > threshold) {
                // E ≈ m(x)·c² = Δ·R(x)
                // (Simplified: ignoring kinetic term for now)
                float mass = Delta * R_field[idx];
                float energy = mass;  // In natural units c=1
                energies.push_back(energy);
            }
        }
    }

    return energies;
}

// Criterion 4: Particle counting
int count_particles(const std::vector<float>& psi_density, int Nx, int Ny) {
    float threshold = 0.1f * (*std::max_element(psi_density.begin(), psi_density.end()));
    int count = 0;

    for (int y = 1; y < Ny - 1; y++) {
        for (int x = 1; x < Nx - 1; x++) {
            int idx = y * Nx + x;
            float rho = psi_density[idx];

            // Check if local maximum
            bool is_peak = true;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int neighbor_idx = (y + dy) * Nx + (x + dx);
                    if (psi_density[neighbor_idx] > rho) {
                        is_peak = false;
                        break;
                    }
                }
                if (!is_peak) break;
            }

            if (is_peak && rho > threshold) {
                count++;
            }
        }
    }

    return count;
}

int main() {
    std::cout << "=== MSFT Dirac Coupling Experiment ===" << std::endl;
    std::cout << "Testing 5 Success Criteria from Dirac-Anomaly.md\n" << std::endl;

    create_output_dirs();

    // Test different coupling strengths
    std::vector<float> lambda_values = {0.1, 1.0, 10.0};

    // Simulation parameters
    const int Nx = 256;
    const int Ny = 256;
    const float dt = 0.01f;
    const float K = 27.21f;
    const float Delta = 2.5f;

    for (float lambda : lambda_values) {
        std::cout << "\n=== Testing λ_D = " << lambda << " ===" << std::endl;

        // Create output directory
        char dirname[256];
        snprintf(dirname, sizeof(dirname), "build/output/dirac_coupling/lambda_%.1f", lambda);
        mkdir(dirname, 0755);

        // Initialize Nova (compute-only mode)
        NovaConfig config = {
            .name = "MSFT Dirac Coupling",
            .screen = {800, 600},
            .debug_level = "info",
            .dimensions = "2D",
            .camera_type = "orthographic",
            .compute = true,
        };
        Nova nova(config);
        nova.initialized = true;

        // Initialize MSFT engine
        MSFTEngine engine(&nova);
        engine.initialize(Nx, Ny, Delta, 0.0f);

        // Phase 1: Vacuum equilibration (500 steps, no Dirac)
        std::cout << "\nPhase 1: Vacuum equilibration..." << std::flush;

        // Set random initial phases
        std::vector<float> theta_init(Nx * Ny);
        srand(42);
        for (int i = 0; i < Nx * Ny; i++) {
            theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
        }
        engine.setInitialPhases(theta_init);

        // Equilibrate vacuum
        for (int step = 0; step < 500; step++) {
            engine.step(dt, K, 0.1f);
            if (step % 100 == 0) std::cout << "." << std::flush;
        }
        std::cout << " done" << std::endl;

        // Get equilibrated vacuum state
        std::vector<float> R_initial = engine.getSyncField();
        DefectStats defects_initial = count_defects(R_initial);

        std::cout << "Initial defects: " << defects_initial.total_defects
                  << " (deep: " << defects_initial.deep_defects << ")" << std::endl;

        // Phase 2: Initialize Dirac field (Gaussian at defect centers)
        std::cout << "\nPhase 2: Initializing Dirac field..." << std::endl;

        // TODO: Implement initializeDiracField() method in MSFTEngine
        // For now, this is a placeholder

        // Phase 3: Coupled evolution (2000 steps)
        std::cout << "\nPhase 3: Coupled evolution..." << std::flush;

        std::vector<float> overlap_series;
        std::vector<float> particle_count_series;
        std::vector<float> total_density_series;
        std::vector<DefectStats> defect_history;

        for (int step = 0; step < 2000; step++) {
            // Step with Dirac coupling
            // TODO: Implement stepWithDirac() method
            engine.step(dt, K, 0.1f);  // Placeholder

            // Measurements every 10 steps
            if (step % 10 == 0) {
                // Get fields
                std::vector<float> R_field = engine.getSyncField();
                // TODO: Get psi_density from engine.getDiracDensity()
                std::vector<float> psi_density(Nx * Ny, 0.0f);  // Placeholder

                // Criterion 1: Localization
                float overlap = compute_overlap(R_field, psi_density);
                overlap_series.push_back(overlap);

                // Criterion 2: Defect survival
                DefectStats defects = count_defects(R_field);
                defect_history.push_back(defects);

                // Criterion 4: Particle count
                int N_particles = count_particles(psi_density, Nx, Ny);
                particle_count_series.push_back(N_particles);

                // Criterion 5: Total density (for stability check)
                float total_density = 0.0f;
                for (float rho : psi_density) total_density += rho;
                total_density_series.push_back(total_density);

                if (step % 500 == 0) {
                    std::cout << "\n  Step " << step << ":"
                              << " O=" << overlap
                              << " N=" << N_particles
                              << " defects=" << defects.total_defects;
                }
            }

            // Save snapshots
            if (step == 0 || step == 500 || step == 1000 || step == 1500 || step == 1999) {
                char stepdir[512];
                snprintf(stepdir, sizeof(stepdir), "%s/step_%04d", dirname, step);
                mkdir(stepdir, 0755);

                std::vector<float> R_field = engine.getSyncField();
                save_field_2d(R_field, Nx, Ny, std::string(stepdir) + "/R_field.dat");

                // TODO: Save psi_density
            }
        }
        std::cout << "\ndone" << std::endl;

        // Final analysis
        std::cout << "\n=== Results for λ_D = " << lambda << " ===" << std::endl;

        // Criterion 1: Localization
        float overlap_mean = 0.0f;
        for (float o : overlap_series) overlap_mean += o;
        overlap_mean /= overlap_series.size();

        std::cout << "\nCriterion 1 (Localization): O = " << overlap_mean;
        if (overlap_mean > 0.7) {
            std::cout << " ✓ PASS (>0.7)" << std::endl;
        } else if (overlap_mean > 0.6) {
            std::cout << " ~ MARGINAL (>0.6)" << std::endl;
        } else {
            std::cout << " ✗ FAIL (<0.6)" << std::endl;
        }

        // Criterion 2: Stabilization
        DefectStats defects_final = defect_history.back();
        float survival_rate = float(defects_final.deep_defects) / float(defects_initial.deep_defects);
        float delta_survival = survival_rate - 0.706f;  // From vacuum-only baseline

        std::cout << "\nCriterion 2 (Stabilization): ΔS = " << (delta_survival * 100) << "%";
        if (delta_survival > 0.10) {
            std::cout << " ✓ PASS (>10%)" << std::endl;
        } else if (delta_survival > 0.05) {
            std::cout << " ~ MARGINAL (>5%)" << std::endl;
        } else {
            std::cout << " ✗ FAIL (<5%)" << std::endl;
        }

        // Criterion 3: Discrete energy levels
        std::vector<float> R_final = engine.getSyncField();
        std::vector<float> psi_density_final(Nx * Ny, 0.0f);  // TODO: Get from engine
        std::vector<float> energies = compute_local_energies(psi_density_final, R_final, Nx, Ny, Delta);

        std::cout << "\nCriterion 3 (Energy Levels): " << energies.size() << " bound states found";
        // TODO: Histogram analysis to count peaks
        std::cout << " (histogram analysis pending)" << std::endl;

        // Criterion 4: Particle count
        int N_particles_final = particle_count_series.empty() ? 0 : particle_count_series.back();
        std::cout << "\nCriterion 4 (Particle Count): N = " << N_particles_final;
        if (N_particles_final >= 20 && N_particles_final <= 150) {
            std::cout << " ✓ PASS ([20,150])" << std::endl;
        } else if (N_particles_final >= 10 && N_particles_final <= 200) {
            std::cout << " ~ MARGINAL ([10,200])" << std::endl;
        } else {
            std::cout << " ✗ FAIL (outside range)" << std::endl;
        }

        // Criterion 5: Stability
        // Check if total density remains constant (no exponential decay)
        if (total_density_series.size() > 100) {
            float density_initial = total_density_series[10];  // After transient
            float density_final = total_density_series.back();
            float decay_ratio = density_final / density_initial;

            std::cout << "\nCriterion 5 (Stability): decay ratio = " << decay_ratio;
            if (decay_ratio > 0.9) {
                std::cout << " ✓ PASS (>0.9, stable)" << std::endl;
            } else if (decay_ratio > 0.7) {
                std::cout << " ~ MARGINAL (>0.7)" << std::endl;
            } else {
                std::cout << " ✗ FAIL (<0.7, unstable)" << std::endl;
            }
        }

        // Save time series
        save_timeseries(overlap_series, std::string(dirname) + "/overlap_timeseries.dat");
        save_timeseries(particle_count_series, std::string(dirname) + "/particle_count_timeseries.dat");
        save_timeseries(total_density_series, std::string(dirname) + "/total_density_timeseries.dat");

        // Save energy histogram
        if (!energies.empty()) {
            std::ofstream f(std::string(dirname) + "/energy_spectrum.dat");
            for (float E : energies) {
                f << E << "\n";
            }
            f.close();
        }
    }

    std::cout << "\n=== Dirac Coupling Experiment Complete ===" << std::endl;
    std::cout << "Results saved to: build/output/dirac_coupling/" << std::endl;
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "1. Analyze energy_spectrum.dat for discrete peaks" << std::endl;
    std::cout << "2. Visualize psi_density overlaid on R_field" << std::endl;
    std::cout << "3. Compare defect survival rates" << std::endl;
    std::cout << "4. Check all 5 criteria → PASS/FAIL decision" << std::endl;

    return 0;
}
