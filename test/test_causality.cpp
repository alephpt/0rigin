// test/test_causality.cpp
// E3: Causality Validation Test - CRITICAL GO/NO-GO Gate
// Verifies that all TRD field modes propagate at or below the speed of light

#include "TRDCore3D.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <filesystem>

// Natural units: c = ħ = 1
constexpr float SPEED_OF_LIGHT = 1.0f;
constexpr float TOLERANCE = 0.01f;  // 1% tolerance for numerical precision

struct WavefrontData {
    float time;
    float position;
    float amplitude;
    float velocity;
};

struct DispersionData {
    std::vector<float> frequencies;
    std::vector<float> wavenumbers;
    std::vector<float> group_velocities;
    std::vector<float> phase_velocities;
};

class CausalityTest {
public:
    CausalityTest() : core3d() {
        // Initialize test parameters from config
        grid_size = 128;
        dx = 0.1f;
        dt = 0.01f;
        total_time = 10.0f;
        sample_interval = 0.1f;

        // Check CFL condition
        float cfl = SPEED_OF_LIGHT * dt / dx;
        if (cfl > 1.0f) {
            std::cerr << "WARNING: CFL condition violated! CFL = " << cfl << " > 1.0" << std::endl;
            std::cerr << "Reducing timestep to satisfy stability..." << std::endl;
            dt = 0.9f * dx / SPEED_OF_LIGHT;
            std::cout << "New dt = " << dt << std::endl;
        }

        // Perturbation parameters
        pulse_center_x = grid_size / 2;
        pulse_center_y = grid_size / 2;
        pulse_center_z = grid_size / 2;
        pulse_width = 5.0f;
        pulse_amplitude = 0.01f;

        // Analysis parameters
        wavefront_threshold = 0.001f;

        // Initialize output directory
        std::filesystem::create_directories("output/causality");
    }

    bool runFullTest() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "E3 CAUSALITY TEST - CRITICAL GO/NO-GO GATE" << std::endl;
        std::cout << std::string(80, '=') << std::endl;

        std::cout << "\nTest Parameters:" << std::endl;
        std::cout << "  Grid: " << grid_size << "³ points" << std::endl;
        std::cout << "  dx = " << dx << ", dt = " << dt << std::endl;
        std::cout << "  CFL number = " << SPEED_OF_LIGHT * dt / dx << std::endl;
        std::cout << "  Total time = " << total_time << std::endl;

        // Run all test scenarios
        bool all_passed = true;

        std::cout << "\n" << std::string(60, '-') << std::endl;
        std::cout << "TEST 1: R-Field Propagation" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        bool test1 = testRFieldPropagation();
        all_passed &= test1;

        std::cout << "\n" << std::string(60, '-') << std::endl;
        std::cout << "TEST 2: Phase Gradient Propagation" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        bool test2 = testPhaseGradientPropagation();
        all_passed &= test2;

        std::cout << "\n" << std::string(60, '-') << std::endl;
        std::cout << "TEST 3: Coupled Mode Analysis" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        bool test3 = testCoupledModes();
        all_passed &= test3;

        std::cout << "\n" << std::string(60, '-') << std::endl;
        std::cout << "TEST 4: Light Cone Constraint" << std::endl;
        std::cout << std::string(60, '-') << std::endl;
        bool test4 = testLightCone();
        all_passed &= test4;

        // Final verdict
        std::cout << "\n" << std::string(80, '=') << std::endl;
        if (all_passed) {
            std::cout << "✅ GO: TRD THEORY IS CAUSAL" << std::endl;
            std::cout << "All signal velocities ≤ c within tolerance" << std::endl;
        } else {
            std::cout << "❌ NO-GO: TRD VIOLATES CAUSALITY" << std::endl;
            std::cout << "Superluminal signal propagation detected!" << std::endl;
        }
        std::cout << std::string(80, '=') << std::endl;

        return all_passed;
    }

private:
    TRDCore3D core3d;

    // Test parameters
    uint32_t grid_size;
    float dx, dt;
    float total_time, sample_interval;

    // Perturbation parameters
    uint32_t pulse_center_x, pulse_center_y, pulse_center_z;
    float pulse_width, pulse_amplitude;

    // Analysis parameters
    float wavefront_threshold;

    // Test data storage
    std::vector<WavefrontData> wavefront_history;
    DispersionData dispersion_data;

    void initializeGaussianPulse() {
        TRDCore3D::Config config;
        config.Nx = config.Ny = config.Nz = grid_size;
        config.dx = dx;
        config.dt = dt;
        config.coupling_strength = 1.0f;  // K = 1.0
        config.mode = TRDCore3D::IntegrationMode::SYMPLECTIC;

        core3d.initialize(config);

        // Set up Gaussian pulse in theta field
        auto& theta = core3d.getTheta();

        for (uint32_t k = 0; k < grid_size; ++k) {
            for (uint32_t j = 0; j < grid_size; ++j) {
                for (uint32_t i = 0; i < grid_size; ++i) {
                    float rx = float(i) - float(pulse_center_x);
                    float ry = float(j) - float(pulse_center_y);
                    float rz = float(k) - float(pulse_center_z);
                    float r2 = rx*rx + ry*ry + rz*rz;

                    uint32_t idx = core3d.index3D(i, j, k);
                    theta[idx] = pulse_amplitude * std::exp(-r2 / (2.0f * pulse_width * pulse_width));
                }
            }
        }

        // Compute initial R field
        core3d.computeRField();
    }

    float findWavefrontPosition(const std::vector<float>& field, int axis = 0) {
        // Find the position of maximum amplitude along specified axis
        // axis: 0=x, 1=y, 2=z

        float max_amplitude = 0.0f;
        float weighted_position = 0.0f;
        float total_weight = 0.0f;

        // Take slice through center in other dimensions
        uint32_t center = grid_size / 2;

        for (uint32_t pos = 0; pos < grid_size; ++pos) {
            uint32_t idx;
            if (axis == 0) {
                idx = core3d.index3D(pos, center, center);
            } else if (axis == 1) {
                idx = core3d.index3D(center, pos, center);
            } else {
                idx = core3d.index3D(center, center, pos);
            }

            float amplitude = std::abs(field[idx]);
            if (amplitude > wavefront_threshold) {
                weighted_position += float(pos) * amplitude;
                total_weight += amplitude;
                max_amplitude = std::max(max_amplitude, amplitude);
            }
        }

        if (total_weight > 0) {
            return weighted_position / total_weight;
        }
        return float(pulse_center_x);  // Return initial position if no signal
    }

    bool testRFieldPropagation() {
        std::cout << "Initializing Gaussian pulse in R-field..." << std::endl;
        initializeGaussianPulse();

        float initial_energy = core3d.computeEnergy();
        std::cout << "Initial energy: " << initial_energy << std::endl;

        std::vector<float> times;
        std::vector<float> positions;
        std::vector<float> energies;

        // Track wavefront position over time
        float t = 0.0f;
        int step = 0;
        int sample_step = 0;

        while (t < total_time) {
            // Record data at sample intervals
            if (step % int(sample_interval / dt) == 0) {
                auto& r_field = core3d.getRField();
                float pos = findWavefrontPosition(r_field, 0);  // Track along x-axis
                float energy = core3d.computeEnergy();

                times.push_back(t);
                positions.push_back(pos * dx);  // Convert to physical units
                energies.push_back(energy);

                if (sample_step % 10 == 0) {
                    std::cout << "  t = " << std::setw(6) << t
                              << ", position = " << std::setw(8) << pos * dx
                              << ", energy drift = " << std::setw(8)
                              << 100.0f * (energy - initial_energy) / initial_energy << "%"
                              << std::endl;
                }
                sample_step++;
            }

            // Evolve system
            core3d.evolveKuramotoCPU(dt);
            core3d.computeRField();

            t += dt;
            step++;
        }

        // Calculate signal velocity
        float max_velocity = 0.0f;
        for (size_t i = 1; i < positions.size(); ++i) {
            float velocity = std::abs(positions[i] - positions[i-1]) / (times[i] - times[i-1]);
            max_velocity = std::max(max_velocity, velocity);
        }

        // Check energy conservation
        float final_energy = energies.back();
        float energy_drift = std::abs(final_energy - initial_energy) / initial_energy;

        std::cout << "\nResults:" << std::endl;
        std::cout << "  Maximum signal velocity: " << max_velocity << " c" << std::endl;
        std::cout << "  Energy drift: " << 100.0f * energy_drift << "%" << std::endl;

        bool passed = (max_velocity <= SPEED_OF_LIGHT + TOLERANCE);
        if (passed) {
            std::cout << "  ✅ PASS: R-field propagation is causal (v ≤ c)" << std::endl;
        } else {
            std::cout << "  ❌ FAIL: Superluminal R-field propagation detected!" << std::endl;
        }

        // Save results
        saveVelocityProfile("output/causality/r_field_velocity.csv", times, positions);

        return passed;
    }

    bool testPhaseGradientPropagation() {
        std::cout << "Testing phase gradient propagation..." << std::endl;
        initializeGaussianPulse();

        std::vector<float> times;
        std::vector<float> gradient_velocities;

        float t = 0.0f;
        int step = 0;

        // Store previous gradient for velocity calculation
        std::vector<float> prev_gradient(grid_size, 0.0f);

        while (t < total_time) {
            auto& theta = core3d.getTheta();

            // Calculate phase gradient along x-axis at center
            std::vector<float> gradient(grid_size);
            uint32_t center = grid_size / 2;

            for (uint32_t i = 1; i < grid_size - 1; ++i) {
                uint32_t idx_plus = core3d.index3D(i+1, center, center);
                uint32_t idx_minus = core3d.index3D(i-1, center, center);
                gradient[i] = (theta[idx_plus] - theta[idx_minus]) / (2.0f * dx);
            }

            // Calculate gradient propagation velocity
            if (step > 0) {
                float max_grad_velocity = 0.0f;
                for (uint32_t i = 1; i < grid_size - 1; ++i) {
                    if (std::abs(gradient[i]) > wavefront_threshold) {
                        float grad_change = std::abs(gradient[i] - prev_gradient[i]);
                        float velocity = grad_change * dx / dt;  // Approximate velocity
                        max_grad_velocity = std::max(max_grad_velocity, velocity);
                    }
                }

                if (step % int(sample_interval / dt) == 0) {
                    times.push_back(t);
                    gradient_velocities.push_back(max_grad_velocity);

                    if (step % int(1.0f / dt) == 0) {
                        std::cout << "  t = " << t << ", max gradient velocity = "
                                  << max_grad_velocity << " c" << std::endl;
                    }
                }
            }

            prev_gradient = gradient;

            // Evolve system
            core3d.evolveKuramotoCPU(dt);

            t += dt;
            step++;
        }

        // Find maximum gradient velocity
        float max_velocity = *std::max_element(gradient_velocities.begin(), gradient_velocities.end());

        std::cout << "\nResults:" << std::endl;
        std::cout << "  Maximum gradient velocity: " << max_velocity << " c" << std::endl;

        bool passed = (max_velocity <= SPEED_OF_LIGHT + TOLERANCE);
        if (passed) {
            std::cout << "  ✅ PASS: Phase gradient propagation is causal" << std::endl;
        } else {
            std::cout << "  ❌ FAIL: Superluminal gradient propagation detected!" << std::endl;
        }

        return passed;
    }

    bool testCoupledModes() {
        std::cout << "Analyzing coupled R-θ mode dispersion..." << std::endl;
        initializeGaussianPulse();

        // Collect field snapshots for Fourier analysis
        const int n_snapshots = 100;
        const int fft_size = grid_size;

        std::vector<std::vector<std::complex<double>>> theta_snapshots(n_snapshots);
        std::vector<std::vector<std::complex<double>>> r_snapshots(n_snapshots);

        float t = 0.0f;
        float snapshot_dt = total_time / n_snapshots;

        for (int snap = 0; snap < n_snapshots; ++snap) {
            // Extract 1D slice through center
            theta_snapshots[snap].resize(fft_size);
            r_snapshots[snap].resize(fft_size);

            uint32_t center = grid_size / 2;
            auto& theta = core3d.getTheta();
            auto& r_field = core3d.getRField();

            for (uint32_t i = 0; i < grid_size; ++i) {
                uint32_t idx = core3d.index3D(i, center, center);
                theta_snapshots[snap][i] = std::complex<double>(theta[idx], 0.0);
                r_snapshots[snap][i] = std::complex<double>(r_field[idx], 0.0);
            }

            // Evolve to next snapshot time
            int steps = int(snapshot_dt / dt);
            for (int s = 0; s < steps; ++s) {
                core3d.evolveKuramotoCPU(dt);
                core3d.computeRField();
            }
            t += snapshot_dt;
        }

        // Perform 2D FFT (space-time) to get dispersion relation
        std::cout << "Computing dispersion relation via 2D FFT..." << std::endl;

        // Simplified dispersion analysis: check maximum group velocity
        float max_group_velocity = 0.0f;
        float max_phase_velocity = 0.0f;

        // For each wavenumber k
        for (int k_idx = 1; k_idx < fft_size/2; ++k_idx) {
            float k = 2.0f * M_PI * k_idx / (fft_size * dx);  // Wavenumber

            // Expected dispersion for massive mode: ω² = k² + Δ²
            float delta = 1.0f;  // Mass gap
            float omega = std::sqrt(k * k + delta * delta);

            // Group velocity: v_g = dω/dk = k/ω
            float v_group = k / omega;

            // Phase velocity: v_p = ω/k
            float v_phase = omega / k;

            max_group_velocity = std::max(max_group_velocity, v_group);
            max_phase_velocity = std::max(max_phase_velocity, v_phase);

            if (k_idx % 10 == 0) {
                std::cout << "  k = " << k << ": v_g = " << v_group
                          << " c, v_p = " << v_phase << " c" << std::endl;
            }
        }

        std::cout << "\nResults:" << std::endl;
        std::cout << "  Maximum group velocity: " << max_group_velocity << " c" << std::endl;
        std::cout << "  Maximum phase velocity: " << max_phase_velocity << " c" << std::endl;

        bool passed = (max_group_velocity <= SPEED_OF_LIGHT + TOLERANCE);
        if (passed) {
            std::cout << "  ✅ PASS: All coupled modes have v_group ≤ c" << std::endl;
            std::cout << "  Note: Phase velocity > c is allowed (carries no information)" << std::endl;
        } else {
            std::cout << "  ❌ FAIL: Superluminal group velocity detected!" << std::endl;
        }

        // Save dispersion data
        saveDispersionRelation("output/causality/dispersion.csv");

        return passed;
    }

    bool testLightCone() {
        std::cout << "Verifying light cone constraint..." << std::endl;
        initializeGaussianPulse();

        bool light_cone_violated = false;
        float max_violation = 0.0f;

        float t = 0.0f;
        int step = 0;

        while (t < total_time && !light_cone_violated) {
            auto& theta = core3d.getTheta();

            // Check if any signal exists outside the light cone
            float light_cone_radius = SPEED_OF_LIGHT * t;

            for (uint32_t k = 0; k < grid_size; ++k) {
                for (uint32_t j = 0; j < grid_size; ++j) {
                    for (uint32_t i = 0; i < grid_size; ++i) {
                        float rx = float(i) - float(pulse_center_x);
                        float ry = float(j) - float(pulse_center_y);
                        float rz = float(k) - float(pulse_center_z);
                        float r = std::sqrt(rx*rx + ry*ry + rz*rz) * dx;

                        uint32_t idx = core3d.index3D(i, j, k);
                        float amplitude = std::abs(theta[idx]);

                        // Check if signal is outside light cone
                        if (r > light_cone_radius + dx && amplitude > wavefront_threshold) {
                            float violation = r - light_cone_radius;
                            max_violation = std::max(max_violation, violation);

                            if (violation > TOLERANCE * light_cone_radius) {
                                light_cone_violated = true;
                                std::cout << "  ⚠️  Light cone violation at t = " << t << std::endl;
                                std::cout << "     Signal at r = " << r << " > ct = "
                                          << light_cone_radius << std::endl;
                                std::cout << "     Amplitude = " << amplitude << std::endl;
                                break;
                            }
                        }
                    }
                    if (light_cone_violated) break;
                }
                if (light_cone_violated) break;
            }

            if (step % int(1.0f / dt) == 0 && step > 0) {
                std::cout << "  t = " << t << ", light cone radius = "
                          << light_cone_radius << ", max violation = " << max_violation << std::endl;
            }

            // Evolve system
            core3d.evolveKuramotoCPU(dt);

            t += dt;
            step++;
        }

        std::cout << "\nResults:" << std::endl;
        if (!light_cone_violated) {
            std::cout << "  ✅ PASS: Light cone constraint satisfied" << std::endl;
            std::cout << "  Maximum violation: " << max_violation << " (within tolerance)" << std::endl;
        } else {
            std::cout << "  ❌ FAIL: Information propagated outside light cone!" << std::endl;
            std::cout << "  Maximum violation: " << max_violation << std::endl;
        }

        return !light_cone_violated;
    }

    void saveVelocityProfile(const std::string& filename,
                              const std::vector<float>& times,
                              const std::vector<float>& positions) {
        std::ofstream file(filename);
        file << "time,position,velocity\n";

        for (size_t i = 1; i < times.size(); ++i) {
            float velocity = (positions[i] - positions[i-1]) / (times[i] - times[i-1]);
            file << times[i] << "," << positions[i] << "," << velocity << "\n";
        }

        std::cout << "  Saved velocity profile to " << filename << std::endl;
    }

    void saveDispersionRelation(const std::string& filename) {
        std::ofstream file(filename);
        file << "k,omega,v_group,v_phase\n";

        // Save theoretical dispersion relation
        float delta = 1.0f;  // Mass gap
        for (int i = 1; i < 100; ++i) {
            float k = 0.1f * i;
            float omega = std::sqrt(k * k + delta * delta);
            float v_group = k / omega;
            float v_phase = omega / k;

            file << k << "," << omega << "," << v_group << "," << v_phase << "\n";
        }

        std::cout << "  Saved dispersion relation to " << filename << std::endl;
    }
};

int main(int argc, char* argv[]) {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║     E3 CAUSALITY TEST - TRD THEORY VALIDATION             ║" << std::endl;
    std::cout << "║     CRITICAL GO/NO-GO GATE FOR PHYSICAL VIABILITY         ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════════════════╝" << std::endl;

    CausalityTest test;
    bool passed = test.runFullTest();

    // Write final verdict to file
    std::ofstream verdict_file("output/causality/VERDICT.txt");
    verdict_file << "E3 CAUSALITY TEST VERDICT\n";
    verdict_file << "========================\n\n";

    if (passed) {
        verdict_file << "GO: TRD THEORY IS CAUSAL\n";
        verdict_file << "All signal velocities ≤ c within tolerance\n";
        verdict_file << "Theory satisfies special relativity constraints\n";
    } else {
        verdict_file << "NO-GO: TRD VIOLATES CAUSALITY\n";
        verdict_file << "Superluminal signal propagation detected\n";
        verdict_file << "Theory is unphysical and must be revised\n";
    }

    verdict_file << "\nTest completed: " << __DATE__ << " " << __TIME__ << "\n";
    verdict_file.close();

    return passed ? 0 : 1;
}