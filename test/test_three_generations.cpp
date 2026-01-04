/**
 * test_three_generations.cpp
 *
 * B3 TEST: Three-Generation Structure from TRD Topology
 *
 * HYPOTHESIS: Topological defect classification in 3D yields exactly 3 fermion families
 *
 * THEORETICAL FRAMEWORK:
 *   - Point defects (0D): Charge quantization Q = 1,2,3 → electron, muon, tau?
 *   - Line defects (1D): Flux tube configurations → quark generations?
 *   - Surface defects (2D): Domain wall structures → neutrino families?
 *   - Volume defects (3D): Phase transitions → generation mixing?
 *
 * APPROACH:
 *   1. Create topological defects with different winding numbers
 *   2. Evolve via TRDCore3D Kuramoto dynamics
 *   3. Classify stable configurations by topological invariants
 *   4. Map invariants to fermion generation structure
 *
 * QUALITY GATE:
 *   - Identify topological classification yielding exactly 3 stable families
 *   - OR document why TRD does NOT naturally predict 3 generations
 *
 * STATUS: Framework implementation - physics hypothesis under investigation
 */

#include "TRDCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <yaml-cpp/yaml.h>

const float PI = 3.14159265358979323846f;

/**
 * Topological charge calculator (winding number)
 * Q = (1/2π) ∮ dθ around closed loop
 */
float computeTopologicalCharge(TRDCore3D& core, uint32_t z_slice = 0) {
    const uint32_t Nx = core.getNx();
    const uint32_t Ny = core.getNy();
    const auto& theta = core.getTheta();

    float total_winding = 0.0f;

    // Compute winding around each plaquette in xy-plane
    for (uint32_t j = 0; j < Ny-1; ++j) {
        for (uint32_t i = 0; i < Nx-1; ++i) {
            // Four corners of plaquette
            float th00 = theta[core.index3D(i, j, z_slice)];
            float th10 = theta[core.index3D(i+1, j, z_slice)];
            float th11 = theta[core.index3D(i+1, j+1, z_slice)];
            float th01 = theta[core.index3D(i, j+1, z_slice)];

            // Phase differences (with proper unwrapping)
            auto wrap = [](float x) {
                while (x > PI) x -= 2*PI;
                while (x < -PI) x += 2*PI;
                return x;
            };

            float d1 = wrap(th10 - th00);
            float d2 = wrap(th11 - th10);
            float d3 = wrap(th01 - th11);
            float d4 = wrap(th00 - th01);

            total_winding += (d1 + d2 + d3 + d4) / (2*PI);
        }
    }

    return total_winding;
}

/**
 * Initialize topological defect with winding number Q
 */
void initTopologicalDefect(TRDCore3D& core, int Q, const std::string& type) {
    const uint32_t Nx = core.getNx();
    const uint32_t Ny = core.getNy();
    const uint32_t Nz = core.getNz();
    const float Lx = Nx / 2.0f;
    const float Ly = Ny / 2.0f;
    const float Lz = Nz / 2.0f;

    auto& theta = core.getTheta();

    if (type == "point") {
        // Point defect: θ = Q·atan2(y,x) constant along z
        for (uint32_t k = 0; k < Nz; ++k) {
            for (uint32_t j = 0; j < Ny; ++j) {
                for (uint32_t i = 0; i < Nx; ++i) {
                    float x = static_cast<float>(i) - Lx;
                    float y = static_cast<float>(j) - Ly;
                    theta[core.index3D(i, j, k)] = Q * std::atan2(y, x);
                }
            }
        }
    } else if (type == "line") {
        // Line defect: Flux tube along z-axis with winding Q
        for (uint32_t k = 0; k < Nz; ++k) {
            for (uint32_t j = 0; j < Ny; ++j) {
                for (uint32_t i = 0; i < Nx; ++i) {
                    float x = static_cast<float>(i) - Lx;
                    float y = static_cast<float>(j) - Ly;
                    float z = static_cast<float>(k) - Lz;

                    // Helical phase structure along z
                    float phi_xy = Q * std::atan2(y, x);
                    float phi_z = Q * z * 0.1f;  // Pitch of helix

                    theta[core.index3D(i, j, k)] = phi_xy + phi_z;
                }
            }
        }
    } else if (type == "surface") {
        // Surface defect: Domain wall at z=Lz with phase jump Q·π
        for (uint32_t k = 0; k < Nz; ++k) {
            for (uint32_t j = 0; j < Ny; ++j) {
                for (uint32_t i = 0; i < Nx; ++i) {
                    float z = static_cast<float>(k) - Lz;
                    theta[core.index3D(i, j, k)] = (z > 0) ? Q * PI : 0.0f;
                }
            }
        }
    }

    std::cout << "  Initialized: " << type << " defect with Q=" << Q << "\n";
}

/**
 * Classify topological stability after evolution
 */
struct TopologicalClassification {
    int winding_number;
    float stability_measure;  // R-field strength
    float energy_level;       // Configuration energy
    std::string defect_type;

    bool is_stable() const { return stability_measure > 0.5f; }
};

TopologicalClassification classifyConfiguration(TRDCore3D& core, const std::string& defect_type) {
    TopologicalClassification result;
    result.defect_type = defect_type;

    // Compute winding number
    result.winding_number = static_cast<int>(std::round(computeTopologicalCharge(core)));

    // Compute R-field (synchronization strength)
    core.computeRField();
    const auto& R_field = core.getRField();
    float R_total = 0.0f;
    const uint32_t total_points = core.getNx() * core.getNy() * core.getNz();
    for (uint32_t i = 0; i < total_points; ++i) {
        R_total += R_field[i];
    }
    result.stability_measure = R_total / total_points;

    // Compute energy (simplified: gradient energy)
    const auto& theta = core.getTheta();
    float energy = 0.0f;
    const uint32_t Nx = core.getNx();
    const uint32_t Ny = core.getNy();
    const uint32_t Nz = core.getNz();

    for (uint32_t k = 0; k < Nz; ++k) {
        for (uint32_t j = 0; j < Ny; ++j) {
            for (uint32_t i = 0; i < Nx; ++i) {
                uint32_t idx = core.index3D(i, j, k);
                auto neighbors = core.getNeighbors(i, j, k);

                float grad_sq = 0.0f;
                grad_sq += std::pow(theta[neighbors.x_plus] - theta[idx], 2);
                grad_sq += std::pow(theta[neighbors.y_plus] - theta[idx], 2);
                grad_sq += std::pow(theta[neighbors.z_plus] - theta[idx], 2);

                energy += grad_sq;
            }
        }
    }
    result.energy_level = energy / total_points;

    return result;
}

int runThreeGenerationsTest() {
    std::cout << "\n===== B3: Three-Generation Structure Test =====\n";
    std::cout << "Hypothesis: Topological defects in 3D yield exactly 3 families\n\n";

    // Load configuration
    YAML::Node config;
    try {
        config = YAML::LoadFile("config/three_generations.yaml");
    } catch (const std::exception& e) {
        std::cerr << "Error loading config: " << e.what() << "\n";
        std::cerr << "Using default parameters\n";
    }

    // Create TRDCore3D instance (CPU mode)
    TRDCore3D core;

    // Configure grid
    TRDCore3D::Config trd_config;
    trd_config.Nx = config["grid"]["Nx"].as<uint32_t>(32);
    trd_config.Ny = config["grid"]["Ny"].as<uint32_t>(32);
    trd_config.Nz = config["grid"]["Nz"].as<uint32_t>(32);
    trd_config.dt = config["physics"]["dt"].as<float>(0.01f);
    trd_config.coupling_strength = config["physics"]["coupling"].as<float>(1.0f);

    core.initialize(trd_config);

    // Test different topological configurations
    std::vector<TopologicalClassification> stable_configs;

    std::cout << "1. Testing Point Defects (0D)\n";
    std::cout << "================================\n";
    for (int Q = 1; Q <= 5; ++Q) {
        std::cout << "\nTesting Q=" << Q << ":\n";

        // Initialize defect
        initTopologicalDefect(core, Q, "point");

        // Evolve system
        const int n_steps = 100;
        for (int step = 0; step < n_steps; ++step) {
            core.evolveKuramotoCPU(trd_config.dt);
        }

        // Classify result
        auto classification = classifyConfiguration(core, "point");

        std::cout << "  Final winding: " << classification.winding_number << "\n";
        std::cout << "  Stability: " << classification.stability_measure << "\n";
        std::cout << "  Energy: " << classification.energy_level << "\n";
        std::cout << "  Stable: " << (classification.is_stable() ? "YES" : "NO") << "\n";

        if (classification.is_stable()) {
            stable_configs.push_back(classification);
        }
    }

    std::cout << "\n2. Testing Line Defects (1D)\n";
    std::cout << "================================\n";
    for (int Q = 1; Q <= 5; ++Q) {
        std::cout << "\nTesting Q=" << Q << ":\n";

        initTopologicalDefect(core, Q, "line");

        const int n_steps = 100;
        for (int step = 0; step < n_steps; ++step) {
            core.evolveKuramotoCPU(trd_config.dt);
        }

        auto classification = classifyConfiguration(core, "line");

        std::cout << "  Final winding: " << classification.winding_number << "\n";
        std::cout << "  Stability: " << classification.stability_measure << "\n";
        std::cout << "  Energy: " << classification.energy_level << "\n";
        std::cout << "  Stable: " << (classification.is_stable() ? "YES" : "NO") << "\n";

        if (classification.is_stable()) {
            stable_configs.push_back(classification);
        }
    }

    std::cout << "\n3. Testing Surface Defects (2D)\n";
    std::cout << "================================\n";
    for (int Q = 1; Q <= 5; ++Q) {
        std::cout << "\nTesting Q=" << Q << ":\n";

        initTopologicalDefect(core, Q, "surface");

        const int n_steps = 100;
        for (int step = 0; step < n_steps; ++step) {
            core.evolveKuramotoCPU(trd_config.dt);
        }

        auto classification = classifyConfiguration(core, "surface");

        std::cout << "  Final winding: " << classification.winding_number << "\n";
        std::cout << "  Stability: " << classification.stability_measure << "\n";
        std::cout << "  Energy: " << classification.energy_level << "\n";
        std::cout << "  Stable: " << (classification.is_stable() ? "YES" : "NO") << "\n";

        if (classification.is_stable()) {
            stable_configs.push_back(classification);
        }
    }

    // Analysis: Do we get exactly 3 families?
    std::cout << "\n===== ANALYSIS: Generation Structure =====\n";
    std::cout << "Total stable configurations: " << stable_configs.size() << "\n\n";

    // Group by defect type
    std::map<std::string, std::vector<TopologicalClassification>> by_type;
    for (const auto& config : stable_configs) {
        by_type[config.defect_type].push_back(config);
    }

    for (const auto& [type, configs] : by_type) {
        std::cout << type << " defects: " << configs.size() << " stable states\n";
        for (const auto& c : configs) {
            std::cout << "  Q=" << c.winding_number
                     << " R=" << c.stability_measure
                     << " E=" << c.energy_level << "\n";
        }
    }

    // CRITICAL ASSESSMENT
    std::cout << "\n===== CRITICAL ASSESSMENT =====\n";

    bool explains_three_generations = false;

    // Check if any defect type yields exactly 3 stable configurations
    for (const auto& [type, configs] : by_type) {
        if (configs.size() == 3) {
            std::cout << "✓ " << type << " defects yield EXACTLY 3 stable states!\n";
            std::cout << "  Possible mapping to fermion generations:\n";
            std::cout << "  Q=" << configs[0].winding_number << " → electron family\n";
            std::cout << "  Q=" << configs[1].winding_number << " → muon family\n";
            std::cout << "  Q=" << configs[2].winding_number << " → tau family\n";
            explains_three_generations = true;
        }
    }

    if (!explains_three_generations) {
        std::cout << "✗ TRD topology does NOT naturally yield 3 generations\n";
        std::cout << "  Additional physics required:\n";
        std::cout << "  - Non-Abelian gauge structure?\n";
        std::cout << "  - Higher-dimensional embedding?\n";
        std::cout << "  - Anthropic selection principle?\n";
    }

    std::cout << "\n===== TEST "
              << (explains_three_generations ? "PASSED" : "FRAMEWORK COMPLETE")
              << " =====\n";

    return explains_three_generations ? 0 : 1;
}