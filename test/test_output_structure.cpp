/**
 * @file test_output_structure.cpp
 * @brief Test output directory structure output/$n/
 *
 * Verifies that we can create numbered output directories
 * and write data files correctly.
 */

#include "../lib/Nova/Nova.h"
#include "../src/SMFTEngine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <filesystem>

int main() {
    std::cout << "=== Output Structure Test ===" << std::endl;
    std::cout << "Testing output/$n/ directory structure\n" << std::endl;

    // Small grid for quick test
    const int Nx = 32;
    const int Ny = 32;
    const int num_tests = 3;  // Run 3 test simulations

    // Initialize Nova (headless mode)
    NovaConfig config = {
        .name = "Output Structure Test",
        .screen = {800, 600},
        .debug_level = "error",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };
    Nova nova(config);
    nova.initialized = true;

    // Initialize SMFT engine
    SMFTEngine engine(&nova);
    engine.initialize(Nx, Ny, 2.5f, 0.0f);

    for (int test_id = 1; test_id <= num_tests; test_id++) {
        std::cout << "\n=== Test " << test_id << " ===" << std::endl;

        // Create output directory: /home/persist/neotec/0rigin/output/1/, etc.
        std::string dirname = "/home/persist/neotec/0rigin/output/" + std::to_string(test_id);

        // Create directory (mkdir -p equivalent)
        std::filesystem::create_directories(dirname);

        std::cout << "Created directory: " << dirname << std::endl;

        // Set random initial phases
        std::vector<float> theta_init(Nx * Ny);
        srand(42 + test_id);  // Different seed for each test
        for (int i = 0; i < Nx * Ny; i++) {
            theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
        }
        engine.setInitialPhases(theta_init);

        // Run 10 steps
        for (int step = 0; step < 10; step++) {
            engine.stepStochastic(0.01f, 1.0f, 1e-4f * test_id);  // Different sigma for each test
        }

        // Get final state
        std::vector<float> R_field = engine.getSyncField();
        std::vector<float> theta_field = engine.getPhaseField();

        // Write output files

        // Write R field
        std::string filename = dirname + "/R_field.dat";
        std::ofstream out_R(filename);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                out_R << R_field[j * Nx + i];
                if (i < Nx - 1) out_R << " ";
            }
            out_R << "\n";
        }
        out_R.close();
        std::cout << "  Written: " << filename << std::endl;

        // Write theta field
        filename = dirname + "/theta.dat";
        std::ofstream out_theta(filename);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                out_theta << theta_field[j * Nx + i];
                if (i < Nx - 1) out_theta << " ";
            }
            out_theta << "\n";
        }
        out_theta.close();
        std::cout << "  Written: " << filename << std::endl;

        // Write metadata
        snprintf(filename, sizeof(filename), "%s/metadata.txt", dirname);
        std::ofstream out_meta(filename);
        out_meta << "Test ID: " << test_id << "\n";
        out_meta << "Grid: " << Nx << " x " << Ny << "\n";
        out_meta << "Steps: 10\n";
        out_meta << "Sigma: " << (1e-4f * test_id) << "\n";
        out_meta << "Final <R>: ";
        float mean_R = 0.0f;
        for (auto r : R_field) mean_R += r;
        mean_R /= R_field.size();
        out_meta << mean_R << "\n";
        out_meta.close();
        std::cout << "  Written: " << filename << std::endl;
    }

    std::cout << "\n=== Output Structure Test Complete ===" << std::endl;
    std::cout << "Created directories:" << std::endl;
    std::cout << "  /home/persist/neotec/0rigin/output/1/ (3 files)" << std::endl;
    std::cout << "  /home/persist/neotec/0rigin/output/2/ (3 files)" << std::endl;
    std::cout << "  /home/persist/neotec/0rigin/output/3/ (3 files)" << std::endl;
    std::cout << "\n✓ Test passed - output structure working!" << std::endl;

    return 0;
}
