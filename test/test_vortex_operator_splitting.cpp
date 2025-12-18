// Vortex operator splitting test - non-uniform R field with phase defect
// Tests N=1 vs N=100 with spatial phase structure to verify coupling

#include "SMFTEngine.h"
#include "DiracEvolution.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void runVortexTest(int N, int total_steps, const std::string& output_dir) {
    std::cout << "\n=== Running Vortex Test N=" << N << " ===" << std::endl;
    std::cout << "Total steps: " << total_steps << std::endl;
    std::cout << "Output: " << output_dir << std::endl;

    // Initialize Nova (minimal config for compute)
    NovaConfig config = {
        .name = "Vortex Operator Splitting Test",
        .screen = {800, 600},
        .debug_level = "info",
        .dimensions = "2D",
        .camera_type = "orthographic",
        .compute = true,
    };

    Nova nova(config);
    nova.initialized = true;

    // Create engine
    SMFTEngine engine(&nova);

    // Initialize 64x64 grid
    uint32_t Nx = 64;
    float Delta = 2.5f;
    float chiral_angle = 0.0f;

    std::cout << "Initializing " << Nx << "x" << Nx << " grid..." << std::endl;
    engine.initialize(Nx, Nx, Delta, chiral_angle);

    // Create VORTEX phase field with winding number 1 at center
    std::vector<float> theta_vortex(Nx * Nx);
    float cx = Nx / 2.0f;
    float cy = Nx / 2.0f;

    for (uint32_t j = 0; j < Nx; j++) {
        for (uint32_t i = 0; i < Nx; i++) {
            uint32_t idx = j * Nx + i;
            float x = static_cast<float>(i);
            float y = static_cast<float>(j);

            // Vortex: θ(x,y) = atan2(y - cy, x - cx)
            // Creates winding number +1 at center
            theta_vortex[idx] = std::atan2(y - cy, x - cx);
        }
    }

    std::cout << "Setting vortex phase field (winding number = 1 at center)..." << std::endl;
    engine.setInitialPhases(theta_vortex);

    // Initialize hybrid system with Dirac field (Gaussian wavepacket)
    float x0 = 32.0f;
    float y0 = 32.0f;
    float sigma = 3.0f;

    std::cout << "Initializing Dirac field..." << std::endl;
    engine.initializeHybrid(x0, y0, sigma);

    // Parameters
    float dt = 0.01f;
    float K = 2.0f;
    float damping = 0.05f;
    float coupling = 0.5f;

    // Create output directory
    system(("mkdir -p " + output_dir).c_str());

    // Open output files
    std::ofstream out(output_dir + "/timeseries.csv");
    out << "step,time,R_avg,R_std,norm,energy,max_density\n";

    std::ofstream R_spatial(output_dir + "/R_field_snapshots.csv");
    R_spatial << "step,time,";
    for (uint32_t j = 0; j < Nx; j++) {
        for (uint32_t i = 0; i < Nx; i++) {
            R_spatial << "R_" << i << "_" << j;
            if (j < Nx - 1 || i < Nx - 1) R_spatial << ",";
        }
    }
    R_spatial << "\n";

    // Run simulation
    for (int step = 0; step < total_steps; ++step) {
        // Evolve with operator splitting
        engine.stepWithDirac(dt, coupling, N, K, damping);

        // Record observables every 10 steps
        if (step % 10 == 0) {
            auto R_field = engine.getSyncField();

            // Compute R statistics
            float R_avg = 0.0f;
            float R_min = 1e10f;
            float R_max = -1e10f;
            for (auto r : R_field) {
                R_avg += r;
                if (r < R_min) R_min = r;
                if (r > R_max) R_max = r;
            }
            R_avg /= R_field.size();

            float R_std = 0.0f;
            for (auto r : R_field) {
                float diff = r - R_avg;
                R_std += diff * diff;
            }
            R_std = std::sqrt(R_std / R_field.size());

            // Access DiracEvolution to get norm and energy
            const DiracEvolution* dirac = engine.getDiracEvolution();
            float norm = dirac ? dirac->getNorm() : 0.0f;

            // Create mass field for energy calculation
            std::vector<float> mass_field(Nx * Nx);
            for (uint32_t i = 0; i < Nx * Nx; i++) {
                mass_field[i] = Delta * R_field[i];
            }
            float KE = 0.0f, PE = 0.0f;
            float energy = dirac ? dirac->getEnergy(mass_field, KE, PE) : 0.0f;

            auto density = engine.getDiracDensity();
            float max_density = 0.0f;
            for (auto d : density) {
                if (d > max_density) max_density = d;
            }

            out << step << "," << (step * dt) << "," << R_avg << ","
                << R_std << "," << norm << "," << energy << "," << max_density << "\n";

            // Save R field snapshots every 1000 steps
            if (step % 1000 == 0) {
                R_spatial << step << "," << (step * dt);
                for (auto r : R_field) {
                    R_spatial << "," << r;
                }
                R_spatial << "\n";
            }

            if (step % 1000 == 0) {
                std::cout << "  Step " << step << "/" << total_steps
                          << " | R_avg=" << R_avg << " ± " << R_std
                          << " | norm=" << norm << std::endl;
            }
        }
    }

    out.close();
    R_spatial.close();
    std::cout << "✓ N=" << N << " complete!" << std::endl;
    std::cout << "  Timeseries: " << output_dir << "/timeseries.csv" << std::endl;
    std::cout << "  R field snapshots: " << output_dir << "/R_field_snapshots.csv" << std::endl;
}

int main() {
    std::cout << "==================================================" << std::endl;
    std::cout << "Vortex Operator Splitting Validation" << std::endl;
    std::cout << "Testing N=1 vs N=100 with phase defect" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Run N=1 and N=100 for comparison
    runVortexTest(1, 10000, "output/vortex_N1");
    runVortexTest(100, 10000, "output/vortex_N100");

    std::cout << "\n=== ALL TESTS COMPLETE ===" << std::endl;
    std::cout << "Results saved to output/vortex_N{1,100}/" << std::endl;
    std::cout << "\nExpected behavior:" << std::endl;
    std::cout << "- R field should have spatial structure (vortex core low R)" << std::endl;
    std::cout << "- N=100 should smooth out R fluctuations vs N=1" << std::endl;
    std::cout << "- Particle trajectory should differ if coupling is N-dependent" << std::endl;

    return 0;
}
