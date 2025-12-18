// Complete 10k timestep operator splitting test - N=1, N=10, N=100
#include "SMFTCore.h"
#include "SMFTEngine.h"
#include "DiracEvolution.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void runTest(int N, int total_steps, const std::string& output_dir) {
    std::cout << "\n=== Running N=" << N << " test ===" << std::endl;
    std::cout << "Total steps: " << total_steps << std::endl;
    std::cout << "Output: " << output_dir << std::endl;

    // Initialize Nova (minimal config for compute)
    NovaConfig config = {
        .name = "Operator Splitting 10k Test",
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

    // Initialize hybrid system with Dirac field (Gaussian wavepacket)
    float x0 = 32.0f;
    float y0 = 32.0f;
    float sigma = 3.0f;

    std::cout << "Initializing hybrid GPU-CPU system..." << std::endl;
    engine.initializeHybrid(x0, y0, sigma);

    // Parameters
    float dt = 0.01f;
    float K = 2.0f;
    float damping = 0.05f;
    float coupling = 0.5f;

    // Create output directory
    system(("mkdir -p " + output_dir).c_str());

    // Open output file
    std::ofstream out(output_dir + "/timeseries.csv");
    out << "step,time,R_avg,norm,energy,max_density\n";

    // Run simulation
    for (int step = 0; step < total_steps; ++step) {
        // Evolve with operator splitting
        engine.stepWithDirac(dt, coupling, N, K, damping);

        // Record observables every 10 steps
        if (step % 10 == 0) {
            auto R_field = engine.getSyncField();
            float R_avg = 0.0f;
            for (auto r : R_field) R_avg += r;
            R_avg /= R_field.size();

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
                << norm << "," << energy << "," << max_density << "\n";

            if (step % 1000 == 0) {
                std::cout << "  Step " << step << "/" << total_steps
                          << " | R_avg=" << R_avg << " | norm=" << norm << std::endl;
            }
        }
    }

    out.close();
    std::cout << "✓ N=" << N << " complete!" << std::endl;
}

int main() {
    std::cout << "==================================================" << std::endl;
    std::cout << "Operator Splitting 10k Timestep Validation" << std::endl;
    std::cout << "Testing N=1, N=10, N=100" << std::endl;
    std::cout << "==================================================" << std::endl;

    runTest(1, 10000, "output/N1_10k");
    runTest(10, 10000, "output/N10_10k");
    runTest(100, 10000, "output/N100_10k");

    std::cout << "\n=== ALL TESTS COMPLETE ===" << std::endl;
    std::cout << "Results saved to output/N{1,10,100}_10k/" << std::endl;

    return 0;
}
