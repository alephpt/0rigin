// test/test_proca_vortex_magnetic_field.cpp
#include "physics/ProcaEM.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

int main() {
    const int nx = 64;
    const int ny = 64;
    const float dx = 1.0f;
    const float photon_mass_coupling = 0.1f;
    const int num_steps = 1000;
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

    // Create Proca EM
    physics::ProcaEM proca(nx, ny, dx, photon_mass_coupling);

    std::cout << "=== PROCA VORTEX→B CRITICAL TEST ===" << std::endl;
    std::cout << "Grid: " << nx << "×" << ny << std::endl;
    std::cout << "Vortex center: (" << vortex_x << ", " << vortex_y << ")" << std::endl;
    std::cout << "Photon mass coupling: " << photon_mass_coupling << std::endl;
    std::cout << "Evolving for " << num_steps << " steps..." << std::endl;

    // Evolve
    for (int step = 0; step < num_steps; ++step) {
        proca.computePotentials(theta.data(), R.data(), nx, ny, dx, dt);

        if (step % 100 == 0) {
            proca.computeFieldStrengths();
            auto F = proca.getFieldAt(nx/2, ny/2);
            std::cout << "Step " << step << ": B_z(center) = " << F.Bz << std::endl;
        }
    }

    // Final measurement
    proca.computeFieldStrengths();

    // Measure B-field at vortex core and surrounding region
    std::cout << "\n=== FINAL B-FIELD MEASUREMENT ===" << std::endl;

    auto F_center = proca.getFieldAt(nx/2, ny/2);
    std::cout << "B_z at vortex core: " << F_center.Bz << std::endl;

    // Measure B in ring around vortex
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

    B_rms = std::sqrt(B_rms / count);

    std::cout << "B_max in ring (r=" << r_min << "-" << r_max << "): " << B_max << std::endl;
    std::cout << "B_rms in ring: " << B_rms << std::endl;

    // Save B-field map
    std::ofstream out("output/proca_vortex_test/B_field.csv");
    out << "i,j,Bz\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            auto F = proca.getFieldAt(i, j);
            out << i << "," << j << "," << F.Bz << "\n";
        }
    }
    out.close();

    std::cout << "B-field map saved to output/proca_vortex_test/B_field.csv" << std::endl;

    // GO/NO-GO DECISION
    std::cout << "\n=== GO/NO-GO VERDICT ===" << std::endl;

    const float threshold = 0.01f;
    bool pass = (B_max > threshold) || (B_rms > threshold * 0.1f);

    if (pass) {
        std::cout << "✅ PASS: B ≠ 0 detected!" << std::endl;
        std::cout << "   Proca mechanism WORKS - EM emerges from SMFT" << std::endl;
        std::cout << "   Proceeding to Week 3 (Stückelberg implementation)" << std::endl;
        return 0;
    } else {
        std::cout << "❌ FAIL: B = 0 (below threshold " << threshold << ")" << std::endl;
        std::cout << "   Proca mechanism DOES NOT generate magnetic fields" << std::endl;
        std::cout << "   Options:" << std::endl;
        std::cout << "     1. Debug Proca implementation" << std::endl;
        std::cout << "     2. Try Stückelberg mechanism" << std::endl;
        std::cout << "     3. Pivot to Option A (drop EM emergence)" << std::endl;
        return 1;
    }
}
