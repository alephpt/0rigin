/**
 * Test if norm drift is from FFT normalization or evolution
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../src/DiracEvolution.h"

int main() {
    const uint32_t Nx = 64;
    const uint32_t Ny = 64;
    const float dt = 0.01f;

    DiracEvolution dirac(Nx, Ny);
    dirac.initialize(Nx/2.0f, Ny/2.0f, 5.0f);
    
    std::vector<float> mass_field(Nx * Ny, 0.0f);  // Zero mass
    
    float initial_norm = dirac.getNorm();
    std::cout << "[INIT] Norm = " << std::setprecision(15) << initial_norm << std::endl;
    
    // Test just kinetic evolution (mass = 0)
    for (int step = 0; step < 100; step++) {
        dirac.step(mass_field, dt);
    }
    
    float norm_after_100 = dirac.getNorm();
    float drift = std::abs(norm_after_100 - initial_norm);
    
    std::cout << "[100 steps, m=0] Norm = " << norm_after_100 
              << ", Drift = " << std::scientific << drift << std::endl;
    
    if (drift < 1e-13) {
        std::cout << "[PASS] Kinetic-only evolution preserves norm" << std::endl;
        return 0;
    } else {
        std::cout << "[FAIL] Kinetic evolution has drift (FFT normalization issue?)" << std::endl;
        return 1;
    }
}
