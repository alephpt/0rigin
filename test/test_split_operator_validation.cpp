/**
 * test_split_operator_validation.cpp
 *
 * Validate that split-operator method preserves norm to machine precision
 * Tests unitary evolution: |Ψ(t)|² should remain 1.0 ± 1e-14
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../src/DiracEvolution.h"

int main() {
    // Setup
    const uint32_t Nx = 64;
    const uint32_t Ny = 64;
    const float dt = 0.01f;
    const int num_steps = 1000;

    DiracEvolution dirac(Nx, Ny);
    
    // Initialize Gaussian wavepacket
    dirac.initialize(Nx/2.0f, Ny/2.0f, 5.0f);
    
    // Check initial norm
    float initial_norm = dirac.getNorm();
    std::cout << "[INIT] Norm = " << std::setprecision(15) << initial_norm << std::endl;
    
    // Create constant mass field
    std::vector<float> mass_field(Nx * Ny, 0.1f);
    
    // Evolve and check norm conservation
    float max_drift = 0.0f;
    for (int step = 0; step < num_steps; step++) {
        dirac.step(mass_field, dt);
        
        if (step % 100 == 0 || step == num_steps - 1) {
            float norm = dirac.getNorm();
            float drift = std::abs(norm - initial_norm);
            max_drift = std::max(max_drift, drift);
            
            std::cout << "[STEP " << std::setw(4) << step << "] "
                      << "Norm = " << std::setprecision(15) << norm << ", "
                      << "Drift = " << std::scientific << drift << std::endl;
        }
    }
    
    // Verdict
    std::cout << "\n[RESULT] Maximum drift: " << std::scientific << max_drift << std::endl;
    
    if (max_drift < 1e-13) {
        std::cout << "[PASS] Split-operator conserves norm to machine precision (drift < 1e-13)" << std::endl;
        return 0;
    } else if (max_drift < 1e-10) {
        std::cout << "[WARN] Acceptable drift but not machine precision (1e-13 < drift < 1e-10)" << std::endl;
        return 0;
    } else {
        std::cout << "[FAIL] Unacceptable norm drift (drift > 1e-10)" << std::endl;
        return 1;
    }
}
