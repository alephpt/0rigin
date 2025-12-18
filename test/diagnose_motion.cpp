/**
 * Diagnose why wavepacket moves when it should be stationary
 */

#include <iostream>
#include <iomanip>
#include "../src/DiracEvolution.h"

int main() {
    DiracEvolution dirac(64, 64);
    dirac.initialize(32.0f, 32.0f, 5.0f);
    
    // Check initial momentum distribution
    std::vector<float> kx, ky, density_k;
    dirac.getMomentumDistribution(kx, ky, density_k);
    
    // Find peak in k-space
    float k_peak_density = 0.0f;
    uint32_t peak_idx = 0;
    for (uint32_t i = 0; i < density_k.size(); i++) {
        if (density_k[i] > k_peak_density) {
            k_peak_density = density_k[i];
            peak_idx = i;
        }
    }
    
    uint32_t Nx = 64;
    uint32_t peak_i = peak_idx % Nx;
    uint32_t peak_j = peak_idx / Nx;
    
    std::cout << "Peak in k-space at: (" << peak_i << ", " << peak_j << ")" << std::endl;
    std::cout << "Corresponding k: (" << kx[peak_i] << ", " << ky[peak_j] << ")" << std::endl;
    
    // Track motion without mass field
    float x_init, y_init;
    dirac.getCenterOfMass(x_init, y_init);
    std::cout << "\nInitial CoM: (" << x_init << ", " << y_init << ")" << std::endl;
    
    std::vector<float> mass_field(64 * 64, 0.0f);  // ZERO mass
    
    for (int step = 0; step < 100; step++) {
        dirac.step(mass_field, 0.01f);
    }
    
    float x_final, y_final;
    dirac.getCenterOfMass(x_final, y_final);
    std::cout << "Final CoM (m=0): (" << x_final << ", " << y_final << ")" << std::endl;
    std::cout << "Displacement: (" << (x_final - x_init) << ", " << (y_final - y_init) << ")" << std::endl;
    
    std::cout << "\nExpected: Zero displacement for symmetric Gaussian at k=0" << std::endl;
    
    return 0;
}
