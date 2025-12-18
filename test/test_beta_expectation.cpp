/**
 * test_beta_expectation.cpp
 *
 * Calculate <Ψ|β|Ψ> to determine expected force direction
 * β = diag(1, 1, -1, -1) for Dirac spinor
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
    
    // Calculate <β> at t=0
    auto comp0 = dirac.getComponent(0);
    auto comp1 = dirac.getComponent(1);
    auto comp2 = dirac.getComponent(2);
    auto comp3 = dirac.getComponent(3);
    
    float upper_norm = 0.0f;  // |ψ₁|² + |ψ₂|²
    float lower_norm = 0.0f;  // |ψ₃|² + |ψ₄|²
    
    for (uint32_t i = 0; i < Nx * Ny; i++) {
        upper_norm += std::norm(comp0[i]) + std::norm(comp1[i]);
        lower_norm += std::norm(comp2[i]) + std::norm(comp3[i]);
    }
    
    float beta_expectation = upper_norm - lower_norm;  // <Ψ|β|Ψ> = Σ(upper) - Σ(lower)
    float total_norm = upper_norm + lower_norm;
    
    std::cout << "=== Initial Spinor Analysis (t=0) ===" << std::endl;
    std::cout << "Upper components |ψ₁|² + |ψ₂|² = " << upper_norm << std::endl;
    std::cout << "Lower components |ψ₃|² + |ψ₄|² = " << lower_norm << std::endl;
    std::cout << "Total norm |Ψ|² = " << total_norm << std::endl;
    std::cout << "<Ψ|β|Ψ> = " << beta_expectation << std::endl;
    std::cout << "<β> (normalized) = " << beta_expectation / total_norm << std::endl;
    
    // Evolve with spatially varying mass
    std::vector<float> mass_field(Nx * Ny);
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            mass_field[idx] = 0.5f * (1.0f + 0.5f * std::sin(x / 10.0f));
        }
    }
    
    // Check mass gradient at initial position
    uint32_t x0 = Nx/2;
    uint32_t y0 = Ny/2;
    uint32_t idx0 = y0 * Nx + x0;
    float m_center = mass_field[idx0];
    float m_right = mass_field[y0 * Nx + (x0 + 1)];
    float dm_dx = m_right - m_center;
    
    std::cout << "\n=== Mass Field at Initial Position ===" << std::endl;
    std::cout << "Position: (" << x0 << ", " << y0 << ")" << std::endl;
    std::cout << "m(x=32) = " << m_center << std::endl;
    std::cout << "m(x=33) = " << m_right << std::endl;
    std::cout << "dm/dx ≈ " << dm_dx << std::endl;
    
    std::cout << "\n=== Expected Force Direction ===" << std::endl;
    std::cout << "F ∝ -∇m(x) · <β>" << std::endl;
    std::cout << "F_x ∝ -(dm/dx) · <β> = " << -dm_dx * beta_expectation / total_norm << std::endl;
    
    if (beta_expectation > 0) {
        if (dm_dx < 0) {
            std::cout << "Expected motion: +x (toward HIGH mass, dm/dx < 0)" << std::endl;
        } else {
            std::cout << "Expected motion: -x (toward HIGH mass, dm/dx > 0)" << std::endl;
        }
    } else {
        if (dm_dx < 0) {
            std::cout << "Expected motion: -x (toward LOW mass, dm/dx < 0)" << std::endl;
        } else {
            std::cout << "Expected motion: +x (toward LOW mass, dm/dx > 0)" << std::endl;
        }
    }
    
    // Now evolve and check
    for (int step = 0; step < 1000; step++) {
        dirac.step(mass_field, dt);
    }
    
    // Calculate <β> at t=final
    comp0 = dirac.getComponent(0);
    comp1 = dirac.getComponent(1);
    comp2 = dirac.getComponent(2);
    comp3 = dirac.getComponent(3);
    
    float upper_norm_final = 0.0f;
    float lower_norm_final = 0.0f;
    
    for (uint32_t i = 0; i < Nx * Ny; i++) {
        upper_norm_final += std::norm(comp0[i]) + std::norm(comp1[i]);
        lower_norm_final += std::norm(comp2[i]) + std::norm(comp3[i]);
    }
    
    float beta_final = upper_norm_final - lower_norm_final;
    float total_final = upper_norm_final + lower_norm_final;
    
    std::cout << "\n=== Final Spinor Analysis (t=10) ===" << std::endl;
    std::cout << "Upper components |ψ₁|² + |ψ₂|² = " << upper_norm_final << std::endl;
    std::cout << "Lower components |ψ₃|² + |ψ₄|² = " << lower_norm_final << std::endl;
    std::cout << "Total norm |Ψ|² = " << total_final << std::endl;
    std::cout << "<Ψ|β|Ψ> = " << beta_final << std::endl;
    std::cout << "<β> (normalized) = " << beta_final / total_final << std::endl;
    
    // Center of mass
    auto density = dirac.getDensity();
    float x_mean = 0.0f, y_mean = 0.0f, total_density = 0.0f;
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            float d = density[idx];
            x_mean += x * d;
            y_mean += y * d;
            total_density += d;
        }
    }
    x_mean /= total_density;
    y_mean /= total_density;
    
    std::cout << "\n=== Observed Motion ===" << std::endl;
    std::cout << "Initial CoM: (32, 32)" << std::endl;
    std::cout << "Final CoM: (" << x_mean << ", " << y_mean << ")" << std::endl;
    std::cout << "Displacement: Δx = " << (x_mean - 32.0f) << std::endl;
    
    float m_final_pos = mass_field[int(y_mean) * Nx + int(x_mean)];
    std::cout << "m(x=32) = " << m_center << std::endl;
    std::cout << "m(x=" << int(x_mean) << ") = " << m_final_pos << std::endl;
    
    if (m_final_pos > m_center) {
        std::cout << "Moved toward HIGHER mass ✓" << std::endl;
    } else {
        std::cout << "Moved toward LOWER mass ✓" << std::endl;
    }
    
    return 0;
}
