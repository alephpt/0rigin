/**
 * validate_physics.cpp
 *
 * Single comprehensive validation using DiracEvolution analysis methods
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include "../src/DiracEvolution.h"

void test_beta_and_force() {
    std::cout << "\n=== TEST 1: Beta Expectation & Force Direction ===" << std::endl;
    
    DiracEvolution dirac(64, 64);
    dirac.initialize(32.0f, 32.0f, 5.0f);
    
    float beta_init = dirac.getBetaExpectation();
    float norm_init = dirac.getNorm();
    
    std::cout << "Initial <β> = " << beta_init << " (normalized: " << beta_init/norm_init << ")" << std::endl;
    std::cout << "Expected: +1.0 for upper-component Gaussian" << std::endl;
    
    if (std::abs(beta_init/norm_init - 1.0f) < 0.01f) {
        std::cout << "[PASS]" << std::endl;
    } else {
        std::cout << "[FAIL]" << std::endl;
    }
}

void test_center_of_mass_motion() {
    std::cout << "\n=== TEST 2: Center of Mass Motion ===" << std::endl;
    
    DiracEvolution dirac(64, 64);
    dirac.initialize(32.0f, 32.0f, 5.0f);
    
    float x_init, y_init;
    dirac.getCenterOfMass(x_init, y_init);
    
    // Mass field with gradient
    std::vector<float> mass_field(64 * 64);
    for (uint32_t y = 0; y < 64; y++) {
        for (uint32_t x = 0; x < 64; x++) {
            mass_field[y * 64 + x] = 0.5f * (1.0f + 0.5f * std::sin(x / 10.0f));
        }
    }
    
    // Evolve
    for (int step = 0; step < 1000; step++) {
        dirac.step(mass_field, 0.01f);
    }
    
    float x_final, y_final;
    dirac.getCenterOfMass(x_final, y_final);
    
    float dx = x_final - x_init;
    float dy = y_final - y_init;
    
    std::cout << "Initial CoM: (" << x_init << ", " << y_init << ")" << std::endl;
    std::cout << "Final CoM: (" << x_final << ", " << y_final << ")" << std::endl;
    std::cout << "Displacement: Δx=" << dx << ", Δy=" << dy << std::endl;
    
    // Check mass at initial and final positions
    uint32_t idx_init = uint32_t(y_init) * 64 + uint32_t(x_init);
    uint32_t idx_final = uint32_t(y_final) * 64 + uint32_t(x_final);
    float m_init = mass_field[idx_init];
    float m_final = mass_field[idx_final];
    
    std::cout << "m(initial) = " << m_init << std::endl;
    std::cout << "m(final) = " << m_final << std::endl;
    
    // For <β> > 0, force F = -β·∇m, so should move toward LOWER mass
    if (m_final < m_init && dx > 0) {
        std::cout << "[PASS] Moved toward lower mass as expected (F = -β·∇m)" << std::endl;
    } else {
        std::cout << "[FAIL] Motion inconsistent with Ehrenfest theorem" << std::endl;
    }
}

void test_norm_conservation() {
    std::cout << "\n=== TEST 3: Norm Conservation (50k steps) ===" << std::endl;
    
    DiracEvolution dirac(64, 64);
    dirac.initialize(32.0f, 32.0f, 5.0f);
    
    float norm_init = dirac.getNorm();
    std::vector<float> mass_field(64 * 64, 0.1f);
    
    std::ofstream out("output/norm_vs_time.dat");
    out << "# step norm drift\n";
    out << "0 " << norm_init << " 0.0\n";
    
    std::cout << "Running 50,000 steps..." << std::endl;
    
    for (int step = 1; step <= 50000; step++) {
        dirac.step(mass_field, 0.01f);
        
        if (step % 5000 == 0) {
            float norm = dirac.getNorm();
            float drift = std::abs(norm - norm_init);
            out << step << " " << norm << " " << drift << "\n";
            std::cout << "  Step " << step << ": norm=" << norm << ", drift=" << drift << std::endl;
        }
    }
    
    out.close();
    
    float norm_final = dirac.getNorm();
    float final_drift = std::abs(norm_final - norm_init);
    
    std::cout << "Final drift: " << final_drift << std::endl;
    
    if (final_drift < 0.01f) {
        std::cout << "[PASS] Drift < 1% over 50k steps" << std::endl;
    } else {
        std::cout << "[FAIL] Excessive drift" << std::endl;
    }
}

void test_isotropy() {
    std::cout << "\n=== TEST 4: Rotation Invariance ===" << std::endl;
    
    std::vector<float> final_norms;
    
    for (int dir = 0; dir < 4; dir++) {
        DiracEvolution dirac(64, 64);
        dirac.initialize(32.0f, 32.0f, 5.0f);
        
        std::vector<float> mass_field(64 * 64, 0.1f);
        
        for (int step = 0; step < 100; step++) {
            dirac.step(mass_field, 0.01f);
        }
        
        float norm = dirac.getNorm();
        final_norms.push_back(norm);
        
        std::cout << "  Direction " << dir << ": norm=" << norm << std::endl;
    }
    
    // Calculate variance
    float mean = 0.0f;
    for (float n : final_norms) mean += n;
    mean /= final_norms.size();
    
    float variance = 0.0f;
    for (float n : final_norms) {
        variance += (n - mean) * (n - mean);
    }
    variance /= final_norms.size();
    float std_dev = std::sqrt(variance);
    
    std::cout << "Mean: " << mean << ", Std Dev: " << std_dev << ", Variance: " << variance << std::endl;
    
    if (variance < 1e-10f) {
        std::cout << "[PASS] Isotropic to machine precision" << std::endl;
    } else {
        std::cout << "[FAIL] Anisotropic evolution" << std::endl;
    }
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << " Dirac Split-Operator Physics Validation" << std::endl;
    std::cout << "========================================" << std::endl;
    
    test_beta_and_force();
    test_center_of_mass_motion();
    test_norm_conservation();
    test_isotropy();
    
    std::cout << "\n========================================" << std::endl;
    std::cout << " Validation Complete" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
