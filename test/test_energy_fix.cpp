/**
 * Quick test to verify energy diagnostic fix
 * Run 5k steps with strong defect and check energy stability
 */

#include "../src/DiracEvolution.h"
#include "../src/SMFTCommon.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <complex>
#include <vector>

const uint32_t NX = 128;
const uint32_t NY = 128;
const uint32_t N_GRID = NX * NY;

const float DELTA = 0.5f;
const float DT = 0.01f;
const float K = 1.0f;
const float DAMPING = 0.1f;

const uint32_t DEFECT_X = 64;
const uint32_t DEFECT_Y = 64;
const float DEFECT_RADIUS = 15.0f;
const float OMEGA_DEFECT = 1.5f;  // Strong mismatch

const float PSI_X0 = 48.0f;
const float PSI_Y0 = 64.0f;
const float PSI_SIGMA = 5.0f;

const int KURAMOTO_WARMUP = 1000;
const int COUPLED_STEPS = 5000;    // Short test
const int OUTPUT_INTERVAL = 100;

using Complex = std::complex<float>;

// Use functions from SMFTCommon instead of duplicating

int main() {
    std::cout << "=== ENERGY DIAGNOSTIC FIX VERIFICATION ===" << std::endl;
    std::cout << "Testing with STRONG defect (ω_defect=1.5)" << std::endl;

    // Initialize
    std::vector<float> theta(N_GRID, 0.0f);
    std::vector<float> omega(N_GRID, 0.0f);

    std::random_device rd;
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);

    for (uint32_t y = 0; y < NY; y++) {
        for (uint32_t x = 0; x < NX; x++) {
            float dx = x - DEFECT_X;
            float dy = y - DEFECT_Y;
            if (std::sqrt(dx*dx + dy*dy) < DEFECT_RADIUS) {
                theta[y * NX + x] = phase_dist(gen);
                omega[y * NX + x] = OMEGA_DEFECT;
            }
        }
    }

    std::cout << "\n[1] Kuramoto warmup..." << std::endl;
    for (int step = 0; step < KURAMOTO_WARMUP; step++) {
        SMFT::stepKuramoto(theta, omega, DT, K, DAMPING, NX, NY);
    }

    auto R_field = SMFT::compute_local_R(theta, NX, NY);
    float defect_R = R_field[DEFECT_Y*NX+DEFECT_X];
    float background_R = R_field[10*NX+10];
    float delta_R = background_R - defect_R;

    std::cout << "  Defect contrast ΔR = " << delta_R << std::endl;
    if (delta_R < 0.5f) {
        std::cout << "  WARNING: Weak defect! Energy test may be inconclusive." << std::endl;
    }

    std::cout << "\n[2] Initialize Dirac..." << std::endl;
    DiracEvolution dirac(NX, NY);
    dirac.initialize(PSI_X0, PSI_Y0, PSI_SIGMA);

    std::ofstream energy_file("output/10/energy_fix_test/energy.dat");
    energy_file << "# step E_total KE PE norm dE_rel\n";

    std::cout << "\n[3] Testing energy stability..." << std::endl;

    float E_initial = 0.0f;
    for (int step = 0; step <= COUPLED_STEPS; step++) {
        if (step > 0) {
            SMFT::stepKuramoto(theta, omega, DT, K, DAMPING, NX, NY);
            R_field = SMFT::compute_local_R(theta, NX, NY);
            std::vector<float> mass_field(N_GRID);
            for (uint32_t i = 0; i < N_GRID; i++) {
                mass_field[i] = DELTA * R_field[i];
            }
            dirac.step(mass_field, DT);
        }

        if (step % OUTPUT_INTERVAL == 0) {
            R_field = SMFT::compute_local_R(theta, NX, NY);
            std::vector<float> mass_field(N_GRID);
            for (uint32_t i = 0; i < N_GRID; i++) {
                mass_field[i] = DELTA * R_field[i];
            }

            float KE, PE;
            float E_total = dirac.getEnergy(mass_field, KE, PE);
            float norm = dirac.getNorm();

            if (step == 0) E_initial = E_total;
            float dE_rel = (E_total - E_initial) / E_initial;

            energy_file << step << " " << E_total << " " << KE << " "
                       << PE << " " << norm << " " << dE_rel << "\n";

            if (step % 1000 == 0) {
                std::cout << "  Step " << std::setw(5) << step
                          << ": E=" << std::setprecision(4) << std::fixed << E_total
                          << ", dE/E=" << std::setprecision(2) << std::scientific << dE_rel
                          << ", norm=" << std::setprecision(6) << std::fixed << norm
                          << std::endl;
            }
        }
    }

    energy_file.close();

    std::cout << "\n=== VERIFICATION COMPLETE ===" << std::endl;
    std::cout << "Check output/10/energy_fix_test/energy.dat" << std::endl;
    std::cout << "Expected: |dE/E| < 0.05 (5% drift acceptable for 5k steps)" << std::endl;
    std::cout << "Old bug: E would explode to 10000+ with strong defect" << std::endl;

    return 0;
}
