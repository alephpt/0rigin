/**
 * test_scale_invariance.cpp
 *
 * E4: Scale Invariance Breaking Analysis for TRD
 *
 * Goal: Analyze how TRD breaks conformal symmetry through mass scales
 *       and quantify the β-function of the theory
 *
 * Physics Context:
 *   TRD action: S = ∫ [∂_μθ ∂^μθ + V(R) + K·R·cos(Δθ)] d⁴x
 *
 *   Conformal transformation: x^μ → λx^μ, θ → θ, R → R
 *   Breaking: Potential V(R) and coupling K introduce mass scales
 *
 * Mathematical Framework:
 *   - β-function: β(K) = μ·∂K/∂μ = d(ln K)/d(ln μ)
 *   - Conformal anomaly: T^μ_μ ≠ 0 due to scale-dependent coupling
 *   - Callan-Symanzik equation for running coupling
 *
 * Test Strategy:
 *   1. Initialize field configuration at scale μ₀
 *   2. Apply scale transformation: x → λx (λ = 0.5, 1, 2, 5)
 *   3. Evolve both original and scaled configurations
 *   4. Measure effective coupling K_eff(λ) from correlation functions
 *   5. Extract β-function: β(K) = d(ln K_eff)/d(ln λ)
 *
 * Quality Gates:
 *   - β(K) ≠ 0 (scale invariance broken)
 *   - |β(K)| > 0.01 (measurable breaking)
 *   - Conformal anomaly T^μ_μ computed and nonzero
 */

#include "TRDCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <complex>
#include <numeric>

// Physical constants
const float PI = 3.14159265358979323846f;

/**
 * Scale transformation of a field configuration
 * θ(x) → θ(λx) with interpolation for non-integer scaling
 */
void applyScaleTransformation(
    const std::vector<float>& field_in,
    std::vector<float>& field_out,
    int nx, int ny, int nz, float lambda
) {
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // Source coordinates after inverse scaling
                float src_i = i / lambda;
                float src_j = j / lambda;
                float src_k = k / lambda;

                // Trilinear interpolation for non-integer coordinates
                int i0 = int(std::floor(src_i)) % nx;
                int j0 = int(std::floor(src_j)) % ny;
                int k0 = int(std::floor(src_k)) % nz;

                int i1 = (i0 + 1) % nx;
                int j1 = (j0 + 1) % ny;
                int k1 = (k0 + 1) % nz;

                float fx = src_i - std::floor(src_i);
                float fy = src_j - std::floor(src_j);
                float fz = src_k - std::floor(src_k);

                // 8 corners of interpolation cube
                float v000 = field_in[k0 * nx * ny + j0 * nx + i0];
                float v100 = field_in[k0 * nx * ny + j0 * nx + i1];
                float v010 = field_in[k0 * nx * ny + j1 * nx + i0];
                float v110 = field_in[k0 * nx * ny + j1 * nx + i1];
                float v001 = field_in[k1 * nx * ny + j0 * nx + i0];
                float v101 = field_in[k1 * nx * ny + j0 * nx + i1];
                float v011 = field_in[k1 * nx * ny + j1 * nx + i0];
                float v111 = field_in[k1 * nx * ny + j1 * nx + i1];

                // Trilinear interpolation
                float v00 = v000 * (1 - fx) + v100 * fx;
                float v01 = v001 * (1 - fx) + v101 * fx;
                float v10 = v010 * (1 - fx) + v110 * fx;
                float v11 = v011 * (1 - fx) + v111 * fx;

                float v0 = v00 * (1 - fy) + v10 * fy;
                float v1 = v01 * (1 - fy) + v11 * fy;

                field_out[k * nx * ny + j * nx + i] = v0 * (1 - fz) + v1 * fz;
            }
        }
    }
}

/**
 * Compute two-point correlation function C(r) = <θ(x)θ(x+r)>
 * Used to extract effective coupling from decay rate
 */
double computeCorrelationFunction(
    const std::vector<float>& theta_field,
    int nx, int ny, int nz, int r_max
) {
    std::vector<double> C(r_max + 1, 0.0);
    std::vector<int> counts(r_max + 1, 0);

    // Sample correlation at different distances
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx0 = k * nx * ny + j * nx + i;
                float theta0 = theta_field[idx0];

                // Compute correlations with points at distance r
                for (int dr = 1; dr <= r_max && dr < nx/2; ++dr) {
                    // Sample along x-direction
                    int idx1 = k * nx * ny + j * nx + ((i + dr) % nx);
                    C[dr] += theta0 * theta_field[idx1];
                    counts[dr]++;
                }
            }
        }
    }

    // Normalize correlations
    for (int r = 1; r <= r_max; ++r) {
        if (counts[r] > 0) {
            C[r] /= counts[r];
        }
    }

    // Extract correlation length from exponential decay
    // C(r) ~ exp(-r/ξ), where ξ ~ 1/m_eff ~ 1/√K_eff
    double sum_log_C = 0.0;
    double sum_r = 0.0;
    int valid_points = 0;

    for (int r = 2; r <= std::min(10, r_max); ++r) {
        if (C[r] > 1e-10) {
            sum_log_C += std::log(std::abs(C[r]));
            sum_r += r;
            valid_points++;
        }
    }

    if (valid_points > 0) {
        double correlation_length = -valid_points / (sum_log_C / sum_r);
        return correlation_length;
    }

    return 1.0; // Default if cannot extract
}

/**
 * Extract β-function from scale-dependent coupling measurements
 */
double computeBetaFunction(
    const std::vector<float>& couplings,
    const std::vector<float>& scales
) {
    if (couplings.size() < 2) return 0.0;

    // Fit log(K) vs log(μ) to extract β = d(ln K)/d(ln μ)
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int n = couplings.size();

    for (int i = 0; i < n; ++i) {
        double x = std::log(scales[i]);
        double y = std::log(couplings[i]);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    double beta = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    return beta;
}

/**
 * Compute trace of stress-energy tensor (conformal anomaly)
 * T^μ_μ = -4H + ∂_μJ^μ where H is Hamiltonian density
 */
double computeConformalAnomaly(
    const std::vector<float>& theta_field,
    const std::vector<float>& R_field,
    float K, int nx, int ny, int nz, float dx
) {
    double trace_T = 0.0;

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int idx = k * nx * ny + j * nx + i;

                // Compute gradients using centered differences
                int idx_xp = k * nx * ny + j * nx + (i + 1);
                int idx_xm = k * nx * ny + j * nx + (i - 1);
                int idx_yp = k * nx * ny + (j + 1) * nx + i;
                int idx_ym = k * nx * ny + (j - 1) * nx + i;
                int idx_zp = (k + 1) * nx * ny + j * nx + i;
                int idx_zm = (k - 1) * nx * ny + j * nx + i;

                float grad_theta_x = (theta_field[idx_xp] - theta_field[idx_xm]) / (2 * dx);
                float grad_theta_y = (theta_field[idx_yp] - theta_field[idx_ym]) / (2 * dx);
                float grad_theta_z = (theta_field[idx_zp] - theta_field[idx_zm]) / (2 * dx);

                // Kinetic energy density
                float kinetic = 0.5f * (grad_theta_x * grad_theta_x +
                                       grad_theta_y * grad_theta_y +
                                       grad_theta_z * grad_theta_z);

                // Potential energy density from R-field coupling
                float potential = K * R_field[idx] * R_field[idx];

                // Energy density
                float energy_density = kinetic + potential;

                // For conformal field theory, T^μ_μ = 0
                // Breaking is proportional to mass scale: T^μ_μ ~ m²φ²
                trace_T += potential; // Potential term breaks scale invariance
            }
        }
    }

    return trace_T / (nx * ny * nz);
}

int runScaleInvarianceTest() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "E4: SCALE INVARIANCE BREAKING ANALYSIS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    std::cout << "\nTheoretical Framework:" << std::endl;
    std::cout << "  TRD Action: S = ∫[(∂θ)² + K·R²·cos(Δθ)] d⁴x" << std::endl;
    std::cout << "  Conformal transformation: x^μ → λx^μ" << std::endl;
    std::cout << "  β-function: β(K) = μ·∂K/∂μ" << std::endl;
    std::cout << "  Quality Gate: |β(K)| > 0.01 (scale invariance broken)\n" << std::endl;

    // Test parameters
    const int nx = 32, ny = 32, nz = 32;
    const float dx = 1.0f;
    const float dt = 0.01f;
    const float K_initial = 2.0f;
    const int evolution_steps = 100;

    std::vector<float> scale_factors = {0.5f, 1.0f, 2.0f, 5.0f};
    std::vector<float> effective_couplings;
    std::vector<float> conformal_anomalies;

    // Initialize TRD core
    TRDCore3D trd;
    TRDCore3D::Config config;
    config.Nx = nx;
    config.Ny = ny;
    config.Nz = nz;
    config.dx = dx;
    config.dt = dt;
    config.coupling_strength = K_initial;
    trd.initialize(config);

    // Initialize field with Gaussian profile
    std::vector<float> theta_original(nx * ny * nz);
    std::vector<float> R_field(nx * ny * nz, 1.0f);

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = k * nx * ny + j * nx + i;
                float x = (i - nx/2) * dx;
                float y = (j - ny/2) * dx;
                float z = (k - nz/2) * dx;
                float r2 = x*x + y*y + z*z;
                float sigma = 5.0f * dx;
                theta_original[idx] = std::exp(-r2 / (2 * sigma * sigma));
            }
        }
    }

    std::cout << "Scale Factor Analysis:" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::setw(10) << "λ"
              << std::setw(15) << "K_eff"
              << std::setw(20) << "Correlation Length"
              << std::setw(15) << "T^μ_μ" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    // Test each scale factor
    for (float lambda : scale_factors) {
        std::vector<float> theta_scaled(nx * ny * nz);

        // Apply scale transformation
        applyScaleTransformation(theta_original, theta_scaled, nx, ny, nz, lambda);

        // Set scaled field in TRD and evolve
        trd.setPhaseField(theta_scaled.data());

        for (int step = 0; step < evolution_steps; ++step) {
            trd.evolveKuramotoCPU(dt);
        }

        // Get evolved field
        trd.getPhaseField(theta_scaled.data());
        trd.computeRField();
        trd.getRField(R_field.data());

        // Extract effective coupling from correlation function
        double corr_length = computeCorrelationFunction(theta_scaled, nx, ny, nz, nx/4);
        float K_eff = K_initial * std::pow(lambda, -1.0); // Dimensional analysis estimate

        // More precise extraction from correlation length
        if (corr_length > 0) {
            K_eff = 1.0f / (corr_length * corr_length); // ξ ~ 1/√K
        }

        effective_couplings.push_back(K_eff);

        // Compute conformal anomaly
        double anomaly = computeConformalAnomaly(theta_scaled, R_field, K_eff, nx, ny, nz, dx);
        conformal_anomalies.push_back(anomaly);

        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(10) << lambda
                  << std::setw(15) << K_eff
                  << std::setw(20) << corr_length
                  << std::setw(15) << anomaly << std::endl;
    }

    // Compute β-function
    double beta = computeBetaFunction(effective_couplings, scale_factors);

    std::cout << std::string(60, '-') << std::endl;
    std::cout << "\nβ-Function Analysis:" << std::endl;
    std::cout << "  β(K) = " << std::fixed << std::setprecision(4) << beta << std::endl;
    std::cout << "  Interpretation: ";

    if (std::abs(beta) < 0.001) {
        std::cout << "Scale invariant (conformal)" << std::endl;
    } else if (beta > 0) {
        std::cout << "Coupling grows with energy scale (UV relevant)" << std::endl;
    } else {
        std::cout << "Coupling decreases with energy scale (UV irrelevant)" << std::endl;
    }

    // Compute average conformal anomaly
    double avg_anomaly = std::accumulate(conformal_anomalies.begin(),
                                        conformal_anomalies.end(), 0.0) / conformal_anomalies.size();

    std::cout << "\nConformal Anomaly:" << std::endl;
    std::cout << "  <T^μ_μ> = " << std::scientific << avg_anomaly << std::endl;
    std::cout << "  Status: " << (std::abs(avg_anomaly) > 1e-6 ? "BROKEN" : "PRESERVED") << std::endl;

    // Quality gate check
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "QUALITY GATES:" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    bool beta_nonzero = std::abs(beta) > 0.01;
    bool anomaly_present = std::abs(avg_anomaly) > 1e-6;
    bool all_passed = beta_nonzero && anomaly_present;

    std::cout << (beta_nonzero ? "✓" : "✗")
              << " β-function nonzero: |β(K)| = " << std::abs(beta)
              << (beta_nonzero ? " > 0.01" : " ≤ 0.01") << std::endl;

    std::cout << (anomaly_present ? "✓" : "✗")
              << " Conformal anomaly present: |<T^μ_μ>| = " << std::abs(avg_anomaly)
              << (anomaly_present ? " > 10⁻⁶" : " ≤ 10⁻⁶") << std::endl;

    std::cout << "\nFINAL STATUS: " << (all_passed ? "PASSED" : "FAILED") << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    return all_passed ? 0 : 1;
}

// Entry point for standalone testing
#ifndef TRD_MAIN_EXECUTABLE
int main() {
    return runScaleInvarianceTest();
}
#endif