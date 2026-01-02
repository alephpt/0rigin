/**
 * test_einstein_field_equations.cpp
 *
 * Einstein Field Equation Derivation Test
 *
 * Goal: Verify G_μν = 8πG·T_μν from SMFT metric
 *
 * Physics:
 *   - SMFT metric: g_μν = R²(x,y,z)·η_μν where η_μν is Minkowski
 *   - Christoffel symbols: Γ^λ_μν = (1/2)g^λρ(∂_μ g_νρ + ∂_ν g_ρμ - ∂_ρ g_μν)
 *   - Riemann tensor: R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
 *   - Ricci tensor: R_μν = R^λ_μλν
 *   - Einstein tensor: G_μν = R_μν - (1/2)g_μν R
 *   - EM stress-energy: T_μν = (1/4π)[F_μα F_ν^α - (1/4)g_μν F²]
 *
 * Quality Gate: |G_μν - 8πG·T_μν| < 10^-12 for all components
 */

#include "Maxwell3D.h"
#include "SMFTCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>

// Physical constants (natural units: c = 1, G = 1)
const float PI = 3.14159265358979323846f;
const float COUPLING_8PI_G = 8.0f * PI;  // 8πG in natural units

/**
 * Metric tensor structure (4x4 spacetime metric)
 */
struct MetricTensor {
    std::array<std::array<float, 4>, 4> g;      // Covariant metric g_μν
    std::array<std::array<float, 4>, 4> g_inv;  // Contravariant metric g^μν

    MetricTensor() {
        // Initialize as flat Minkowski metric
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                g[mu][nu] = 0.0f;
                g_inv[mu][nu] = 0.0f;
            }
        }
        g[0][0] = -1.0f;  // g_00 = -1 (time component)
        g[1][1] = 1.0f;   // g_11 = 1 (space component)
        g[2][2] = 1.0f;   // g_22 = 1
        g[3][3] = 1.0f;   // g_33 = 1

        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
        g_inv[2][2] = 1.0f;
        g_inv[3][3] = 1.0f;
    }

    // Scale by R² field
    void scaleByRField(float R) {
        float R2 = R * R;
        g[0][0] = -R2;
        g[1][1] = R2;
        g[2][2] = R2;
        g[3][3] = R2;

        // Inverse scales by 1/R²
        float inv_R2 = 1.0f / R2;
        g_inv[0][0] = -inv_R2;
        g_inv[1][1] = inv_R2;
        g_inv[2][2] = inv_R2;
        g_inv[3][3] = inv_R2;
    }
};

/**
 * Christoffel symbols structure
 * Γ^λ_μν with indices: lambda (upper), mu (lower), nu (lower)
 */
struct ChristoffelSymbols {
    std::array<std::array<std::array<float, 4>, 4>, 4> gamma;  // gamma[lambda][mu][nu]

    ChristoffelSymbols() {
        for (int lambda = 0; lambda < 4; ++lambda) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    gamma[lambda][mu][nu] = 0.0f;
                }
            }
        }
    }
};

/**
 * Calculate Christoffel symbols from metric and its derivatives
 * Γ^λ_μν = (1/2)g^λρ(∂_μ g_νρ + ∂_ν g_ρμ - ∂_ρ g_μν)
 */
ChristoffelSymbols calculateChristoffel(
    const MetricTensor& metric,
    const std::array<std::array<std::array<float, 4>, 4>, 4>& dg  // dg[alpha][mu][nu] = ∂_α g_μν
) {
    ChristoffelSymbols christoffel;

    for (int lambda = 0; lambda < 4; ++lambda) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                float sum = 0.0f;
                for (int rho = 0; rho < 4; ++rho) {
                    sum += metric.g_inv[lambda][rho] *
                           (dg[mu][nu][rho] + dg[nu][rho][mu] - dg[rho][mu][nu]);
                }
                christoffel.gamma[lambda][mu][nu] = 0.5f * sum;
            }
        }
    }

    return christoffel;
}

/**
 * Calculate Ricci tensor from Christoffel symbols and their derivatives
 * R_μν = ∂_λ Γ^λ_μν - ∂_ν Γ^λ_μλ + Γ^λ_λρ Γ^ρ_μν - Γ^λ_νρ Γ^ρ_μλ
 */
std::array<std::array<float, 4>, 4> calculateRicciTensor(
    const ChristoffelSymbols& christoffel,
    const std::array<std::array<std::array<std::array<float, 4>, 4>, 4>, 4>& d_gamma  // d_gamma[alpha][lambda][mu][nu] = ∂_α Γ^λ_μν
) {
    std::array<std::array<float, 4>, 4> ricci;

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            float sum = 0.0f;

            // ∂_λ Γ^λ_μν term
            for (int lambda = 0; lambda < 4; ++lambda) {
                sum += d_gamma[lambda][lambda][mu][nu];
            }

            // -∂_ν Γ^λ_μλ term
            for (int lambda = 0; lambda < 4; ++lambda) {
                sum -= d_gamma[nu][lambda][mu][lambda];
            }

            // Γ^λ_λρ Γ^ρ_μν term
            for (int lambda = 0; lambda < 4; ++lambda) {
                for (int rho = 0; rho < 4; ++rho) {
                    sum += christoffel.gamma[lambda][lambda][rho] *
                           christoffel.gamma[rho][mu][nu];
                }
            }

            // -Γ^λ_νρ Γ^ρ_μλ term
            for (int lambda = 0; lambda < 4; ++lambda) {
                for (int rho = 0; rho < 4; ++rho) {
                    sum -= christoffel.gamma[lambda][nu][rho] *
                           christoffel.gamma[rho][mu][lambda];
                }
            }

            ricci[mu][nu] = sum;
        }
    }

    return ricci;
}

/**
 * Calculate Einstein tensor from Ricci tensor and scalar curvature
 * G_μν = R_μν - (1/2)g_μν R
 */
std::array<std::array<float, 4>, 4> calculateEinsteinTensor(
    const std::array<std::array<float, 4>, 4>& ricci,
    const MetricTensor& metric
) {
    // Calculate Ricci scalar: R = g^μν R_μν
    float ricci_scalar = 0.0f;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            ricci_scalar += metric.g_inv[mu][nu] * ricci[mu][nu];
        }
    }

    // Calculate Einstein tensor: G_μν = R_μν - (1/2)g_μν R
    std::array<std::array<float, 4>, 4> einstein;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            einstein[mu][nu] = ricci[mu][nu] - 0.5f * metric.g[mu][nu] * ricci_scalar;
        }
    }

    return einstein;
}

/**
 * Calculate electromagnetic stress-energy tensor
 * T_μν = (1/4π)[F_μα F_ν^α - (1/4)g_μν F²]
 * where F² = F_μν F^μν
 */
std::array<std::array<float, 4>, 4> calculateEMStressEnergy(
    float Ex, float Ey, float Ez,
    float Bx, float By, float Bz,
    const MetricTensor& metric
) {
    // Build electromagnetic field tensor F_μν
    // F_0i = -E_i, F_ij = ε_ijk B_k
    std::array<std::array<float, 4>, 4> F;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            F[mu][nu] = 0.0f;
        }
    }

    // Electric field components
    F[0][1] = -Ex; F[1][0] = Ex;
    F[0][2] = -Ey; F[2][0] = Ey;
    F[0][3] = -Ez; F[3][0] = Ez;

    // Magnetic field components
    F[1][2] = -Bz; F[2][1] = Bz;
    F[2][3] = -Bx; F[3][2] = Bx;
    F[3][1] = -By; F[1][3] = By;

    // Raise one index: F^μ_ν = g^μρ F_ρν
    std::array<std::array<float, 4>, 4> F_up;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            float sum = 0.0f;
            for (int rho = 0; rho < 4; ++rho) {
                sum += metric.g_inv[mu][rho] * F[rho][nu];
            }
            F_up[mu][nu] = sum;
        }
    }

    // Calculate F²: F² = F_μν F^μν
    float F_squared = 0.0f;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            F_squared += F[mu][nu] * F_up[mu][nu];
        }
    }

    // Calculate stress-energy tensor
    // T_μν = (1/4π)[F_μα F^α_ν - (1/4)g_μν F²]
    std::array<std::array<float, 4>, 4> T;
    const float prefactor = 1.0f / (4.0f * PI);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            float sum = 0.0f;
            // F_μα F^α_ν term
            for (int alpha = 0; alpha < 4; ++alpha) {
                sum += F[mu][alpha] * F_up[alpha][nu];
            }
            // Subtract (1/4)g_μν F² term
            T[mu][nu] = prefactor * (sum - 0.25f * metric.g[mu][nu] * F_squared);
        }
    }

    return T;
}

/**
 * Compute spatial derivatives using central differences
 */
float computeDerivative(
    const std::vector<float>& field,
    uint32_t i, uint32_t j, uint32_t k,
    uint32_t Nx, uint32_t Ny, uint32_t Nz,
    int direction,  // 0=x, 1=y, 2=z
    float dx
) {
    auto index3D = [Nx, Ny](uint32_t i, uint32_t j, uint32_t k) {
        return k * (Nx * Ny) + j * Nx + i;
    };

    auto wrapX = [Nx](int32_t x) { return (x + Nx) % Nx; };
    auto wrapY = [Ny](int32_t y) { return (y + Ny) % Ny; };
    auto wrapZ = [Nz](int32_t z) { return (z + Nz) % Nz; };

    float derivative = 0.0f;

    if (direction == 0) {  // ∂/∂x
        uint32_t ip = wrapX(i + 1);
        uint32_t im = wrapX(i - 1);
        derivative = (field[index3D(ip, j, k)] - field[index3D(im, j, k)]) / (2.0f * dx);
    } else if (direction == 1) {  // ∂/∂y
        uint32_t jp = wrapY(j + 1);
        uint32_t jm = wrapY(j - 1);
        derivative = (field[index3D(i, jp, k)] - field[index3D(i, jm, k)]) / (2.0f * dx);
    } else if (direction == 2) {  // ∂/∂z
        uint32_t kp = wrapZ(k + 1);
        uint32_t km = wrapZ(k - 1);
        derivative = (field[index3D(i, j, kp)] - field[index3D(i, j, km)]) / (2.0f * dx);
    }

    return derivative;
}

/**
 * Einstein Field Equations Test Entry Point
 * Called from main.cpp when: ./smft --test config/einstein_field_equations.yaml
 */
int runEinsteinFieldEquationsTest() {
    std::cout << "========================================\n";
    std::cout << "  Einstein Field Equations Validation\n";
    std::cout << "========================================\n\n";

    // Grid parameters
    const uint32_t Nx = 32;
    const uint32_t Ny = 32;
    const uint32_t Nz = 32;
    const float dx = 1.0f;
    const float dt = 0.01f;
    const int evolution_steps = 100;
    const int sample_points = 10;

    std::cout << "Grid: " << Nx << "x" << Ny << "x" << Nz << "\n";
    std::cout << "Evolution steps: " << evolution_steps << "\n";
    std::cout << "Sample points: " << sample_points << "\n\n";

    // Initialize 3D SMFT and Maxwell
    SMFTCore3D core;
    SMFTCore3D::Config config;
    config.Nx = Nx;
    config.Ny = Ny;
    config.Nz = Nz;
    config.dx = dx;
    config.dt = dt;
    config.coupling_strength = 1.0f;
    core.initialize(config);

    Maxwell3D maxwell(Nx, Ny, Nz);

    // Initialize with uniform magnetic field in z-direction
    std::cout << "Initializing uniform B_z = 1.0 field...\n";
    uint32_t N_total = Nx * Ny * Nz;
    std::vector<float> Ex(N_total, 0.0f);
    std::vector<float> Ey(N_total, 0.0f);
    std::vector<float> Ez(N_total, 0.0f);
    std::vector<float> Bx(N_total, 0.0f);
    std::vector<float> By(N_total, 0.0f);
    std::vector<float> Bz(N_total, 1.0f);  // Uniform B_z field

    maxwell.initialize(Ex, Ey, Ez, Bx, By, Bz);

    // Initialize SMFT with random phase to generate non-trivial R-field
    core.initializeRandom(42);

    // Evolve system to steady state
    std::cout << "Evolving to steady state...\n";
    for (int step = 0; step < evolution_steps; ++step) {
        // Evolve SMFT Kuramoto dynamics
        core.evolveKuramotoCPU(dt);
        core.computeRField();

        // Evolve Maxwell fields
        maxwell.step(dt);

        if (step % 20 == 0) {
            float avg_R = core.getAverageR();
            float em_energy = maxwell.getTotalEnergy();
            std::cout << "  Step " << step << ": <R> = " << avg_R
                      << ", EM energy = " << em_energy << "\n";
        }
    }

    std::cout << "\n=== Validating Einstein Field Equations ===\n\n";

    // Sample grid points for validation
    std::vector<float> max_residuals(10, 0.0f);  // Track max residual for each component
    std::vector<std::string> component_names = {
        "G_00", "G_01", "G_02", "G_03",
        "G_11", "G_12", "G_13",
        "G_22", "G_23",
        "G_33"
    };

    // Generate sample points evenly distributed in grid
    for (int sample = 0; sample < sample_points; ++sample) {
        // Sample point coordinates
        uint32_t i = (sample * Nx / sample_points + Nx/4) % Nx;
        uint32_t j = (sample * Ny / sample_points + Ny/3) % Ny;
        uint32_t k = (sample * Nz / sample_points + Nz/5) % Nz;
        uint32_t idx = core.index3D(i, j, k);

        std::cout << "Sample point " << (sample + 1) << " at ("
                  << i << ", " << j << ", " << k << "):\n";

        // Get R-field value at this point
        float R = core.getRField()[idx];
        std::cout << "  R-field value: " << R << "\n";

        // Construct metric: g_μν = R²·η_μν
        MetricTensor metric;
        metric.scaleByRField(R);

        // Calculate metric derivatives
        // For SMFT metric, only spatial derivatives of R matter
        // ∂_μ g_νρ = 2R·∂_μR·η_νρ for spatial components
        std::array<std::array<std::array<float, 4>, 4>, 4> dg;
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    dg[alpha][mu][nu] = 0.0f;
                }
            }
        }

        // Compute spatial derivatives of R
        float dR_dx = computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 0, dx);
        float dR_dy = computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 1, dx);
        float dR_dz = computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 2, dx);

        // Set metric derivatives (2R·∂_i R for diagonal spatial components)
        dg[1][0][0] = 2.0f * R * dR_dx * (-1.0f);  // ∂_x g_00
        dg[1][1][1] = 2.0f * R * dR_dx;             // ∂_x g_11
        dg[1][2][2] = 2.0f * R * dR_dx;             // ∂_x g_22
        dg[1][3][3] = 2.0f * R * dR_dx;             // ∂_x g_33

        dg[2][0][0] = 2.0f * R * dR_dy * (-1.0f);  // ∂_y g_00
        dg[2][1][1] = 2.0f * R * dR_dy;             // ∂_y g_11
        dg[2][2][2] = 2.0f * R * dR_dy;             // ∂_y g_22
        dg[2][3][3] = 2.0f * R * dR_dy;             // ∂_y g_33

        dg[3][0][0] = 2.0f * R * dR_dz * (-1.0f);  // ∂_z g_00
        dg[3][1][1] = 2.0f * R * dR_dz;             // ∂_z g_11
        dg[3][2][2] = 2.0f * R * dR_dz;             // ∂_z g_22
        dg[3][3][3] = 2.0f * R * dR_dz;             // ∂_z g_33

        // Calculate Christoffel symbols
        ChristoffelSymbols christoffel = calculateChristoffel(metric, dg);

        // For simplicity in this first implementation, we'll use numerical derivatives
        // of Christoffel symbols. In production, this should be analytical.
        // For now, we'll compute a simplified Ricci tensor

        // Simplified Ricci calculation for diagonal metric
        // For g_μν = R²·η_μν, the Ricci tensor has a known form
        std::array<std::array<float, 4>, 4> ricci;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                ricci[mu][nu] = 0.0f;
            }
        }

        // For conformally flat metric g_μν = R²·η_μν in 3+1D:
        // R_μν = -2∇²(ln R)·η_μν - 2∂_μ(ln R)∂_ν(ln R)
        float laplacian_ln_R = (1.0f / R) * (
            computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 0, dx) * dR_dx +
            computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 1, dx) * dR_dy +
            computeDerivative(core.getRField(), i, j, k, Nx, Ny, Nz, 2, dx) * dR_dz
        ) / R;

        // Diagonal components (simplified)
        ricci[0][0] = -2.0f * laplacian_ln_R * metric.g[0][0];
        ricci[1][1] = -2.0f * laplacian_ln_R * metric.g[1][1];
        ricci[2][2] = -2.0f * laplacian_ln_R * metric.g[2][2];
        ricci[3][3] = -2.0f * laplacian_ln_R * metric.g[3][3];

        // Calculate Einstein tensor
        std::array<std::array<float, 4>, 4> einstein = calculateEinsteinTensor(ricci, metric);

        // Get EM fields at this point
        float Ex_val = maxwell.getEx()[idx];
        float Ey_val = maxwell.getEy()[idx];
        float Ez_val = maxwell.getEz()[idx];
        float Bx_val = maxwell.getBx()[idx];
        float By_val = maxwell.getBy()[idx];
        float Bz_val = maxwell.getBz()[idx];

        // Calculate EM stress-energy tensor
        std::array<std::array<float, 4>, 4> T_em = calculateEMStressEnergy(
            Ex_val, Ey_val, Ez_val,
            Bx_val, By_val, Bz_val,
            metric
        );

        // Verify Einstein field equations: G_μν = 8πG·T_μν
        std::cout << "  Component residuals |G_μν - 8πG·T_μν|:\n";
        int comp_idx = 0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu; nu < 4; ++nu) {  // Only upper triangular (symmetric)
                float G_component = einstein[mu][nu];
                float T_component = COUPLING_8PI_G * T_em[mu][nu];
                float residual = std::abs(G_component - T_component);

                max_residuals[comp_idx] = std::max(max_residuals[comp_idx], residual);

                std::cout << "    " << component_names[comp_idx] << ": "
                          << std::scientific << std::setprecision(3)
                          << "G=" << G_component
                          << ", 8πG·T=" << T_component
                          << ", residual=" << residual << "\n";

                comp_idx++;
            }
        }
        std::cout << "\n";
    }

    // Report overall results
    std::cout << "=== Maximum Residuals Across All Sample Points ===\n";
    float max_overall_residual = 0.0f;
    for (int i = 0; i < 10; ++i) {
        std::cout << "  " << component_names[i] << ": "
                  << std::scientific << std::setprecision(3)
                  << max_residuals[i] << "\n";
        max_overall_residual = std::max(max_overall_residual, max_residuals[i]);
    }

    std::cout << "\n=== Quality Gate ===\n";
    const float threshold = 1.0e-12f;
    bool pass = max_overall_residual < threshold;

    std::cout << "Maximum residual: " << std::scientific << std::setprecision(3)
              << max_overall_residual << "\n";
    std::cout << "Threshold: " << threshold << "\n";
    std::cout << "Status: " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";

    // Save residual field for visualization (optional)
    std::cout << "\nSaving residual field to output/einstein_residuals.dat...\n";
    std::ofstream outfile("output/einstein_residuals.dat");
    if (outfile.is_open()) {
        outfile << "# x y z residual\n";
        for (uint32_t k = 0; k < Nz; k += 4) {  // Subsample for manageable file size
            for (uint32_t j = 0; j < Ny; j += 4) {
                for (uint32_t i = 0; i < Nx; i += 4) {
                    // Compute max residual at this point (simplified)
                    outfile << i << " " << j << " " << k << " "
                            << max_overall_residual << "\n";  // Placeholder
                }
            }
        }
        outfile.close();
    }

    std::cout << "\n========================================\n";
    std::cout << "FINAL VERDICT: " << (pass ? "PASS ✓" : "FAIL ✗") << "\n";
    std::cout << "========================================\n";

    return pass ? 0 : 1;
}