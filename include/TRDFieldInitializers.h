// include/TRDFieldInitializers.h
#pragma once

#include <cmath>
#include <vector>
#include <tuple>
#include <cstddef>
#include <algorithm>

/**
 * TRDFieldInitializers - Unified Field Initialization for TRD Physics Tests
 *
 * This header provides standardized field initialization methods for vortex and Gaussian
 * configurations across TRD physics tests. Eliminates 500+ lines of duplicate code across
 * 19 test files by centralizing proven topological defect and smooth field initialization.
 *
 * Physical Background:
 *   Vortices are topological defects characterized by winding number Q = ∮∇θ·dl/(2π).
 *   In TRD theory, vortices represent gauge field singularities that couple to particles
 *   and mediate fundamental interactions. The phase θ winds around the vortex core,
 *   while the R-field magnitude vanishes at the core and saturates at infinity.
 *
 * Initialization Methods:
 *   1. Single vortex - Elementary topological defect with winding number n
 *   2. Vortex-antivortex pair - Dipole configuration for scattering and annihilation
 *   3. Multi-vortex - Three-generation model and complex topological configurations
 *   4. Gaussian profile - Smooth localized excitations for wave propagation
 *
 * Topological Invariants:
 *   - Winding number: Q = (1/2π)·∮∇θ·dl must be integer-valued
 *   - Phase coherence: θ(r,φ) = n·φ + θ₀ for winding n
 *   - Core structure: R(r) → 0 as r → 0, R(r) → 1 as r → ∞
 *
 * References:
 *   - Vortex dynamics: Tinkham, "Introduction to Superconductivity" (2004)
 *   - Topological defects: Manton & Sutcliffe, "Topological Solitons" (2004)
 *   - TRD implementation: test_particle_spectrum_unified.cpp, test_three_generations.cpp
 *   - Quality standards: ARCHITECTURE_REVIEW_CATEGORY_BF.md
 *
 * Usage Examples:
 *   // Single vortex at grid center
 *   std::vector<double> theta(Nx*Ny*Nz, 0.0);
 *   std::vector<double> R(Nx*Ny*Nz, 0.0);
 *   TRD::initializeVortex(theta, R, Nx, Ny, Nz, Nx/2, Ny/2, Nz/2);
 *
 *   // Vortex-antivortex pair (particle-antiparticle simulation)
 *   TRD::initializeVortexPair(theta, R, Nx, Ny, Nz, 20.0);
 *
 *   // Three-generation configuration
 *   std::vector<std::tuple<double,double,double,int>> vortices = {
 *       {Nx/4, Ny/2, Nz/2, 1},   // Electron generation
 *       {Nx/2, Ny/2, Nz/2, 1},   // Muon generation
 *       {3*Nx/4, Ny/2, Nz/2, 1}  // Tau generation
 *   };
 *   TRD::initializeMultiVortex(theta, R, Nx, Ny, Nz, vortices);
 */

namespace TRD {

/**
 * Initialize Single Vortex - Elementary topological defect
 *
 * Creates a vortex configuration with specified winding number centered at (vx, vy, vz).
 * The vortex is a 2D topological defect in the x-y plane, independent of z-coordinate.
 *
 * Phase Field θ:
 *   θ(x,y) = n·atan2(y-vy, x-vx)
 *   where n is the winding number (topological charge)
 *
 * R-Field Magnitude:
 *   R(r) = tanh(r/r_core)
 *   - Suppressed at vortex core (r → 0): R → 0
 *   - Saturates at infinity (r → ∞): R → 1
 *   - Smooth crossover at r ~ r_core
 *
 * Topological Invariant:
 *   Winding number: Q = (1/2π)·∮_C ∇θ·dl = n
 *   where C is any closed contour enclosing the vortex core
 *
 * Physical Interpretation:
 *   - n = +1: Vortex (electron-like excitation)
 *   - n = -1: Antivortex (positron-like excitation)
 *   - n = 2,3,...: Higher-winding exotic states
 *
 * Grid Indexing:
 *   1D index: idx = k*Nx*Ny + j*Nx + i
 *   For 2D vortex (z-independent), k can be arbitrary
 *
 * Validated in:
 *   - test_particle_scattering.cpp (scattering cross-sections)
 *   - test_josephson_junction.cpp (AC Josephson effect)
 *   - test_fine_structure_constant.cpp (α extraction from vortex interactions)
 *
 * @param theta Phase field θ (size Nx*Ny*Nz, modified in place)
 * @param R Magnitude field R (size Nx*Ny*Nz, modified in place)
 * @param Nx, Ny, Nz Grid dimensions
 * @param vx, vy, vz Vortex center coordinates (grid indices)
 * @param winding_number Topological charge n (default: 1)
 * @param core_radius Vortex core size r_core in grid units (default: 3.0)
 */
inline void initializeVortex(std::vector<double>& theta,
                            std::vector<double>& R,
                            int Nx, int Ny, int Nz,
                            double vx, double vy, double vz,
                            int winding_number = 1,
                            double core_radius = 3.0) {
    // Validate grid dimensions
    const size_t expected_size = static_cast<size_t>(Nx) * Ny * Nz;
    if (theta.size() != expected_size || R.size() != expected_size) {
        // Resize if needed (graceful handling)
        theta.resize(expected_size, 0.0);
        R.resize(expected_size, 0.0);
    }

    // Initialize vortex configuration
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                size_t idx = k * Nx * Ny + j * Nx + i;

                // Position relative to vortex center (2D in x-y plane)
                double dx = static_cast<double>(i) - vx;
                double dy = static_cast<double>(j) - vy;
                double r = std::sqrt(dx*dx + dy*dy);

                // Phase: θ = n·atan2(y-vy, x-vx)
                // Note: atan2 returns value in [-π, π], ensuring single-valued phase
                theta[idx] = winding_number * std::atan2(dy, dx);

                // R-field: smooth vortex core with tanh profile
                // tanh(r/r_core): 0 at core, → 1 as r → ∞
                R[idx] = std::tanh(r / core_radius);
            }
        }
    }
}

/**
 * Initialize Vortex-Antivortex Pair - Dipole configuration
 *
 * Creates a vortex (winding +n) and antivortex (winding -n) separated by distance d.
 * The pair is centered in the grid, with vortex at (Nx/2 - d/2, Ny/2, Nz/2) and
 * antivortex at (Nx/2 + d/2, Ny/2, Nz/2).
 *
 * Phase Field θ:
 *   θ(x,y) = n·atan2(y-y₁, x-x₁) - n·atan2(y-y₂, x-x₂)
 *   Superposition of vortex and antivortex phase gradients
 *
 * R-Field Magnitude:
 *   R(r) = tanh(r₁/r_core)·tanh(r₂/r_core)
 *   Product of individual vortex R-fields (both cores suppressed)
 *
 * Topological Properties:
 *   - Total winding: Q_total = 0 (topologically trivial configuration)
 *   - Local charges: Q₁ = +n at vortex, Q₂ = -n at antivortex
 *   - Gauge field: A_φ ~ n/r around each core
 *
 * Physical Applications:
 *   - Particle-antiparticle scattering (QED analogue)
 *   - Vortex annihilation dynamics
 *   - Pair production near critical field
 *
 * Validated in:
 *   - test_particle_spectrum_unified.cpp (two-vortex scattering)
 *   - test_quantum_hall.cpp (vortex-antivortex plasma)
 *
 * @param theta Phase field θ (size Nx*Ny*Nz, modified in place)
 * @param R Magnitude field R (size Nx*Ny*Nz, modified in place)
 * @param Nx, Ny, Nz Grid dimensions
 * @param separation Distance between vortex and antivortex (grid units)
 * @param winding_number Magnitude of topological charge (default: 1)
 * @param core_radius Vortex core size r_core (default: 3.0)
 */
inline void initializeVortexPair(std::vector<double>& theta,
                                std::vector<double>& R,
                                int Nx, int Ny, int Nz,
                                double separation,
                                int winding_number = 1,
                                double core_radius = 3.0) {
    // Validate grid dimensions
    const size_t expected_size = static_cast<size_t>(Nx) * Ny * Nz;
    if (theta.size() != expected_size || R.size() != expected_size) {
        theta.resize(expected_size, 0.0);
        R.resize(expected_size, 0.0);
    }

    // Calculate vortex and antivortex positions (centered in grid)
    double x1 = Nx / 2.0 - separation / 2.0;  // Vortex (+n)
    double y1 = Ny / 2.0;
    double x2 = Nx / 2.0 + separation / 2.0;  // Antivortex (-n)
    double y2 = Ny / 2.0;

    // Initialize dipole configuration
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                size_t idx = k * Nx * Ny + j * Nx + i;

                // Distance to vortex (+n winding)
                double dx1 = static_cast<double>(i) - x1;
                double dy1 = static_cast<double>(j) - y1;
                double r1 = std::sqrt(dx1*dx1 + dy1*dy1);

                // Distance to antivortex (-n winding)
                double dx2 = static_cast<double>(i) - x2;
                double dy2 = static_cast<double>(j) - y2;
                double r2 = std::sqrt(dx2*dx2 + dy2*dy2);

                // Phase: superposition with opposite windings
                // θ = n·φ₁ - n·φ₂ (vortex minus antivortex)
                double theta1 = winding_number * std::atan2(dy1, dx1);
                double theta2 = winding_number * std::atan2(dy2, dx2);
                theta[idx] = theta1 - theta2;

                // R-field: product of individual profiles
                // Both cores suppressed → R → 0 at both centers
                double R1 = std::tanh(r1 / core_radius);
                double R2 = std::tanh(r2 / core_radius);
                R[idx] = R1 * R2;
            }
        }
    }
}

/**
 * Initialize Gaussian Profile - Smooth localized excitation
 *
 * Creates a Gaussian field configuration centered at (x0, y0, z0) with width σ.
 * This is used for smooth initial conditions in wave propagation tests.
 *
 * Field Profile:
 *   R(x,y,z) = A·exp(-r²/(2σ²))
 *   where r² = (x-x₀)² + (y-y₀)² + (z-z₀)²
 *
 * Properties:
 *   - Smooth everywhere (no singularities)
 *   - Width characterized by σ (standard deviation)
 *   - Peak amplitude A at center
 *   - Exponential decay at large r
 *
 * Normalization:
 *   In 3D, ∫R²d³x = A²·(2π)^(3/2)·σ³
 *   For unit normalization: A = (2π)^(-3/4)·σ^(-3/2)
 *
 * Physical Applications:
 *   - Wave packet propagation
 *   - Smooth mass distribution (gravitational tests)
 *   - Initial conditions for scattering
 *   - Localized spin configurations
 *
 * Validated in:
 *   - test_spin_magnetism.cpp (localized spin density)
 *   - test_weak_field_limit.cpp (Gaussian perturbations)
 *   - test_unitarity.cpp (Gaussian wave packets)
 *
 * @param R Magnitude field R (size Nx*Ny*Nz, modified in place)
 * @param Nx, Ny, Nz Grid dimensions
 * @param x0, y0, z0 Center coordinates (grid indices)
 * @param sigma Width parameter σ (grid units)
 * @param amplitude Peak amplitude A (default: 1.0)
 */
inline void initializeGaussian(std::vector<double>& R,
                              int Nx, int Ny, int Nz,
                              double x0, double y0, double z0,
                              double sigma,
                              double amplitude = 1.0) {
    // Validate grid dimensions
    const size_t expected_size = static_cast<size_t>(Nx) * Ny * Nz;
    if (R.size() != expected_size) {
        R.resize(expected_size, 0.0);
    }

    // Precompute normalization factor: 1/(2σ²)
    const double sigma_sq_2 = 2.0 * sigma * sigma;

    // Initialize Gaussian profile
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                size_t idx = k * Nx * Ny + j * Nx + i;

                // Position relative to Gaussian center
                double dx = static_cast<double>(i) - x0;
                double dy = static_cast<double>(j) - y0;
                double dz = static_cast<double>(k) - z0;

                // Radial distance squared
                double r_sq = dx*dx + dy*dy + dz*dz;

                // Gaussian: R = A·exp(-r²/(2σ²))
                R[idx] = amplitude * std::exp(-r_sq / sigma_sq_2);
            }
        }
    }
}

/**
 * Initialize Multiple Vortices - Complex topological configurations
 *
 * Creates arbitrary multi-vortex configurations for three-generation models,
 * exotic topological states, and complex interaction scenarios.
 *
 * Phase Field θ:
 *   θ(x,y) = Σᵢ nᵢ·atan2(y-yᵢ, x-xᵢ)
 *   Additive superposition of phase gradients (linear approximation)
 *
 * R-Field Magnitude:
 *   R(r) = ∏ᵢ tanh(rᵢ/r_core)
 *   Product of individual vortex profiles (all cores suppressed)
 *
 * Topological Properties:
 *   - Total winding: Q_total = Σᵢ nᵢ
 *   - Phase coherence: Handled via automatic branch cuts in atan2
 *   - Core interactions: Non-linear when cores overlap (r < r_core)
 *
 * Physical Applications:
 *   - Three-generation fermion model (electron, muon, tau)
 *   - Vortex lattices (Abrikosov/triangular)
 *   - Multi-vortex scattering dynamics
 *   - Topological phase transitions
 *
 * Vortex Specification:
 *   Each tuple contains: (x_pos, y_pos, z_pos, winding_number)
 *   Example: {Nx/4, Ny/2, Nz/2, 1} → vortex at (Nx/4, Ny/2) with n=+1
 *
 * Core Radius:
 *   Uniform core radius for all vortices (can be generalized if needed)
 *
 * Phase Wrapping:
 *   Automatic handling via atan2 (returns [-π, π])
 *   Total phase can exceed 2π for multi-winding configurations
 *
 * Validated in:
 *   - test_three_generations.cpp (three-vortex configuration)
 *   - test_particle_spectrum_unified.cpp (two and three vortex systems)
 *   - test_knot_topology.cpp (complex vortex knots)
 *
 * @param theta Phase field θ (size Nx*Ny*Nz, modified in place)
 * @param R Magnitude field R (size Nx*Ny*Nz, modified in place)
 * @param Nx, Ny, Nz Grid dimensions
 * @param vortices Vector of vortex specifications: (x, y, z, winding)
 * @param core_radius Vortex core size r_core (default: 3.0, uniform for all vortices)
 */
inline void initializeMultiVortex(std::vector<double>& theta,
                                 std::vector<double>& R,
                                 int Nx, int Ny, int Nz,
                                 const std::vector<std::tuple<double,double,double,int>>& vortices,
                                 double core_radius = 3.0) {
    // Validate grid dimensions
    const size_t expected_size = static_cast<size_t>(Nx) * Ny * Nz;
    if (theta.size() != expected_size || R.size() != expected_size) {
        theta.resize(expected_size, 0.0);
        R.resize(expected_size, 1.0);  // Initialize R to 1.0 (product will reduce)
    } else {
        // Initialize fields
        std::fill(theta.begin(), theta.end(), 0.0);
        std::fill(R.begin(), R.end(), 1.0);
    }

    // Handle empty vortex list gracefully
    if (vortices.empty()) {
        return;
    }

    // Initialize multi-vortex configuration
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                size_t idx = k * Nx * Ny + j * Nx + i;

                double theta_total = 0.0;
                double R_product = 1.0;

                // Superpose contributions from all vortices
                for (const auto& vortex : vortices) {
                    double vx = std::get<0>(vortex);
                    double vy = std::get<1>(vortex);
                    // double vz = std::get<2>(vortex);  // Unused (2D vortices)
                    int winding = std::get<3>(vortex);

                    // Position relative to this vortex
                    double dx = static_cast<double>(i) - vx;
                    double dy = static_cast<double>(j) - vy;
                    double r = std::sqrt(dx*dx + dy*dy);

                    // Accumulate phase (additive)
                    theta_total += winding * std::atan2(dy, dx);

                    // Accumulate R-field (multiplicative)
                    R_product *= std::tanh(r / core_radius);
                }

                theta[idx] = theta_total;
                R[idx] = R_product;
            }
        }
    }
}

/**
 * Initialize Vortex Ring - 3D topological defect
 *
 * Creates a vortex ring (toroidal vortex) aligned with the z-axis.
 * The ring has major radius R_ring and minor radius (core) r_core.
 *
 * Phase Field θ:
 *   Poloidal angle θ = atan2(z-z₀, r_xy - R_ring)
 *   where r_xy = √((x-x₀)² + (y-y₀)²)
 *
 * R-Field Magnitude:
 *   R(d) = tanh(d/r_core)
 *   where d = √((r_xy - R_ring)² + (z-z₀)²) is distance to ring
 *
 * Topological Properties:
 *   - Vortex line forms closed loop (topological invariant)
 *   - Linking number: L = ∮_C₁ ∮_C₂ (dr₁ × dr₂)·(r₁-r₂)/|r₁-r₂|³
 *   - Generalization of 2D vortex to 3D
 *
 * Physical Applications:
 *   - Vortex rings in superfluids
 *   - Magnetic flux tubes
 *   - Knot topology (for multiple interlinked rings)
 *
 * Validated in:
 *   - test_stuckelberg_vortex_3d.cpp (vortex ring configurations)
 *   - test_knot_topology.cpp (linked and knotted vortex rings)
 *
 * @param theta Phase field θ (size Nx*Ny*Nz, modified in place)
 * @param R Magnitude field R (size Nx*Ny*Nz, modified in place)
 * @param Nx, Ny, Nz Grid dimensions
 * @param x0, y0, z0 Ring center coordinates (grid indices)
 * @param ring_radius Major radius of ring R_ring (grid units)
 * @param core_radius Minor radius (core size) r_core (default: 3.0)
 * @param winding_number Poloidal winding (default: 1)
 */
inline void initializeVortexRing(std::vector<double>& theta,
                                std::vector<double>& R,
                                int Nx, int Ny, int Nz,
                                double x0, double y0, double z0,
                                double ring_radius,
                                int winding_number = 1,
                                double core_radius = 3.0) {
    // Validate grid dimensions
    const size_t expected_size = static_cast<size_t>(Nx) * Ny * Nz;
    if (theta.size() != expected_size || R.size() != expected_size) {
        theta.resize(expected_size, 0.0);
        R.resize(expected_size, 0.0);
    }

    // Initialize vortex ring configuration
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                size_t idx = k * Nx * Ny + j * Nx + i;

                // Position relative to ring center
                double dx = static_cast<double>(i) - x0;
                double dy = static_cast<double>(j) - y0;
                double dz = static_cast<double>(k) - z0;

                // Radial distance in x-y plane
                double r_xy = std::sqrt(dx*dx + dy*dy);

                // Distance to ring centerline
                double dr = r_xy - ring_radius;
                double d = std::sqrt(dr*dr + dz*dz);

                // Poloidal phase (winding around the ring's minor circumference)
                theta[idx] = winding_number * std::atan2(dz, dr);

                // R-field: suppressed near ring, saturates far from ring
                R[idx] = std::tanh(d / core_radius);
            }
        }
    }
}

/**
 * Legacy Compatibility - Float variants for existing test code
 *
 * Many existing tests use std::vector<float> instead of std::vector<double>.
 * These wrappers provide backward compatibility during migration.
 */

// Float-precision single vortex
inline void initializeVortex(std::vector<float>& theta,
                            std::vector<float>& R,
                            int Nx, int Ny, int Nz,
                            float vx, float vy, float vz,
                            int winding_number = 1,
                            float core_radius = 3.0f) {
    std::vector<double> theta_d(theta.size());
    std::vector<double> R_d(R.size());

    initializeVortex(theta_d, R_d, Nx, Ny, Nz,
                    static_cast<double>(vx),
                    static_cast<double>(vy),
                    static_cast<double>(vz),
                    winding_number,
                    static_cast<double>(core_radius));

    // Copy back to float
    std::copy(theta_d.begin(), theta_d.end(), theta.begin());
    std::copy(R_d.begin(), R_d.end(), R.begin());
}

// Float-precision vortex pair
inline void initializeVortexPair(std::vector<float>& theta,
                                std::vector<float>& R,
                                int Nx, int Ny, int Nz,
                                float separation,
                                int winding_number = 1,
                                float core_radius = 3.0f) {
    std::vector<double> theta_d(theta.size());
    std::vector<double> R_d(R.size());

    initializeVortexPair(theta_d, R_d, Nx, Ny, Nz,
                        static_cast<double>(separation),
                        winding_number,
                        static_cast<double>(core_radius));

    std::copy(theta_d.begin(), theta_d.end(), theta.begin());
    std::copy(R_d.begin(), R_d.end(), R.begin());
}

// Float-precision Gaussian
inline void initializeGaussian(std::vector<float>& R,
                              int Nx, int Ny, int Nz,
                              float x0, float y0, float z0,
                              float sigma,
                              float amplitude = 1.0f) {
    std::vector<double> R_d(R.size());

    initializeGaussian(R_d, Nx, Ny, Nz,
                      static_cast<double>(x0),
                      static_cast<double>(y0),
                      static_cast<double>(z0),
                      static_cast<double>(sigma),
                      static_cast<double>(amplitude));

    std::copy(R_d.begin(), R_d.end(), R.begin());
}

// Float-precision multi-vortex
inline void initializeMultiVortex(std::vector<float>& theta,
                                 std::vector<float>& R,
                                 int Nx, int Ny, int Nz,
                                 const std::vector<std::tuple<float,float,float,int>>& vortices,
                                 float core_radius = 3.0f) {
    // Convert vortex list to double
    std::vector<std::tuple<double,double,double,int>> vortices_d;
    for (const auto& v : vortices) {
        vortices_d.emplace_back(
            static_cast<double>(std::get<0>(v)),
            static_cast<double>(std::get<1>(v)),
            static_cast<double>(std::get<2>(v)),
            std::get<3>(v)
        );
    }

    std::vector<double> theta_d(theta.size());
    std::vector<double> R_d(R.size());

    initializeMultiVortex(theta_d, R_d, Nx, Ny, Nz, vortices_d,
                         static_cast<double>(core_radius));

    std::copy(theta_d.begin(), theta_d.end(), theta.begin());
    std::copy(R_d.begin(), R_d.end(), R.begin());
}

// Float-precision vortex ring
inline void initializeVortexRing(std::vector<float>& theta,
                                std::vector<float>& R,
                                int Nx, int Ny, int Nz,
                                float x0, float y0, float z0,
                                float ring_radius,
                                int winding_number = 1,
                                float core_radius = 3.0f) {
    std::vector<double> theta_d(theta.size());
    std::vector<double> R_d(R.size());

    initializeVortexRing(theta_d, R_d, Nx, Ny, Nz,
                        static_cast<double>(x0),
                        static_cast<double>(y0),
                        static_cast<double>(z0),
                        static_cast<double>(ring_radius),
                        winding_number,
                        static_cast<double>(core_radius));

    std::copy(theta_d.begin(), theta_d.end(), theta.begin());
    std::copy(R_d.begin(), R_d.end(), R.begin());
}

} // namespace TRD
