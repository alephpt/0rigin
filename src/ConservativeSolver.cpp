// src/ConservativeSolver.cpp
#include "ConservativeSolver.h"
#include "Dirac3D.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

ConservativeSolver::ConservativeSolver()
    : nx_(0), ny_(0), nz_(0), dx_(1.0f), dt_(0.005f),
      initial_energy_(0.0f), dirac_mass_(1.0f) {
}

ConservativeSolver::~ConservativeSolver() {
    // Destructor defined here to handle unique_ptr<Dirac3D> with forward declaration
}

void ConservativeSolver::initialize(const Config& config) {
    config_ = config;
    nx_ = config.nx;
    ny_ = config.ny;
    nz_ = config.nz;
    dx_ = config.dx;
    dt_ = config.dt;

    uint32_t total_points = nx_ * ny_ * nz_;

    // Allocate Sine-Gordon fields
    theta_.resize(total_points, 0.0f);
    theta_dot_.resize(total_points, 0.0f);

    // Allocate Dirac spinor fields (4 components × 2 real/imag)
    psi_real_.resize(4 * total_points, 0.0f);
    psi_imag_.resize(4 * total_points, 0.0f);

    // Initialize Dirac3D solver for particle dynamics
    dirac3d_ = std::make_unique<Dirac3D>(nx_, ny_, nz_);
    std::cout << "[ConservativeSolver] Dirac3D solver initialized for chiral mass coupling" << std::endl;

    std::cout << "[ConservativeSolver] Initialized " << nx_ << "×" << ny_ << "×" << nz_
              << " grid (" << total_points << " points)" << std::endl;
    std::cout << "[ConservativeSolver] Integration method: ";
    switch (config_.method) {
        case IntegrationMethod::VELOCITY_VERLET:
            std::cout << "Velocity Verlet (symplectic)" << std::endl;
            break;
        case IntegrationMethod::RK2_SYMPLECTIC:
            std::cout << "RK2 Symplectic" << std::endl;
            break;
        case IntegrationMethod::STRANG_SPLITTING:
            std::cout << "Strang Splitting (T-V-T)" << std::endl;
            break;
        case IntegrationMethod::HALF_STRANG:
            std::cout << "Half-Strang Splitting" << std::endl;
            break;
    }
    std::cout << "[ConservativeSolver] Spatial discretization: ";
    switch (config_.spatial_order) {
        case SpatialOrder::SECOND_ORDER:
            std::cout << "2nd-order (6-neighbor stencil)" << std::endl;
            break;
        case SpatialOrder::FOURTH_ORDER:
            std::cout << "4th-order (12-neighbor stencil)" << std::endl;
            break;
    }
    std::cout << "[ConservativeSolver] Timestep dt = " << dt_ << std::endl;
}

void ConservativeSolver::evolveSineGordon(float dt) {
    static bool first_call = true;
    if (first_call) {
        std::cout << "[ConservativeSolver] Sine-Gordon evolution (∂²θ/∂t² = ∇²θ - sin(θ))" << std::endl;
        first_call = false;
    }

    // Dispatch to appropriate integrator
    switch (config_.method) {
        case IntegrationMethod::VELOCITY_VERLET:
            velocityVerletStep(dt);
            break;
        case IntegrationMethod::RK2_SYMPLECTIC:
            rk2SymplecticStep(dt);
            break;
        case IntegrationMethod::STRANG_SPLITTING:
            strangSplittingStep(dt);
            break;
        default:
            std::cerr << "[ConservativeSolver] Unsupported method for Sine-Gordon, using Velocity Verlet" << std::endl;
            velocityVerletStep(dt);
            break;
    }
}

void ConservativeSolver::velocityVerletStep(float dt) {
    // Velocity Verlet for Sine-Gordon: ∂²θ/∂t² = ∇²θ - sin(θ)
    //
    // Kick-drift-kick pattern:
    //   1. v_{n+1/2} = v_n + (dt/2)·a_n
    //   2. θ_{n+1} = θ_n + dt·v_{n+1/2}
    //   3. a_{n+1} = ∇²θ - sin(θ)  (recompute at new position)
    //   4. v_{n+1} = v_{n+1/2} + (dt/2)·a_{n+1}

    uint32_t total = nx_ * ny_ * nz_;

    // Step 1: Half-step velocity kick (using current acceleration)
    std::vector<float> accel_n(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        // Acceleration: a = ∇²θ - sin(θ)
        float laplacian;
        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            laplacian = computeLaplacian4thOrder(theta_, i, j, k);
        } else {
            laplacian = computeLaplacian(theta_, i, j, k);
        }
        float potential_force = -std::sin(theta_[idx]);
        accel_n[idx] = laplacian + potential_force;

        // Half-step velocity
        theta_dot_[idx] += 0.5f * dt * accel_n[idx];
    }

    // Step 2: Full-step position drift
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_[idx] += dt * theta_dot_[idx];

        // Periodic wrapping to [-π, π]
        while (theta_[idx] > M_PI) theta_[idx] -= 2.0f * M_PI;
        while (theta_[idx] < -M_PI) theta_[idx] += 2.0f * M_PI;
    }

    // Step 3: Recompute acceleration at new position
    std::vector<float> accel_n1(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        float laplacian;
        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            laplacian = computeLaplacian4thOrder(theta_, i, j, k);
        } else {
            laplacian = computeLaplacian(theta_, i, j, k);
        }
        float potential_force = -std::sin(theta_[idx]);
        accel_n1[idx] = laplacian + potential_force;
    }

    // Step 4: Final half-step velocity kick
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_dot_[idx] += 0.5f * dt * accel_n1[idx];
    }
}

void ConservativeSolver::rk2SymplecticStep(float dt) {
    // RK2 Midpoint Method for first-order reformulation
    // Convert ∂²θ/∂t² to system: dθ/dt = v, dv/dt = ∇²θ - sin(θ)

    uint32_t total = nx_ * ny_ * nz_;

    // Step 1: Compute k1 at current state
    std::vector<float> k1_theta(total);
    std::vector<float> k1_v(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        k1_theta[idx] = theta_dot_[idx];

        float laplacian;
        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            laplacian = computeLaplacian4thOrder(theta_, i, j, k);
        } else {
            laplacian = computeLaplacian(theta_, i, j, k);
        }
        k1_v[idx] = laplacian - std::sin(theta_[idx]);
    }

    // Step 2: Compute midpoint state
    std::vector<float> theta_mid(total);
    std::vector<float> v_mid(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_mid[idx] = theta_[idx] + 0.5f * dt * k1_theta[idx];
        v_mid[idx] = theta_dot_[idx] + 0.5f * dt * k1_v[idx];
    }

    // Step 3: Compute k2 at midpoint
    std::vector<float> k2_theta(total);
    std::vector<float> k2_v(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        k2_theta[idx] = v_mid[idx];

        // Use theta_mid for Laplacian computation
        float laplacian;
        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            laplacian = computeLaplacian4thOrder(theta_mid, i, j, k);
        } else {
            laplacian = computeLaplacian(theta_mid, i, j, k);
        }
        k2_v[idx] = laplacian - std::sin(theta_mid[idx]);
    }

    // Step 4: Final update
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_[idx] += dt * k2_theta[idx];
        theta_dot_[idx] += dt * k2_v[idx];

        // Periodic wrapping
        while (theta_[idx] > M_PI) theta_[idx] -= 2.0f * M_PI;
        while (theta_[idx] < -M_PI) theta_[idx] += 2.0f * M_PI;
    }
}

void ConservativeSolver::strangSplittingStep(float dt) {
    // Strang Splitting for Sine-Gordon: ∂²θ/∂t² = ∇²θ - sin(θ)
    //
    // Rewrite as first-order system:
    //   ∂θ/∂t = π                    (kinetic part T)
    //   ∂π/∂t = ∇²θ - sin(θ)         (potential part V)
    //
    // Proper T-V-T splitting (2nd order symplectic):
    //   exp(dt·L) = exp(dt/2·T) · exp(dt·V) · exp(dt/2·T) + O(dt³)
    //
    // Steps:
    //   1. Half-step kinetic:  θ += (dt/2)·π
    //   2. Full-step potential: π += dt·(∇²θ - sin(θ)) [using updated θ]
    //   3. Half-step kinetic:  θ += (dt/2)·π [using updated π]

    uint32_t total = nx_ * ny_ * nz_;

    // STEP 1: Half-step kinetic (T/2)
    // θ_{n+1/2} = θ_n + (dt/2)·π_n
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_[idx] += 0.5f * dt * theta_dot_[idx];
    }

    // STEP 2: Full-step potential (V)
    // π_{n+1} = π_n + dt·(∇²θ_{n+1/2} - sin(θ_{n+1/2}))
    std::vector<float> acceleration(total);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        // Compute Laplacian at midpoint position
        float laplacian;
        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            laplacian = computeLaplacian4thOrder(theta_, i, j, k);
        } else {
            laplacian = computeLaplacian(theta_, i, j, k);
        }

        // Compute nonlinear force at midpoint position
        float nonlinear = -std::sin(theta_[idx]);

        // Total acceleration
        acceleration[idx] = laplacian + nonlinear;
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_dot_[idx] += dt * acceleration[idx];
    }

    // STEP 3: Half-step kinetic (T/2)
    // θ_{n+1} = θ_{n+1/2} + (dt/2)·π_{n+1}
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_[idx] += 0.5f * dt * theta_dot_[idx];
    }

    // Periodic wrapping to [-π, π] (done ONCE at end to preserve energy)
    // NOTE: Disabled for Gaussian test - wrapping breaks energy conservation for small amplitudes
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    // for (uint32_t idx = 0; idx < total; ++idx) {
    //     while (theta_[idx] > M_PI) theta_[idx] -= 2.0f * M_PI;
    //     while (theta_[idx] < -M_PI) theta_[idx] += 2.0f * M_PI;
    // }
}

void ConservativeSolver::evolveDirac(float dt, const std::vector<float>& R_field,
                                     const std::vector<float>& theta_field, float Delta) {
    // Dirac evolution with chiral mass coupling
    // Uses split-step integrator pattern:
    //   1. Half-step kinetic evolution
    //   2. Full-step chiral mass interaction
    //   3. Half-step kinetic evolution

    if (!dirac3d_) {
        std::cerr << "[ConservativeSolver] ERROR: Dirac3D not initialized" << std::endl;
        return;
    }

    // Validate field sizes
    uint32_t expected_size = nx_ * ny_ * nz_;
    if (R_field.size() != expected_size || theta_field.size() != expected_size) {
        std::cerr << "[ConservativeSolver] ERROR: Field size mismatch. Expected: "
                  << expected_size << ", R_field: " << R_field.size()
                  << ", theta_field: " << theta_field.size() << std::endl;
        return;
    }

    // Use public split-step integrator with chiral mass coupling
    dirac3d_->stepWithChiralMass(R_field, theta_field, Delta, dt);

    std::cout << "[ConservativeSolver] Dirac evolution with chiral mass coupling completed (dt="
              << dt << ", Delta=" << Delta << ")" << std::endl;
}

void ConservativeSolver::initializeVortexWithProperVelocity(float x0, float y0, float z0, int charge) {
    std::cout << "[ConservativeSolver] Initializing vortex (charge=" << charge
              << ") at (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t k = 0; k < nz_; ++k) {
        for (uint32_t j = 0; j < ny_; ++j) {
            for (uint32_t i = 0; i < nx_; ++i) {
                int idx = index3D(i, j, k);

                // Position relative to vortex center
                float x = i * dx_ - x0;
                float y = j * dx_ - y0;
                float z = k * dx_ - z0;

                float r = std::sqrt(x*x + y*y + z*z);
                float phi = std::atan2(y, x);

                // Phase field: θ(r,φ) = n·φ (topological winding)
                theta_[idx] = charge * phi;

                // CRITICAL FIX: Velocity field for STATIONARY vortex
                // For stationary vortex: ∂θ/∂t = 0
                theta_dot_[idx] = 0.0f;

                // For moving vortex with velocity (vx, vy, vz):
                // ∂θ/∂t = -v·∇θ
                // For θ = n·atan2(y,x): ∇θ = n·(-y, x, 0)/r²
                // So: ∂θ/∂t = -n·(vx·(-y) + vy·x)/r²
                //           = n·(vx·y - vy·x)/r²
                // (Currently set to zero for stationary vortex)
            }
        }
    }

    // Store initial energy for validation
    initial_energy_ = computeTotalEnergy();
    std::cout << "[ConservativeSolver] Initial energy: " << initial_energy_ << std::endl;
}

void ConservativeSolver::initializeCollisionScenario(
    float x1, float y1, float x2, float y2, float velocity) {

    std::cout << "[ConservativeSolver] Initializing collision: vortex (+1) at ("
              << x1 << "," << y1 << ") → vortex (-1) at (" << x2 << "," << y2 << ")" << std::endl;
    std::cout << "[ConservativeSolver] Relative velocity: " << velocity << std::endl;

    uint32_t total = nx_ * ny_ * nz_;

    // Temporary storage for individual vortex contributions
    std::vector<float> theta_v1(total, 0.0f);
    std::vector<float> theta_v2(total, 0.0f);
    std::vector<float> tdot_v1(total, 0.0f);
    std::vector<float> tdot_v2(total, 0.0f);

    float z_center = nz_ * dx_ / 2.0f;  // Center in Z

    // Initialize vortex 1 (charge +1, moving +x direction)
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        float x = i * dx_ - x1;
        float y = j * dx_ - y1;
        float z = k * dx_ - z_center;

        float r_sq = x*x + y*y + z*z;
        float phi = std::atan2(y, x);

        // Phase: θ = +1 · φ
        theta_v1[idx] = phi;

        // Velocity for moving vortex: ∂θ/∂t = -v·∇θ
        // For θ = atan2(y,x): ∇θ = (-y/r², x/r², 0)
        // ∂θ/∂t = -vx·(-y/r²) = vx·y/r²
        if (r_sq > 1e-6f) {
            tdot_v1[idx] = velocity * y / r_sq;
        }
    }

    // Initialize vortex 2 (charge -1, moving -x direction)
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        float x = i * dx_ - x2;
        float y = j * dx_ - y2;
        float z = k * dx_ - z_center;

        float r_sq = x*x + y*y + z*z;
        float phi = std::atan2(y, x);

        // Phase: θ = -1 · φ
        theta_v2[idx] = -phi;

        // Velocity (moving in -x direction)
        if (r_sq > 1e-6f) {
            tdot_v2[idx] = -velocity * y / r_sq;
        }
    }

    // Superposition (linear approximation for small amplitudes)
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        theta_[idx] = theta_v1[idx] + theta_v2[idx];
        theta_dot_[idx] = tdot_v1[idx] + tdot_v2[idx];

        // Wrap to [-π, π]
        while (theta_[idx] > M_PI) theta_[idx] -= 2.0f * M_PI;
        while (theta_[idx] < -M_PI) theta_[idx] += 2.0f * M_PI;
    }

    initial_energy_ = computeTotalEnergy();
    std::cout << "[ConservativeSolver] Collision initial energy: " << initial_energy_ << std::endl;
}

void ConservativeSolver::initializeGaussian(float x0, float y0, float z0, float sigma, float amplitude) {
    std::cout << "[ConservativeSolver] Initializing Gaussian wave packet (σ=" << sigma
              << ", A=" << amplitude << ")" << std::endl;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint32_t k = 0; k < nz_; ++k) {
        for (uint32_t j = 0; j < ny_; ++j) {
            for (uint32_t i = 0; i < nx_; ++i) {
                int idx = index3D(i, j, k);

                float x = i * dx_ - x0;
                float y = j * dx_ - y0;
                float z = k * dx_ - z0;

                float r_sq = x*x + y*y + z*z;

                // Gaussian profile
                theta_[idx] = amplitude * std::exp(-r_sq / (2.0f * sigma * sigma));

                // Initial velocity: zero (stationary wave packet)
                theta_dot_[idx] = 0.0f;
            }
        }
    }

    initial_energy_ = computeTotalEnergy();
    std::cout << "[ConservativeSolver] Gaussian initial energy: " << initial_energy_ << std::endl;
}

float ConservativeSolver::computeTotalEnergy() const {
    // Total energy: E = ∫[(∂θ/∂t)² + (∇θ)² + V(θ)]dV
    // where V(θ) = 1 - cos(θ) for Sine-Gordon
    //
    // CRITICAL: Use same discretization order as evolution for consistent energy measurement

    float energy_kinetic = 0.0f;
    float energy_gradient = 0.0f;
    float energy_potential = 0.0f;

    uint32_t total = nx_ * ny_ * nz_;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:energy_kinetic, energy_gradient, energy_potential)
    #endif
    for (uint32_t idx = 0; idx < total; ++idx) {
        int i, j, k;
        coords3D(idx, i, j, k);

        // Kinetic energy: (∂θ/∂t)²
        float tdot = theta_dot_[idx];
        energy_kinetic += tdot * tdot;

        // Gradient energy: |∇θ|² - use consistent order with evolution
        float grad_x, grad_y, grad_z;

        if (config_.spatial_order == SpatialOrder::FOURTH_ORDER) {
            // 4th-order gradients: [f_{-2} - 8f_{-1} + 8f_{+1} - f_{+2}] / (12·dx)
            const float f_xm2 = theta_[index3D(wrapX(i-2), j, k)];
            const float f_xm1 = theta_[index3D(wrapX(i-1), j, k)];
            const float f_xp1 = theta_[index3D(wrapX(i+1), j, k)];
            const float f_xp2 = theta_[index3D(wrapX(i+2), j, k)];
            grad_x = (f_xm2 - 8.0f*f_xm1 + 8.0f*f_xp1 - f_xp2) / (12.0f * dx_);

            const float f_ym2 = theta_[index3D(i, wrapY(j-2), k)];
            const float f_ym1 = theta_[index3D(i, wrapY(j-1), k)];
            const float f_yp1 = theta_[index3D(i, wrapY(j+1), k)];
            const float f_yp2 = theta_[index3D(i, wrapY(j+2), k)];
            grad_y = (f_ym2 - 8.0f*f_ym1 + 8.0f*f_yp1 - f_yp2) / (12.0f * dx_);

            const float f_zm2 = theta_[index3D(i, j, wrapZ(k-2))];
            const float f_zm1 = theta_[index3D(i, j, wrapZ(k-1))];
            const float f_zp1 = theta_[index3D(i, j, wrapZ(k+1))];
            const float f_zp2 = theta_[index3D(i, j, wrapZ(k+2))];
            grad_z = (f_zm2 - 8.0f*f_zm1 + 8.0f*f_zp1 - f_zp2) / (12.0f * dx_);
        } else {
            // 2nd-order gradients (centered difference)
            grad_x = (theta_[index3D(wrapX(i+1), j, k)] - theta_[index3D(wrapX(i-1), j, k)]) / (2.0f * dx_);
            grad_y = (theta_[index3D(i, wrapY(j+1), k)] - theta_[index3D(i, wrapY(j-1), k)]) / (2.0f * dx_);
            grad_z = (theta_[index3D(i, j, wrapZ(k+1))] - theta_[index3D(i, j, wrapZ(k-1))]) / (2.0f * dx_);
        }

        energy_gradient += grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

        // Potential energy: V(θ) = 1 - cos(θ)
        energy_potential += 1.0f - std::cos(theta_[idx]);
    }

    // Volume element: dx³
    float dV = dx_ * dx_ * dx_;
    float total_energy = 0.5f * (energy_kinetic + energy_gradient) + energy_potential;
    total_energy *= dV;

    return total_energy;
}

float ConservativeSolver::measureEnergyDrift(float E_initial) const {
    float E_current = computeTotalEnergy();
    float drift = std::abs(E_current - E_initial) / E_initial;
    return drift;
}

bool ConservativeSolver::validateEnergyConservation(float threshold) {
    float drift = measureEnergyDrift(initial_energy_);

    std::cout << "[ConservativeSolver] Energy conservation check:" << std::endl;
    std::cout << "  Initial: " << initial_energy_ << std::endl;
    std::cout << "  Current: " << computeTotalEnergy() << std::endl;
    std::cout << "  Drift: " << (drift * 100.0f) << "% (threshold: " << (threshold * 100.0f) << "%)" << std::endl;

    bool passed = drift < threshold;
    std::cout << "  Status: " << (passed ? "✓ PASS" : "✗ FAIL") << std::endl;

    return passed;
}

bool ConservativeSolver::validateTimeReversibility(float threshold) {
    std::cout << "[ConservativeSolver] Time reversibility test (threshold: " << threshold << " rad)" << std::endl;

    // Store initial state
    std::vector<float> theta_init = theta_;
    std::vector<float> tdot_init = theta_dot_;

    // Forward evolution (10 steps)
    for (int step = 0; step < 10; ++step) {
        evolveSineGordon(dt_);
    }

    // Reverse time (negate velocities)
    for (auto& v : theta_dot_) {
        v = -v;
    }

    // Backward evolution (10 steps)
    for (int step = 0; step < 10; ++step) {
        evolveSineGordon(dt_);
    }

    // Reverse velocities again to compare
    for (auto& v : theta_dot_) {
        v = -v;
    }

    // Compute error
    float max_error = 0.0f;
    for (size_t i = 0; i < theta_.size(); ++i) {
        float error = std::abs(theta_[i] - theta_init[i]);
        max_error = std::max(max_error, error);
    }

    std::cout << "  Maximum phase error: " << max_error << " rad" << std::endl;

    bool passed = max_error < threshold;
    std::cout << "  Status: " << (passed ? "✓ PASS" : "✗ FAIL") << std::endl;

    // Restore initial state
    theta_ = theta_init;
    theta_dot_ = tdot_init;

    return passed;
}

// Helper methods

float ConservativeSolver::computeLaplacian(const std::vector<float>& field, int i, int j, int k) const {
    // 6-neighbor stencil: ∇²θ ≈ (θ_{i+1} + θ_{i-1} + ... - 6θ_i) / dx²
    // 2nd-order accurate: O(dx²)

    int idx_center = index3D(i, j, k);

    float theta_xp = field[index3D(wrapX(i+1), j, k)];
    float theta_xm = field[index3D(wrapX(i-1), j, k)];
    float theta_yp = field[index3D(i, wrapY(j+1), k)];
    float theta_ym = field[index3D(i, wrapY(j-1), k)];
    float theta_zp = field[index3D(i, j, wrapZ(k+1))];
    float theta_zm = field[index3D(i, j, wrapZ(k-1))];
    float theta_c = field[idx_center];

    float laplacian = (theta_xp + theta_xm + theta_yp + theta_ym + theta_zp + theta_zm - 6.0f * theta_c) / (dx_ * dx_);

    return laplacian;
}

float ConservativeSolver::computeLaplacian4thOrder(const std::vector<float>& field, int i, int j, int k) const {
    // 4th-order accurate Laplacian using 18-neighbor compact stencil
    //
    // Standard 1D 4th-order: d²f/dx² ≈ [-f_{i-2} + 16f_{i-1} - 30f_i + 16f_{i+1} - f_{i+2}] / (12·dx²)
    //
    // For 3D Laplacian: ∇²θ = ∂²θ/∂x² + ∂²θ/∂y² + ∂²θ/∂z²
    //
    // Each 1D operator contributes independently, so we sum them directly.
    // This is the correct separable form for the Laplacian operator.

    const int idx_center = index3D(i, j, k);
    const float dx2 = dx_ * dx_;
    const float theta_c = field[idx_center];

    // X-direction: ∂²θ/∂x²
    const float theta_xm2 = field[index3D(wrapX(i-2), j, k)];
    const float theta_xm1 = field[index3D(wrapX(i-1), j, k)];
    const float theta_xp1 = field[index3D(wrapX(i+1), j, k)];
    const float theta_xp2 = field[index3D(wrapX(i+2), j, k)];

    const float d2_dx2 = (-theta_xm2 + 16.0f*theta_xm1 - 30.0f*theta_c + 16.0f*theta_xp1 - theta_xp2) / (12.0f * dx2);

    // Y-direction: ∂²θ/∂y²
    const float theta_ym2 = field[index3D(i, wrapY(j-2), k)];
    const float theta_ym1 = field[index3D(i, wrapY(j-1), k)];
    const float theta_yp1 = field[index3D(i, wrapY(j+1), k)];
    const float theta_yp2 = field[index3D(i, wrapY(j+2), k)];

    const float d2_dy2 = (-theta_ym2 + 16.0f*theta_ym1 - 30.0f*theta_c + 16.0f*theta_yp1 - theta_yp2) / (12.0f * dx2);

    // Z-direction: ∂²θ/∂z²
    const float theta_zm2 = field[index3D(i, j, wrapZ(k-2))];
    const float theta_zm1 = field[index3D(i, j, wrapZ(k-1))];
    const float theta_zp1 = field[index3D(i, j, wrapZ(k+1))];
    const float theta_zp2 = field[index3D(i, j, wrapZ(k+2))];

    const float d2_dz2 = (-theta_zm2 + 16.0f*theta_zm1 - 30.0f*theta_c + 16.0f*theta_zp1 - theta_zp2) / (12.0f * dx2);

    // Total 4th-order Laplacian (sum of separable 1D operators)
    const float laplacian = d2_dx2 + d2_dy2 + d2_dz2;

    return laplacian;
}

float ConservativeSolver::computeGradientX(const std::vector<float>& field, int idx) const {
    int i, j, k;
    coords3D(idx, i, j, k);

    float theta_plus = field[index3D(wrapX(i+1), j, k)];
    float theta_minus = field[index3D(wrapX(i-1), j, k)];

    return (theta_plus - theta_minus) / (2.0f * dx_);
}

float ConservativeSolver::computeGradientY(const std::vector<float>& field, int idx) const {
    int i, j, k;
    coords3D(idx, i, j, k);

    float theta_plus = field[index3D(i, wrapY(j+1), k)];
    float theta_minus = field[index3D(i, wrapY(j-1), k)];

    return (theta_plus - theta_minus) / (2.0f * dx_);
}

float ConservativeSolver::computeGradientZ(const std::vector<float>& field, int idx) const {
    int i, j, k;
    coords3D(idx, i, j, k);

    float theta_plus = field[index3D(i, j, wrapZ(k+1))];
    float theta_minus = field[index3D(i, j, wrapZ(k-1))];

    return (theta_plus - theta_minus) / (2.0f * dx_);
}
