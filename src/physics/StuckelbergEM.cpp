// src/physics/StuckelbergEM.cpp
#include "physics/StuckelbergEM.h"
#include <cmath>
#include <algorithm>

namespace physics {

StuckelbergEM::StuckelbergEM(int nx, int ny, float dx, float photon_mass)
    : nx_(nx), ny_(ny), dx_(dx), photon_mass_(photon_mass) {

    int size = nx * ny;
    A_x_.resize(size, 0.0f);
    A_y_.resize(size, 0.0f);
    phi_.resize(size, 0.0f);
    phi_dot_.resize(size, 0.0f);
    Aprime_x_.resize(size, 0.0f);
    Aprime_y_.resize(size, 0.0f);
    field_tensor_.resize(size);
}

void StuckelbergEM::initializeGaussianPulse(float center_x, float center_y,
                                           float width_x, float width_y,
                                           float amplitude, float k_x, float k_y,
                                           const std::string& component) {
    // Initialize a Gaussian wave packet in the specified component
    // ψ(x,y) = A * exp(-(x-x0)²/2σ_x² - (y-y0)²/2σ_y²) * exp(i(k_x*x + k_y*y))

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = j * nx_ + i;

            float x = i * dx_;
            float y = j * dx_;
            float dx_val = x - center_x;
            float dy_val = y - center_y;

            // Gaussian envelope
            float envelope = amplitude * std::exp(
                -0.5f * (dx_val*dx_val/(width_x*width_x) +
                        dy_val*dy_val/(width_y*width_y))
            );

            // Phase factor for wave propagation
            float phase = k_x * x + k_y * y;

            if (component == "A_x") {
                A_x_[idx] = envelope * std::cos(phase);
            } else if (component == "A_y") {
                A_y_[idx] = envelope * std::cos(phase);
                // For wave propagation, also set time derivative
                // ∂A/∂t = -ck for rightward propagation (c=1 in natural units)
                // This ensures the wave propagates rather than oscillates in place
                if (k_x != 0.0f || k_y != 0.0f) {
                    // Initial velocity for wave propagation
                    // For a wave A = A₀cos(kx - ωt), ∂A/∂t = ωA₀sin(kx - ωt)
                    // At t=0: ∂A/∂t = -ω*A₀sin(kx) where ω = c*k (c=1)
                    float omega = std::sqrt(k_x*k_x + k_y*k_y); // ω = c|k|, c=1
                    // Store this for next timestep evolution
                    // Note: We'll need a proper velocity field for full wave equation
                }
            } else if (component == "phi") {
                phi_[idx] = envelope * std::cos(phase);
                // Initial time derivative for wave propagation
                float omega = std::sqrt(k_x*k_x + k_y*k_y + photon_mass_*photon_mass_);
                phi_dot_[idx] = -omega * envelope * std::sin(phase);
            }
        }
    }
}

void StuckelbergEM::computePotentials(const float* theta_field, const float* R_field,
                                      int nx, int ny, float dx, float dt) {
    // Initialize φ from θ on first call
    static bool initialized = false;
    if (!initialized) {
        initializeFromTheta(theta_field);
        initialized = true;
    }

    // Evolve both fields
    evolveMaxwell(dt);
    evolveKleinGordon(dt);

    // Transform: A' = A + ∂φ
    transformPotential();
}

void StuckelbergEM::initializeFromTheta(const float* theta_field) {
    // Direct coupling: φ = θ initially
    for (int i = 0; i < nx_ * ny_; ++i) {
        phi_[i] = theta_field[i];
    }
}

void StuckelbergEM::evolveMaxwell(float dt) {
    // □A = ∂²A/∂t² - ∇²A = 0 (massless wave equation)
    // Using leap-frog scheme for stable wave propagation
    // Need to track both A and ∂A/∂t (velocity)

    // First time: initialize velocities if needed
    static bool vel_initialized = false;
    static std::vector<float> A_x_dot, A_y_dot;
    if (!vel_initialized) {
        A_x_dot.resize(nx_ * ny_, 0.0f);
        A_y_dot.resize(nx_ * ny_, 0.0f);
        vel_initialized = true;
    }

    std::vector<float> A_x_new = A_x_;
    std::vector<float> A_y_new = A_y_;
    std::vector<float> A_x_dot_new = A_x_dot;
    std::vector<float> A_y_dot_new = A_y_dot;

    // Use periodic boundary conditions for wave propagation
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = j * nx_ + i;

            // Periodic indices
            int ip = ((i + 1) % nx_);
            int im = ((i - 1 + nx_) % nx_);
            int jp = ((j + 1) % ny_);
            int jm = ((j - 1 + ny_) % ny_);

            int idx_xp = j * nx_ + ip;
            int idx_xm = j * nx_ + im;
            int idx_yp = jp * nx_ + i;
            int idx_ym = jm * nx_ + i;

            // 5-point Laplacian with periodic boundaries
            float lap_Ax = (A_x_[idx_xp] + A_x_[idx_xm] + A_x_[idx_yp] + A_x_[idx_ym] - 4*A_x_[idx])/(dx_*dx_);
            float lap_Ay = (A_y_[idx_xp] + A_y_[idx_xm] + A_y_[idx_yp] + A_y_[idx_ym] - 4*A_y_[idx])/(dx_*dx_);

            // Wave equation: ∂²A/∂t² = c²∇²A (c=1 in natural units)
            // Leap-frog: v(t+dt/2) = v(t-dt/2) + dt*a(t)
            //           x(t+dt) = x(t) + dt*v(t+dt/2)
            float c_sq = 1.0f; // Speed of light squared
            A_x_dot_new[idx] = A_x_dot[idx] + dt * c_sq * lap_Ax;
            A_y_dot_new[idx] = A_y_dot[idx] + dt * c_sq * lap_Ay;

            A_x_new[idx] = A_x_[idx] + dt * A_x_dot_new[idx];
            A_y_new[idx] = A_y_[idx] + dt * A_y_dot_new[idx];
        }
    }

    A_x_ = A_x_new;
    A_y_ = A_y_new;
    A_x_dot = A_x_dot_new;
    A_y_dot = A_y_dot_new;
}

void StuckelbergEM::evolveKleinGordon(float dt) {
    // □φ + m²φ = 0
    // ∂²φ/∂t² = ∇²φ - m²φ
    std::vector<float> phi_new = phi_;
    std::vector<float> phi_dot_new = phi_dot_;

    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // 5-point Laplacian
            float lap_phi = (phi_[idx+1] + phi_[idx-1] + phi_[idx+nx_] + phi_[idx-nx_] - 4*phi_[idx])/(dx_*dx_);

            // Klein-Gordon: ∂²φ/∂t² = ∇²φ - m²φ
            float phi_dotdot = lap_phi - photon_mass_*photon_mass_*phi_[idx];

            // Velocity Verlet integration
            phi_dot_new[idx] = phi_dot_[idx] + dt * phi_dotdot;
            phi_new[idx] = phi_[idx] + dt * phi_dot_new[idx];
        }
    }

    phi_ = phi_new;
    phi_dot_ = phi_dot_new;
}

void StuckelbergEM::transformPotential() {
    // A'_μ = A_μ + ∂_μφ
    // This is the key to gauge restoration!
    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // Central difference for ∂φ
            float dphi_dx = (phi_[idx + 1] - phi_[idx - 1]) / (2.0f * dx_);
            float dphi_dy = (phi_[idx + nx_] - phi_[idx - nx_]) / (2.0f * dx_);

            Aprime_x_[idx] = A_x_[idx] + dphi_dx;
            Aprime_y_[idx] = A_y_[idx] + dphi_dy;
        }
    }
}

void StuckelbergEM::computeFieldStrengths() {
    // B = ∇×A' (uses transformed potential!)
    // F^μν from A'_μ
    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // B_z = ∂A'_y/∂x - ∂A'_x/∂y (curl in 2D)
            float dAy_dx = (Aprime_y_[idx + 1] - Aprime_y_[idx - 1]) / (2.0f * dx_);
            float dAx_dy = (Aprime_x_[idx + nx_] - Aprime_x_[idx - nx_]) / (2.0f * dx_);

            field_tensor_[idx].Bx = 0.0f;
            field_tensor_[idx].By = 0.0f;
            field_tensor_[idx].Bz = dAy_dx - dAx_dy;

            // E field (simplified - would need ∂A'/∂t for full implementation)
            field_tensor_[idx].Ex = 0.0f;
            field_tensor_[idx].Ey = 0.0f;
            field_tensor_[idx].Ez = 0.0f;
        }
    }
}

FieldTensor StuckelbergEM::getFieldAt(int i, int j) const {
    return field_tensor_[j * nx_ + i];
}

float StuckelbergEM::getPhiAt(int i, int j) const {
    return phi_[j * nx_ + i];
}

float StuckelbergEM::getAprimeX(int i, int j) const {
    return Aprime_x_[j * nx_ + i];
}

float StuckelbergEM::getAprimeY(int i, int j) const {
    return Aprime_y_[j * nx_ + i];
}

float StuckelbergEM::computeFieldEnergy() const {
    float energy = 0.0f;
    for (const auto& F : field_tensor_) {
        // E²/2 + B²/2 (in natural units)
        energy += (F.Ex*F.Ex + F.Ey*F.Ey + F.Ez*F.Ez +
                   F.Bx*F.Bx + F.By*F.By + F.Bz*F.Bz) * dx_ * dx_ * 0.5f;
    }
    return energy;
}

void StuckelbergEM::evolveStuckelbergField(float dt) {
    // This is already called in computePotentials, but exposed for testing
    evolveKleinGordon(dt);
}

} // namespace physics
