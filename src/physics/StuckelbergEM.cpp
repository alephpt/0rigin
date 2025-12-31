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
    // Using simple explicit scheme: A(t+dt) = A(t) + dt² * ∇²A
    std::vector<float> A_x_new = A_x_;
    std::vector<float> A_y_new = A_y_;

    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // 5-point Laplacian
            float lap_Ax = (A_x_[idx+1] + A_x_[idx-1] + A_x_[idx+nx_] + A_x_[idx-nx_] - 4*A_x_[idx])/(dx_*dx_);
            float lap_Ay = (A_y_[idx+1] + A_y_[idx-1] + A_y_[idx+nx_] + A_y_[idx-nx_] - 4*A_y_[idx])/(dx_*dx_);

            // Wave equation: ∂²A/∂t² = ∇²A
            A_x_new[idx] = A_x_[idx] + dt*dt * lap_Ax;
            A_y_new[idx] = A_y_[idx] + dt*dt * lap_Ay;
        }
    }

    A_x_ = A_x_new;
    A_y_ = A_y_new;
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
