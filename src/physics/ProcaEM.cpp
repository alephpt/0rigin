// src/physics/ProcaEM.cpp
#include "physics/ProcaEM.h"
#include <cmath>
#include <algorithm>

namespace physics {

ProcaEM::ProcaEM(int nx, int ny, float dx, float photon_mass_coupling, float alpha_coupling)
    : nx_(nx), ny_(ny), dx_(dx), photon_mass_coupling_(photon_mass_coupling), alpha_coupling_(alpha_coupling) {

    int size = nx * ny;
    phi_.resize(size, 0.0f);
    A_x_.resize(size, 0.0f);
    A_y_.resize(size, 0.0f);
    A_z_.resize(size, 0.0f);

    field_tensor_.resize(size);

    j_t_.resize(size, 0.0f);
    j_x_.resize(size, 0.0f);
    j_y_.resize(size, 0.0f);
}

void ProcaEM::computePotentials(const float* theta_field, const float* R_field,
                                int nx, int ny, float dx, float dt) {
    (void)nx; (void)ny; (void)dx; // Used internally by member variables

    // Compute current from SMFT fields
    computeCurrent(theta_field, R_field);

    // Evolve Proca equation: (□ + m²)A_μ = j_μ
    evolveProcaField(dt);

    // Enforce Lorenz gauge
    applyLorenzGauge();
}

void ProcaEM::computeCurrent(const float* theta_field, const float* R_field) {
    // Noether current: j_μ = α ∂_μθ
    // alpha_coupling_ is now configurable

    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // ∂_x θ
            float dtheta_dx = (theta_field[idx + 1] - theta_field[idx - 1]) / (2.0f * dx_);
            // ∂_y θ
            float dtheta_dy = (theta_field[idx + nx_] - theta_field[idx - nx_]) / (2.0f * dx_);

            j_x_[idx] = alpha_coupling_ * dtheta_dx * R_field[idx];
            j_y_[idx] = alpha_coupling_ * dtheta_dy * R_field[idx];
            j_t_[idx] = 0.0f; // No charge density (current conserved)
        }
    }
}

void ProcaEM::evolveProcaField(float dt) {
    // Simplified: Proca equation (□ + m²)A_μ = j_μ
    // Using finite differences and forward Euler (placeholder - use better integrator)

    std::vector<float> phi_new = phi_;
    std::vector<float> A_x_new = A_x_;
    std::vector<float> A_y_new = A_y_;

    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            float m_gamma_sq = getPhotonMass(i, j);
            m_gamma_sq = m_gamma_sq * m_gamma_sq;

            // Laplacian of A_x
            float lap_Ax = (A_x_[idx+1] + A_x_[idx-1] + A_x_[idx+nx_] + A_x_[idx-nx_] - 4*A_x_[idx])/(dx_*dx_);

            // Proca: ∂²A/∂t² - ∇²A + m²A = j
            // Simplified to: A_new = A_old + dt²(∇²A - m²A + j)
            A_x_new[idx] = A_x_[idx] + dt * dt * (lap_Ax - m_gamma_sq * A_x_[idx] + j_x_[idx]);

            // Same for A_y
            float lap_Ay = (A_y_[idx+1] + A_y_[idx-1] + A_y_[idx+nx_] + A_y_[idx-nx_] - 4*A_y_[idx])/(dx_*dx_);
            A_y_new[idx] = A_y_[idx] + dt * dt * (lap_Ay - m_gamma_sq * A_y_[idx] + j_y_[idx]);
        }
    }

    phi_ = phi_new;
    A_x_ = A_x_new;
    A_y_ = A_y_new;
}

void ProcaEM::computeFieldStrengths() {
    // E = -∇φ - ∂A/∂t (∂A/∂t ≈ 0 in steady state)
    // B = ∇×A

    for (int j = 1; j < ny_ - 1; ++j) {
        for (int i = 1; i < nx_ - 1; ++i) {
            int idx = j * nx_ + i;

            // Electric field
            float dphi_dx = (phi_[idx + 1] - phi_[idx - 1]) / (2.0f * dx_);
            float dphi_dy = (phi_[idx + nx_] - phi_[idx - nx_]) / (2.0f * dx_);

            field_tensor_[idx].Ex = -dphi_dx;
            field_tensor_[idx].Ey = -dphi_dy;
            field_tensor_[idx].Ez = 0.0f;

            // Magnetic field: B_z = ∂A_y/∂x - ∂A_x/∂y
            float dAy_dx = (A_y_[idx + 1] - A_y_[idx - 1]) / (2.0f * dx_);
            float dAx_dy = (A_x_[idx + nx_] - A_x_[idx - nx_]) / (2.0f * dx_);

            field_tensor_[idx].Bx = 0.0f;
            field_tensor_[idx].By = 0.0f;
            field_tensor_[idx].Bz = dAy_dx - dAx_dy;
        }
    }
}

FieldTensor ProcaEM::getFieldAt(int i, int j) const {
    return field_tensor_[j * nx_ + i];
}

float ProcaEM::getPhotonMass(int i, int j) const {
    (void)i; (void)j; // Would use R_field[j*nx_+i] for spatially-varying mass
    // Placeholder: m_γ = g (would use R_field here)
    return photon_mass_coupling_;
}

float ProcaEM::computeFieldEnergy() const {
    float energy = 0.0f;
    for (const auto& F : field_tensor_) {
        float E_sq = F.Ex*F.Ex + F.Ey*F.Ey + F.Ez*F.Ez;
        float B_sq = F.Bx*F.Bx + F.By*F.By + F.Bz*F.Bz;
        energy += (E_sq + B_sq) * dx_ * dx_; // ∫(E²+B²)/2 dV
    }
    return energy * 0.5f;
}

void ProcaEM::applyLorenzGauge() {
    // Lorenz gauge: ∂_μ A^μ = 0
    // Simplified projection (placeholder)
}

} // namespace physics
