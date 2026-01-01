// include/physics/StuckelbergEM.h
#pragma once
#include "GaugeTheory.h"
#include <vector>

namespace physics {

/**
 * Stückelberg Gauge-Restored EM Implementation
 *
 * Key difference from Proca:
 * - Proca: A_μ → B via current j_μ = ∂_μθ (FAILED - too weak, B~10⁻⁹)
 * - Stückelberg: A'_μ = A_μ + ∂_μφ/e, where φ = θ (DIRECT coupling)
 *
 * Gauge transformation:
 *   Under θ → θ + eα:
 *     A_μ → A_μ + ∂_μα
 *     φ → φ - eα
 *     A'_μ unchanged (gauge invariant!)
 *
 * Evolution equations:
 *   □A_μ = 0 (massless Maxwell)
 *   □φ = -m²φ (Klein-Gordon with mass)
 *   B = ∇×(A + ∂φ) (field from BOTH A and φ)
 */
class StuckelbergEM : public GaugeTheory {
public:
    StuckelbergEM(int nx, int ny, float dx, float photon_mass);

    void computePotentials(const float* theta_field, const float* R_field,
                          int nx, int ny, float dx, float dt) override;
    void computeFieldStrengths() override;
    FieldTensor getFieldAt(int i, int j) const override;

    bool isGaugeInvariant() const override { return true; } // Restored!
    Mechanism getMechanism() const override { return Mechanism::STUCKELBERG; }
    std::string getName() const override { return "Stückelberg Gauge-Restored"; }

    float computeFieldEnergy() const override;

    // Stückelberg-specific accessors
    void evolveStuckelbergField(float dt);
    float getPhiAt(int i, int j) const;
    float getAprimeX(int i, int j) const;
    float getAprimeY(int i, int j) const;

    // Wave propagation initialization
    void initializeGaussianPulse(float center_x, float center_y,
                                float width_x, float width_y,
                                float amplitude, float k_x, float k_y,
                                const std::string& component);

private:
    int nx_, ny_;
    float dx_;
    float photon_mass_;

    // A_μ (massless Maxwell)
    std::vector<float> A_x_, A_y_;

    // φ (Stückelberg scalar) - couples to θ
    std::vector<float> phi_;
    std::vector<float> phi_dot_; // Time derivative

    // Transformed potential: A'_μ = A_μ + ∂_μφ
    std::vector<float> Aprime_x_, Aprime_y_;

    // Field strengths from A'
    std::vector<FieldTensor> field_tensor_;

    void initializeFromTheta(const float* theta_field);
    void evolveMaxwell(float dt); // □A = 0
    void evolveKleinGordon(float dt); // □φ + m²φ = 0
    void transformPotential(); // A' = A + ∂φ
};

} // namespace physics
