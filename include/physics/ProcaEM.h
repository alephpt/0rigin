// include/physics/ProcaEM.h
#pragma once
#include "GaugeTheory.h"
#include <vector>

namespace physics {

class ProcaEM : public GaugeTheory {
public:
    ProcaEM(int nx, int ny, float dx, float photon_mass_coupling, float alpha_coupling = 0.1f);
    ~ProcaEM() override = default;

    // GaugeTheory interface
    void computePotentials(const float* theta_field, const float* R_field,
                          int nx, int ny, float dx, float dt) override;
    void computeFieldStrengths() override;
    FieldTensor getFieldAt(int i, int j) const override;

    bool isGaugeInvariant() const override { return false; } // Proca breaks U(1)
    Mechanism getMechanism() const override { return Mechanism::PROCA; }
    std::string getName() const override { return "Proca Massive Photon"; }

    float computeFieldEnergy() const override;

    // Proca-specific
    void evolveProcaField(float dt);
    float getPhotonMass(int i, int j) const;

private:
    int nx_, ny_;
    float dx_;
    float photon_mass_coupling_;  // g in m_γ = g(1-R)
    float alpha_coupling_;        // α in j_μ = α ∂_μθ

    // Field storage: A_μ = (φ, A_x, A_y, A_z)
    std::vector<float> phi_;      // Scalar potential
    std::vector<float> A_x_;      // Vector potential x
    std::vector<float> A_y_;      // Vector potential y
    std::vector<float> A_z_;      // Vector potential z (=0 in 2D)

    // Field strengths
    std::vector<FieldTensor> field_tensor_;

    // Current source j_μ from SMFT
    std::vector<float> j_t_;      // Charge density
    std::vector<float> j_x_;      // Current density x
    std::vector<float> j_y_;      // Current density y

    void computeCurrent(const float* theta_field, const float* R_field);
    void applyLorenzGauge();
};

} // namespace physics
