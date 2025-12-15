#pragma once

#include <cmath>
#include <algorithm>

namespace Nova {
namespace MSFT {

struct MSFTParams {
    float dt = 0.01f;
    float K = 1.0f;
    float damping = 0.1f;
    float Delta = 1.0f;
    float chiral_angle = 0.0f;
    uint32_t Nx = 64;
    uint32_t Ny = 64;
    uint32_t N_total = 4096;
    uint32_t enable_feedback = 1;
    uint32_t enable_normalization = 1;
    uint32_t neighborhood_radius = 1;
    uint32_t normalize_density = 1;
    uint32_t dirac_substeps = 1;

    float compute_safe_dt(float dx = 1.0f, float R_max = 1.0f) const {
        // Dirac CFL condition (accounts for mass gap)
        float m_max = Delta * R_max;
        float dt_max_dirac = (m_max > 1e-6f) ? (0.5f * dx / m_max) : 0.1f;

        // Kuramoto stability condition
        float dt_max_kuramoto = (K > 1e-6f) ? (0.1f / K) : 0.1f;

        // Return most restrictive constraint
        return std::min(dt_max_dirac, dt_max_kuramoto);
    }

    bool validate_dt(float dx = 1.0f, float R_max = 1.0f) const {
        float dt_safe = compute_safe_dt(dx, R_max);
        return dt <= dt_safe;
    }

    void auto_adjust_dt(float dx = 1.0f, float R_max = 1.0f) {
        float dt_safe = compute_safe_dt(dx, R_max);
        if (dt > dt_safe) {
            dt = dt_safe * 0.9f;  // Add 10% safety margin
        }
    }

    uint32_t compute_optimal_substeps(float dt_kuramoto, float dx = 1.0f, float R_max = 1.0f) const {
        float m_max = Delta * R_max;
        if (m_max < 1e-6f) return 1;  // No mass, no sub-stepping needed

        // Required Dirac timestep for stability
        float dt_dirac_safe = 0.5f * dx / m_max;

        // How many sub-steps needed?
        uint32_t N_substeps = static_cast<uint32_t>(std::ceil(dt_kuramoto / dt_dirac_safe));

        // Clamp to reasonable range
        return std::max(1u, std::min(N_substeps, 100u));
    }

    void auto_configure_substeps(float dx = 1.0f, float R_max = 1.0f) {
        dirac_substeps = compute_optimal_substeps(dt, dx, R_max);
    }
};

} // namespace MSFT
} // namespace Nova
