#include "observables_cpu.hpp"
#include "Dirac3D.h"

#include <complex>

namespace observables_cpu {

std::vector<float> scalarDensity(const Dirac3D& dirac) {
    const auto& p0 = dirac.getComponent(0);
    const auto& p1 = dirac.getComponent(1);
    const auto& p2 = dirac.getComponent(2);
    const auto& p3 = dirac.getComponent(3);
    const uint32_t N = dirac.getTotalPoints();
    std::vector<float> out(N);
    for (uint32_t i = 0; i < N; ++i) {
        out[i] = std::norm(p0[i]) + std::norm(p1[i])
               - std::norm(p2[i]) - std::norm(p3[i]);
    }
    return out;
}

std::vector<float> pseudoscalarDensity(const Dirac3D& dirac) {
    const auto& p0 = dirac.getComponent(0);
    const auto& p1 = dirac.getComponent(1);
    const auto& p2 = dirac.getComponent(2);
    const auto& p3 = dirac.getComponent(3);
    const uint32_t N = dirac.getTotalPoints();
    std::vector<float> out(N);
    for (uint32_t i = 0; i < N; ++i) {
        // βγ⁵ψ swaps upper/lower then flips lower-block sign:
        //   (βγ⁵ψ)_0 = ψ₂, (βγ⁵ψ)_1 = ψ₃,
        //   (βγ⁵ψ)_2 = -ψ₀, (βγ⁵ψ)_3 = -ψ₁
        // ψ†(βγ⁵)ψ = ψ̄₀ψ₂ + ψ̄₁ψ₃ - ψ̄₂ψ₀ - ψ̄₃ψ₁
        //          = (ψ̄₀ψ₂ + ψ̄₁ψ₃) - conj(ψ̄₀ψ₂ + ψ̄₁ψ₃)
        //          = 2i·Im(ψ̄₀ψ₂ + ψ̄₁ψ₃)         [purely imaginary]
        // i·(ψ†βγ⁵ψ) = -2·Im(ψ̄₀ψ₂ + ψ̄₁ψ₃)   [purely real]
        const std::complex<float> bg5 = std::conj(p0[i]) * p2[i]
                                        + std::conj(p1[i]) * p3[i];
        out[i] = -2.0f * bg5.imag();
    }
    return out;
}

float scalarCondensate(const Dirac3D& dirac) {
    auto s = scalarDensity(dirac);
    double total = 0.0;
    for (float v : s) total += v;
    return static_cast<float>(total);
}

float pseudoscalarCondensate(const Dirac3D& dirac) {
    auto p = pseudoscalarDensity(dirac);
    double total = 0.0;
    for (float v : p) total += v;
    return static_cast<float>(total);
}

}  // namespace observables_cpu
