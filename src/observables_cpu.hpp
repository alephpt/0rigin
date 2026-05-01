#pragma once

#include <vector>
#include <cstdint>

class Dirac3D;

// CPU-side observables for the Dirac field.
//
// The Vulkan ObservablesEngine in src/observables.hpp/.cpp is currently a
// stub (declared but the engine bodies are empty). The functions here
// provide the same physical quantities on the CPU side — they read the
// spinor components directly out of Dirac3D, which keeps the rest of the
// test framework portable and avoids reviving the GPU pipeline solely for
// observable reduction.
//
// Conventions match the corrected Dirac-basis γ⁵ in src/Dirac3D.cpp:
//   β  = diag(1, 1, -1, -1)
//   γ⁵ = anti-diagonal block [[0, I], [I, 0]]
// In particular, βγ⁵ is anti-Hermitian, so ⟨ψ̄γ⁵ψ⟩ is purely imaginary
// and ⟨ψ̄iγ⁵ψ⟩ is purely real (the standard pseudoscalar density).
namespace observables_cpu {

/**
 * Local scalar density: σ(x) = ψ̄ψ(x) = ψ†βψ(x)
 *                            = |ψ₀|² + |ψ₁|² − |ψ₂|² − |ψ₃|²
 * Returns one float per lattice site.
 */
std::vector<float> scalarDensity(const Dirac3D& dirac);

/**
 * Local pseudoscalar density: π(x) = ⟨ψ̄ iγ⁵ ψ⟩(x) = i·ψ†βγ⁵ψ(x)
 *                                  = −2·Im( ψ̄₀ψ₂ + ψ̄₁ψ₃ )
 * Returns one float per lattice site.
 */
std::vector<float> pseudoscalarDensity(const Dirac3D& dirac);

/**
 * Volume-integrated condensates. Both return the lattice sum (no dx³ factor;
 * tests can multiply by the lattice volume element if they want a continuum
 * density).
 */
float scalarCondensate(const Dirac3D& dirac);
float pseudoscalarCondensate(const Dirac3D& dirac);

}  // namespace observables_cpu
