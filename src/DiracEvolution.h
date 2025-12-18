/**
 * DiracEvolution.h
 *
 * Split-operator method for 2D Dirac equation evolution
 * Implements unitary time propagation via Strang splitting
 */

#pragma once

#include <complex>
#include <vector>
#include <cstdint>

class DiracEvolution {
public:
    DiracEvolution(uint32_t Nx, uint32_t Ny);
    ~DiracEvolution();

    // Initialize 4-component spinor field
    void initialize(float x0, float y0, float sigma);

    // Split-operator evolution step
    void step(const std::vector<float>& mass_field, float dt);

    // Get spinor density |Ψ|² at each point
    std::vector<float> getDensity() const;

    // Get individual spinor components (for analysis)
    const std::vector<std::complex<float>>& getComponent(int c) const;

    // Norm check (should always be 1.0 to machine precision)
    float getNorm() const;

    // === Physics Analysis Methods ===

    // Get beta expectation value: <Ψ|β|Ψ> = (upper) - (lower)
    // Returns: +1 for pure upper, -1 for pure lower, 0 for equal mix
    float getBetaExpectation() const;

    // Get center of mass of wavepacket
    void getCenterOfMass(float& x_mean, float& y_mean) const;

    // Get momentum space density |Ψ̃(k)|² for dispersion analysis
    // Returns kx, ky grids and momentum-space density
    void getMomentumDistribution(std::vector<float>& kx_out,
                                  std::vector<float>& ky_out,
                                  std::vector<float>& density_k) const;

    // === CRITICAL DIAGNOSTIC: Energy ===

    // Compute total energy E = <Ψ|H|Ψ> where H = -iα·∇ + βm(x)
    // Returns: Total energy (kinetic + potential)
    //
    // For bound state: E < 0
    // For unbound: E > 0
    //
    // Also returns KE and PE separately via output parameters
    float getEnergy(const std::vector<float>& mass_field,
                   float& KE_out, float& PE_out) const;

private:
    uint32_t _Nx, _Ny;
    uint32_t _N_points;

    // 4-component spinor: psi[0..3][y*Nx + x]
    std::vector<std::complex<float>> _psi[4];

    // FFT plans and buffers (FFTW)
    void* _fft_forward[4];
    void* _fft_backward[4];
    std::vector<std::complex<float>> _psi_k[4];

    // Momentum grid
    std::vector<float> _kx, _ky;

    // K-space cache validity flag (for lazy energy computation)
    mutable bool _psi_k_valid;

    // Split-operator sub-steps
    void applyPotentialStep(const std::vector<float>& mass_field, float dt);
    void applyKineticHalfStep(float dt_half);

    // Helper: apply Dirac kinetic matrix at single k-point
    void applyDiracKineticMatrix(std::complex<float> psi_k[4],
                                  float k_mag, float kx, float ky,
                                  float dt);

    // Setup FFT plans and momentum grid
    void setupFFT();
    void setupMomentumGrid();
    void cleanupFFT();
};
