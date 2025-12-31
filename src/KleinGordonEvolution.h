/**
 * KleinGordonEvolution.h
 *
 * Split-operator method for 2D Klein-Gordon equation evolution
 * Implements unitary time propagation via Strang splitting
 *
 * Physics:
 *   Klein-Gordon equation: (∂²_t - ∇² + m²)φ = 0
 *
 *   where:
 *     φ(x,y,t): Complex scalar field
 *     m(x,y) = Δ·R(x,y): Position-dependent mass from Kuramoto synchronization
 *     Δ: Mass gap parameter (constant, typically 1.0)
 *     R(x,y): Synchronization order parameter
 *
 * Comparison with Dirac:
 *   - Dirac: 4-component spinor (particle + antiparticle, spin up/down)
 *   - Klein-Gordon: Single complex scalar (no spin structure)
 *   - Both are relativistic, but different particle types
 *   - Question: Does scalar field perform better at v>0.5c?
 */

#pragma once

#include <complex>
#include <vector>
#include <cstdint>

class KleinGordonEvolution {
public:
    KleinGordonEvolution(uint32_t Nx, uint32_t Ny);
    ~KleinGordonEvolution();

    // Initialize scalar field as Gaussian wavepacket
    void initialize(float x0, float y0, float sigma);

    // Initialize as plane wave with momentum (kx, ky)
    // Used for dispersion relation analysis
    void initializePlaneWave(float kx, float ky);

    // Split-operator evolution step
    // Evolves both φ and ∂_tφ (Klein-Gordon is second-order in time)
    void step(const std::vector<float>& mass_field, float dt);

    // Get scalar field density |φ|² at each point
    std::vector<float> getDensity() const;

    // Get individual field components (for analysis)
    const std::vector<std::complex<float>>& getField() const { return _phi; }
    const std::vector<std::complex<float>>& getFieldDot() const { return _phi_dot; }

    // Get full scalar field as flat array [Nx * Ny]
    std::vector<std::complex<double>> getScalarField() const;

    // Get grid parameters
    int getNx() const { return _Nx; }
    int getNy() const { return _Ny; }
    double getDx() const { return 1.0; } // Unit grid spacing

    // Norm check (should always be conserved to machine precision)
    float getNorm() const;

    // === Physics Analysis Methods ===

    // Get center of mass of wavepacket
    void getCenterOfMass(float& x_mean, float& y_mean) const;

    // Get momentum space density |φ̃(k)|² for dispersion analysis
    // Returns kx, ky grids and momentum-space density
    void getMomentumDistribution(std::vector<float>& kx_out,
                                  std::vector<float>& ky_out,
                                  std::vector<float>& density_k) const;

    // === CRITICAL DIAGNOSTIC: Energy ===

    // Compute total energy E = <H> where H is the Klein-Gordon Hamiltonian
    // For Klein-Gordon: E = ∫[|∂_tφ|² + |∇φ|² + m²|φ|²] dx
    //
    // Returns: Total energy (kinetic + gradient + mass)
    //
    // Also returns KE and PE separately via output parameters
    float getEnergy(const std::vector<float>& mass_field,
                   float& KE_out, float& PE_out) const;

private:
    uint32_t _Nx, _Ny;
    uint32_t _N_points;

    // Scalar field and time derivative: φ(x,y,t) and ∂_tφ(x,y,t)
    std::vector<std::complex<float>> _phi;      // Scalar field
    std::vector<std::complex<float>> _phi_dot;  // Time derivative

    // FFT plans and buffers (FFTW)
    void* _fft_forward_phi;
    void* _fft_backward_phi;
    void* _fft_forward_phi_dot;
    void* _fft_backward_phi_dot;

    std::vector<std::complex<float>> _phi_k;      // φ in momentum space
    std::vector<std::complex<float>> _phi_dot_k;  // ∂_tφ in momentum space

    // Momentum grid
    std::vector<float> _kx, _ky;

    // K-space cache validity flag (for lazy energy computation)
    mutable bool _phi_k_valid;

    // Split-operator sub-steps
    void applyMassStep(const std::vector<float>& mass_field, float dt);
    void applyKineticHalfStep(float dt_half);

    // Helper: apply Klein-Gordon kinetic evolution in momentum space
    void applyKGKineticMatrix(std::complex<float>& phi_k,
                              std::complex<float>& phi_dot_k,
                              float k_mag, float dt);

    // Setup FFT plans and momentum grid
    void setupFFT();
    void setupMomentumGrid();
    void cleanupFFT();
};
