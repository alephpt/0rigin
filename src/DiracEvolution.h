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
    DiracEvolution(uint32_t Nx, uint32_t Ny, float beta_sign = +1.0f);
    ~DiracEvolution();

    // Initialize 4-component spinor field
    void initialize(float x0, float y0, float sigma);

    // Initialize as plane wave with momentum (kx, ky)
    // Used for dispersion relation analysis
    void initializePlaneWave(float kx, float ky);

    // Split-operator evolution step
    // If time_dilation_mode enabled and R_field provided, uses dt_eff = R_center·dt
    void step(const std::vector<float>& mass_field, float dt,
              const std::vector<float>* R_field = nullptr);

    // Get spinor density |Ψ|² at each point
    std::vector<float> getDensity() const;

    // Get individual spinor components (for analysis)
    const std::vector<std::complex<float>>& getComponent(int c) const;

    // Get full spinor field as flat array [4 * Nx * Ny]
    // Element [i*4 + alpha] = spinor component alpha at grid point i
    std::vector<std::complex<double>> getSpinorField() const;

    // Get grid parameters
    int getNx() const { return _Nx; }
    int getNy() const { return _Ny; }
    double getDx() const { return 1.0; } // Unit grid spacing

    // Get beta sign (particle/antiparticle charge)
    float getBetaSign() const { return _beta_sign; }

    // Norm check (should always be 1.0 to machine precision)
    float getNorm() const;

    // === Time Dilation Mode (Phase 4 Test 4.1) ===

    // Enable/disable time dilation: dτ = R(x)·dt
    // When enabled, effective timestep becomes dt_eff = R_local·dt
    void setTimeDilationMode(bool enable) { _time_dilation_mode = enable; }
    bool getTimeDilationMode() const { return _time_dilation_mode; }

    // === Electromagnetic Coupling Mode (Phase 5 Test 5.1) ===

    // Enable/disable EM coupling: Minimal coupling ∇ → ∇ - iqA
    // When enabled, gauge potential A from Kuramoto phase affects evolution
    void setEMCouplingMode(bool enable, float coupling_strength = 1.0f) {
        _em_coupling_enabled = enable;
        _em_coupling_strength = coupling_strength;
    }
    bool getEMCouplingMode() const { return _em_coupling_enabled; }
    float getEMCouplingStrength() const { return _em_coupling_strength; }

    // Get R-field value at wavepacket center via bilinear interpolation
    // Used for time-dilation evolution: dt_eff = R_center·dt
    float getRFieldAtPosition(const std::vector<float>& R_field, float x, float y) const;

    // === Physics Analysis Methods ===

    // Get beta expectation value: <Ψ|β|Ψ> = (upper) - (lower)
    // Returns: +1 for pure upper, -1 for pure lower, 0 for equal mix
    float getBetaExpectation() const;

    // Get center of mass of wavepacket
    void getCenterOfMass(float& x_mean, float& y_mean) const;

    // === Electromagnetic Coupling Evolution ===

    /**
     * Apply electromagnetic potential step
     *
     * Physics: U_EM = exp(-i q φ dt)
     * Where φ = ∂_t θ is scalar potential from Kuramoto phase
     *
     * This is a diagonal operator in position space
     *
     * @param phi_field: Scalar potential φ(x,y) [Nx × Ny as flat vector]
     * @param dt: Timestep
     */
    void applyEMPotentialStep(const std::vector<float>& phi_field, float dt);

    /**
     * Apply minimal coupling to kinetic step
     *
     * Physics: Replace ∇ → ∇ - iq A in momentum operator
     * Modified Hamiltonian: H = α·(p - qA) + βm
     *
     * Implementation: Peierls substitution in real space
     * ∂_x ψ → (ψ[i+1] - ψ[i-1])/(2dx) - iq A_x ψ
     *
     * Note: This must be called BEFORE the kinetic half-step
     *
     * @param A_x_field: Vector potential x-component [Nx × Ny as flat vector]
     * @param A_y_field: Vector potential y-component [Nx × Ny as flat vector]
     */
    void applyMinimalCoupling(const std::vector<float>& A_x_field,
                              const std::vector<float>& A_y_field);

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
    float _beta_sign;  // +1 for particle, -1 for antiparticle

    // Time dilation mode flag (Phase 4)
    bool _time_dilation_mode;

    // Electromagnetic coupling mode flags (Phase 5)
    bool _em_coupling_enabled;
    float _em_coupling_strength;  // Effective charge q

    // Vector potential storage (for Peierls substitution)
    std::vector<float> _A_x_field;
    std::vector<float> _A_y_field;

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
