/**
 * Dirac3D.h
 *
 * 3D+1 Dirac Equation Solver
 *
 * Implements 4-component Dirac spinor evolution in 3D space:
 *   iℏ ∂ψ/∂t = H ψ
 *   H = -iℏc α·∇ + βmc²
 *
 * where:
 *   - α = (α₁, α₂, α₃) are 4×4 Dirac matrices
 *   - β is the 4×4 mass matrix
 *   - ψ = (ψ₁, ψ₂, ψ₃, ψ₄) is the 4-component spinor
 *
 * Split-operator method:
 *   1. Kinetic step: exp(-iα·p dt/2) in momentum space
 *   2. Mass step: exp(-iβm dt) in position space
 *   3. Kinetic step: exp(-iα·p dt/2) in momentum space
 */

#pragma once

#include <complex>
#include <vector>
#include <array>
#include <cstdint>

class Dirac3D {
public:
    /**
     * Constructor
     * @param Nx Grid size in X dimension
     * @param Ny Grid size in Y dimension
     * @param Nz Grid size in Z dimension
     */
    Dirac3D(uint32_t Nx, uint32_t Ny, uint32_t Nz);

    /**
     * Destructor (cleanup FFTW resources)
     */
    ~Dirac3D();

    /**
     * Initialize 4-component spinor field
     * @param psi_init Initial spinor field (4 × N_total complex values)
     */
    void initialize(const std::vector<std::complex<float>>& psi_init);

    /**
     * Initialize Gaussian wavepacket (for testing)
     * @param x0 Center X coordinate
     * @param y0 Center Y coordinate
     * @param z0 Center Z coordinate
     * @param sigma Wavepacket width
     */
    void initializeGaussian(float x0, float y0, float z0, float sigma);

    /**
     * Split-operator evolution step
     * @param mass_field Scalar mass field m(x,y,z)
     * @param dt Time step
     */
    void step(const std::vector<float>& mass_field, float dt);

    /**
     * Split-operator evolution step with chiral mass coupling
     * @param R_field Vacuum R(x) field
     * @param theta_field Vacuum θ(x) field
     * @param Delta Coupling strength
     * @param dt Time step
     */
    void stepWithChiralMass(const std::vector<float>& R_field,
                           const std::vector<float>& theta_field,
                           float Delta, float dt);

    /**
     * Get spinor density |ψ|² at each point
     * ρ = Σ_α |ψ_α|²
     */
    std::vector<float> getDensity() const;

    /**
     * Get probability current j^i = ψ† α^i ψ
     * @param component Which spatial component (0=x, 1=y, 2=z)
     * @return Current density at each grid point
     */
    std::vector<float> getCurrent(int component) const;

    /**
     * Get individual spinor component
     * @param component Which spinor component (0-3)
     */
    const std::vector<std::complex<float>>& getComponent(int component) const;

    /**
     * Get total norm (should be conserved)
     * ∫ ψ†ψ d³x
     */
    float getNorm() const;

    /**
     * Get grid dimensions
     */
    uint32_t getNx() const { return _Nx; }
    uint32_t getNy() const { return _Ny; }
    uint32_t getNz() const { return _Nz; }
    uint32_t getTotalPoints() const { return _N_total; }

    /**
     * Map 3D coordinates to linear index
     */
    uint32_t index3D(uint32_t i, uint32_t j, uint32_t k) const {
        return k * (_Nx * _Ny) + j * _Nx + i;
    }

private:
    // Grid dimensions
    uint32_t _Nx, _Ny, _Nz;
    uint32_t _N_total;

    // 4-component spinor field (position space)
    std::vector<std::complex<float>> _psi[4];

    // Momentum space buffers
    std::vector<std::complex<float>> _psi_k[4];

    // FFT plans (FFTW)
    void* _fft_forward[4];   // Position → Momentum
    void* _fft_backward[4];  // Momentum → Position

    // Momentum grid (kx, ky, kz)
    std::vector<float> _kx, _ky, _kz;

    // Grid spacing
    float _dx;

    // Setup FFT plans and momentum grid
    void setupFFT();
    void setupMomentumGrid();
    void cleanupFFT();

    // Split-operator sub-steps
    void applyKineticHalfStep(float dt_half);
    void applyMassStep(const std::vector<float>& mass_field, float dt);
    // NOTE: applyChiralMassStep() removed - obsolete scalar approximation
    // Use stepWithChiralMass() for production with correct eigenvalue decomposition

    // Apply Dirac kinetic operator in momentum space
    void applyDiracKineticMatrix(std::complex<float> psi_k[4],
                                 float kx, float ky, float kz,
                                 float dt);

    // Velocity Verlet for mass evolution with full chiral coupling
    void applyMassVelocityVerlet(
        const std::vector<float>& R_field,
        const std::vector<float>& theta_field,
        float Delta, float dt);

    // Compute mass derivative: dΨ/dt = -i·β·M·Ψ
    void computeMassDerivative(
        const std::vector<float>& R_field,
        const std::vector<float>& theta_field,
        float Delta,
        const std::vector<std::complex<float>> psi_in[4],
        std::vector<std::complex<float>> dpsi_dt[4]);

    // Dirac alpha matrices (4×4 complex, stored as flat arrays)
    // α^i are the spatial Dirac matrices
    static const std::array<std::complex<float>, 16> alpha_x;
    static const std::array<std::complex<float>, 16> alpha_y;
    static const std::array<std::complex<float>, 16> alpha_z;

    // Dirac beta matrix (4×4 complex)
    static const std::array<std::complex<float>, 16> beta;

    // Dirac gamma5 matrix (4×4 complex)
    static const std::array<std::complex<float>, 16> gamma5;
};
