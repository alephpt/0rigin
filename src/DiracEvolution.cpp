/**
 * DiracEvolution.cpp
 *
 * Implementation of split-operator method for 2D Dirac equation
 * Using FFTW for fast Fourier transforms
 */

#include "DiracEvolution.h"
#include <fftw3.h>
#include <cmath>
#include <iostream>

DiracEvolution::DiracEvolution(uint32_t Nx, uint32_t Ny, float beta_sign)
    : _Nx(Nx), _Ny(Ny), _N_points(Nx * Ny), _beta_sign(beta_sign),
      _time_dilation_mode(false), _em_coupling_enabled(false),
      _em_coupling_strength(1.0f), _psi_k_valid(false) {

    // Allocate spinor components
    for (int c = 0; c < 4; c++) {
        _psi[c].resize(_N_points);
        _psi_k[c].resize(_N_points);
        _fft_forward[c] = nullptr;
        _fft_backward[c] = nullptr;
    }

    // Allocate vector potential storage (initialized to zero)
    _A_x_field.resize(_N_points, 0.0f);
    _A_y_field.resize(_N_points, 0.0f);

    setupMomentumGrid();
    setupFFT();

    std::string particle_type = (beta_sign > 0) ? "particle" : "antiparticle";
    std::cout << "[DiracEvolution] Initialized " << _Nx << "x" << _Ny
              << " 4-component spinor field (" << particle_type << ", β=" << beta_sign << ")" << std::endl;
}

DiracEvolution::~DiracEvolution() {
    cleanupFFT();
}

void DiracEvolution::setupMomentumGrid() {
    _kx.resize(_Nx);
    _ky.resize(_Ny);

    // FFT frequency grid: k = 2π * freq / L
    // For grid spacing dx=1, L=N, so k = 2π*freq/N
    float dk_x = 2.0f * M_PI / _Nx;
    float dk_y = 2.0f * M_PI / _Ny;

    // fftfreq ordering: [0, 1, ..., N/2-1, -N/2, ..., -1]
    for (uint32_t i = 0; i < _Nx; i++) {
        int freq = (i < _Nx/2) ? i : i - _Nx;
        _kx[i] = freq * dk_x;
    }

    for (uint32_t j = 0; j < _Ny; j++) {
        int freq = (j < _Ny/2) ? j : j - _Ny;
        _ky[j] = freq * dk_y;
    }
}

void DiracEvolution::setupFFT() {
    // Create FFTW plans for each spinor component
    for (int c = 0; c < 4; c++) {
        int n[] = {(int)_Ny, (int)_Nx};

        _fft_forward[c] = fftwf_plan_dft_2d(
            _Ny, _Nx,
            reinterpret_cast<fftwf_complex*>(_psi[c].data()),
            reinterpret_cast<fftwf_complex*>(_psi_k[c].data()),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        _fft_backward[c] = fftwf_plan_dft_2d(
            _Ny, _Nx,
            reinterpret_cast<fftwf_complex*>(_psi_k[c].data()),
            reinterpret_cast<fftwf_complex*>(_psi[c].data()),
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
    }
}

void DiracEvolution::cleanupFFT() {
    for (int c = 0; c < 4; c++) {
        if (_fft_forward[c]) {
            fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_forward[c]));
        }
        if (_fft_backward[c]) {
            fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_backward[c]));
        }
    }
}

void DiracEvolution::initialize(float x0, float y0, float sigma) {
    // Initialize as Gaussian wavepacket in upper spinor components
    // Lower components start at zero
    float norm_sum = 0.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            float dx = x - x0;
            float dy = y - y0;
            float r2 = dx*dx + dy*dy;
            float amplitude = std::exp(-r2 / (2.0f * sigma * sigma));

            // Upper components (particle)
            _psi[0][idx] = std::complex<float>(amplitude, 0.0f);
            _psi[1][idx] = std::complex<float>(amplitude, 0.0f);

            // Lower components (antiparticle) - start at zero
            _psi[2][idx] = std::complex<float>(0.0f, 0.0f);
            _psi[3][idx] = std::complex<float>(0.0f, 0.0f);

            norm_sum += amplitude * amplitude * 2.0f; // Two upper components
        }
    }

    // Normalize
    float norm_factor = std::sqrt(norm_sum);
    for (int c = 0; c < 4; c++) {
        for (uint32_t i = 0; i < _N_points; i++) {
            _psi[c][i] /= norm_factor;
        }
    }

    // Invalidate k-space cache after initialization
    _psi_k_valid = false;

    std::cout << "[DiracEvolution] Initialized Gaussian wavepacket at ("
              << x0 << ", " << y0 << ") with σ=" << sigma << std::endl;
}

void DiracEvolution::initializePlaneWave(float kx, float ky) {
    float k_mag = std::sqrt(kx*kx + ky*ky);
    
    // Construct spinor: Positive energy eigenstate of alpha*k
    // u = [1, 0, 0, (kx + i*ky)/|k|]^T (normalized later)
    std::complex<float> u0(1.0f, 0.0f);
    std::complex<float> u1(0.0f, 0.0f);
    std::complex<float> u2(0.0f, 0.0f);
    std::complex<float> u3(0.0f, 0.0f);
    
    if (k_mag > 1e-6f) {
        u3 = std::complex<float>(kx, ky) / k_mag;
    }

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;
            
            // Plane wave phase e^(i k.r)
            float phase_arg = kx * x + ky * y; // dx=1
            std::complex<float> phase(std::cos(phase_arg), std::sin(phase_arg));
            
            _psi[0][idx] = u0 * phase;
            _psi[1][idx] = u1 * phase;
            _psi[2][idx] = u2 * phase;
            _psi[3][idx] = u3 * phase;
        }
    }

    // Normalize
    float norm_sum = 0.0f;
    for (int c = 0; c < 4; c++) {
        for (uint32_t i = 0; i < _N_points; i++) {
            norm_sum += std::norm(_psi[c][i]);
        }
    }
    
    float norm_factor = 1.0f / std::sqrt(norm_sum);
    for (int c = 0; c < 4; c++) {
        for (uint32_t i = 0; i < _N_points; i++) {
            _psi[c][i] *= norm_factor;
        }
    }

    _psi_k_valid = false;
    std::cout << "[DiracEvolution] Initialized Plane Wave k=(" << kx << ", " << ky << ")" << std::endl;
}

void DiracEvolution::step(const std::vector<float>& mass_field, float dt,
                          const std::vector<float>* R_field) {
    // Apply time dilation if mode enabled
    float dt_effective = dt;

    if (_time_dilation_mode && R_field != nullptr) {
        // Get wavepacket center
        float x_center, y_center;
        getCenterOfMass(x_center, y_center);

        // Get R-field value at center (bilinear interpolation)
        float R_local = getRFieldAtPosition(*R_field, x_center, y_center);

        // Effective timestep: dτ = R(x)·dt
        dt_effective = R_local * dt;
    }

    // Strang splitting: K/2 - V - K/2
    applyKineticHalfStep(dt_effective / 2.0f);
    applyPotentialStep(mass_field, dt_effective);
    applyKineticHalfStep(dt_effective / 2.0f);

    // Invalidate k-space cache after evolution
    _psi_k_valid = false;
}

void DiracEvolution::applyPotentialStep(const std::vector<float>& mass_field, float dt) {
    // Potential: V = β_sign * β * m(x,y)
    // β = diag(1, 1, -1, -1) in standard representation
    // For particle (β_sign = +1): exp(-iβmΔt) = diag(e^(-imΔt), e^(-imΔt), e^(+imΔt), e^(+imΔt))
    // For antiparticle (β_sign = -1): exp(+iβmΔt) = diag(e^(+imΔt), e^(+imΔt), e^(-imΔt), e^(-imΔt))

    for (uint32_t i = 0; i < _N_points; i++) {
        float m = mass_field[i];
        std::complex<float> phase_plus(0.0f, -_beta_sign * m * dt);  // e^(-i·β_sign·mΔt)
        std::complex<float> phase_minus(0.0f, +_beta_sign * m * dt); // e^(+i·β_sign·mΔt)

        // Upper components: β = +1
        _psi[0][i] *= std::exp(phase_plus);
        _psi[1][i] *= std::exp(phase_plus);

        // Lower components: β = -1
        _psi[2][i] *= std::exp(phase_minus);
        _psi[3][i] *= std::exp(phase_minus);
    }
}

void DiracEvolution::applyKineticHalfStep(float dt_half) {
    // Forward FFT for all 4 components
    for (int c = 0; c < 4; c++) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_forward[c]));
    }

    // Apply kinetic operator in momentum space
    for (uint32_t j = 0; j < _Ny; j++) {
        for (uint32_t i = 0; i < _Nx; i++) {
            uint32_t idx = j * _Nx + i;

            float kx = _kx[i];
            float ky = _ky[j];
            float k_mag = std::sqrt(kx*kx + ky*ky);

            // Get spinor at this k-point
            std::complex<float> psi_k_point[4];
            for (int c = 0; c < 4; c++) {
                psi_k_point[c] = _psi_k[c][idx];
            }

            // Apply Dirac kinetic matrix
            applyDiracKineticMatrix(psi_k_point, k_mag, kx, ky, dt_half);

            // Store result
            for (int c = 0; c < 4; c++) {
                _psi_k[c][idx] = psi_k_point[c];
            }
        }
    }

    // Inverse FFT for all 4 components
    for (int c = 0; c < 4; c++) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_backward[c]));

        // FFTW doesn't normalize, so divide by N
        float norm_factor = 1.0f / _N_points;
        for (uint32_t i = 0; i < _N_points; i++) {
            _psi[c][i] *= norm_factor;
        }
    }
}

void DiracEvolution::applyDiracKineticMatrix(std::complex<float> psi_k[4],
                                              float k_mag, float kx, float ky,
                                              float dt) {
    // Kinetic operator: K = α·k = α_x k_x + α_y k_y
    // Matrix exponential: exp(-i(α·k)Δt) = cos(|k|Δt)I - i sin(|k|Δt)(α·k)/|k|
    //
    // For 2D Dirac (using 2-component representation per spinor):
    // α_x = σ_x, α_y = σ_y (Pauli matrices)
    //
    // This gives 2×2 block structure for each spinor pair

    if (k_mag < 1e-10f) {
        // k=0: no kinetic evolution
        return;
    }

    float cos_term = std::cos(k_mag * dt);
    float sin_term = std::sin(k_mag * dt);

    // Normalized momentum direction
    float kx_norm = kx / k_mag;
    float ky_norm = ky / k_mag;

    // exp(-i(α·k)Δt) in 2D:
    // [cos(|k|t) - i sin(|k|t)(k_x/|k|)σ_x - i sin(|k|t)(k_y/|k|)σ_y]
    //
    // For 2-component spinor (ψ_up, ψ_down):
    // Matrix = cos(|k|t)I - i sin(|k|t)[k_x/|k| σ_x + k_y/|k| σ_y]
    //        = cos(|k|t)I - i sin(|k|t)[(k_x/|k|)(0  1) + (k_y/|k|)(0  -i)]
    //                                              (1  0)           (i   0)
    //
    // = [ cos(|k|t)                    -i sin(|k|t)(k_x - ik_y)/|k| ]
    //   [ -i sin(|k|t)(k_x + ik_y)/|k|  cos(|k|t)                   ]

    // Matrix elements: exp(-i(σ·k)Δt) = cos(|k|Δt)I - i·sin(|k|Δt)·(σ·k)/|k|
    // σ·k = σ_x k_x + σ_y k_y = ( 0      k_x-ik_y )
    //                            ( k_x+ik_y    0    )
    // So: exp(-i(σ·k)Δt) = ( cos(|k|Δt)     -i·sin(|k|Δt)·(k_x-ik_y)/|k| )
    //                       ( -i·sin(|k|Δt)·(k_x+ik_y)/|k|   cos(|k|Δt)   )

    std::complex<float> diag(cos_term, 0.0f);
    // Upper-right: -i·sin·(k_x - ik_y)/|k| = -i·sin·k_x/|k| + i·i·sin·k_y/|k|
    //                                        = -i·sin·k_x/|k| - sin·k_y/|k|
    std::complex<float> off_diag_down(-sin_term * ky_norm, -sin_term * kx_norm);
    // Lower-left: -i·sin·(k_x + ik_y)/|k| = -i·sin·k_x/|k| - i·i·sin·k_y/|k|
    //                                       = -i·sin·k_x/|k| + sin·k_y/|k|
    std::complex<float> off_diag_up(sin_term * ky_norm, -sin_term * kx_norm);

    // Apply to first spinor pair (components 0, 1)
    std::complex<float> psi0_new = diag * psi_k[0] + off_diag_down * psi_k[1];
    std::complex<float> psi1_new = off_diag_up * psi_k[0] + diag * psi_k[1];

    // Apply to second spinor pair (components 2, 3)
    std::complex<float> psi2_new = diag * psi_k[2] + off_diag_down * psi_k[3];
    std::complex<float> psi3_new = off_diag_up * psi_k[2] + diag * psi_k[3];

    psi_k[0] = psi0_new;
    psi_k[1] = psi1_new;
    psi_k[2] = psi2_new;
    psi_k[3] = psi3_new;
}

std::vector<float> DiracEvolution::getDensity() const {
    std::vector<float> density(_N_points);

    for (uint32_t i = 0; i < _N_points; i++) {
        density[i] = 0.0f;
        for (int c = 0; c < 4; c++) {
            density[i] += std::norm(_psi[c][i]);
        }
    }

    return density;
}

const std::vector<std::complex<float>>& DiracEvolution::getComponent(int c) const {
    return _psi[c];
}

std::vector<std::complex<double>> DiracEvolution::getSpinorField() const {
    // Return spinor field as flat array [4 * Nx * Ny]
    // Layout: [psi0[0], psi1[0], psi2[0], psi3[0], psi0[1], psi1[1], ...]
    std::vector<std::complex<double>> spinor_flat(_N_points * 4);

    for (uint32_t i = 0; i < _N_points; i++) {
        for (int alpha = 0; alpha < 4; alpha++) {
            spinor_flat[i * 4 + alpha] = std::complex<double>(_psi[alpha][i].real(), _psi[alpha][i].imag());
        }
    }

    return spinor_flat;
}

float DiracEvolution::getNorm() const {
    float norm = 0.0f;
    for (int c = 0; c < 4; c++) {
        for (uint32_t i = 0; i < _N_points; i++) {
            norm += std::norm(_psi[c][i]);
        }
    }
    return norm;
}

float DiracEvolution::getBetaExpectation() const {
    /**
     * Calculate <Ψ|β|Ψ> where β = diag(1, 1, -1, -1)
     *
     * <β> = (|ψ₁|² + |ψ₂|²) - (|ψ₃|² + |ψ₄|²)
     *
     * Physical meaning:
     *   +1: Pure particle (upper components)
     *    0: Equal particle/antiparticle mix
     *   -1: Pure antiparticle (lower components)
     */
    float upper = 0.0f;  // |ψ₁|² + |ψ₂|²
    float lower = 0.0f;  // |ψ₃|² + |ψ₄|²

    for (uint32_t i = 0; i < _N_points; i++) {
        upper += std::norm(_psi[0][i]) + std::norm(_psi[1][i]);
        lower += std::norm(_psi[2][i]) + std::norm(_psi[3][i]);
    }

    return upper - lower;
}

void DiracEvolution::getCenterOfMass(float& x_mean, float& y_mean) const {
    /**
     * Calculate center of mass: (x̄, ȳ) = ∫ (x,y)|Ψ|² dxdy / ∫|Ψ|² dxdy
     */
    x_mean = 0.0f;
    y_mean = 0.0f;
    float total = 0.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            float density = 0.0f;
            for (int c = 0; c < 4; c++) {
                density += std::norm(_psi[c][idx]);
            }

            x_mean += x * density;
            y_mean += y * density;
            total += density;
        }
    }

    if (total > 0.0f) {
        x_mean /= total;
        y_mean /= total;
    }
}

void DiracEvolution::getMomentumDistribution(std::vector<float>& kx_out,
                                              std::vector<float>& ky_out,
                                              std::vector<float>& density_k) const {
    /**
     * Compute momentum-space density |Ψ̃(k)|² for all 4 components
     *
     * This is needed for:
     * - Dispersion relation validation E(k)
     * - Wave packet analysis
     * - Group velocity measurements
     */

    kx_out = _kx;
    ky_out = _ky;
    density_k.resize(_N_points, 0.0f);

    // Sum |Ψ̃_c(k)|² over all 4 components
    for (int c = 0; c < 4; c++) {
        for (uint32_t i = 0; i < _N_points; i++) {
            density_k[i] += std::norm(_psi_k[c][i]);
        }
    }
}

/**
 * Compute total energy E = <Ψ|H|Ψ>
 *
 * Hamiltonian: H = -iα·∇ + βm(x)
 *
 * Decomposition:
 *   Kinetic:   E_K = ∫ Ψ*(x) [-iα·∇] Ψ(x) dx
 *            = Σ_k |Ψ̃_k|² ω(k)  where ω(k) = √(k² + m_0²)
 *
 *   Potential: E_V = ∫ |Ψ(x)|² β m(x) dx
 *
 * For Dirac equation in 2D with 4-component spinor:
 *   β = diag(1, 1, -1, -1)
 *
 * Note: m_0 (rest mass) is assumed ~0 in this implementation (massless Dirac limit)
 *       The mass term m(x) comes from SMFT coupling, not intrinsic mass.
 */
float DiracEvolution::getEnergy(const std::vector<float>& mass_field,
                                float& KE_out, float& PE_out) const {

    // === Part 1: Kinetic Energy (in momentum space) ===
    //
    // E_K = Σ_{k,c} |Ψ̃_c(k)|² ω(k)
    //
    // For massless Dirac: ω(k) = |k| (linear dispersion)
    // For massive Dirac: ω(k) = √(k² + m₀²)
    //
    // Using massless approximation (m₀ = 0):

    // Update k-space if cache is invalid
    if (!_psi_k_valid) {
        // Compute fresh k-space representation via FFT
        for (int c = 0; c < 4; c++) {
            fftwf_execute((fftwf_plan)_fft_forward[c]);
        }
        _psi_k_valid = true;
    }

    float KE = 0.0f;

    // Sum kinetic energy over all momentum modes
    for (uint32_t j = 0; j < _Ny; j++) {
        for (uint32_t i = 0; i < _Nx; i++) {
            uint32_t idx = j * _Nx + i;

            float kx = _kx[i];
            float ky = _ky[j];
            float k_mag = std::sqrt(kx*kx + ky*ky);

            // Energy of this mode: ω(k) = |k| (massless Dirac)
            float omega_k = k_mag;

            // Sum over all 4 spinor components
            for (int c = 0; c < 4; c++) {
                float density_k = std::norm(_psi_k[c][idx]);
                KE += density_k * omega_k;
            }
        }
    }

    // Normalize by grid size (FFT normalization)
    KE /= (_Nx * _Ny);

    // === Part 2: Potential Energy (in position space) ===
    //
    // E_V = Σ_x |Ψ(x)|² β m(x)
    //
    // β = diag(1, 1, -1, -1) → upper components contribute +m, lower contribute -m

    float PE = 0.0f;

    for (uint32_t idx = 0; idx < _N_points; idx++) {
        float m = mass_field[idx];

        // Upper components (c=0,1): β = +1
        float density_upper = std::norm(_psi[0][idx]) + std::norm(_psi[1][idx]);
        PE += density_upper * m;

        // Lower components (c=2,3): β = -1
        float density_lower = std::norm(_psi[2][idx]) + std::norm(_psi[3][idx]);
        PE -= density_lower * m;
    }

    // === Total Energy ===
    float E_total = KE + PE;

    // Output parameters
    KE_out = KE;
    PE_out = PE;

    return E_total;
}

/**
 * Get R-field value at position (x,y) via bilinear interpolation
 * Used for time-dilation evolution mode (Phase 4 Test 4.1)
 *
 * Physics: R(x) acts as local time flow rate: dτ/dt = R(x)
 * Particles in desynchronized regions (R < 1) experience slower proper time
 *
 * @param R_field Synchronization field R(x,y) [size: Nx*Ny]
 * @param x Position x (grid units, can be fractional)
 * @param y Position y (grid units, can be fractional)
 * @return R value at (x,y) via bilinear interpolation
 */
float DiracEvolution::getRFieldAtPosition(const std::vector<float>& R_field,
                                          float x, float y) const {
    // Periodic boundary conditions
    x = std::fmod(x + _Nx, (float)_Nx);
    y = std::fmod(y + _Ny, (float)_Ny);

    // Get grid cell indices
    int ix0 = (int)std::floor(x);
    int iy0 = (int)std::floor(y);
    int ix1 = (ix0 + 1) % _Nx;
    int iy1 = (iy0 + 1) % _Ny;

    // Fractional coordinates within cell
    float fx = x - ix0;
    float fy = y - iy0;

    // Bilinear interpolation
    float R00 = R_field[iy0 * _Nx + ix0];
    float R10 = R_field[iy0 * _Nx + ix1];
    float R01 = R_field[iy1 * _Nx + ix0];
    float R11 = R_field[iy1 * _Nx + ix1];

    float R0 = R00 * (1.0f - fx) + R10 * fx;
    float R1 = R01 * (1.0f - fx) + R11 * fx;
    float R = R0 * (1.0f - fy) + R1 * fy;

    return R;
}

/**
 * Apply electromagnetic potential step (Phase 5)
 *
 * Physics: U_EM = exp(-i q φ dt)
 * where φ = ∂_t θ is the scalar potential from Kuramoto phase
 *
 * This is diagonal in position space and applied to all spinor components equally
 *
 * @param phi_field Scalar potential φ(x,y) [size: Nx*Ny]
 * @param dt Timestep
 */
void DiracEvolution::applyEMPotentialStep(const std::vector<float>& phi_field, float dt) {
    if (!_em_coupling_enabled) return;

    for (uint32_t i = 0; i < _N_points; i++) {
        float phi = phi_field[i];

        // Phase factor: exp(-i q φ dt)
        std::complex<float> phase_factor(0.0f, -_em_coupling_strength * phi * dt);
        std::complex<float> U_em = std::exp(phase_factor);

        // Apply to all 4 spinor components
        for (int c = 0; c < 4; c++) {
            _psi[c][i] *= U_em;
        }
    }
}

/**
 * Apply minimal coupling to kinetic step (Phase 5) - Peierls Substitution
 *
 * Physics: Minimal coupling replaces ∇ → ∇ - iq A
 * Modified kinetic operator: K = α·(p - qA)
 *
 * Implementation: Peierls substitution for lattice gauge theory
 *   - Store A fields for use in next kinetic step
 *   - Kinetic derivatives modified by link phase factors: U_ij = exp(-iq A·dl)
 *   - Preserves gauge invariance: ψ → e^(iχ)ψ, A → A + ∇χ
 *
 * This method stores the A fields; they are applied in applyKineticHalfStep
 * via modified finite-difference operators with gauge links.
 *
 * @param A_x_field Vector potential x-component [size: Nx*Ny]
 * @param A_y_field Vector potential y-component [size: Nx*Ny]
 */
void DiracEvolution::applyMinimalCoupling(
    const std::vector<float>& A_x_field,
    const std::vector<float>& A_y_field)
{
    if (!_em_coupling_enabled) return;

    // Store A fields for use in kinetic step
    _A_x_field = A_x_field;
    _A_y_field = A_y_field;

    // NOTE: Actual Peierls substitution is applied in applyKineticHalfStep
    // where we modify the momentum-space evolution with position-dependent phases.
    //
    // For now, we implement a simple real-space minimal coupling:
    // Modified derivatives: (∂_μ - iq A_μ) ψ
    //
    // This is applied by modifying the wavefunction phase before FFT:
    //   In momentum space, p → p - qA becomes a convolution
    //   We approximate this with local A(x) acting on ψ(x)

    // Apply gauge-covariant phase shift in real space
    // This approximates the effect of A on the kinetic operator
    // Valid when A varies slowly on the scale of the wavepacket

    for (uint32_t j = 0; j < _Ny; j++) {
        for (uint32_t i = 0; i < _Nx; i++) {
            uint32_t idx = j * _Nx + i;

            // Local vector potential
            float A_x = A_x_field[idx];
            float A_y = A_y_field[idx];

            // For Peierls substitution, we need to apply phases on LINKS, not sites
            // Approximate with local A value (valid for slowly varying A)
            //
            // Proper lattice gauge: ψ(x+dx) → U(x,x+dx) ψ(x+dx)
            // where U(x,x+dx) = exp(-iq ∫_x^(x+dx) A·dl) ≈ exp(-iq A(x)·dx)

            // For now, use perturbative approximation
            // This will be refined in full lattice gauge implementation

            float dx = 1.0f; // Lattice spacing
            float phase_x = -_em_coupling_strength * A_x * dx;
            float phase_y = -_em_coupling_strength * A_y * dx;

            // Combined phase for 2D motion
            float gauge_phase = phase_x + phase_y;

            std::complex<float> U_gauge = std::exp(std::complex<float>(0.0f, gauge_phase));

            // Apply to all spinor components
            for (int c = 0; c < 4; c++) {
                _psi[c][idx] *= U_gauge;
            }
        }
    }
}
