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

DiracEvolution::DiracEvolution(uint32_t Nx, uint32_t Ny)
    : _Nx(Nx), _Ny(Ny), _N_points(Nx * Ny), _psi_k_valid(false) {

    // Allocate spinor components
    for (int c = 0; c < 4; c++) {
        _psi[c].resize(_N_points);
        _psi_k[c].resize(_N_points);
        _fft_forward[c] = nullptr;
        _fft_backward[c] = nullptr;
    }

    setupMomentumGrid();
    setupFFT();

    std::cout << "[DiracEvolution] Initialized " << _Nx << "x" << _Ny
              << " 4-component spinor field" << std::endl;
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

void DiracEvolution::step(const std::vector<float>& mass_field, float dt,
                         const std::vector<float>& A_x,
                         const std::vector<float>& A_y,
                         const std::vector<float>& A_z,
                         const std::vector<float>& phi) {
    // Strang splitting: K/2 - V - K/2
    // With EM fields: kinetic operator includes minimal coupling -iα·(∇ - ieA)
    // Potential step includes both mass term βm(x) and scalar potential eφ

    applyKineticHalfStep(dt / 2.0f, A_x, A_y, A_z);
    applyPotentialStep(mass_field, dt, phi);
    applyKineticHalfStep(dt / 2.0f, A_x, A_y, A_z);

    // Invalidate k-space cache after evolution
    _psi_k_valid = false;
}

void DiracEvolution::applyPotentialStep(const std::vector<float>& mass_field, float dt,
                                        const std::vector<float>& phi) {
    /**
     * Potential step: V = β·m(x,y) + e·φ(x,y)
     *
     * β = diag(1, 1, -1, -1) in standard representation
     * The scalar potential φ contributes equally to all components (charge coupling)
     *
     * Evolution operator: exp(-i(βm + eφ)Δt)
     *
     * For upper components (β=+1): exp(-i(m + eφ)Δt)
     * For lower components (β=-1): exp(-i(-m + eφ)Δt) = exp(i(m - eφ)Δt)
     *
     * Note: We use e=1 (natural units) for the elementary charge
     */

    // Check if scalar potential is provided
    bool has_phi = !phi.empty();
    const float e = 1.0f;  // Elementary charge in natural units

    for (uint32_t i = 0; i < _N_points; i++) {
        float m = mass_field[i];
        float phi_val = has_phi ? phi[i] : 0.0f;

        // Upper components: β = +1
        std::complex<float> phase_upper(0.0f, -(m + e * phi_val) * dt);
        _psi[0][i] *= std::exp(phase_upper);
        _psi[1][i] *= std::exp(phase_upper);

        // Lower components: β = -1
        std::complex<float> phase_lower(0.0f, (m - e * phi_val) * dt);
        _psi[2][i] *= std::exp(phase_lower);
        _psi[3][i] *= std::exp(phase_lower);
    }
}

void DiracEvolution::applyKineticHalfStep(float dt_half,
                                          const std::vector<float>& A_x,
                                          const std::vector<float>& A_y,
                                          const std::vector<float>& A_z) {
    /**
     * Kinetic operator with minimal coupling: K = α·π where π = -i∇ - eA
     *
     * Split-operator approach for gauge-invariant evolution:
     * 1. Apply gauge phase exp(ieA·x) in position space (Peierls substitution)
     * 2. FFT to momentum space
     * 3. Apply free kinetic evolution exp(-iα·k·t)
     * 4. Inverse FFT to position space
     * 5. Apply inverse gauge phase exp(-ieA·x)
     *
     * This maintains gauge covariance and numerical stability.
     *
     * Note: For 2D system, we only use A_x and A_y. A_z is ignored if provided.
     * Note: We use e=1 (natural units) for the elementary charge
     */

    const float e = 1.0f;  // Elementary charge in natural units
    bool has_A = !A_x.empty() || !A_y.empty();

    // Step 1: Apply gauge coupling via local phase transformation
    // The proper implementation for minimal coupling in split-operator method
    // uses exp(ieA·r) as a gauge transformation
    if (has_A) {
        // Apply gauge-covariant derivative transformation
        // This implements π = p - eA via local U(1) gauge transformation
        for (uint32_t j = 0; j < _Ny; j++) {
            for (uint32_t i = 0; i < _Nx; i++) {
                uint32_t idx = j * _Nx + i;

                // Get vector potential at this point
                float Ax = !A_x.empty() ? A_x[idx] : 0.0f;
                float Ay = !A_y.empty() ? A_y[idx] : 0.0f;

                // Apply gauge transformation that converts ∇ → ∇ - ieA
                // The phase factor exp(ieA·r·dt) implements minimal coupling
                // We use position relative to center to avoid boundary issues
                float x_rel = (float)i - _Nx/2.0f;
                float y_rel = (float)j - _Ny/2.0f;

                float phase_arg = e * (Ax * x_rel + Ay * y_rel) * dt_half;
                std::complex<float> gauge_phase = std::exp(std::complex<float>(0.0f, phase_arg));

                // Apply to all spinor components
                for (int c = 0; c < 4; c++) {
                    _psi[c][idx] *= gauge_phase;
                }
            }
        }
    }

    // Step 2: Forward FFT for all 4 components
    for (int c = 0; c < 4; c++) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_forward[c]));
    }

    // Step 3: Apply free kinetic operator in momentum space
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

            // Apply Dirac kinetic matrix (free evolution)
            applyDiracKineticMatrix(psi_k_point, k_mag, kx, ky, dt_half);

            // Store result
            for (int c = 0; c < 4; c++) {
                _psi_k[c][idx] = psi_k_point[c];
            }
        }
    }

    // Step 4: Inverse FFT for all 4 components
    for (int c = 0; c < 4; c++) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_backward[c]));

        // FFTW doesn't normalize, so divide by N
        float norm_factor = 1.0f / _N_points;
        for (uint32_t i = 0; i < _N_points; i++) {
            _psi[c][i] *= norm_factor;
        }
    }

    // Step 5: Apply inverse gauge transformation
    if (has_A) {
        for (uint32_t j = 0; j < _Ny; j++) {
            for (uint32_t i = 0; i < _Nx; i++) {
                uint32_t idx = j * _Nx + i;

                // Get vector potential at this point
                float Ax = !A_x.empty() ? A_x[idx] : 0.0f;
                float Ay = !A_y.empty() ? A_y[idx] : 0.0f;

                // Apply inverse gauge transformation
                float x_rel = (float)i - _Nx/2.0f;
                float y_rel = (float)j - _Ny/2.0f;

                float phase_arg = -e * (Ax * x_rel + Ay * y_rel) * dt_half;
                std::complex<float> inv_gauge_phase = std::exp(std::complex<float>(0.0f, phase_arg));

                // Apply to all spinor components
                for (int c = 0; c < 4; c++) {
                    _psi[c][idx] *= inv_gauge_phase;
                }
            }
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
 * Hamiltonian with EM fields: H = -iα·(∇ - ieA) + βm(x) + eφ
 * Without EM fields: H = -iα·∇ + βm(x)
 *
 * Decomposition:
 *   Kinetic:   E_K = ∫ Ψ*(x) [-iα·(∇ - ieA)] Ψ(x) dx
 *            = Σ_k |Ψ̃_k|² ω(k)  (in gauge where A=0, or after gauge transformation)
 *
 *   Potential: E_V = ∫ |Ψ(x)|² [β m(x) + eφ(x)] dx
 *
 * For Dirac equation in 2D with 4-component spinor:
 *   β = diag(1, 1, -1, -1)
 *
 * Note: m_0 (rest mass) is assumed ~0 in this implementation (massless Dirac limit)
 *       The mass term m(x) comes from SMFT coupling, not intrinsic mass.
 *       The EM interaction energy includes both A·j (current coupling) and φ·ρ (charge coupling).
 */
float DiracEvolution::getEnergy(const std::vector<float>& mass_field,
                                float& KE_out, float& PE_out,
                                const std::vector<float>& A_x,
                                const std::vector<float>& A_y,
                                const std::vector<float>& A_z,
                                const std::vector<float>& phi) const {

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
    // E_V = Σ_x |Ψ(x)|² [β m(x) + eφ(x)]
    //
    // β = diag(1, 1, -1, -1) → upper components contribute +m, lower contribute -m
    // The scalar potential φ contributes to all components with same sign (charge density)

    float PE = 0.0f;
    const float e = 1.0f;  // Elementary charge in natural units
    bool has_phi = !phi.empty();

    for (uint32_t idx = 0; idx < _N_points; idx++) {
        float m = mass_field[idx];
        float phi_val = has_phi ? phi[idx] : 0.0f;

        // Upper components (c=0,1): β = +1
        float density_upper = std::norm(_psi[0][idx]) + std::norm(_psi[1][idx]);
        PE += density_upper * (m + e * phi_val);

        // Lower components (c=2,3): β = -1
        float density_lower = std::norm(_psi[2][idx]) + std::norm(_psi[3][idx]);
        PE += density_lower * (-m + e * phi_val);
    }

    // === Part 2b: Magnetic Energy (from vector potential) ===
    //
    // The interaction energy with the vector potential A contributes through
    // the current density: E_mag = ∫ j·A dx where j = eΨ†αΨ
    //
    // For numerical stability, we compute this in position space directly
    // using the Dirac current j^μ = Ψ†γ^μΨ where γ^i = α^i in our convention

    if (!A_x.empty() || !A_y.empty()) {
        for (uint32_t idx = 0; idx < _N_points; idx++) {
            // Compute current density components j_x, j_y
            // j_x = e·Ψ†α_xΨ where α_x = σ_x ⊗ I (acts on pairs)
            // j_y = e·Ψ†α_yΨ where α_y = σ_y ⊗ I

            // For α_x (σ_x couples components 0↔1 and 2↔3)
            std::complex<float> j_x = e * (
                std::conj(_psi[0][idx]) * _psi[1][idx] +
                std::conj(_psi[1][idx]) * _psi[0][idx] +
                std::conj(_psi[2][idx]) * _psi[3][idx] +
                std::conj(_psi[3][idx]) * _psi[2][idx]
            );

            // For α_y (σ_y couples components 0↔1 and 2↔3 with i factor)
            std::complex<float> j_y = e * std::complex<float>(0.0f, 1.0f) * (
                -std::conj(_psi[0][idx]) * _psi[1][idx] +
                std::conj(_psi[1][idx]) * _psi[0][idx] +
                -std::conj(_psi[2][idx]) * _psi[3][idx] +
                std::conj(_psi[3][idx]) * _psi[2][idx]
            );

            // Interaction energy: j·A
            float Ax = !A_x.empty() ? A_x[idx] : 0.0f;
            float Ay = !A_y.empty() ? A_y[idx] : 0.0f;

            PE += j_x.real() * Ax + j_y.real() * Ay;
        }
    }

    // === Total Energy ===
    float E_total = KE + PE;

    // Output parameters
    KE_out = KE;
    PE_out = PE;

    return E_total;
}
