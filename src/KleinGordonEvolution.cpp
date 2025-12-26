/**
 * KleinGordonEvolution.cpp
 *
 * Implementation of split-operator method for 2D Klein-Gordon equation
 * Using FFTW for fast Fourier transforms
 *
 * Klein-Gordon equation (second-order in time):
 *   (∂²_t - ∇² + m²)φ = 0
 *
 * Converted to first-order system:
 *   ∂_tφ = φ_dot
 *   ∂_tφ_dot = ∇²φ - m²φ
 *
 * Split-operator method (Strang splitting):
 *   1. Half-step kinetic: evolve gradient term ∇²
 *   2. Full-step mass: evolve mass term m²
 *   3. Half-step kinetic: evolve gradient term ∇²
 */

#include "KleinGordonEvolution.h"
#include <fftw3.h>
#include <cmath>
#include <iostream>

KleinGordonEvolution::KleinGordonEvolution(uint32_t Nx, uint32_t Ny)
    : _Nx(Nx), _Ny(Ny), _N_points(Nx * Ny), _phi_k_valid(false) {

    // Allocate scalar field and time derivative
    _phi.resize(_N_points);
    _phi_dot.resize(_N_points);
    _phi_k.resize(_N_points);
    _phi_dot_k.resize(_N_points);

    _fft_forward_phi = nullptr;
    _fft_backward_phi = nullptr;
    _fft_forward_phi_dot = nullptr;
    _fft_backward_phi_dot = nullptr;

    setupMomentumGrid();
    setupFFT();

    std::cout << "[KleinGordonEvolution] Initialized " << _Nx << "x" << _Ny
              << " scalar field" << std::endl;
}

KleinGordonEvolution::~KleinGordonEvolution() {
    cleanupFFT();
}

void KleinGordonEvolution::setupMomentumGrid() {
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

void KleinGordonEvolution::setupFFT() {
    // Create FFTW plans for φ and φ_dot
    _fft_forward_phi = fftwf_plan_dft_2d(
        _Ny, _Nx,
        reinterpret_cast<fftwf_complex*>(_phi.data()),
        reinterpret_cast<fftwf_complex*>(_phi_k.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    _fft_backward_phi = fftwf_plan_dft_2d(
        _Ny, _Nx,
        reinterpret_cast<fftwf_complex*>(_phi_k.data()),
        reinterpret_cast<fftwf_complex*>(_phi.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    _fft_forward_phi_dot = fftwf_plan_dft_2d(
        _Ny, _Nx,
        reinterpret_cast<fftwf_complex*>(_phi_dot.data()),
        reinterpret_cast<fftwf_complex*>(_phi_dot_k.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    _fft_backward_phi_dot = fftwf_plan_dft_2d(
        _Ny, _Nx,
        reinterpret_cast<fftwf_complex*>(_phi_dot_k.data()),
        reinterpret_cast<fftwf_complex*>(_phi_dot.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );
}

void KleinGordonEvolution::cleanupFFT() {
    if (_fft_forward_phi) {
        fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_forward_phi));
    }
    if (_fft_backward_phi) {
        fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_backward_phi));
    }
    if (_fft_forward_phi_dot) {
        fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_forward_phi_dot));
    }
    if (_fft_backward_phi_dot) {
        fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_backward_phi_dot));
    }
}

void KleinGordonEvolution::initialize(float x0, float y0, float sigma) {
    // Initialize as Gaussian wavepacket
    // Start with φ_dot = 0 (initially at rest)
    float norm_sum = 0.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            float dx = x - x0;
            float dy = y - y0;
            float r2 = dx*dx + dy*dy;
            float amplitude = std::exp(-r2 / (2.0f * sigma * sigma));

            _phi[idx] = std::complex<float>(amplitude, 0.0f);
            _phi_dot[idx] = std::complex<float>(0.0f, 0.0f);  // Initially at rest

            norm_sum += amplitude * amplitude;
        }
    }

    // Normalize
    float norm_factor = std::sqrt(norm_sum);
    for (uint32_t i = 0; i < _N_points; i++) {
        _phi[i] /= norm_factor;
    }

    // Invalidate k-space cache after initialization
    _phi_k_valid = false;

    std::cout << "[KleinGordonEvolution] Initialized Gaussian wavepacket at ("
              << x0 << ", " << y0 << ") with σ=" << sigma << std::endl;
}

void KleinGordonEvolution::initializePlaneWave(float kx, float ky) {
    float k_mag = std::sqrt(kx*kx + ky*ky);

    // For Klein-Gordon plane wave: φ = A·e^(i(k·r - ωt))
    // where ω = √(k² + m²)
    // At t=0, just spatial part: φ = A·e^(ik·r)
    // φ_dot = -iω·φ at t=0

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            // Plane wave phase e^(i k.r)
            float phase_arg = kx * x + ky * y; // dx=1
            std::complex<float> phase(std::cos(phase_arg), std::sin(phase_arg));

            _phi[idx] = phase;

            // For free particle (m=0): ω = |k|
            // φ_dot = -iω·φ
            float omega = k_mag;
            _phi_dot[idx] = std::complex<float>(0.0f, -omega) * phase;
        }
    }

    // Normalize
    float norm_sum = 0.0f;
    for (uint32_t i = 0; i < _N_points; i++) {
        norm_sum += std::norm(_phi[i]);
    }

    float norm_factor = 1.0f / std::sqrt(norm_sum);
    for (uint32_t i = 0; i < _N_points; i++) {
        _phi[i] *= norm_factor;
        _phi_dot[i] *= norm_factor;
    }

    _phi_k_valid = false;
    std::cout << "[KleinGordonEvolution] Initialized Plane Wave k=(" << kx << ", " << ky << ")" << std::endl;
}

void KleinGordonEvolution::step(const std::vector<float>& mass_field, float dt) {
    // Strang splitting: K/2 - M - K/2
    // where K = kinetic (gradient) operator, M = mass operator

    applyKineticHalfStep(dt / 2.0f);
    applyMassStep(mass_field, dt);
    applyKineticHalfStep(dt / 2.0f);

    // Invalidate k-space cache after evolution
    _phi_k_valid = false;
}

void KleinGordonEvolution::applyMassStep(const std::vector<float>& mass_field, float dt) {
    /**
     * Mass operator: V = m²(x,y)·φ
     *
     * From Klein-Gordon: ∂²_tφ = ∇²φ - m²φ
     * As first-order system:
     *   ∂_tφ = φ_dot
     *   ∂_tφ_dot = ∇²φ - m²φ
     *
     * Mass step evolution (exact for m² term):
     *   φ(t+dt) = φ(t)·cos(m·dt) + (φ_dot(t)/m)·sin(m·dt)
     *   φ_dot(t+dt) = φ_dot(t)·cos(m·dt) - φ(t)·m·sin(m·dt)
     *
     * This is rotation in (φ, φ_dot/m) phase space with frequency m
     */

    for (uint32_t i = 0; i < _N_points; i++) {
        float m = mass_field[i];

        if (m < 1e-8f) {
            // Zero mass: no evolution
            continue;
        }

        float m_dt = m * dt;
        float cos_m_dt = std::cos(m_dt);
        float sin_m_dt = std::sin(m_dt);

        std::complex<float> phi_old = _phi[i];
        std::complex<float> phi_dot_old = _phi_dot[i];

        // Coupled evolution
        _phi[i] = phi_old * cos_m_dt + (phi_dot_old / m) * sin_m_dt;
        _phi_dot[i] = phi_dot_old * cos_m_dt - phi_old * m * sin_m_dt;
    }
}

void KleinGordonEvolution::applyKineticHalfStep(float dt_half) {
    /**
     * Kinetic operator: K = ∇² (Laplacian)
     *
     * From Klein-Gordon first-order system:
     *   ∂_tφ = φ_dot
     *   ∂_tφ_dot = ∇²φ - m²φ
     *
     * Kinetic step (gradient term only, m=0):
     *   ∂_tφ = φ_dot
     *   ∂_tφ_dot = ∇²φ
     *
     * In momentum space: ∇² → -k²
     *   ∂_tφ_k = φ_dot_k
     *   ∂_tφ_dot_k = -k²·φ_k
     *
     * This is a harmonic oscillator with frequency ω_k = |k|
     * Exact solution:
     *   φ_k(t+dt) = φ_k(t)·cos(|k|dt) + (φ_dot_k(t)/|k|)·sin(|k|dt)
     *   φ_dot_k(t+dt) = φ_dot_k(t)·cos(|k|dt) - φ_k(t)·|k|·sin(|k|dt)
     */

    // Forward FFT for both φ and φ_dot
    fftwf_execute(static_cast<fftwf_plan>(_fft_forward_phi));
    fftwf_execute(static_cast<fftwf_plan>(_fft_forward_phi_dot));

    // Apply kinetic operator in momentum space
    for (uint32_t j = 0; j < _Ny; j++) {
        for (uint32_t i = 0; i < _Nx; i++) {
            uint32_t idx = j * _Nx + i;

            float kx = _kx[i];
            float ky = _ky[j];
            float k_mag = std::sqrt(kx*kx + ky*ky);

            // Apply Klein-Gordon kinetic evolution
            applyKGKineticMatrix(_phi_k[idx], _phi_dot_k[idx], k_mag, dt_half);
        }
    }

    // Inverse FFT for both φ and φ_dot
    fftwf_execute(static_cast<fftwf_plan>(_fft_backward_phi));
    fftwf_execute(static_cast<fftwf_plan>(_fft_backward_phi_dot));

    // FFTW doesn't normalize, so divide by N
    float norm_factor = 1.0f / _N_points;
    for (uint32_t i = 0; i < _N_points; i++) {
        _phi[i] *= norm_factor;
        _phi_dot[i] *= norm_factor;
    }
}

void KleinGordonEvolution::applyKGKineticMatrix(std::complex<float>& phi_k,
                                                 std::complex<float>& phi_dot_k,
                                                 float k_mag, float dt) {
    /**
     * Exact evolution under kinetic operator in momentum space
     *
     * Hamiltonian in k-space: H_k = |k|² (free wave dispersion)
     * Equations of motion:
     *   dφ_k/dt = φ_dot_k
     *   dφ_dot_k/dt = -|k|²·φ_k
     *
     * This is a 2D harmonic oscillator with frequency ω = |k|
     * Exact solution (rotation in phase space):
     *   [φ_k(t)]     = [ cos(ωt)    sin(ωt)/ω ] [φ_k(0)    ]
     *   [φ_dot_k(t)]   [-ω·sin(ωt)  cos(ωt)   ] [φ_dot_k(0)]
     */

    if (k_mag < 1e-10f) {
        // k=0: no kinetic evolution (constant mode)
        return;
    }

    float omega = k_mag;
    float cos_omega_dt = std::cos(omega * dt);
    float sin_omega_dt = std::sin(omega * dt);

    std::complex<float> phi_old = phi_k;
    std::complex<float> phi_dot_old = phi_dot_k;

    // Coupled evolution (exact for free Klein-Gordon)
    phi_k = phi_old * cos_omega_dt + (phi_dot_old / omega) * sin_omega_dt;
    phi_dot_k = phi_dot_old * cos_omega_dt - phi_old * omega * sin_omega_dt;
}

std::vector<float> KleinGordonEvolution::getDensity() const {
    std::vector<float> density(_N_points);

    for (uint32_t i = 0; i < _N_points; i++) {
        density[i] = std::norm(_phi[i]);
    }

    return density;
}

std::vector<std::complex<double>> KleinGordonEvolution::getScalarField() const {
    // Return scalar field as flat array [Nx * Ny]
    std::vector<std::complex<double>> scalar_flat(_N_points);

    for (uint32_t i = 0; i < _N_points; i++) {
        scalar_flat[i] = std::complex<double>(_phi[i].real(), _phi[i].imag());
    }

    return scalar_flat;
}

float KleinGordonEvolution::getNorm() const {
    float norm = 0.0f;
    for (uint32_t i = 0; i < _N_points; i++) {
        norm += std::norm(_phi[i]);
    }
    return norm;
}

void KleinGordonEvolution::getCenterOfMass(float& x_mean, float& y_mean) const {
    /**
     * Calculate center of mass: (x̄, ȳ) = ∫ (x,y)|φ|² dxdy / ∫|φ|² dxdy
     */
    x_mean = 0.0f;
    y_mean = 0.0f;
    float total = 0.0f;

    for (uint32_t y = 0; y < _Ny; y++) {
        for (uint32_t x = 0; x < _Nx; x++) {
            uint32_t idx = y * _Nx + x;

            float density = std::norm(_phi[idx]);

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

void KleinGordonEvolution::getMomentumDistribution(std::vector<float>& kx_out,
                                                     std::vector<float>& ky_out,
                                                     std::vector<float>& density_k) const {
    /**
     * Compute momentum-space density |φ̃(k)|²
     *
     * This is needed for:
     * - Dispersion relation validation E(k)
     * - Wave packet analysis
     * - Group velocity measurements
     */

    kx_out = _kx;
    ky_out = _ky;
    density_k.resize(_N_points, 0.0f);

    // Compute |φ̃(k)|²
    for (uint32_t i = 0; i < _N_points; i++) {
        density_k[i] = std::norm(_phi_k[i]);
    }
}

/**
 * Compute total energy E = <H>
 *
 * Klein-Gordon Hamiltonian density:
 *   ℋ = |∂_tφ|² + |∇φ|² + m²|φ|²
 *
 * Total energy:
 *   E = ∫ ℋ dx = ∫ [|φ_dot|² + |∇φ|² + m²|φ|²] dx
 *
 * Decomposition:
 *   Kinetic:  E_K = ∫ |φ_dot|² dx
 *   Gradient: E_G = ∫ |∇φ|² dx = Σ_k k²|φ̃_k|² (in momentum space)
 *   Mass:     E_M = ∫ m²|φ|² dx
 *
 * Output:
 *   KE_out = E_K + E_G (total kinetic + gradient energy)
 *   PE_out = E_M (mass potential energy)
 */
float KleinGordonEvolution::getEnergy(const std::vector<float>& mass_field,
                                       float& KE_out, float& PE_out) const {

    // === Part 1: Kinetic Energy (time derivative term) ===
    float E_time_derivative = 0.0f;
    for (uint32_t i = 0; i < _N_points; i++) {
        E_time_derivative += std::norm(_phi_dot[i]);
    }

    // === Part 2: Gradient Energy (in momentum space) ===
    // E_gradient = Σ_k k²|φ̃_k|²

    // Update k-space if cache is invalid
    if (!_phi_k_valid) {
        // Compute fresh k-space representation via FFT
        fftwf_execute((fftwf_plan)_fft_forward_phi);
        _phi_k_valid = true;
    }

    float E_gradient = 0.0f;
    for (uint32_t j = 0; j < _Ny; j++) {
        for (uint32_t i = 0; i < _Nx; i++) {
            uint32_t idx = j * _Nx + i;

            float kx = _kx[i];
            float ky = _ky[j];
            float k2 = kx*kx + ky*ky;

            float density_k = std::norm(_phi_k[idx]);
            E_gradient += k2 * density_k;
        }
    }

    // Normalize by grid size (FFT normalization)
    E_gradient /= (_Nx * _Ny);

    // === Part 3: Mass Potential Energy ===
    float E_mass = 0.0f;
    for (uint32_t i = 0; i < _N_points; i++) {
        float m = mass_field[i];
        float density = std::norm(_phi[i]);
        E_mass += m * m * density;
    }

    // === Total Energy ===
    float KE = E_time_derivative + E_gradient;  // Kinetic + gradient
    float PE = E_mass;                           // Mass potential

    float E_total = KE + PE;

    // Output parameters
    KE_out = KE;
    PE_out = PE;

    return E_total;
}
