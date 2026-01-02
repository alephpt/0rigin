/**
 * Dirac3D.cpp
 *
 * Implementation of 3D Dirac equation solver
 */

#include "Dirac3D.h"
#include <fftw3.h>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// Dirac alpha matrices in standard (Dirac) representation
// α^x, α^y, α^z are 4×4 matrices
// Stored in row-major order: M[i,j] = data[i*4 + j]

// α^x = [ 0   0   0   1 ]
//       [ 0   0   1   0 ]
//       [ 0   1   0   0 ]
//       [ 1   0   0   0 ]
const std::array<std::complex<float>, 16> Dirac3D::alpha_x = {{
    {0,0}, {0,0}, {0,0}, {1,0},
    {0,0}, {0,0}, {1,0}, {0,0},
    {0,0}, {1,0}, {0,0}, {0,0},
    {1,0}, {0,0}, {0,0}, {0,0}
}};

// α^y = [ 0   0   0  -i ]
//       [ 0   0   i   0 ]
//       [ 0  -i   0   0 ]
//       [ i   0   0   0 ]
const std::array<std::complex<float>, 16> Dirac3D::alpha_y = {{
    {0,0}, {0,0}, {0,0}, {0,-1},
    {0,0}, {0,0}, {0,1}, {0,0},
    {0,0}, {0,-1}, {0,0}, {0,0},
    {0,1}, {0,0}, {0,0}, {0,0}
}};

// α^z = [ 0   0   1   0 ]
//       [ 0   0   0  -1 ]
//       [ 1   0   0   0 ]
//       [ 0  -1   0   0 ]
const std::array<std::complex<float>, 16> Dirac3D::alpha_z = {{
    {0,0}, {0,0}, {1,0}, {0,0},
    {0,0}, {0,0}, {0,0}, {-1,0},
    {1,0}, {0,0}, {0,0}, {0,0},
    {0,0}, {-1,0}, {0,0}, {0,0}
}};

// β = [ 1   0   0   0 ]
//     [ 0   1   0   0 ]
//     [ 0   0  -1   0 ]
//     [ 0   0   0  -1 ]
const std::array<std::complex<float>, 16> Dirac3D::beta = {{
    {1,0}, {0,0}, {0,0}, {0,0},
    {0,0}, {1,0}, {0,0}, {0,0},
    {0,0}, {0,0}, {-1,0}, {0,0},
    {0,0}, {0,0}, {0,0}, {-1,0}
}};

Dirac3D::Dirac3D(uint32_t Nx, uint32_t Ny, uint32_t Nz)
    : _Nx(Nx), _Ny(Ny), _Nz(Nz), _N_total(Nx * Ny * Nz), _dx(1.0f)
{
    // Allocate spinor storage (4 components)
    for (int c = 0; c < 4; ++c) {
        _psi[c].resize(_N_total, {0.0f, 0.0f});
        _psi_k[c].resize(_N_total, {0.0f, 0.0f});
        _fft_forward[c] = nullptr;
        _fft_backward[c] = nullptr;
    }

    setupMomentumGrid();
    setupFFT();
}

Dirac3D::~Dirac3D() {
    cleanupFFT();
}

void Dirac3D::setupMomentumGrid() {
    // Create momentum grid for FFT
    // k_i = 2π n_i / (N_i * dx)
    // where n_i ∈ [-N_i/2, N_i/2)

    _kx.resize(_Nx);
    _ky.resize(_Ny);
    _kz.resize(_Nz);

    const float dk_x = 2.0f * M_PI / (_Nx * _dx);
    const float dk_y = 2.0f * M_PI / (_Ny * _dx);
    const float dk_z = 2.0f * M_PI / (_Nz * _dx);

    for (uint32_t i = 0; i < _Nx; ++i) {
        int32_t n = (i < _Nx/2) ? i : i - _Nx;
        _kx[i] = n * dk_x;
    }

    for (uint32_t i = 0; i < _Ny; ++i) {
        int32_t n = (i < _Ny/2) ? i : i - _Ny;
        _ky[i] = n * dk_y;
    }

    for (uint32_t i = 0; i < _Nz; ++i) {
        int32_t n = (i < _Nz/2) ? i : i - _Nz;
        _kz[i] = n * dk_z;
    }
}

void Dirac3D::setupFFT() {
    // Create FFTW plans for 3D FFT (all 4 spinor components)
    const int N[3] = {static_cast<int>(_Nz), static_cast<int>(_Ny), static_cast<int>(_Nx)};

    for (int c = 0; c < 4; ++c) {
        _fft_forward[c] = fftwf_plan_dft_3d(
            _Nz, _Ny, _Nx,
            reinterpret_cast<fftwf_complex*>(_psi[c].data()),
            reinterpret_cast<fftwf_complex*>(_psi_k[c].data()),
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );

        _fft_backward[c] = fftwf_plan_dft_3d(
            _Nz, _Ny, _Nx,
            reinterpret_cast<fftwf_complex*>(_psi_k[c].data()),
            reinterpret_cast<fftwf_complex*>(_psi[c].data()),
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
    }
}

void Dirac3D::cleanupFFT() {
    for (int c = 0; c < 4; ++c) {
        if (_fft_forward[c]) {
            fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_forward[c]));
            _fft_forward[c] = nullptr;
        }
        if (_fft_backward[c]) {
            fftwf_destroy_plan(static_cast<fftwf_plan>(_fft_backward[c]));
            _fft_backward[c] = nullptr;
        }
    }
}

void Dirac3D::initialize(const std::vector<std::complex<float>>& psi_init) {
    if (psi_init.size() != 4 * _N_total) {
        throw std::runtime_error("Dirac3D::initialize: spinor must have 4 × N_total components");
    }

    for (int c = 0; c < 4; ++c) {
        for (uint32_t i = 0; i < _N_total; ++i) {
            _psi[c][i] = psi_init[c * _N_total + i];
        }
    }
}

void Dirac3D::initializeGaussian(float x0, float y0, float z0, float sigma) {
    // Create Gaussian wavepacket centered at (x0, y0, z0)
    // Initialize in upper components (ψ₁, ψ₂) with equal amplitudes

    const float cx = _Nx / 2.0f;
    const float cy = _Ny / 2.0f;
    const float cz = _Nz / 2.0f;

    float norm_sum = 0.0f;

    for (uint32_t iz = 0; iz < _Nz; ++iz) {
        for (uint32_t iy = 0; iy < _Ny; ++iy) {
            for (uint32_t ix = 0; ix < _Nx; ++ix) {
                const uint32_t idx = index3D(ix, iy, iz);

                const float x = ix - cx + x0;
                const float y = iy - cy + y0;
                const float z = iz - cz + z0;
                const float r2 = x*x + y*y + z*z;

                const float psi_val = std::exp(-r2 / (2.0f * sigma * sigma));
                norm_sum += psi_val * psi_val;

                // Upper components (particle-like)
                _psi[0][idx] = {psi_val / std::sqrt(2.0f), 0.0f};
                _psi[1][idx] = {psi_val / std::sqrt(2.0f), 0.0f};

                // Lower components (antiparticle-like) - initially zero
                _psi[2][idx] = {0.0f, 0.0f};
                _psi[3][idx] = {0.0f, 0.0f};
            }
        }
    }

    // Normalize
    const float norm = std::sqrt(norm_sum * _dx * _dx * _dx);
    for (int c = 0; c < 4; ++c) {
        for (uint32_t i = 0; i < _N_total; ++i) {
            _psi[c][i] /= norm;
        }
    }
}

void Dirac3D::applyDiracKineticMatrix(std::complex<float> psi_k[4],
                                     float kx, float ky, float kz,
                                     float dt)
{
    // Apply exp(-i α·k dt) to spinor in momentum space
    // H_kinetic = α·k has eigenvalues ±|k|
    // Use exact exponential: exp(-i α·k dt) = cos(|k| dt) - i (α·k/|k|) sin(|k| dt)

    const float k_mag = std::sqrt(kx*kx + ky*ky + kz*kz);

    // Handle zero momentum case
    if (k_mag < 1e-10f) {
        return; // No evolution at k=0
    }

    const float k_inv = 1.0f / k_mag;
    const float phase = k_mag * dt;
    const float cos_phase = std::cos(phase);
    const float sin_phase = std::sin(phase);

    // Normalized direction: k̂ = k/|k|
    const float kx_norm = kx * k_inv;
    const float ky_norm = ky * k_inv;
    const float kz_norm = kz * k_inv;

    // Compute (α·k̂) × psi
    // α·k̂ = α_x k̂_x + α_y k̂_y + α_z k̂_z
    std::complex<float> alpha_k_psi[4] = {{0,0}, {0,0}, {0,0}, {0,0}};

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const std::complex<float> alpha_k_ij =
                alpha_x[i*4 + j] * kx_norm +
                alpha_y[i*4 + j] * ky_norm +
                alpha_z[i*4 + j] * kz_norm;
            alpha_k_psi[i] += alpha_k_ij * psi_k[j];
        }
    }

    // Apply: psi_new = cos(phase) * psi - i * sin(phase) * (α·k̂) * psi
    std::complex<float> psi_new[4];
    const std::complex<float> i_unit(0.0f, 1.0f);

    for (int i = 0; i < 4; ++i) {
        psi_new[i] = cos_phase * psi_k[i] - i_unit * sin_phase * alpha_k_psi[i];
    }

    for (int i = 0; i < 4; ++i) {
        psi_k[i] = psi_new[i];
    }
}

void Dirac3D::applyKineticHalfStep(float dt_half) {
    // Transform to momentum space
    for (int c = 0; c < 4; ++c) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_forward[c]));
    }

    // Apply kinetic operator at each k-point
    for (uint32_t iz = 0; iz < _Nz; ++iz) {
        for (uint32_t iy = 0; iy < _Ny; ++iy) {
            for (uint32_t ix = 0; ix < _Nx; ++ix) {
                const uint32_t idx = index3D(ix, iy, iz);

                std::complex<float> psi_k_local[4];
                for (int c = 0; c < 4; ++c) {
                    psi_k_local[c] = _psi_k[c][idx];
                }

                applyDiracKineticMatrix(psi_k_local, _kx[ix], _ky[iy], _kz[iz], dt_half);

                for (int c = 0; c < 4; ++c) {
                    _psi_k[c][idx] = psi_k_local[c];
                }
            }
        }
    }

    // Transform back to position space
    for (int c = 0; c < 4; ++c) {
        fftwf_execute(static_cast<fftwf_plan>(_fft_backward[c]));
    }

    // Normalize after inverse FFT
    const float norm_factor = 1.0f / _N_total;
    for (int c = 0; c < 4; ++c) {
        for (uint32_t i = 0; i < _N_total; ++i) {
            _psi[c][i] *= norm_factor;
        }
    }
}

void Dirac3D::applyMassStep(const std::vector<float>& mass_field, float dt) {
    if (mass_field.size() != _N_total) {
        throw std::runtime_error("Dirac3D::applyMassStep: mass_field size mismatch");
    }

    // Apply exp(-i β m dt) at each spatial point
    // β is diagonal, so this is simple

    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        const float m = mass_field[idx];
        const float phase = m * dt;

        const std::complex<float> exp_plus = std::exp(std::complex<float>(0.0f, -phase));
        const std::complex<float> exp_minus = std::exp(std::complex<float>(0.0f, phase));

        // β = diag(1, 1, -1, -1)
        // Upper components: multiply by exp(-i m dt)
        _psi[0][idx] *= exp_plus;
        _psi[1][idx] *= exp_plus;

        // Lower components: multiply by exp(+i m dt)
        _psi[2][idx] *= exp_minus;
        _psi[3][idx] *= exp_minus;
    }
}

void Dirac3D::step(const std::vector<float>& mass_field, float dt) {
    // Strang splitting:
    // 1. Apply kinetic half-step
    // 2. Apply mass full-step
    // 3. Apply kinetic half-step

    applyKineticHalfStep(dt / 2.0f);
    applyMassStep(mass_field, dt);
    applyKineticHalfStep(dt / 2.0f);
}

std::vector<float> Dirac3D::getDensity() const {
    std::vector<float> density(_N_total);

    for (uint32_t i = 0; i < _N_total; ++i) {
        float rho = 0.0f;
        for (int c = 0; c < 4; ++c) {
            rho += std::norm(_psi[c][i]);
        }
        density[i] = rho;
    }

    return density;
}

std::vector<float> Dirac3D::getCurrent(int component) const {
    std::vector<float> current(_N_total);

    // j^i = ψ† α^i ψ
    const std::array<std::complex<float>, 16>* alpha = nullptr;

    if (component == 0) alpha = &alpha_x;
    else if (component == 1) alpha = &alpha_y;
    else if (component == 2) alpha = &alpha_z;
    else throw std::runtime_error("Invalid component for getCurrent");

    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        std::complex<float> psi[4];
        for (int c = 0; c < 4; ++c) {
            psi[c] = _psi[c][idx];
        }

        // Compute ψ† α^i ψ
        std::complex<float> j_val = {0.0f, 0.0f};
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                j_val += std::conj(psi[i]) * (*alpha)[i*4 + j] * psi[j];
            }
        }

        current[idx] = j_val.real();
    }

    return current;
}

const std::vector<std::complex<float>>& Dirac3D::getComponent(int component) const {
    if (component < 0 || component >= 4) {
        throw std::runtime_error("Invalid spinor component");
    }
    return _psi[component];
}

float Dirac3D::getNorm() const {
    float norm = 0.0f;

    for (uint32_t i = 0; i < _N_total; ++i) {
        for (int c = 0; c < 4; ++c) {
            norm += std::norm(_psi[c][i]);
        }
    }

    return norm * (_dx * _dx * _dx);
}
