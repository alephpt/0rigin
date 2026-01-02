/**
 * Maxwell3D.cpp
 *
 * Implementation of 3D Maxwell electromagnetic field solver
 */

#include "Maxwell3D.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

Maxwell3D::Maxwell3D(uint32_t Nx, uint32_t Ny, uint32_t Nz)
    : _Nx(Nx), _Ny(Ny), _Nz(Nz), _N_total(Nx * Ny * Nz), _dx(1.0f)
{
    // Allocate field storage
    _Ex.resize(_N_total, 0.0f);
    _Ey.resize(_N_total, 0.0f);
    _Ez.resize(_N_total, 0.0f);
    _Bx.resize(_N_total, 0.0f);
    _By.resize(_N_total, 0.0f);
    _Bz.resize(_N_total, 0.0f);
}

void Maxwell3D::initialize(const std::vector<float>& Ex_init,
                          const std::vector<float>& Ey_init,
                          const std::vector<float>& Ez_init,
                          const std::vector<float>& Bx_init,
                          const std::vector<float>& By_init,
                          const std::vector<float>& Bz_init)
{
    if (Ex_init.size() != _N_total || Ey_init.size() != _N_total || Ez_init.size() != _N_total ||
        Bx_init.size() != _N_total || By_init.size() != _N_total || Bz_init.size() != _N_total) {
        throw std::runtime_error("Maxwell3D::initialize: field sizes must match grid size");
    }

    _Ex = Ex_init;
    _Ey = Ey_init;
    _Ez = Ez_init;
    _Bx = Bx_init;
    _By = By_init;
    _Bz = Bz_init;
}

void Maxwell3D::initializeSphericalWave(float wavelength, float amplitude) {
    // Create spherical electromagnetic wave
    // For simplicity: linearly polarized wave with E in x-direction, B in y-direction
    // Propagating radially outward from center

    const float k = 2.0f * M_PI / wavelength;
    const float cx = _Nx / 2.0f;
    const float cy = _Ny / 2.0f;
    const float cz = _Nz / 2.0f;

    for (uint32_t iz = 0; iz < _Nz; ++iz) {
        for (uint32_t iy = 0; iy < _Ny; ++iy) {
            for (uint32_t ix = 0; ix < _Nx; ++ix) {
                const uint32_t idx = index3D(ix, iy, iz);

                // Compute radius from center
                const float dx = ix - cx;
                const float dy = iy - cy;
                const float dz = iz - cz;
                const float r = std::sqrt(dx*dx + dy*dy + dz*dz);

                // Spherical wave profile
                const float phase = k * r;
                const float envelope = (r > 0.1f) ? (1.0f / r) : 10.0f; // Avoid singularity at center

                // E field (radial polarization for simplicity)
                const float E_mag = amplitude * envelope * std::cos(phase);
                if (r > 0.1f) {
                    _Ex[idx] = E_mag * dx / r;
                    _Ey[idx] = E_mag * dy / r;
                    _Ez[idx] = E_mag * dz / r;
                } else {
                    _Ex[idx] = 0.0f;
                    _Ey[idx] = 0.0f;
                    _Ez[idx] = 0.0f;
                }

                // B field (perpendicular to E and propagation direction)
                // For radial E, B is tangential
                const float B_mag = amplitude * envelope * std::cos(phase);
                if (r > 0.1f) {
                    // Tangential component (perpendicular to radial direction)
                    _Bx[idx] = B_mag * (-dy) / r;
                    _By[idx] = B_mag * dx / r;
                    _Bz[idx] = 0.0f;
                } else {
                    _Bx[idx] = 0.0f;
                    _By[idx] = 0.0f;
                    _Bz[idx] = 0.0f;
                }
            }
        }
    }
}

std::vector<float> Maxwell3D::curl(const std::vector<float>& Fx,
                                  const std::vector<float>& Fy,
                                  const std::vector<float>& Fz,
                                  int component) const
{
    std::vector<float> result(_N_total, 0.0f);

    for (uint32_t iz = 0; iz < _Nz; ++iz) {
        for (uint32_t iy = 0; iy < _Ny; ++iy) {
            for (uint32_t ix = 0; ix < _Nx; ++ix) {
                const uint32_t idx = index3D(ix, iy, iz);

                // Get neighbor indices with periodic boundaries
                const uint32_t ix_plus = wrapX(ix + 1);
                const uint32_t ix_minus = wrapX(ix - 1);
                const uint32_t iy_plus = wrapY(iy + 1);
                const uint32_t iy_minus = wrapY(iy - 1);
                const uint32_t iz_plus = wrapZ(iz + 1);
                const uint32_t iz_minus = wrapZ(iz - 1);

                // Central differences for derivatives
                // ∂/∂x: [f(x+dx) - f(x-dx)] / (2*dx)
                // ∂/∂y: [f(y+dy) - f(y-dy)] / (2*dy)
                // ∂/∂z: [f(z+dz) - f(z-dz)] / (2*dz)

                if (component == 0) {
                    // (∇×F)_x = ∂F_z/∂y - ∂F_y/∂z
                    const float dFz_dy = (Fz[index3D(ix, iy_plus, iz)] - Fz[index3D(ix, iy_minus, iz)]) / (2.0f * _dx);
                    const float dFy_dz = (Fy[index3D(ix, iy, iz_plus)] - Fy[index3D(ix, iy, iz_minus)]) / (2.0f * _dx);
                    result[idx] = dFz_dy - dFy_dz;
                }
                else if (component == 1) {
                    // (∇×F)_y = ∂F_x/∂z - ∂F_z/∂x
                    const float dFx_dz = (Fx[index3D(ix, iy, iz_plus)] - Fx[index3D(ix, iy, iz_minus)]) / (2.0f * _dx);
                    const float dFz_dx = (Fz[index3D(ix_plus, iy, iz)] - Fz[index3D(ix_minus, iy, iz)]) / (2.0f * _dx);
                    result[idx] = dFx_dz - dFz_dx;
                }
                else if (component == 2) {
                    // (∇×F)_z = ∂F_y/∂x - ∂F_x/∂y
                    const float dFy_dx = (Fy[index3D(ix_plus, iy, iz)] - Fy[index3D(ix_minus, iy, iz)]) / (2.0f * _dx);
                    const float dFx_dy = (Fx[index3D(ix, iy_plus, iz)] - Fx[index3D(ix, iy_minus, iz)]) / (2.0f * _dx);
                    result[idx] = dFy_dx - dFx_dy;
                }
            }
        }
    }

    return result;
}

void Maxwell3D::evolveElectricField(float dt) {
    // ∂E/∂t = ∇×B
    // E(t + dt) = E(t) + dt * (∇×B)

    const auto curl_Bx = curl(_Bx, _By, _Bz, 0);
    const auto curl_By = curl(_Bx, _By, _Bz, 1);
    const auto curl_Bz = curl(_Bx, _By, _Bz, 2);

    for (uint32_t i = 0; i < _N_total; ++i) {
        _Ex[i] += dt * curl_Bx[i];
        _Ey[i] += dt * curl_By[i];
        _Ez[i] += dt * curl_Bz[i];
    }
}

void Maxwell3D::evolveMagneticField(float dt) {
    // ∂B/∂t = -∇×E
    // B(t + dt) = B(t) - dt * (∇×E)

    const auto curl_Ex = curl(_Ex, _Ey, _Ez, 0);
    const auto curl_Ey = curl(_Ex, _Ey, _Ez, 1);
    const auto curl_Ez = curl(_Ex, _Ey, _Ez, 2);

    for (uint32_t i = 0; i < _N_total; ++i) {
        _Bx[i] -= dt * curl_Ex[i];
        _By[i] -= dt * curl_Ey[i];
        _Bz[i] -= dt * curl_Ez[i];
    }
}

void Maxwell3D::step(float dt) {
    // Second-order Strang splitting:
    // 1. Evolve B by dt/2
    // 2. Evolve E by dt
    // 3. Evolve B by dt/2

    evolveMagneticField(dt / 2.0f);
    evolveElectricField(dt);
    evolveMagneticField(dt / 2.0f);
}

std::vector<float> Maxwell3D::getEnergyDensity() const {
    std::vector<float> energy(_N_total);

    for (uint32_t i = 0; i < _N_total; ++i) {
        const float E2 = _Ex[i]*_Ex[i] + _Ey[i]*_Ey[i] + _Ez[i]*_Ez[i];
        const float B2 = _Bx[i]*_Bx[i] + _By[i]*_By[i] + _Bz[i]*_Bz[i];
        energy[i] = 0.5f * (E2 + B2);
    }

    return energy;
}

float Maxwell3D::getTotalEnergy() const {
    float total_energy = 0.0f;

    for (uint32_t i = 0; i < _N_total; ++i) {
        const float E2 = _Ex[i]*_Ex[i] + _Ey[i]*_Ey[i] + _Ez[i]*_Ez[i];
        const float B2 = _Bx[i]*_Bx[i] + _By[i]*_By[i] + _Bz[i]*_Bz[i];
        total_energy += 0.5f * (E2 + B2);
    }

    return total_energy * (_dx * _dx * _dx); // Multiply by volume element
}
