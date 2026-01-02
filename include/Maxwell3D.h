/**
 * Maxwell3D.h
 *
 * 3D Maxwell Electromagnetic Field Solver
 *
 * Implements Maxwell's equations in 3D:
 *   ∂E/∂t = ∇×B   (Faraday's law)
 *   ∂B/∂t = -∇×E  (Ampère's law, no currents/charges)
 *
 * Features:
 *   - 6 field components: Ex, Ey, Ez, Bx, By, Bz
 *   - Periodic boundary conditions
 *   - Central difference curl operator
 *   - CPU implementation for Week 5-6
 */

#pragma once

#include <vector>
#include <cstdint>

class Maxwell3D {
public:
    /**
     * Constructor
     * @param Nx Grid size in X dimension
     * @param Ny Grid size in Y dimension
     * @param Nz Grid size in Z dimension
     */
    Maxwell3D(uint32_t Nx, uint32_t Ny, uint32_t Nz);

    /**
     * Initialize electromagnetic fields
     * @param Ex_init Initial electric field X component
     * @param Ey_init Initial electric field Y component
     * @param Ez_init Initial electric field Z component
     * @param Bx_init Initial magnetic field X component
     * @param By_init Initial magnetic field Y component
     * @param Bz_init Initial magnetic field Z component
     */
    void initialize(const std::vector<float>& Ex_init,
                   const std::vector<float>& Ey_init,
                   const std::vector<float>& Ez_init,
                   const std::vector<float>& Bx_init,
                   const std::vector<float>& By_init,
                   const std::vector<float>& Bz_init);

    /**
     * Initialize with spherical EM wave (for testing)
     * @param wavelength Wavelength λ = 2π/k
     * @param amplitude Field amplitude
     */
    void initializeSphericalWave(float wavelength, float amplitude);

    /**
     * Evolve electric field: ∂E/∂t = ∇×B
     * @param dt Time step
     */
    void evolveElectricField(float dt);

    /**
     * Evolve magnetic field: ∂B/∂t = -∇×E
     * @param dt Time step
     */
    void evolveMagneticField(float dt);

    /**
     * Full Maxwell step (Strang splitting)
     * Evolves fields by dt using second-order accurate scheme
     * @param dt Time step
     */
    void step(float dt);

    /**
     * Compute curl of vector field (∇×F)
     * @param Fx X component of input field
     * @param Fy Y component of input field
     * @param Fz Z component of input field
     * @param component Which component of curl to compute (0=x, 1=y, 2=z)
     * @return Curl component at all grid points
     */
    std::vector<float> curl(const std::vector<float>& Fx,
                           const std::vector<float>& Fy,
                           const std::vector<float>& Fz,
                           int component) const;

    /**
     * Get electric field components (read-only)
     */
    const std::vector<float>& getEx() const { return _Ex; }
    const std::vector<float>& getEy() const { return _Ey; }
    const std::vector<float>& getEz() const { return _Ez; }

    /**
     * Get magnetic field components (read-only)
     */
    const std::vector<float>& getBx() const { return _Bx; }
    const std::vector<float>& getBy() const { return _By; }
    const std::vector<float>& getBz() const { return _Bz; }

    /**
     * Get grid dimensions
     */
    uint32_t getNx() const { return _Nx; }
    uint32_t getNy() const { return _Ny; }
    uint32_t getNz() const { return _Nz; }
    uint32_t getTotalPoints() const { return _N_total; }

    /**
     * Compute electromagnetic energy density
     * u = (E² + B²) / 2
     * @return Energy density at each grid point
     */
    std::vector<float> getEnergyDensity() const;

    /**
     * Compute total electromagnetic energy (integral of energy density)
     * @return Total EM energy
     */
    float getTotalEnergy() const;

    /**
     * Map 3D coordinates to linear index
     */
    uint32_t index3D(uint32_t i, uint32_t j, uint32_t k) const {
        return k * (_Nx * _Ny) + j * _Nx + i;
    }

    /**
     * Periodic boundary wrapping
     */
    uint32_t wrapX(int32_t x) const { return (x + _Nx) % _Nx; }
    uint32_t wrapY(int32_t y) const { return (y + _Ny) % _Ny; }
    uint32_t wrapZ(int32_t z) const { return (z + _Nz) % _Nz; }

private:
    // Grid dimensions
    uint32_t _Nx, _Ny, _Nz;
    uint32_t _N_total;

    // Electric field components
    std::vector<float> _Ex, _Ey, _Ez;

    // Magnetic field components
    std::vector<float> _Bx, _By, _Bz;

    // Grid spacing (unit spacing for now)
    float _dx;
};
