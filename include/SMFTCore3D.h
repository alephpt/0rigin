// include/SMFTCore3D.h
#pragma once

#include <vector>
#include <cstdint>
#include <vulkan/vulkan.h>

/**
 * SMFTCore3D - 3D Grid Infrastructure for SMFT Simulation
 *
 * This class provides the foundational 3D grid structure for the
 * Schwinger Mean Field Theory simulation migration from 2D to 3D.
 *
 * Sprint 1 Goals (Weeks 1-2):
 *   - Week 1: Core 3D grid allocation, index mapping, neighbor access
 *   - Week 2: GPU integration with Vulkan compute shaders
 *
 * Grid Layout:
 *   - Row-major ordering: idx = k * (Nx * Ny) + j * Nx + i
 *   - Periodic boundary conditions in all dimensions
 *   - 6-neighbor stencil (±x, ±y, ±z)
 *
 * Memory Requirements:
 *   - 32³ grid: 32,768 points, ~128KB per field (development)
 *   - 64³ grid: 262,144 points, ~1MB per field (production)
 */
class SMFTCore3D {
public:
    /**
     * Configuration for 3D grid
     */
    struct Config {
        uint32_t Nx = 32;  // X dimension
        uint32_t Ny = 32;  // Y dimension
        uint32_t Nz = 32;  // Z dimension
        float dx = 1.0f;   // Spatial discretization
        float dt = 0.01f;  // Time step
        float coupling_strength = 1.0f;  // Kuramoto coupling K
    };

    /**
     * 6-neighbor structure for efficient access
     */
    struct Neighbors3D {
        uint32_t x_plus;   // (i+1, j, k)
        uint32_t x_minus;  // (i-1, j, k)
        uint32_t y_plus;   // (i, j+1, k)
        uint32_t y_minus;  // (i, j-1, k)
        uint32_t z_plus;   // (i, j, k+1)
        uint32_t z_minus;  // (i, j, k-1)
    };

    /**
     * Constructor
     * @param device Vulkan logical device (for GPU integration)
     * @param physicalDevice Vulkan physical device (for memory queries)
     */
    SMFTCore3D(VkDevice device, VkPhysicalDevice physicalDevice);

    /**
     * Alternative constructor for CPU-only mode (Week 1 development)
     */
    SMFTCore3D();

    /**
     * Destructor
     */
    ~SMFTCore3D() = default;

    /**
     * Initialize 3D grid with configuration
     * @param config Grid configuration parameters
     */
    void initialize(const Config& config);

    /**
     * Map 3D coordinates to linear index
     * @param i X coordinate [0, Nx)
     * @param j Y coordinate [0, Ny)
     * @param k Z coordinate [0, Nz)
     * @return Linear index in [0, N_total)
     */
    uint32_t index3D(uint32_t i, uint32_t j, uint32_t k) const {
        return k * (_Nx * _Ny) + j * _Nx + i;
    }

    /**
     * Map linear index to 3D coordinates
     * @param idx Linear index
     * @param i Output X coordinate
     * @param j Output Y coordinate
     * @param k Output Z coordinate
     */
    void coords3D(uint32_t idx, uint32_t& i, uint32_t& j, uint32_t& k) const {
        k = idx / (_Nx * _Ny);
        uint32_t rem = idx % (_Nx * _Ny);
        j = rem / _Nx;
        i = rem % _Nx;
    }

    /**
     * Get 6 neighbors with periodic boundary conditions
     * @param i X coordinate
     * @param j Y coordinate
     * @param k Z coordinate
     * @return Structure with 6 neighbor indices
     */
    Neighbors3D getNeighbors(uint32_t i, uint32_t j, uint32_t k) const;

    /**
     * Wrap X coordinate with periodic boundary
     */
    uint32_t wrapX(int32_t x) const {
        return (x + _Nx) % _Nx;
    }

    /**
     * Wrap Y coordinate with periodic boundary
     */
    uint32_t wrapY(int32_t y) const {
        return (y + _Ny) % _Ny;
    }

    /**
     * Wrap Z coordinate with periodic boundary
     */
    uint32_t wrapZ(int32_t z) const {
        return (z + _Nz) % _Nz;
    }

    /**
     * Evolve phase field using Kuramoto dynamics (CPU version for Week 1)
     * @param dt Time step
     */
    void evolveKuramotoCPU(float dt);

    /**
     * Compute synchronization order parameter R
     * R = |<e^(i*theta)>| = |sum(cos(theta) + i*sin(theta))| / N
     */
    void computeRField();

    /**
     * Get average R value over entire grid
     */
    float getAverageR() const;

    /**
     * Initialize uniform phase field (for testing)
     * @param phase Initial phase value for all points
     */
    void initializeUniform(float phase = 0.0f);

    /**
     * Initialize random phase field
     * @param seed Random seed for reproducibility
     */
    void initializeRandom(uint32_t seed = 42);

    /**
     * Get total number of grid points
     */
    uint32_t getTotalPoints() const { return _N_total; }

    /**
     * Get grid dimensions
     */
    uint32_t getNx() const { return _Nx; }
    uint32_t getNy() const { return _Ny; }
    uint32_t getNz() const { return _Nz; }

    /**
     * Direct field access (for GPU upload/download)
     */
    std::vector<float>& getTheta() { return _theta_data; }
    const std::vector<float>& getTheta() const { return _theta_data; }

    std::vector<float>& getOmega() { return _omega_data; }
    const std::vector<float>& getOmega() const { return _omega_data; }

    std::vector<float>& getRField() { return _R_field_data; }
    const std::vector<float>& getRField() const { return _R_field_data; }

    /**
     * Enable GPU compute (Week 2)
     */
    void enableGPU();

    /**
     * Check if GPU is enabled
     */
    bool isGPUEnabled() const { return _gpu_enabled; }

private:
    // Grid dimensions
    uint32_t _Nx;
    uint32_t _Ny;
    uint32_t _Nz;
    uint32_t _N_total;

    // Configuration
    Config _config;

    // 3D field data (stored as linear arrays)
    std::vector<float> _theta_data;    // Phase field θ(x,y,z)
    std::vector<float> _omega_data;    // Intrinsic frequencies ω(x,y,z)
    std::vector<float> _R_field_data;  // Synchronization field R(x,y,z)

    // Vulkan resources (Week 2)
    VkDevice _device;
    VkPhysicalDevice _physical_device;
    bool _gpu_enabled;

    // Helper for Kuramoto coupling computation
    float computeKuramotoCoupling(uint32_t idx) const;
};