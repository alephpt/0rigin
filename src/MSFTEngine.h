#pragma once
#include "Nova.h"
#include <vulkan/vulkan.h>
#include <vector>
#include <complex>

/**
 * MSFTEngine - Mass Synchronization Field Theory Physics Compute Engine
 *
 * Implements the core MSFT equation: m(x) = Δ · R(x)
 * where:
 *   - m(x): Effective mass field at position x
 *   - Δ: Mass gap parameter (control parameter)
 *   - R(x): Local synchronization order parameter from Kuramoto dynamics
 *
 * The engine manages GPU computation of:
 *   1. Kuramoto phase evolution: dθ/dt = ω + K·Σsin(θⱼ - θᵢ)
 *   2. Synchronization field: R(x) = |⟨e^(iθ)⟩|
 *   3. Mass field: m(x) = Δ · R(x)
 */
class MSFTEngine {
public:
    /**
     * Constructor - stores Nova instance pointer for Vulkan access
     * @param nova Pointer to initialized Nova graphics engine
     */
    MSFTEngine(Nova* nova);

    /**
     * Destructor - cleans up all Vulkan resources
     */
    ~MSFTEngine();

    /**
     * Initialize the simulation grid and parameters
     * @param Nx Grid width (number of oscillators in x)
     * @param Ny Grid height (number of oscillators in y)
     * @param Delta Mass gap parameter (MSFT control parameter)
     * @param chiral_angle Chiral mass angle for spinor coupling
     */
    void initialize(uint32_t Nx, uint32_t Ny, float Delta, float chiral_angle);

    /**
     * Set initial phase distribution θ(x,y,t=0)
     * @param theta Vector of initial phases (size = Nx * Ny)
     */
    void setInitialPhases(const std::vector<float>& theta);

    /**
     * Set natural frequency distribution ω(x,y)
     * @param omega Vector of natural frequencies (size = Nx * Ny)
     */
    void setNaturalFrequencies(const std::vector<float>& omega);

    /**
     * Execute one time step of the simulation
     * @param dt Time step size
     * @param K Kuramoto coupling strength
     * @param damping Phase damping coefficient
     */
    void step(float dt, float K, float damping);

    /**
     * Get the current synchronization field R(x,y)
     * @return Vector of R values (size = Nx * Ny)
     */
    std::vector<float> getSyncField() const;

    /**
     * Get the current mass field m(x,y) = Δ · R(x,y)
     * @return Vector of mass values (size = Nx * Ny)
     */
    std::vector<float> getMassField() const;

    /**
     * Get the current phase field θ(x,y)
     * @return Vector of phase values (size = Nx * Ny)
     */
    std::vector<float> getPhaseField() const;

    /**
     * Get the gravitational field g(x,y) = -Δ · ∇R(x,y)
     * Returns 2D vector field as interleaved (gx, gy) pairs
     *
     * Physical Interpretation (0.md Step 8 - Bekenstein-Hawking):
     * - Gravity is surface tension of synchronization field R(x)
     * - Spatial gradients ∇R create "pull" toward high-R (massive) regions
     * - This IS gravitational attraction - no separate gravity force
     * - g points from low-R to high-R regions (attracts to mass)
     *
     * Example Usage:
     *   auto g = getGravitationalField();
     *   float gx_at_ij = g[2*(j*Nx + i) + 0];  // x-component at (i,j)
     *   float gy_at_ij = g[2*(j*Nx + i) + 1];  // y-component at (i,j)
     *
     * @return Vector of interleaved (gx, gy) values (size = 2 * Nx * Ny)
     */
    std::vector<float> getGravitationalField() const;

private:
    // Nova graphics engine instance
    Nova* _nova;

    // Simulation parameters
    uint32_t _Nx, _Ny;  // Grid dimensions
    float _Delta;        // Mass gap parameter
    float _chiral_angle; // Chiral mass angle

    // Vulkan GPU resources (to be implemented in Phase 2)
    VkBuffer _theta_buffer;        // Current phase field θ(x,y)
    VkBuffer _theta_out_buffer;    // Updated phase field after step
    VkBuffer _omega_buffer;        // Natural frequencies ω(x,y)
    VkBuffer _R_field_buffer;      // Synchronization field R(x,y)
    VkBuffer _gravity_x_buffer;    // Gravitational field x-component gx(x,y)
    VkBuffer _gravity_y_buffer;    // Gravitational field y-component gy(x,y)

    VkDeviceMemory _theta_memory;       // Memory for theta buffer
    VkDeviceMemory _theta_out_memory;   // Memory for theta_out buffer
    VkDeviceMemory _omega_memory;       // Memory for omega buffer
    VkDeviceMemory _R_field_memory;     // Memory for R_field buffer
    VkDeviceMemory _gravity_x_memory;   // Memory for gravity_x buffer
    VkDeviceMemory _gravity_y_memory;   // Memory for gravity_y buffer

    VkPipeline _kuramoto_pipeline;     // Compute pipeline for phase evolution
    VkPipeline _sync_pipeline;         // Compute pipeline for R field calculation
    VkPipeline _gravity_pipeline;      // Compute pipeline for gravity field ∇R(x,y)
    VkDescriptorSet _descriptor_set;   // Descriptor set for buffer bindings
    VkDescriptorSetLayout _descriptor_layout;  // Layout for descriptors
    VkPipelineLayout _pipeline_layout;  // Pipeline layout

    // CPU-side data mirrors for host access
    std::vector<float> _theta_data;     // Host mirror of phase field
    std::vector<float> _omega_data;     // Host mirror of frequencies
    std::vector<float> _R_field_data;   // Host mirror of sync field
    std::vector<float> _gravity_x_data; // Host mirror of gravity field x-component
    std::vector<float> _gravity_y_data; // Host mirror of gravity field y-component

    // Spinor field data (4-component complex Dirac spinor)
    // Each point (x,y) has 4 complex components for the Dirac spinor
    std::vector<std::complex<float>> _spinor_field;  // 4 * Nx * Ny components

    // Additional buffers for spinor evolution (Phase 2)
    VkBuffer _spinor_buffer;        // Current spinor field Ψ(x,y)
    VkDeviceMemory _spinor_memory;  // Memory for spinor buffer

    // Dirac evolution pipeline (Phase 3)
    VkPipeline _dirac_pipeline;     // Compute pipeline for Dirac evolution

    // Helper to access spinor component at grid point (x,y)
    // component: 0-3 for the 4 Dirac components
    std::complex<float> getSpinorComponent(uint32_t x, uint32_t y, uint32_t component) const;

    // Internal helper methods (stubs for Phase 1)
    void createBuffers();       // Phase 2: Allocate Vulkan buffers
    void createPipelines();     // Phase 3: Load and compile shaders
    void destroyResources();    // Cleanup all Vulkan resources
};