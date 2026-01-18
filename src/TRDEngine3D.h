#pragma once
#include "Nova.h"
#include "TRDPipelineFactory.h"
#include "TRDBufferManager.h"
#include "TRDCompute.h"
#include "TRDDescriptorManager.h"
#include "TRDCore3D.h"
#include "ConservativeSolver.h"
#include <vulkan/vulkan.h>
#include <vector>
#include <complex>
#include <memory>
#include <string>

/**
 * TRDEngine3D - 3D Topological Resonance Dynamics GPU Compute Engine
 *
 * Weeks 3-4: GPU Integration + 3D Kuramoto
 * Weeks 5-6: 3D Electromagnetic Fields
 * Weeks 7-8: Dirac 3D+1 Spinor
 * Weeks 9-10: Testing + Validation
 *
 * GPU Compute Pipeline:
 *   1. kuramoto3d.comp - 3D phase evolution
 *   2. sync_field3d.comp - R-field computation
 *   3. maxwell_evolve_E.comp - Electric field evolution
 *   4. maxwell_evolve_B.comp - Magnetic field evolution
 *   5. dirac3d.comp - 3D+1 Dirac spinor evolution
 */
class TRDEngine3D {
public:
    /**
     * Constructor
     * @param nova Pointer to initialized Nova graphics engine
     */
    TRDEngine3D(Nova* nova);

    /**
     * Destructor - cleans up all Vulkan resources
     */
    ~TRDEngine3D();

    /**
     * Initialize 3D simulation grid
     * @param Nx Grid size in X dimension
     * @param Ny Grid size in Y dimension
     * @param Nz Grid size in Z dimension
     * @param Delta Mass gap parameter (vacuum potential)
     */
    void initialize(uint32_t Nx, uint32_t Ny, uint32_t Nz, float Delta);

    /**
     * Set initial 3D phase distribution θ(x,y,z,t=0)
     * @param theta Vector of initial phases (size = Nx * Ny * Nz)
     */
    void setInitialPhases(const std::vector<float>& theta);

    /**
     * Set 3D natural frequency distribution ω(x,y,z)
     * @param omega Vector of natural frequencies (size = Nx * Ny * Nz)
     */
    void setNaturalFrequencies(const std::vector<float>& omega);

    /**
     * Execute one GPU time step of 3D Kuramoto dynamics
     * @param dt Time step size
     * @param K Kuramoto coupling strength
     * @param damping Phase damping coefficient
     */
    void stepKuramoto3D(float dt, float K, float damping);

    /**
     * Compute 3D synchronization field R(x,y,z)
     * Dispatches sync_field3d.comp shader
     */
    void computeSyncField3D();

    /**
     * Get current 3D phase field θ(x,y,z)
     * @return Vector of phase values (size = Nx * Ny * Nz)
     */
    std::vector<float> getPhaseField3D() const;

    /**
     * Get current 3D synchronization field R(x,y,z)
     * @return Vector of R values (size = Nx * Ny * Nz)
     */
    std::vector<float> getSyncField3D() const;

    /**
     * Get current 3D mass field m(x,y,z) = Δ · R(x,y,z)
     * @return Vector of mass values (size = Nx * Ny * Nz)
     */
    std::vector<float> getMassField3D() const;

    /**
     * Initialize 3D electromagnetic fields
     * Allocates buffers for Ex, Ey, Ez, Bx, By, Bz
     */
    void initializeEM3D();

    /**
     * Step 3D Maxwell equations
     * ∂E/∂t = ∇×B - J
     * ∂B/∂t = -∇×E
     *
     * @param dt Time step size
     */
    void stepMaxwell3D(float dt);

    /**
     * Initialize 3D Stückelberg gauge field
     * Scalar field φ(x,y,z) for gauge restoration
     */
    void initializeStuckelberg3D();

    /**
     * Apply 3D Stückelberg gauge transformation
     * A'μ = Aμ + ∂μφ/e
     */
    void applyGaugeTransform3D();

    /**
     * Initialize 3D vortex line (closed loop)
     * @param center_x X center of vortex loop
     * @param center_y Y center of vortex loop
     * @param center_z Z center of vortex loop
     * @param radius Radius of vortex loop
     * @param axis Axis of rotation (0=x, 1=y, 2=z)
     */
    void initializeVortexLine(float center_x, float center_y, float center_z,
                             float radius, int axis = 2);

    /**
     * Initialize 4-component 3D+1 Dirac spinor field
     * Allocates GPU buffers for spinor components
     */
    void initializeDirac3D();

    /**
     * Step 3D+1 Dirac evolution with EM coupling
     * i∂ψ/∂t = (-iγⁱ∂ᵢ + m)ψ with Dμ = ∂μ - ieA'μ
     *
     * @param dt Time step size
     */
    void stepDirac3D(float dt);

    /**
     * Get grid dimensions
     */
    uint32_t getNx() const { return _Nx; }
    uint32_t getNy() const { return _Ny; }
    uint32_t getNz() const { return _Nz; }
    uint32_t getTotalPoints() const { return _Nx * _Ny * _Nz; }

    /**
     * Check if GPU compute is ready
     */
    bool isGPUReady() const { return _gpu_ready; }

    /**
     * Set physics model (dual-solver routing)
     * @param model "vacuum_kuramoto" | "particle_sine_gordon" | "particle_dirac" | "coupled_vacuum_particle"
     */
    void setPhysicsModel(const std::string& model);

    /**
     * Configure ConservativeSolver integration method
     * @param method "velocity_verlet" | "rk2_symplectic" | "strang_splitting"
     */
    void setIntegrationMethod(const std::string& method);

    /**
     * Get current physics model
     */
    std::string getPhysicsModel() const { return _physics_model; }

    /**
     * Unified simulation step (routes to appropriate solver)
     * @param dt Time step
     */
    void runSimulation(float dt);

private:
    // Nova graphics engine
    Nova* _nova;

    // DUAL SOLVER ARCHITECTURE
    // - TRDCore3D: Vacuum (dissipative/thermodynamic Kuramoto)
    // - ConservativeSolver: Particle (conservative/unitary Sine-Gordon/Dirac)
    std::unique_ptr<TRDCore3D> _core3d;
    std::unique_ptr<ConservativeSolver> _conservative_solver;

    // Physics model selection
    std::string _physics_model;  // "vacuum_kuramoto" | "particle_sine_gordon" | etc.

    // Vulkan resource managers
    std::unique_ptr<TRDPipelineFactory> _pipelineFactory;
    std::unique_ptr<TRDBufferManager> _bufferManager;
    std::unique_ptr<TRDCompute> _compute;
    std::unique_ptr<TRDDescriptorManager> _descriptorManager;

    // Grid dimensions
    uint32_t _Nx, _Ny, _Nz;
    uint32_t _N_total;

    // Physics parameters
    float _Delta;  // Vacuum potential
    uint32_t _time_step;

    // GPU state
    bool _gpu_ready;

    // Vulkan buffers - Kuramoto fields
    VkBuffer _theta_buffer;
    VkBuffer _theta_out_buffer;
    VkBuffer _omega_buffer;
    VkBuffer _R_field_buffer;

    VkDeviceMemory _theta_memory;
    VkDeviceMemory _theta_out_memory;
    VkDeviceMemory _omega_memory;
    VkDeviceMemory _R_field_memory;

    // Vulkan buffers - EM fields (6 components)
    VkBuffer _Ex_buffer, _Ey_buffer, _Ez_buffer;
    VkBuffer _Bx_buffer, _By_buffer, _Bz_buffer;
    VkDeviceMemory _Ex_memory, _Ey_memory, _Ez_memory;
    VkDeviceMemory _Bx_memory, _By_memory, _Bz_memory;

    // Vulkan buffers - Stückelberg gauge
    VkBuffer _phi_buffer;  // Scalar gauge field
    VkBuffer _Ax_buffer, _Ay_buffer, _Az_buffer, _A0_buffer;  // 4-potential
    VkDeviceMemory _phi_memory;
    VkDeviceMemory _Ax_memory, _Ay_memory, _Az_memory, _A0_memory;

    // Vulkan buffers - Dirac spinor (4 complex components = 8 floats)
    VkBuffer _spinor_buffer;
    VkDeviceMemory _spinor_memory;

    // Compute pipelines
    VkPipeline _kuramoto3d_pipeline;
    VkPipeline _sync_field3d_pipeline;
    VkPipeline _maxwell_E_pipeline;
    VkPipeline _maxwell_B_pipeline;
    VkPipeline _dirac3d_pipeline;

    // Pipeline layouts
    VkPipelineLayout _kuramoto3d_layout;
    VkPipelineLayout _sync_field3d_layout;
    VkPipelineLayout _maxwell_E_layout;
    VkPipelineLayout _maxwell_B_layout;
    VkPipelineLayout _dirac3d_layout;

    // Descriptor sets and layouts
    VkDescriptorSet _kuramoto3d_descriptor_set;
    VkDescriptorSet _sync_field3d_descriptor_set;
    VkDescriptorSet _maxwell_descriptor_set;
    VkDescriptorSet _dirac3d_descriptor_set;

    VkDescriptorSetLayout _kuramoto3d_descriptor_layout;
    VkDescriptorSetLayout _sync_field3d_descriptor_layout;
    VkDescriptorSetLayout _maxwell_descriptor_layout;
    VkDescriptorSetLayout _dirac3d_descriptor_layout;

    VkDescriptorPool _descriptor_pool;

    // CPU-side data for initialization
    std::vector<float> _theta_data;
    std::vector<float> _omega_data;
    std::vector<float> _R_field_data;

    // Helper functions
    void createBuffers();
    void createPipelines();
    void createDescriptors();
    void destroyResources();

    void uploadFieldData(VkBuffer buffer, const std::vector<float>& data);
    void downloadFieldData(VkBuffer buffer, std::vector<float>& data) const;
};
