#pragma once
#include "Nova.h"
#include "SMFTPipelineFactory.h"
#include "SMFTBufferManager.h"
#include "SMFTCompute.h"
#include "SMFTDescriptorManager.h"
#include <vulkan/vulkan.h>
#include <vector>
#include <complex>
#include <memory>
#include <random>

// Forward declarations
class DiracEvolution;
class KleinGordonEvolution;

/**
 * SMFTEngine - Synchronization Mass Field Theory Physics Compute Engine
 *
 * Implements the core SMFT equation: m(x) = Δ · R(x)
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
class SMFTEngine {
public:
    /**
     * Constructor - stores Nova instance pointer for Vulkan access
     * @param nova Pointer to initialized Nova graphics engine
     */
    SMFTEngine(Nova* nova);

    /**
     * Destructor - cleans up all Vulkan resources
     */
    ~SMFTEngine();

    /**
     * Initialize the simulation grid and parameters
     * @param Nx Grid width (number of oscillators in x)
     * @param Ny Grid height (number of oscillators in y)
     * @param Delta Mass gap parameter (SMFT control parameter)
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
     * Execute stochastic time step with MSR noise formalism (GPU)
     * @param dt Time step size (0.01 baseline)
     * @param K Kuramoto coupling strength (1.0 baseline)
     * @param damping Phase damping coefficient (0.1 baseline)
     * @param sigma_theta Phase noise amplitude (0.05 baseline)
     * @param sigma_psi Spinor noise amplitude (0.05 baseline)
     */
    void stepStochastic(float dt, float K, float damping,
                       float sigma_theta, float sigma_psi);

    /**
     * Execute CPU-only stochastic Kuramoto step (for Phase Transition tests)
     * Uses Langevin noise: dθ/dt = ω + K·∇²θ + σ·ξ(t) where ξ ~ N(0,1)
     * @param dt Time step size
     * @param K Kuramoto coupling strength
     * @param damping Phase damping coefficient
     * @param sigma Noise amplitude (σ)
     */
    void stepKuramotoCPUStochastic(float dt, float K, float damping, float sigma);

    /**
     * Get the current synchronization field R(x,y)
     * @return Vector of R values (size = Nx * Ny)
     */
    std::vector<float> getSyncField() const;

    /**
     * Get R-field time derivative ∂R/∂t (Phase 4 Test 4.2)
     * Computed via centered finite differences from R-field history
     * @return Vector of ∂R/∂t values (size = Nx * Ny)
     */
    std::vector<float> getRFieldDerivative() const;

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

    /**
     * Initialize Dirac spinor field with Gaussian wavepacket
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * @param x0 Center x-coordinate of Gaussian
     * @param y0 Center y-coordinate of Gaussian
     * @param sigma Gaussian width parameter
     * @param amplitude Initial amplitude (normalized automatically)
     */
    void initializeDiracField(float x0, float y0, float sigma, float amplitude);

    /**
     * Initialize Dirac spinor field with plane wave for dispersion analysis
     * @param kx Momentum x component
     * @param ky Momentum y component
     */
    void initializeDiracPlaneWave(float kx, float ky);

    /**
     * Initialize Klein-Gordon scalar field with Gaussian wavepacket
     * CPU-only implementation for Phase 2.5A comparison
     *
     * @param x0 Center x-coordinate of Gaussian
     * @param y0 Center y-coordinate of Gaussian
     * @param sigma Gaussian width parameter
     * @param amplitude Initial amplitude (normalized automatically)
     */
    void initializeKleinGordonField(float x0, float y0, float sigma, float amplitude);

    /**
     * Initialize Klein-Gordon scalar field with boosted Gaussian wavepacket (Phase 2.5A)
     * CPU-only implementation with relativistic momentum boost
     *
     * @param x0 Center x-coordinate
     * @param y0 Center y-coordinate
     * @param sigma Gaussian width
     * @param vx Boost velocity x (c=1 in Planck units)
     * @param vy Boost velocity y
     * @param R_bg Background synchronization parameter for mass calculation
     */
    void initializeBoostedKleinGordonField(float x0, float y0, float sigma,
                                           float vx, float vy, float R_bg);

    /**
     * Initialize Klein-Gordon scalar field with plane wave for dispersion analysis
     * @param kx Momentum x component
     * @param ky Momentum y component
     */
    void initializeKleinGordonPlaneWave(float kx, float ky);

    /**
     * Initialize Dirac spinor field with boosted Gaussian wavepacket (Scenario 2.3)
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * Creates a Gaussian with relativistic momentum boost p = γ·m·v
     * @param x0 Center x coordinate (grid units)
     * @param y0 Center y coordinate (grid units)
     * @param sigma Gaussian width parameter (grid units)
     * @param vx Boost velocity in x direction (c = 1 in Planck units)
     * @param vy Boost velocity in y direction (c = 1 in Planck units)
     * @param R_bg Background synchronization parameter R for mass calculation
     */
    void initializeBoostedDiracField(float x0, float y0, float sigma,
                                    float vx, float vy, float R_bg);

    /**
     * Step coupled Kuramoto-Dirac evolution with mass coupling
     * CPU-only implementation - GPU shaders exceed timeout budget
     *
     * Physics:
     * - Dirac evolution: i·dΨ/dt = [-iα·∇ + β·m(x)]Ψ
     * - Mass coupling: m(x,y) = Δ·R(x,y) from synchronization field
     * - Feedback: Ψ density influences phase field via coupling λ
     *
     * @param dt Time step size
     * @param lambda_coupling Coupling strength for Ψ→θ feedback
     */
    void stepWithDirac(float dt, float lambda_coupling, int substep_ratio = 1, float K = 1.0f, float damping = 0.1f);

    /**
     * Get the current spinor density |Ψ|²(x,y)
     * @return Vector of density values (size = Nx * Ny)
     */
    std::vector<float> getDiracDensity() const;

    /**
     * Get the full spinor field data
     * @return Vector of complex values (size = 4 * Nx * Ny)
     */
    std::vector<std::complex<double>> getSpinorField() const;

    /**
     * Get internal DiracEvolution object for observable computation
     * @return Pointer to DiracEvolution (nullptr if not initialized)
     */
    const DiracEvolution* getDiracEvolution() const;

    /**
     * Get internal DiracEvolution object for modification (e.g., for specific analysis)
     * @return Non-const pointer to DiracEvolution (nullptr if not initialized)
     */
    DiracEvolution* getDiracEvolutionNonConst();

    /**
     * Initialize two-particle system (particle + antiparticle) for Test 3.4
     * Creates two DiracEvolution instances with opposite beta signs
     *
     * @param x1 Particle center x
     * @param y1 Particle center y
     * @param x2 Antiparticle center x
     * @param y2 Antiparticle center y
     * @param sigma Gaussian width for both
     */
    void initializeTwoParticleSystem(float x1, float y1, float x2, float y2, float sigma);

    /**
     * Get antiparticle DiracEvolution object (Test 3.4)
     * @return Pointer to antiparticle field (nullptr if not initialized)
     */
    const DiracEvolution* getAntiparticleEvolution() const;

    /**
     * Get antiparticle DiracEvolution object for modification (Test 3.4)
     * @return Non-const pointer to antiparticle field (nullptr if not initialized)
     */
    DiracEvolution* getAntiparticleEvolutionNonConst();

    /**
     * Get internal KleinGordonEvolution object for observable computation
     * @return Pointer to KleinGordonEvolution (nullptr if not initialized)
     */
    const KleinGordonEvolution* getKleinGordonEvolution() const;

    /**
     * Get internal KleinGordonEvolution object for modification
     * @return Non-const pointer to KleinGordonEvolution (nullptr if not initialized)
     */
    KleinGordonEvolution* getKleinGordonEvolutionNonConst();

    /**
     * Set the substep ratio N for operator splitting adiabatic approximation
     * N = ratio of fast (Kuramoto) to slow (Dirac) timescales
     * Typical values: 10 (testing), 100 (production)
     *
     * @param N Number of Kuramoto substeps per Dirac step
     */
    void setSubstepRatio(int N);

    // ========================================================================
    // Electromagnetic Field Access (Phase 5 - EM Coupling)
    // ========================================================================

    /**
     * Get vector potential x-component A_x = ∂_x θ
     * @return Vector of A_x values [Nx × Ny]
     */
    const std::vector<float>& getVectorPotentialX() const { return _A_x_data; }

    /**
     * Get vector potential y-component A_y = ∂_y θ
     * @return Vector of A_y values [Nx × Ny]
     */
    const std::vector<float>& getVectorPotentialY() const { return _A_y_data; }

    /**
     * Get scalar potential φ = ∂_t θ
     * @return Vector of φ values [Nx × Ny]
     */
    const std::vector<float>& getScalarPotential() const { return _phi_data; }

    /**
     * Get electric field x-component E_x
     * @return Vector of E_x values [Nx × Ny]
     */
    const std::vector<float>& getElectricFieldX() const { return _E_x_data; }

    /**
     * Get electric field y-component E_y
     * @return Vector of E_y values [Nx × Ny]
     */
    const std::vector<float>& getElectricFieldY() const { return _E_y_data; }

    /**
     * Get magnetic field z-component B_z (2D only)
     * @return Vector of B_z values [Nx × Ny]
     */
    const std::vector<float>& getMagneticFieldZ() const { return _B_z_data; }

    /**
     * Enable/disable electromagnetic coupling
     * @param enable True to enable EM field computation from phase gradients
     * @param coupling_strength EM coupling strength parameter (default: 1.0)
     */
    void setEMCoupling(bool enable, float coupling_strength = 1.0f) {
        _em_coupling_enabled = enable;
        _em_coupling_strength = coupling_strength;
    }

    /**
     * Check if EM coupling is enabled
     * @return True if EM coupling is active
     */
    bool getEMCouplingEnabled() const { return _em_coupling_enabled; }

    /**
     * Initialize hybrid GPU-CPU system with operator splitting
     * Sets up Kuramoto (GPU), Dirac (CPU), and accumulator buffers
     *
     * @param x0 Dirac wavepacket center x
     * @param y0 Dirac wavepacket center y
     * @param sigma Gaussian width
     */
    void initializeHybrid(float x0, float y0, float sigma);

    /**
     * Update time-averaged fields for Dirac evolution
     * Called internally after accumulation completes
     *
     * @param theta_avg Time-averaged phase field
     * @param R_avg Time-averaged sync field
     */
    void updateAveragedFields(const std::vector<float>& theta_avg,
                              const std::vector<float>& R_avg);

private:
    // Nova graphics engine instance
    Nova* _nova;

    // Simulation parameters
    uint32_t _Nx, _Ny;  // Grid dimensions
    float _Delta;        // Mass gap parameter
    float _chiral_angle; // Chiral mass angle
    uint32_t _time_step; // Current timestep (for PRNG seeding)

    // CPU stochastic evolution state
    std::mt19937 _cpu_rng;  // Random number generator for CPU stochastic evolution

    // Vulkan GPU resources (to be implemented in Phase 2)
    VkBuffer _theta_buffer;        // Current phase field θ(x,y)
    VkBuffer _theta_out_buffer;    // Updated phase field after step
    VkBuffer _omega_buffer;        // Natural frequencies ω(x,y)
    VkBuffer _R_field_buffer;      // Synchronization field R(x,y)
    VkBuffer _gravity_x_buffer;    // Gravitational field x-component gx(x,y)
    VkBuffer _gravity_y_buffer;    // Gravitational field y-component gy(x,y)
    VkBuffer _spinor_density_buffer; // Spinor density |Ψ|² for quantum-classical feedback

    VkDeviceMemory _theta_memory;       // Memory for theta buffer
    VkDeviceMemory _theta_out_memory;   // Memory for theta_out buffer
    VkDeviceMemory _omega_memory;       // Memory for omega buffer
    VkDeviceMemory _R_field_memory;     // Memory for R_field buffer
    VkDeviceMemory _gravity_x_memory;   // Memory for gravity_x buffer
    VkDeviceMemory _gravity_y_memory;   // Memory for gravity_y buffer
    VkDeviceMemory _spinor_density_memory; // Memory for spinor density buffer

    VkPipeline _kuramoto_pipeline;     // Compute pipeline for phase evolution
    VkPipeline _sync_pipeline;         // Compute pipeline for R field calculation
    VkPipeline _gravity_pipeline;      // Compute pipeline for gravity field ∇R(x,y)

    // Separate descriptor sets for each pipeline (each shader has different bindings)
    VkDescriptorSet _kuramoto_descriptor_set;   // Descriptor set for kuramoto shader
    VkDescriptorSet _sync_descriptor_set;       // Descriptor set for sync shader
    VkDescriptorSet _gravity_descriptor_set;    // Descriptor set for gravity shader

    // Separate descriptor layouts for each pipeline
    VkDescriptorSetLayout _kuramoto_descriptor_layout;  // Layout for kuramoto descriptors
    VkDescriptorSetLayout _sync_descriptor_layout;      // Layout for sync descriptors
    VkDescriptorSetLayout _gravity_descriptor_layout;   // Layout for gravity descriptors

    // Separate pipeline layouts
    VkPipelineLayout _kuramoto_pipeline_layout;  // Pipeline layout for kuramoto
    VkPipelineLayout _sync_pipeline_layout;      // Pipeline layout for sync
    VkPipelineLayout _gravity_pipeline_layout;   // Pipeline layout for gravity

    VkDescriptorPool _descriptor_pool;  // Descriptor pool for all descriptor sets

    // CPU-side data mirrors for host access
    std::vector<float> _theta_data;     // Host mirror of phase field
    std::vector<float> _omega_data;     // Host mirror of frequencies
    std::vector<float> _R_field_data;   // Host mirror of sync field
    std::vector<float> _gravity_x_data; // Host mirror of gravity field x-component
    std::vector<float> _gravity_y_data; // Host mirror of gravity field y-component

    // R-field history for temporal derivative computation (Phase 4 Test 4.2)
    // Ring buffer: [0] = t-dt, [1] = t, [2] = t+dt
    std::vector<float> _R_history[3];
    int _R_history_index;               // Current position in ring buffer
    float _last_dt;                     // Last timestep size for derivative computation

    // ========================================================================
    // Electromagnetic Field Storage (Phase 5 - EM Coupling)
    // ========================================================================
    // Fields computed from Kuramoto phase gradients: A_μ = ∂_μ θ
    // Used for minimal coupling in Dirac evolution

    std::vector<float> _A_x_data;       // Vector potential A_x = ∂_x θ [Nx × Ny]
    std::vector<float> _A_y_data;       // Vector potential A_y = ∂_y θ [Nx × Ny]
    std::vector<float> _phi_data;       // Scalar potential φ = ∂_t θ [Nx × Ny]
    std::vector<float> _E_x_data;       // Electric field E_x [Nx × Ny]
    std::vector<float> _E_y_data;       // Electric field E_y [Nx × Ny]
    std::vector<float> _B_z_data;       // Magnetic field B_z (2D only) [Nx × Ny]

    // EM field history for temporal derivatives
    std::vector<float> _theta_previous; // Previous timestep phase for ∂_t computation

    // EM coupling configuration
    bool _em_coupling_enabled;          // Flag to enable/disable EM field computation
    float _em_coupling_strength;        // Coupling strength parameter α

    // Spinor field data (4-component complex Dirac spinor)
    // Each point (x,y) has 4 complex components for the Dirac spinor
    std::vector<std::complex<float>> _spinor_field;  // 4 * Nx * Ny components

    // Relativistic field state (CPU-side, avoiding GPU timeouts)
    class DiracEvolution* _dirac_evolution;           // Split-operator Dirac evolution (spinor) - particle
    class DiracEvolution* _dirac_antiparticle;        // Second Dirac field for antiparticle (Test 3.4)
    class KleinGordonEvolution* _kg_evolution;        // Split-operator Klein-Gordon evolution (scalar)
    bool _dirac_initialized;
    bool _kg_initialized;
    bool _two_particle_mode;                          // True if simulating particle + antiparticle

    // Operator splitting state for GPU-CPU hybrid
    int _substep_count;              // Current substep counter
    int _substep_ratio;              // N = Kuramoto steps per Dirac step
    std::vector<float> _theta_avg;   // Time-averaged theta for Dirac
    std::vector<float> _R_avg;       // Time-averaged R for Dirac
    float _lambda_coupling;          // Coupling strength for operator splitting

    // GPU accumulator buffers for time averaging
    VkBuffer _theta_sum_buffer;      // Accumulator for theta averaging
    VkBuffer _R_sum_buffer;          // Accumulator for R averaging
    VkDeviceMemory _theta_sum_memory;
    VkDeviceMemory _R_sum_memory;

    // Accumulation pipeline and descriptors
    VkPipeline _accumulation_pipeline;
    VkDescriptorSet _accumulation_descriptor_set;
    VkDescriptorSetLayout _accumulation_descriptor_layout;
    VkPipelineLayout _accumulation_pipeline_layout;

    // Additional buffers for spinor evolution (Phase 2)
    VkBuffer _spinor_buffer;        // Current spinor field Ψ(x,y)
    VkDeviceMemory _spinor_memory;  // Memory for spinor buffer

    // Dirac evolution pipeline (Phase 3)
    VkPipeline _dirac_pipeline;     // Compute pipeline for Dirac evolution

    // Stochastic pipelines for MSR formalism
    VkPipeline _kuramoto_stochastic_pipeline;  // Stochastic phase evolution
    VkPipeline _dirac_stochastic_pipeline;     // Stochastic Dirac evolution
    VkDescriptorSet _dirac_descriptor_set;     // Descriptor set for Dirac shader
    VkDescriptorSetLayout _dirac_descriptor_layout; // Layout for Dirac descriptors
    VkPipelineLayout _dirac_pipeline_layout;   // Pipeline layout for Dirac

    // Helper to access spinor component at grid point (x,y)
    // component: 0-3 for the 4 Dirac components
    std::complex<float> getSpinorComponent(uint32_t x, uint32_t y, uint32_t component) const;

    // Internal helper methods
    void createBuffers();       // Phase 2: Allocate Vulkan buffers
    void createPipelines();     // Phase 3: Load and compile shaders
    void uploadToGPU();         // Upload CPU data to GPU buffers
    void downloadFromGPU();     // Download GPU results to CPU arrays
    void destroyResources();    // Cleanup all Vulkan resources

    // Pipeline factory for managing shader compilation
    std::unique_ptr<SMFTPipelineFactory> _pipelineFactory;

    // Compute dispatcher for GPU operations
    std::unique_ptr<SMFTCompute> _compute;

    // Buffer manager for handling Vulkan buffer operations
    std::unique_ptr<SMFTBufferManager> _bufferManager;

    // Descriptor manager for handling Vulkan descriptor operations
    std::unique_ptr<SMFTDescriptorManager> _descriptorManager;
};