# Sprint 3: GPU Acceleration of EM Field Computation - Technical Plan

**Date**: 2025-12-28  
**Phase**: Phase 5 (Scenario 2.6) - Electromagnetic Coupling  
**Goal**: Accelerate EM field extraction from Kuramoto phase using Nova's Vulkan compute pipeline

---

## Executive Summary

**Current State**: CPU-based EM field computation via `EMFieldComputer` (Eigen matrices, O(Nx·Ny) per kernel)  
**Target State**: GPU-accelerated compute shaders for A_μ, E, B extraction (5-10× speedup)  
**Feasibility**: HIGH - Nova has production compute infrastructure (SMFTCompute, SMFTBufferManager, existing shaders)  
**Timeline**: 2 weeks (Steps 1-7), feasible if no major blockers  

**Key Discovery**: Nova already supports full compute pipeline infrastructure:
- ✅ Compute pipeline creation (`vkCreateComputePipelines`)
- ✅ Buffer management (`SMFTBufferManager` with VMA allocator)
- ✅ Compute dispatch (`SMFTCompute` command buffer wrapper)
- ✅ Shader ecosystem (8 production shaders in `shaders/smft/`)
- ✅ Descriptor management (`SMFTDescriptorManager`)

**Strategy**: Reuse existing Nova compute infrastructure, add 3 new compute shaders for EM field extraction.

---

## 1. Nova Compute API Summary

### 1.1 Compute Pipeline Creation

**Location**: `src/SMFTPipelineFactory.cpp:84-95`, `src/smft_pipeline.cpp:344-352`

**Pattern** (from existing code):
```cpp
// 1. Load SPIR-V shader
VkShaderModule shader = loadShader("shaders/smft/em_potentials.comp.spv");

// 2. Define shader stage
VkPipelineShaderStageCreateInfo stageInfo{};
stageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
stageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
stageInfo.module = shader;
stageInfo.pName = "main";

// 3. Create compute pipeline
VkComputePipelineCreateInfo pipelineInfo{};
pipelineInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
pipelineInfo.stage = stageInfo;
pipelineInfo.layout = pipelineLayout;  // From descriptor bindings

VkPipeline pipeline;
vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &pipeline);
```

**Nova Abstraction**: `SMFTPipelineFactory` handles shader loading, SPIR-V compilation, pipeline creation

### 1.2 GPU Buffer Allocation

**Location**: `src/SMFTBufferManager.h:45`, `src/SMFTEngine.cpp:564-632`

**API**:
```cpp
SMFTBufferManager* bufferMgr = new SMFTBufferManager(device, physicalDevice);

// Allocate device buffer (GPU-side, HOST_VISIBLE for CPU access)
auto [buffer, memory] = bufferMgr->createStorageBuffer(size_in_bytes);

// Upload data CPU→GPU
bufferMgr->uploadData(memory, data_ptr, size);

// Download data GPU→CPU
bufferMgr->downloadData(memory, data_ptr, size);

// Direct memory mapping (for frequent access)
void* mapped = bufferMgr->mapMemory(memory, size);
memcpy(mapped, data, size);
bufferMgr->unmapMemory(memory);
```

**Buffer Properties**: All buffers use `VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT`  
**Memory Management**: Tracked in `_managedBuffers` vector, automatic cleanup in destructor

### 1.3 Compute Shader Dispatch

**Location**: `src/SMFTCompute.cpp:76-104`

**API**:
```cpp
SMFTCompute* compute = new SMFTCompute(device, queue, queueFamilyIndex);
compute->initialize();

// Begin command buffer recording
compute->beginBatch();

// Dispatch kernel
struct PushConstants { float dt; uint Nx; uint Ny; } params;
uint32_t workgroupsX, workgroupsY;
SMFTCompute::calculateWorkgroups(Nx, Ny, workgroupsX, workgroupsY, 16, 16);

compute->dispatchCompute(
    pipeline, pipelineLayout, descriptorSet,
    &params, sizeof(params),
    workgroupsX, workgroupsY
);

// Insert memory barrier (GPU synchronization between kernels)
compute->insertMemoryBarrier();

// Submit and wait for completion
compute->submitBatch(true);  // waitForCompletion=true
```

**Key Features**:
- Command buffer reset (`vkResetCommandBuffer`) before each batch
- Pipeline binding, descriptor binding, push constants
- Automatic workgroup calculation (16×16 local size default)
- Memory barriers for read-after-write dependencies
- Fence-based synchronization

### 1.4 Descriptor Set Layout

**Location**: `src/SMFTEngine.cpp:660-663`

**Pattern** (from kuramoto shader):
```cpp
std::vector<VkDescriptorSetLayoutBinding> bindings = {
    {0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},  // theta_in
    {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},  // theta_out
    {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},  // omega
};

VkDescriptorSetLayoutCreateInfo layoutInfo{};
layoutInfo.bindingCount = bindings.size();
layoutInfo.pBindings = bindings.data();

VkDescriptorSetLayout layout;
vkCreateDescriptorSetLayout(device, &layoutInfo, nullptr, &layout);
```

**Descriptor Pool**: Created once, supports up to 32 storage buffers, 10 descriptor sets (src/SMFTEngine.cpp:641-655)

---

## 2. GPU Kernel Design

### 2.1 Kernel 1: `computeEMPotentials.comp`

**Purpose**: Extract gauge potential A_μ = ∂_μ θ from phase field

**Physics**:
```
φ(x,y) = ∂_t θ = (θ_current - θ_previous) / dt
A_x(x,y) = ∂_x θ = (θ[i+1,j] - θ[i-1,j]) / (2dx)
A_y(x,y) = ∂_y θ = (θ[i,j+1] - θ[i,j-1]) / (2dy)
```

**Shader Interface**:
```glsl
#version 450
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;

layout(push_constant) uniform PushConstants {
    float dx;
    float dy;
    float dt;
    uint Nx;
    uint Ny;
} params;

layout(set = 0, binding = 0) readonly buffer ThetaCurrent {
    float theta_current[];
};
layout(set = 0, binding = 1) readonly buffer ThetaPrevious {
    float theta_previous[];
};
layout(set = 0, binding = 2) writeonly buffer Phi {
    float phi[];          // Output: scalar potential
};
layout(set = 0, binding = 3) writeonly buffer A_x {
    float A_x_out[];      // Output: vector potential x
};
layout(set = 0, binding = 4) writeonly buffer A_y {
    float A_y_out[];      // Output: vector potential y
};

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    
    if (i >= params.Nx || j >= params.Ny) return;
    
    uint idx = j * params.Nx + i;
    
    // Temporal derivative: φ = ∂_t θ
    phi[idx] = (theta_current[idx] - theta_previous[idx]) / params.dt;
    
    // Spatial derivatives with periodic BC
    uint i_plus = (i + 1) % params.Nx;
    uint i_minus = (i - 1 + params.Nx) % params.Nx;
    uint j_plus = (j + 1) % params.Ny;
    uint j_minus = (j - 1 + params.Ny) % params.Ny;
    
    uint idx_xp = j * params.Nx + i_plus;
    uint idx_xm = j * params.Nx + i_minus;
    uint idx_yp = j_plus * params.Nx + i;
    uint idx_ym = j_minus * params.Nx + i;
    
    // A_x = ∂_x θ
    A_x_out[idx] = (theta_current[idx_xp] - theta_current[idx_xm]) / (2.0 * params.dx);
    
    // A_y = ∂_y θ
    A_y_out[idx] = (theta_current[idx_yp] - theta_current[idx_ym]) / (2.0 * params.dy);
}
```

**Thread Grid**: `Nx × Ny` invocations (one thread per grid point)  
**Workgroup Size**: 16×16 threads (256 threads per workgroup)  
**Workgroups**: `ceil(Nx/16) × ceil(Ny/16)`  
**Shared Memory**: None (no inter-thread communication needed)  
**Complexity**: O(1) per thread, O(Nx·Ny) total

### 2.2 Kernel 2: `computeFieldStrengths.comp`

**Purpose**: Compute electric and magnetic field strengths from A_μ

**Physics**:
```
E_x = -∂_x φ - ∂_t A_x  (simplified: E_x ≈ -∂_x φ for first implementation)
E_y = -∂_y φ - ∂_t A_y
B_z = ∂_x A_y - ∂_y A_x  (curl of vector potential)
```

**Shader Interface**:
```glsl
layout(set = 0, binding = 0) readonly buffer Phi { float phi[]; };
layout(set = 0, binding = 1) readonly buffer A_x { float A_x_in[]; };
layout(set = 0, binding = 2) readonly buffer A_y { float A_y_in[]; };
layout(set = 0, binding = 3) writeonly buffer E_x { float E_x_out[]; };
layout(set = 0, binding = 4) writeonly buffer E_y { float E_y_out[]; };
layout(set = 0, binding = 5) writeonly buffer B_z { float B_z_out[]; };

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    uint idx = j * params.Nx + i;
    
    // Periodic neighbors
    uint idx_xp = j * params.Nx + ((i+1) % params.Nx);
    uint idx_xm = j * params.Nx + ((i-1+params.Nx) % params.Nx);
    uint idx_yp = ((j+1) % params.Ny) * params.Nx + i;
    uint idx_ym = ((j-1+params.Ny) % params.Ny) * params.Nx + i;
    
    // E = -∇φ (neglecting ∂_t A for first implementation)
    E_x_out[idx] = -(phi[idx_xp] - phi[idx_xm]) / (2.0 * params.dx);
    E_y_out[idx] = -(phi[idx_yp] - phi[idx_ym]) / (2.0 * params.dy);
    
    // B_z = ∂_x A_y - ∂_y A_x
    float dAy_dx = (A_y_in[idx_xp] - A_y_in[idx_xm]) / (2.0 * params.dx);
    float dAx_dy = (A_x_in[idx_yp] - A_x_in[idx_ym]) / (2.0 * params.dy);
    B_z_out[idx] = dAy_dx - dAx_dy;
}
```

**Thread Grid**: `Nx × Ny`  
**Workgroup Size**: 16×16  
**Shared Memory**: None  
**Dependencies**: Reads from A_μ computed by Kernel 1 → requires memory barrier

### 2.3 Kernel 3: `reduceFieldEnergy.comp`

**Purpose**: Compute total field energy U = ∫(E² + B²)/(8π) dV via parallel reduction

**Physics**:
```
u(i,j) = (E_x² + E_y² + B_z²) / (8π)  // Energy density at grid point
U = Σ u(i,j) · dx·dy                   // Total energy (trapezoid integration)
```

**Algorithm**: Tree-based parallel reduction (logarithmic depth)

**Shader Interface** (simplified):
```glsl
layout(local_size_x = 256) in;  // 1D workgroup for reduction

layout(set = 0, binding = 0) readonly buffer E_x { float E_x_in[]; };
layout(set = 0, binding = 1) readonly buffer E_y { float E_y_in[]; };
layout(set = 0, binding = 2) readonly buffer B_z { float B_z_in[]; };
layout(set = 0, binding = 3) buffer EnergySum { float energy_sum[]; };

shared float s_energy[256];

void main() {
    uint tid = gl_LocalInvocationID.x;
    uint gid = gl_GlobalInvocationID.x;
    uint N = params.Nx * params.Ny;
    
    // Step 1: Compute local energy density
    float local_energy = 0.0;
    if (gid < N) {
        float E2 = E_x_in[gid]*E_x_in[gid] + E_y_in[gid]*E_y_in[gid];
        float B2 = B_z_in[gid] * B_z_in[gid];
        local_energy = (E2 + B2) / (8.0 * 3.14159265359);
    }
    s_energy[tid] = local_energy;
    barrier();
    
    // Step 2: Tree reduction in shared memory
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) {
            s_energy[tid] += s_energy[tid + s];
        }
        barrier();
    }
    
    // Step 3: Write workgroup sum to global memory
    if (tid == 0) {
        energy_sum[gl_WorkGroupID.x] = s_energy[0];
    }
}
```

**Thread Grid**: `ceil(Nx·Ny / 256)` threads (1D)  
**Workgroup Size**: 256 threads (single dimension)  
**Shared Memory**: 256 floats (1 KB per workgroup)  
**Multi-Pass**: If N > 256, requires second reduction pass on workgroup sums

**Alternative** (simpler): Use atomic addition (slower but easier):
```glsl
atomicAdd(energy_sum[0], local_energy);
```

---

## 3. Memory Management Plan

### 3.1 Buffer Allocation

**Required Buffers** (per timestep):

| Buffer Name | Size (bytes) | Usage | Notes |
|-------------|--------------|-------|-------|
| `theta_current` | `4·Nx·Ny` | Input (readonly) | Already exists in SMFTEngine |
| `theta_previous` | `4·Nx·Ny` | Input (readonly) | Already exists in SMFTEngine |
| `phi` | `4·Nx·Ny` | Intermediate | New, GPU-only |
| `A_x` | `4·Nx·Ny` | Intermediate | New, GPU-only |
| `A_y` | `4·Nx·Ny` | Intermediate | New, GPU-only |
| `E_x` | `4·Nx·Ny` | Output | New, download if needed |
| `E_y` | `4·Nx·Ny` | Output | New, download if needed |
| `B_z` | `4·Nx·Ny` | Output | New, download if needed |
| `energy_sum` | `4·ceil(N/256)` | Reduction temp | New, GPU-only |
| `energy_total` | `4` | Final output | New, download always (4 bytes!) |

**Total GPU Memory** (for 512×512 grid):
- Existing: 2 × 4MB = 8 MB (`theta_current`, `theta_previous`)
- New persistent: 6 × 4MB = 24 MB (`phi`, `A_x`, `A_y`, `E_x`, `E_y`, `B_z`)
- Reduction temp: ~10 KB (`energy_sum`)
- **Total**: 32 MB (negligible for modern GPUs)

**For 1024×1024 grid**: 128 MB total (still very reasonable)

### 3.2 Buffer Lifecycle

**Initialization** (once at engine start):
```cpp
// Allocate persistent EM field buffers
auto [phi_buf, phi_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [A_x_buf, A_x_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [A_y_buf, A_y_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [E_x_buf, E_x_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [E_y_buf, E_y_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [B_z_buf, B_z_mem] = bufferMgr->createStorageBuffer(Nx * Ny * sizeof(float));
auto [energy_buf, energy_mem] = bufferMgr->createStorageBuffer(sizeof(float));
```

**Per Timestep** (in `SMFTEngine::computeEMFields()`):
```cpp
// 1. theta_current, theta_previous already on GPU (from Kuramoto step)
// 2. Dispatch Kernel 1: theta → A_μ
compute->beginBatch();
compute->dispatchCompute(em_potentials_pipeline, ...);
compute->insertMemoryBarrier();  // Ensure A_μ writes visible

// 3. Dispatch Kernel 2: A_μ → E, B
compute->dispatchCompute(field_strengths_pipeline, ...);
compute->insertMemoryBarrier();  // Ensure E, B writes visible

// 4. Dispatch Kernel 3: E, B → U
compute->dispatchCompute(reduce_energy_pipeline, ...);
compute->submitBatch(true);  // Wait for completion

// 5. Download only energy scalar (4 bytes!)
float field_energy;
bufferMgr->downloadData(energy_mem, &field_energy, sizeof(float));

// 6. (Optional) Download E, B fields for visualization
if (needVisualization) {
    std::vector<float> E_x_cpu(Nx * Ny);
    bufferMgr->downloadData(E_x_mem, E_x_cpu.data(), Nx * Ny * sizeof(float));
    // ... (same for E_y, B_z)
}
```

### 3.3 Transfer Schedule & Bandwidth

**Upload per Timestep**: 0 bytes (theta already on GPU!)  
**Download per Timestep**: 4 bytes (energy scalar only)  
**Download for Visualization** (optional): `3 × Nx × Ny × 4` bytes (E_x, E_y, B_z)

**Bandwidth Analysis** (512×512 grid):
- Energy download: 4 bytes (~0.001 ms over PCIe Gen3)
- Visualization download: 3 MB (~0.3 ms over PCIe Gen3)

**Optimization**: Skip visualization downloads during performance benchmarks → effectively zero CPU↔GPU traffic!

---

## 4. Performance Estimates

### 4.1 CPU Baseline (Current Implementation)

**Source**: `src/physics/EMFieldComputer.cpp`

**Operations per Timestep**:
1. `computeFromPhase()`:
   - Temporal derivative (Nx·Ny): ~1 FLOP/element
   - Spatial derivatives (4 calls, 3 FLOP/element each): ~12 FLOP/element
   - Total: ~13 FLOP/element
2. `computeFieldStrengths()`:
   - 4 spatial derivatives: ~12 FLOP/element
3. `computeFieldEnergy()`:
   - E², B², sum: ~5 FLOP/element
4. `computeDiagnostics()`:
   - Max/avg: ~4 FLOP/element

**Total**: ~34 FLOP/element × Nx·Ny elements

**Measured Performance** (estimated, CPU single-threaded):
- 512×512 grid: 262,144 elements × 34 FLOP = 8.9 MFLOP
- At 3 GHz CPU, ~1 cycle/FLOP → ~3 ms per call

### 4.2 GPU Performance (Projected)

**GPU Specs** (NVIDIA RTX 3060 Ti):
- 4864 CUDA cores @ 1.67 GHz
- 16 TFLOPS FP32 peak
- Memory bandwidth: 448 GB/s

**Kernel 1 (computeEMPotentials)**:
- Compute: 13 FLOP/element × 262,144 = 3.4 MFLOP
- Memory: Read 2 × 1 MB, Write 3 × 1 MB = 5 MB
- Compute time: 3.4 MFLOP / 16 TFLOPS = 0.0002 ms (negligible)
- Memory time: 5 MB / 448 GB/s = 0.011 ms
- **Total**: ~0.02 ms (memory-bound)

**Kernel 2 (computeFieldStrengths)**:
- Compute: 12 FLOP/element × 262,144 = 3.1 MFLOP
- Memory: Read 3 × 1 MB, Write 3 × 1 MB = 6 MB
- Memory time: 6 MB / 448 GB/s = 0.013 ms
- **Total**: ~0.02 ms

**Kernel 3 (reduceFieldEnergy)**:
- Compute: 5 FLOP/element × 262,144 = 1.3 MFLOP
- Reduction passes: log₂(262144/256) ≈ 10 passes
- **Total**: ~0.01 ms

**GPU Total**: 0.02 + 0.02 + 0.01 = **0.05 ms**

**Speedup**: 3 ms (CPU) / 0.05 ms (GPU) = **60×** (theoretical max)

**Realistic Speedup** (accounting for kernel launch overhead, synchronization):
- Overhead: ~0.1 ms per dispatch × 3 kernels = 0.3 ms
- Effective GPU time: 0.05 + 0.3 = 0.35 ms
- **Realistic speedup**: 3 ms / 0.35 ms = **8-10×**

### 4.3 Grid Size Scalability

| Grid Size | Elements | CPU Time (est) | GPU Time (est) | Speedup |
|-----------|----------|----------------|----------------|---------|
| 128×128 | 16,384 | 0.5 ms | 0.3 ms | 1.7× |
| 256×256 | 65,536 | 2 ms | 0.32 ms | 6× |
| 512×512 | 262,144 | 8 ms | 0.35 ms | 23× |
| 1024×1024 | 1,048,576 | 32 ms | 0.50 ms | 64× |
| 2048×2048 | 4,194,304 | 128 ms | 1.5 ms | 85× |

**Key Insight**: GPU advantage grows with grid size (amortizes kernel launch overhead)

### 4.4 Enabled Capabilities

**Current CPU Limits** (single-threaded Eigen):
- Practical max: 512×512 (8 ms per call, 125 Hz update rate)
- Marginal: 1024×1024 (32 ms, 31 Hz)

**GPU Enabled** (with 8-10× speedup):
- Comfortable: 1024×1024 (0.5 ms, 2000 Hz!)
- Feasible: 2048×2048 (1.5 ms, 666 Hz)
- Stretch: 4096×4096 (6 ms, 166 Hz)

**Impact**: Enables high-resolution EM field analysis for vortex defect localization

---

## 5. Implementation Roadmap (Steps 2-7)

### Step 2: Definition (2 days)

**Deliverables**:
- [ ] Compute shader interfaces specification (binding layouts, push constants)
- [ ] Buffer allocation plan (sizes, lifetimes, transfer schedule)
- [ ] Integration points with SMFTEngine (new method `computeEMFieldsGPU()`)
- [ ] Performance validation criteria (correctness thresholds, speedup targets)

**Tasks**:
1. Define descriptor set layouts for 3 kernels
2. Design push constant structures
3. Specify buffer creation in `SMFTEngine::createBuffers()`
4. Define new pipeline variables in `SMFTEngine.h`

### Step 3: Design (2 days)

**Deliverables**:
- [ ] `computeEMPotentials.comp` shader code (GLSL)
- [ ] `computeFieldStrengths.comp` shader code (GLSL)
- [ ] `reduceFieldEnergy.comp` shader code (GLSL)
- [ ] SPIR-V compilation workflow (CMake integration)

**Tasks**:
1. Write 3 compute shaders in `shaders/smft/`
2. Add compilation rules to `CMakeLists.txt`
3. Test SPIR-V generation with `glslangValidator`

### Step 4: Development (4 days)

**Deliverables**:
- [ ] `SMFTEngine::computeEMFieldsGPU()` implementation
- [ ] Pipeline creation in `SMFTEngine::createPipelines()`
- [ ] Descriptor set creation and binding
- [ ] Compute dispatch calls via `SMFTCompute`
- [ ] Buffer allocation for EM fields

**Tasks**:
1. Add EM field buffers to `SMFTEngine::createBuffers()`
2. Create 3 compute pipelines in `createPipelines()`
3. Implement `computeEMFieldsGPU()`:
   - Dispatch Kernel 1 (A_μ)
   - Barrier
   - Dispatch Kernel 2 (E, B)
   - Barrier
   - Dispatch Kernel 3 (energy)
   - Download energy scalar
4. Integrate with `EMFieldComputer` (fallback to CPU if GPU fails)

**Code Structure**:
```cpp
// In SMFTEngine.h
class SMFTEngine {
    // EM field buffers
    VkBuffer _phi_buffer, _A_x_buffer, _A_y_buffer;
    VkBuffer _E_x_buffer, _E_y_buffer, _B_z_buffer;
    VkBuffer _energy_buffer;
    VkDeviceMemory _phi_memory, _A_x_memory, _A_y_memory;
    VkDeviceMemory _E_x_memory, _E_y_memory, _B_z_memory;
    VkDeviceMemory _energy_memory;
    
    // EM field pipelines
    VkPipeline _em_potentials_pipeline;
    VkPipeline _field_strengths_pipeline;
    VkPipeline _reduce_energy_pipeline;
    VkDescriptorSet _em_potentials_desc_set;
    VkDescriptorSet _field_strengths_desc_set;
    VkDescriptorSet _reduce_energy_desc_set;
    
public:
    EMFieldComputer::EMFields computeEMFieldsGPU(double dx, double dy, double dt);
};
```

### Step 5: Testing (3 days)

**Deliverables**:
- [ ] Correctness validation (GPU vs CPU comparison)
- [ ] Performance benchmarks (grid size sweep)
- [ ] Regression tests (existing EM tests pass with GPU backend)

**Test Plan**:
1. **Unit Test**: Single grid point validation
   - Input: Known phase θ(x,y,t)
   - Expected: Hand-calculated A_μ, E, B
   - Tolerance: 10⁻⁵ (FP32 precision)

2. **Integration Test**: Full grid comparison
   - Input: 64×64 traveling wave phase field
   - Compare: GPU vs CPU `EMFieldComputer` results
   - Metrics: L2 norm difference, max pointwise error
   - Tolerance: L2 < 10⁻⁴, max < 10⁻³

3. **Performance Test**: Grid size sweep
   - Grids: 64×64, 128×128, 256×256, 512×512, 1024×1024
   - Measure: CPU time, GPU time, speedup factor
   - Target: 5-10× speedup for 512×512

4. **Regression Test**: Phase 2.6 scenarios
   - Run existing EM coupling tests with GPU backend
   - Verify: Same physical results (energy conservation, field localization)

**Validation Config** (`config/em_gpu_validation.yaml`):
```yaml
validation:
  name: "EM GPU Correctness"
  grid: { size_x: 64, size_y: 64 }
  initial_conditions:
    kuramoto:
      type: "traveling_wave"
      velocity: 0.3
  comparison:
    cpu_backend: true
    gpu_backend: true
    tolerance: 1e-4
  output:
    save_fields: true
    save_comparison: true
```

### Step 6: Launch (2 days)

**Deliverables**:
- [ ] Production deployment of GPU EM computation
- [ ] Documentation update (usage guide, performance notes)
- [ ] Integration with Phase 2.6 scenarios

**Tasks**:
1. Enable GPU backend by default in `SMFTEngine`
2. Add fallback to CPU if GPU unavailable
3. Update `EMObservables` to use `computeEMFieldsGPU()`
4. Update documentation (`docs/GPU_COMPUTE.md`)
5. Run full Phase 2.6 test suite with GPU enabled

**Deployment Checklist**:
- [ ] GPU buffers allocated without leaks
- [ ] Compute shaders compile on CI
- [ ] Performance meets 5× speedup target
- [ ] All Phase 2.6 tests pass
- [ ] No regressions in CPU fallback mode

### Step 7: Growth (1 day)

**Deliverables**:
- [ ] Performance optimization iteration
- [ ] Analysis of bottlenecks
- [ ] Recommendations for future work

**Optimization Ideas**:
1. **Kernel Fusion**: Combine Kernel 1 + Kernel 2 (save one dispatch + barrier)
2. **Shared Memory**: Use shared memory for neighbor access (reduce global reads)
3. **Atomic Energy Reduction**: Replace tree reduction with atomic add (simpler, maybe faster for small grids)
4. **FP64 Support**: Use double precision for energy accumulation (if GPU supports)
5. **Multi-Pass Reduction**: Implement full tree reduction for large grids (>256 workgroups)

**Measurement**:
- Profile with `VK_LAYER_LUNARG_api_dump`
- Identify memory-bound vs compute-bound kernels
- Optimize hot path (likely Kernel 1)

---

## 6. Risk Assessment & Mitigation

### Risk 1: GPU Hardware Availability
**Probability**: Low  
**Impact**: High (blocks entire sprint)  
**Mitigation**: 
- Verify GPU access before Step 4
- Maintain CPU fallback path (already implemented)
- Can develop/test on integrated GPU (slower, but functional)

### Risk 2: Vulkan Compute Pipeline Complexity
**Probability**: Medium  
**Impact**: Medium (delays development)  
**Mitigation**: 
- Nova already has 8 working compute shaders → established patterns
- SMFTCompute abstraction hides most Vulkan boilerplate
- Extensive examples in `src/SMFTEngine.cpp`

### Risk 3: Numerical Accuracy (FP32 vs FP64)
**Probability**: Medium  
**Impact**: Low (affects high-precision scenarios)  
**Mitigation**: 
- Use FP32 initially (matches CPU Eigen default)
- Validate against CPU reference with tolerance 10⁻⁴
- If needed: add FP64 support (check `GL_ARB_gpu_shader_fp64`)

### Risk 4: Memory Bandwidth Bottleneck
**Probability**: High  
**Impact**: Low (still faster than CPU)  
**Mitigation**: 
- Kernels are memory-bound by design (simple stencil operations)
- GPU memory bandwidth (448 GB/s) >> CPU (50 GB/s) → still 8× faster
- Optimization: kernel fusion reduces memory traffic

### Risk 5: Kernel Launch Overhead
**Probability**: High (for small grids)  
**Impact**: Medium (reduces speedup for <256×256)  
**Mitigation**: 
- Target grids: 512×512 and larger (where GPU shines)
- Kernel fusion reduces dispatch count (3→2 or 3→1)
- CPU fallback for small grids (<128×128)

---

## 7. Success Criteria

### Correctness (PASS/FAIL)
- [ ] GPU results match CPU within tolerance: `||E_gpu - E_cpu||₂ / ||E_cpu||₂ < 10⁻⁴`
- [ ] Energy conservation: `|U_final - U_initial| / U_initial < 10⁻³` over 1000 steps
- [ ] Field localization: Vortex defects produce E, B peaks at expected locations
- [ ] Regression: All Phase 2.6 tests pass with GPU backend

### Performance (QUANTITATIVE)
- [ ] 512×512 grid: ≥5× speedup vs CPU
- [ ] 1024×1024 grid: ≥10× speedup vs CPU
- [ ] Zero CPU↔GPU transfers per timestep (except 4-byte energy download)
- [ ] GPU memory footprint <100 MB for 1024×1024 grid

### Integration (QUALITATIVE)
- [ ] Seamless integration with `SMFTEngine`
- [ ] Zero code changes required in Phase 2.6 scenarios
- [ ] CPU fallback works (no GPU dependency for CI)
- [ ] Documentation complete (usage, performance notes)

---

## 8. Timeline Feasibility

**Total Effort**: 14 person-days (2 weeks @ 1 FTE)

| Step | Duration | Dependencies | Risk |
|------|----------|--------------|------|
| 1. Discovery | 1 day | None | ✅ COMPLETE |
| 2. Definition | 2 days | Step 1 | Low |
| 3. Design | 2 days | Step 2 | Low |
| 4. Development | 4 days | Step 3 | Medium |
| 5. Testing | 3 days | Step 4 | Medium |
| 6. Launch | 2 days | Step 5 | Low |
| 7. Growth | 1 day | Step 6 | Low |

**Critical Path**: Discovery → Definition → Design → Development → Testing → Launch

**Parallelization Opportunities**: None (sequential by nature)

**Buffer**: 1 day (included in 2-week estimate)

**Assessment**: **FEASIBLE** if:
- No major blockers in Vulkan pipeline creation
- Compute shaders work first try (have good examples)
- Performance meets expectations (no unexpected GPU issues)

**Fallback Plan**: If Step 4 takes >6 days → reduce scope:
- Skip Kernel 3 (energy reduction) → compute on CPU
- Skip optimization (Step 7) → defer to future sprint

---

## 9. References

### Nova Compute Infrastructure
- `lib/Nova/Core/core.h:51` - `constructComputePipeline()` declaration
- `src/SMFTCompute.h` - Compute dispatch abstraction
- `src/SMFTBufferManager.h` - Buffer allocation API
- `src/SMFTPipelineFactory.cpp:84-95` - Compute pipeline creation example
- `src/SMFTEngine.cpp:564-632` - Buffer creation pattern
- `src/SMFTEngine.cpp:634-663` - Pipeline + descriptor creation

### Existing Compute Shaders (Reference Implementation)
- `shaders/smft/kuramoto_step.comp` - 2D grid stencil pattern
- `shaders/smft/sync_field.comp` - Field computation from phases
- `shaders/smft/accumulate.comp` - Reduction operation (simpler than tree)
- `shaders/include/precision.glsl` - FP32/FP64 handling
- `shaders/include/complex.glsl` - Complex arithmetic (not needed for EM)

### CPU Reference Implementation
- `src/physics/EMFieldComputer.h:245-248` - Main entry point
- `src/physics/EMFieldComputer.cpp:11-40` - Potential computation
- `src/physics/EMFieldComputer.cpp:67-90` - Field strength computation
- `src/physics/EMFieldComputer.cpp:112-138` - Energy integration
- `src/physics/EMFieldComputer.cpp:248-279` - Spatial derivative (reference numerics)

### Validation Framework
- `src/simulations/EMObservables.h` - High-level EM field validation
- `config/scenario_2.6B_electromagnetic_coupling.yaml` - Target scenario

---

## 10. Next Steps (Immediate Actions)

1. **Verify GPU Access** (Step 2 prerequisite):
   ```bash
   vulkaninfo | grep -A5 "VkPhysicalDeviceProperties"
   nvidia-smi  # Check GPU availability
   ```

2. **Create Sprint 3 Branch**:
   ```bash
   git checkout -b sprint3-em-gpu-acceleration
   ```

3. **Define Descriptor Layouts** (Step 2):
   - Write binding specifications for 3 kernels
   - Document push constant structures

4. **Write Compute Shaders** (Step 3):
   - Start with `computeEMPotentials.comp` (simplest)
   - Compile to SPIR-V, verify with `spirv-val`

5. **Implement `createBuffers()` extension** (Step 4):
   - Add EM field buffer creation
   - Test memory allocation before pipeline creation

---

**End of Technical Plan**
