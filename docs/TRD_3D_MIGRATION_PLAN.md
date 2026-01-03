# TRD 2D → 3D Migration Plan
**Version**: 1.0
**Date**: January 2026
**Status**: Architecture Planning

## Executive Summary

The current TRD implementation operates on a 2D grid (64×64) with simplified physics. Full 3D+1 spacetime simulation is required for physically accurate modeling of:
- Vortex lines and rings (topological defects)
- Full electromagnetic field dynamics (6 components)
- Proper Dirac spinor evolution in 3D+1 spacetime
- Realistic phase synchronization in 3D lattices

**Recommendation**: Option C - New TRDCore3D class with gradual migration path.

## 1. Current Architecture Analysis

### 1.1 Grid Topology
**Current (2D)**:
- Grid: 64×64 = 4,096 points
- Index mapping: `idx = j * Nx + i`
- Memory per field: 16 KB (float32)
- Total memory (10 fields): ~160 KB

**Required (3D)**:
- Grid: 64×64×64 = 262,144 points
- Index mapping: `idx = k * Nx * Ny + j * Nx + i`
- Memory per field: 1 MB (float32)
- Total memory (10 fields): ~10 MB
- **64× memory increase**

### 1.2 Component Changes Required

| Component | Current (2D) | Required (3D) | Impact |
|-----------|-------------|---------------|--------|
| TRDCore | nx×ny grid | nx×ny×nz grid | Major refactor |
| Kuramoto shader | 4-8 neighbors | 6-26 neighbors | Shader rewrite |
| EM Fields | Ex, Ey, Bz | Ex, Ey, Ez, Bx, By, Bz | Double fields |
| Dirac Evolution | 2D+1 (suppressed z) | 3D+1 full | New gamma matrices |
| Buffers | 4K elements | 262K elements | 64× size |
| Workgroups | 16×16×1 | 8×8×8 or 4×4×4 | GPU optimization |

## 2. Physics Implementation Details

### 2.1 Kuramoto Oscillator (3D)

**Current (2D)**:
```glsl
// 4-neighbor coupling (von Neumann)
coupling = K * (sin(θ_left - θ_i) + sin(θ_right - θ_i) +
                sin(θ_up - θ_i) + sin(θ_down - θ_i))
```

**3D Implementation**:
```glsl
// 6-neighbor coupling (3D von Neumann)
coupling = K * (sin(θ_left - θ_i) + sin(θ_right - θ_i) +
                sin(θ_up - θ_i) + sin(θ_down - θ_i) +
                sin(θ_front - θ_i) + sin(θ_back - θ_i))

// Or 26-neighbor (3D Moore neighborhood)
for (int dz = -1; dz <= 1; dz++)
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++)
            if (!(dx==0 && dy==0 && dz==0))
                coupling += sin(θ_neighbor - θ_i)
```

**Vortex Topology**:
- 2D: Point vortices (winding number ∈ ℤ)
- 3D: Vortex lines, rings, knots (linking number)
- New topological invariants needed

### 2.2 Electromagnetic Fields (3D)

**Maxwell Equations (3D)**:
```
∇·E = ρ/ε₀        (Gauss)
∇·B = 0           (No monopoles)
∇×E = -∂B/∂t      (Faraday)
∇×B = μ₀J + μ₀ε₀∂E/∂t  (Ampère-Maxwell)
```

**Field Components**:
```cpp
// Current (2D+1)
struct EMField2D {
    float Ex, Ey;  // In-plane E
    float Bz;      // Out-of-plane B
};

// Required (3D+1)
struct EMField3D {
    float Ex, Ey, Ez;  // Full E vector
    float Bx, By, Bz;  // Full B vector
};
```

**Stückelberg Mechanism (3D)**:
```cpp
// Gauge-invariant potential
A'_μ = A_μ + ∂_μφ/e

// 3D implementation
Aprime_x = A_x + dx_phi/e;
Aprime_y = A_y + dy_phi/e;
Aprime_z = A_z + dz_phi/e;  // New component
```

### 2.3 Dirac Spinor (3D+1)

**Gamma Matrices**:
```cpp
// Current (2+1D) - 2×2 Pauli matrices
γ⁰ = σ_z, γ¹ = iσ_y, γ² = -iσ_x

// Required (3+1D) - 4×4 Dirac matrices
γ⁰ = [[I,  0],     γ¹ = [[0,  σ_x],
      [0, -I]]           [-σ_x, 0]]

γ² = [[0,  σ_y],    γ³ = [[0,  σ_z],
      [-σ_y, 0]]          [-σ_z, 0]]
```

**Dirac Equation (3D+1)**:
```
(iγ^μ(∂_μ - ieA_μ) - m)ψ = 0

// Expanded:
iγ⁰∂_t ψ + iγ¹(∂_x - ieA_x)ψ + iγ²(∂_y - ieA_y)ψ +
iγ³(∂_z - ieA_z)ψ - mψ = 0
```

**Spinor Components**:
- 2D+1: 2-component spinor (reduced from 4)
- 3D+1: 4-component spinor (full Dirac)

## 3. Implementation Phases

### Phase 1: Core 3D Infrastructure (Week 1-2)
**Deliverables**:
- [ ] TRDCore3D class with nx×ny×nz grid
- [ ] 3D index mapping utilities
- [ ] 3D buffer allocation (1MB per field)
- [ ] Boundary condition handling (periodic/absorbing)
- [ ] 3D visualization slice extraction

**Key Files**:
- `include/TRDCore3D.h`
- `src/TRDCore3D.cpp`
- `shaders/include/index3d.glsl`

### Phase 2: Kuramoto 3D Evolution (Week 3)
**Deliverables**:
- [ ] 3D Kuramoto compute shader
- [ ] 6-neighbor and 26-neighbor coupling
- [ ] 3D synchronization field computation
- [ ] Vortex line detection algorithm

**Key Files**:
- `shaders/trd3d/kuramoto_step_3d.comp`
- `shaders/trd3d/sync_field_3d.comp`
- `src/topology/VortexLineDetector.cpp`

### Phase 3: EM Fields 3D (Week 4-5)
**Deliverables**:
- [ ] 6-component EM field storage
- [ ] 3D Maxwell evolution shaders
- [ ] 3D Stückelberg implementation
- [ ] Field energy computation

**Key Files**:
- `include/physics/StuckelbergEM3D.h`
- `shaders/em3d/maxwell_evolve_3d.comp`
- `shaders/em3d/stuckelberg_3d.comp`

### Phase 4: Dirac 3D+1 (Week 6-7)
**Deliverables**:
- [ ] 4×4 gamma matrix implementation
- [ ] 4-component spinor evolution
- [ ] 3D split-operator method
- [ ] EM coupling in 3D+1

**Key Files**:
- `src/DiracEvolution3D.cpp`
- `shaders/dirac3d/split_operator_3d.comp`

### Phase 5: Testing & Validation (Week 8-10)
**Deliverables**:
- [ ] Port all 2D tests to 3D
- [ ] 3D-specific tests (vortex lines, field divergence)
- [ ] Performance benchmarks
- [ ] Memory usage validation
- [ ] Convergence tests

## 4. Risk Assessment

### 4.1 Memory Constraints

| Grid Size | Points | Memory/Field | Total (10 fields) | GPU VRAM Required |
|-----------|--------|--------------|-------------------|-------------------|
| 32³ | 32K | 128 KB | 1.3 MB | 4 MB |
| 64³ | 262K | 1 MB | 10 MB | 32 MB |
| 128³ | 2M | 8 MB | 80 MB | 256 MB |
| 256³ | 16M | 64 MB | 640 MB | 2 GB |

**Risk**: 64³ requires 32MB VRAM minimum (comfortable on modern GPUs)
**Mitigation**: Start with 32³ for development, scale up for production

### 4.2 Compute Performance

**Operation Scaling**:
- Kuramoto coupling: O(N) → O(N³/²) for 3D neighbors
- FFT: O(N log N) → O(N³ log N) for 3D
- Field evolution: 3× more PDEs to solve

**Risk**: 10-100× slower per timestep
**Mitigation**:
- Optimize workgroup sizes (8×8×8 typical)
- Use shared memory for neighbor access
- Consider multi-GPU for large grids

### 4.3 Breaking Changes

**Incompatibilities**:
- All existing test configs invalid (2D→3D)
- Shader bindings change (more fields)
- Output formats incompatible
- Visualization needs 3D→2D slicing

**Risk**: Cannot run existing experiments
**Mitigation**: Maintain 2D branch until 3D validated

## 5. Migration Strategy

### Option A: In-Place Migration
**Approach**: Replace all 2D code with 3D directly in main branch
- ✅ Clean, no duplication
- ❌ Breaks all existing work
- ❌ No fallback if 3D fails
- **Timeline**: 4 weeks
- **Risk**: HIGH

### Option B: Parallel 3D Branch
**Approach**: Create `feature-3d` branch, maintain both
- ✅ Safe, can switch back
- ✅ A/B testing possible
- ❌ Double maintenance
- ❌ Divergence over time
- **Timeline**: 6 weeks
- **Risk**: MEDIUM

### Option C: TRDCore3D Class (RECOMMENDED)
**Approach**: New `TRDCore3D` alongside `TRDCore`
- ✅ Gradual migration path
- ✅ Can run 2D and 3D tests
- ✅ Shared infrastructure (buffers, engine)
- ✅ Production can choose 2D/3D
- ❌ Some code duplication initially
- **Timeline**: 8 weeks
- **Risk**: LOW

**Implementation Plan (Option C)**:
1. Create `TRDCore3D` inheriting from base class
2. Add `--dimension=3` flag to executable
3. Port shaders incrementally (can test each)
4. Migrate tests one by one
5. Deprecate 2D after 3D validated
6. Eventually merge to single 3D implementation

## 6. Technical Specifications

### 6.1 Data Layout

**Memory Layout (Row-Major)**:
```cpp
// 2D: [j][i] → idx = j * Nx + i
// 3D: [k][j][i] → idx = k * Ny * Nx + j * Nx + i

inline uint32_t index3D(uint32_t i, uint32_t j, uint32_t k,
                       uint32_t Nx, uint32_t Ny) {
    return k * Ny * Nx + j * Nx + i;
}
```

### 6.2 Boundary Conditions

**Periodic (3D)**:
```cpp
i = (i + Nx) % Nx;
j = (j + Ny) % Ny;
k = (k + Nz) % Nz;
```

**Absorbing (3D PML)**:
```cpp
// Perfectly Matched Layer in 3D
float sigma_x = PML_strength * pow(x/PML_width, 2);
float sigma_y = PML_strength * pow(y/PML_width, 2);
float sigma_z = PML_strength * pow(z/PML_width, 2);
```

### 6.3 Workgroup Optimization

**2D Current**:
```glsl
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;
// 256 threads per workgroup
```

**3D Options**:
```glsl
// Option 1: 8×8×8 = 512 threads
layout(local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

// Option 2: 4×4×4 = 64 threads (better for complex kernels)
layout(local_size_x = 4, local_size_y = 4, local_size_z = 4) in;

// Option 3: 16×16×1 slices (for z-direction processing)
layout(local_size_x = 16, local_size_y = 16, local_size_z = 1) in;
```

## 7. Validation Requirements

### 7.1 Correctness Tests
- [ ] Energy conservation in 3D
- [ ] Gauge invariance preserved
- [ ] Norm conservation (unitary evolution)
- [ ] Momentum conservation
- [ ] Charge conservation

### 7.2 Physics Tests
- [ ] 3D vortex line stability
- [ ] EM wave propagation (all 6 components)
- [ ] Spinor dispersion relation
- [ ] Phase synchronization in 3D lattice
- [ ] Topological charge conservation

### 7.3 Performance Benchmarks
- [ ] Time per step vs grid size
- [ ] Memory bandwidth utilization
- [ ] GPU occupancy
- [ ] Weak scaling (fixed work per GPU)
- [ ] Strong scaling (fixed total work)

## 8. Success Criteria

**Minimum Viable 3D (MVP)**:
- ✅ 32³ grid running at 10 fps
- ✅ Kuramoto 6-neighbor coupling working
- ✅ Basic EM field evolution (no Stückelberg)
- ✅ Energy conservation < 1% drift over 1000 steps
- ✅ One working 3D test case

**Production Ready**:
- ✅ 64³ grid at 30+ fps
- ✅ Full Stückelberg EM coupling
- ✅ 4-component Dirac evolution
- ✅ All 2D tests ported and passing
- ✅ Vortex line tracking working
- ✅ Energy conservation < 0.1% over 10,000 steps
- ✅ Documentation complete

## 9. Timeline Summary

| Week | Phase | Deliverables | Risk Points |
|------|-------|--------------|-------------|
| 1-2 | Core 3D | TRDCore3D, buffers | Memory allocation |
| 3 | Kuramoto | 3D coupling, sync | Neighbor access patterns |
| 4-5 | EM Fields | 6-component fields | Shader complexity |
| 6-7 | Dirac | 4×4 matrices, 4-spinor | Numerical stability |
| 8-10 | Testing | Validation, benchmarks | Performance targets |

**Total Duration**: 8-10 weeks for full 3D implementation

## 10. Recommendations

1. **Start with Option C** (TRDCore3D class) for safest migration path
2. **Begin with 32³ grids** for development, scale to 64³ for production
3. **Implement in phases** with validation at each step
4. **Maintain 2D branch** until 3D fully validated
5. **Focus on correctness first**, optimize performance second
6. **Document all breaking changes** for users
7. **Create migration guide** for existing experiments

## Appendix A: Memory Calculations

```python
# 3D Grid Memory Requirements
def calculate_memory(Nx, Ny, Nz, num_fields=10, bytes_per_element=4):
    total_points = Nx * Ny * Nz
    memory_per_field = total_points * bytes_per_element
    total_memory = memory_per_field * num_fields

    return {
        'points': total_points,
        'per_field_MB': memory_per_field / (1024**2),
        'total_MB': total_memory / (1024**2),
        'gpu_vram_min_MB': total_memory * 3 / (1024**2)  # 3x for buffers
    }

# Examples:
32³:  {'points': 32768, 'per_field_MB': 0.125, 'total_MB': 1.25, 'gpu_vram_min_MB': 3.75}
64³:  {'points': 262144, 'per_field_MB': 1.0, 'total_MB': 10.0, 'gpu_vram_min_MB': 30.0}
128³: {'points': 2097152, 'per_field_MB': 8.0, 'total_MB': 80.0, 'gpu_vram_min_MB': 240.0}
```

## Appendix B: Code Migration Examples

### B.1 Index Mapping
```cpp
// 2D → 3D Index Conversion
// OLD (2D):
uint idx = j * Nx + i;

// NEW (3D):
uint idx = k * Ny * Nx + j * Nx + i;

// With bounds checking:
uint index3D_safe(uint i, uint j, uint k) {
    assert(i < Nx && j < Ny && k < Nz);
    return k * Ny * Nx + j * Nx + i;
}
```

### B.2 Neighbor Access
```glsl
// 2D → 3D Neighbor Access
// OLD (2D - 4 neighbors):
float left  = texture(field, ivec2(i-1, j));
float right = texture(field, ivec2(i+1, j));
float up    = texture(field, ivec2(i, j-1));
float down  = texture(field, ivec2(i, j+1));

// NEW (3D - 6 neighbors):
float left  = texture(field, ivec3(i-1, j, k));
float right = texture(field, ivec3(i+1, j, k));
float up    = texture(field, ivec3(i, j-1, k));
float down  = texture(field, ivec3(i, j+1, k));
float front = texture(field, ivec3(i, j, k-1));
float back  = texture(field, ivec3(i, j, k+1));
```

## Document Status
- **Created**: January 2026
- **Status**: APPROVED FOR IMPLEMENTATION
- **Next Steps**: Begin Phase 1 with TRDCore3D class