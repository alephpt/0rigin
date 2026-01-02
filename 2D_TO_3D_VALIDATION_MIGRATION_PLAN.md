# 2D to 3D Validation Migration Plan

**Date**: 2026-01-02
**Status**: Ready to Execute
**Mission**: Port ALL existing 2D physics validation tests to 3D framework

---

## Executive Summary

**Context**:
- 3D infrastructure COMPLETE (Maxwell3D, Dirac3D, SMFTCore3D) ✅
- Multiple 2D validation tests exist but NOT verified in 3D
- User requires Option A: Full 2D→3D porting of ALL validation tests

**Current 2D Validation Tests Identified**:
1. ✅ `test_stuckelberg_vortex_bfield.cpp` - Stückelberg vortex → B-field (2D: x-y plane)
2. ✅ `test_geodesic_verification.cpp` - Geodesic equation (2D: x-y plane)
3. ✅ `test_weak_field_limit.cpp` - Weak field → Newtonian limit (2D: x-y plane)
4. ✅ `test_maxwell3d_wave.cpp` - EM wave propagation (ALREADY 3D) ✅
5. ⚠️ `test_em_gravity_coupling.cpp` - Incomplete (placeholder)

**Missing 2D Tests** (from TODO.md requirements):
- Lorentz force / cyclotron motion (not found in codebase)
- Three-body EM dynamics (G2 - not implemented)

**3D Infrastructure Available**:
- `Maxwell3D` (6-component EM: Ex, Ey, Ez, Bx, By, Bz)
- `Dirac3D` (4-spinor evolution with FFTW)
- `SMFTCore3D` (Kuramoto synchronization in 3D)

---

## Phase 1: Electromagnetic Validation (Priority 1)

### 1.1 Lorentz Force in 3D ⭐ **HIGH PRIORITY**

**Current Status**: No 2D test exists - CREATE FROM SCRATCH

**Physical Requirements**:
- Charged particle in 3D electromagnetic field
- Boris integrator (symplectic, energy-conserving)
- Test scenarios:
  1. Pure B_z field → cyclotron in x-y plane (validates 2D equivalent)
  2. Pure B_x field → cyclotron in y-z plane (NEW 3D)
  3. Pure B_y field → cyclotron in x-z plane (NEW 3D)
  4. General B=(Bx,By,Bz) → helical 3D motion
  5. E×B drift in 3D

**Implementation**:
```cpp
// File: test/test_lorentz_force_3d.cpp
// - Initialize particle (x, y, z, vx, vy, vz)
// - Boris integrator (velocity Verlet with magnetic field)
// - Evolve in 3D Maxwell3D fields
// - Measure: cyclotron frequency ω = qB/m
// - Verify: energy conservation < 0.01%
```

**Quality Gates**:
- ✅ Cyclotron frequency: ω within 3% of qB/m (all field orientations)
- ✅ Energy conservation: ΔE/E < 0.01%
- ✅ Helical pitch: Consistent with parallel velocity component

**Deliverables**:
- `test/test_lorentz_force_3d.cpp`
- `config/lorentz_force_3d.yaml`

---

### 1.2 Stückelberg Vortex in 3D ⭐ **HIGH PRIORITY**

**Current Status**: ✅ 2D test exists (`test_stuckelberg_vortex_bfield.cpp`)

**3D Extension Requirements**:
- **Vortex lines** (not points) - topological structures in 3D
- Phase field: θ(x,y,z) for vortex line along axis
- B-field from 3D curl: B = ∇×A where A = ∇θ
- Vortex ring stability (toroidal geometry)

**Test Scenarios**:
1. Vortex line along z-axis → B_z (reproduces 2D result)
2. Vortex line along x-axis → B_x (NEW 3D)
3. Vortex line along y-axis → B_y (NEW 3D)
4. Vortex ring (toroidal) → Full 3D B-field
5. Flux quantization: Φ = ∮ A·dl = 2πn (topological invariant)

**Implementation**:
```cpp
// File: test/test_stuckelberg_vortex_3d.cpp
// - Initialize vortex line: θ = atan2(y-y0, x-x0) for z-axis
// - Compute A = ∇θ via StuckelbergEM extended to 3D
// - Compute B = ∇×A using Maxwell3D curl operator
// - Measure B_max, flux quantization
// - Verify topological stability over 1000 steps
```

**Quality Gates**:
- ✅ B_max ≈ 1.567 (same as 2D for equivalent vortex)
- ✅ Flux quantization: Φ/2π = integer ± 0.01
- ✅ Vortex line stability: B_max drift < 1% over 1000 steps

**Deliverables**:
- `test/test_stuckelberg_vortex_3d.cpp`
- `config/stuckelberg_vortex_3d.yaml`
- Extension to `physics/StuckelbergEM.h` for 3D (if needed)

---

### 1.3 EM Wave Propagation (G1)

**Current Status**: ✅ ALREADY 3D (`test_maxwell3d_wave.cpp`)

**Enhancement Requirements**:
- Current test: Spherical wave initialization ✅
- Add: Poynting vector verification S = E × B
- Add: Wave packet dispersion in 3D
- Add: Phase velocity verification v_phase = ω/k

**Implementation**:
```cpp
// File: test/test_maxwell3d_wave_enhanced.cpp (or update existing)
// - Spherical wave packet (Gaussian envelope)
// - Track wave front position vs time
// - Measure: phase velocity = c (speed of light)
// - Compute: Poynting vector S = E × B
// - Verify: Energy flux direction = wave propagation direction
```

**Quality Gates**:
- ✅ Wave speed: v = 1.0 ± 1% (natural units, c=1)
- ✅ Energy conservation: ΔE/E < 5% (already passing)
- ✅ Poynting vector: S ⊥ E, S ⊥ B, |S| = |E||B|

**Deliverables**:
- Update `test/test_maxwell3d_wave.cpp` OR create enhanced version
- `config/maxwell3d_wave_enhanced.yaml`

---

## Phase 2: General Relativity (Category A) in 3D

### 2.1 Geodesic Deviation (A3)

**Current Status**: ✅ 2D test exists (`test_geodesic_verification.cpp`)

**3D Extension Requirements**:
- Full 3D metric: g_μν(x,y,z)
- Christoffel symbols: Γ^μ_αβ in 3D
- Geodesic equation: d²x^μ/dτ² = -Γ^μ_αβ (dx^α/dτ)(dx^β/dτ)
- Test particle evolution in 3D curved spacetime

**Implementation**:
```cpp
// File: test/test_geodesic_3d.cpp
// - Generate 3D R-field: R(x,y,z) = 1 + A·exp(-r²/σ²)
// - Compute Christoffel symbols in 3D
// - Evolve test particle (Dirac wavepacket center of mass)
// - Compare to analytical geodesic solution
// - Measure trajectory deviation
```

**Quality Gates**:
- ✅ Trajectory deviation < 1% from analytical (weak field)
- ✅ Energy conservation along geodesic < 0.1%
- ✅ Norm conservation < 0.01%

**Deliverables**:
- `test/test_geodesic_3d.cpp`
- Extension to `GeodesicIntegrator.cpp` for 3D
- `config/geodesic_3d.yaml`

---

### 2.2 Weak Field Limit (A2)

**Current Status**: ✅ 2D test exists (`test_weak_field_limit.cpp`)

**3D Extension Requirements**:
- 3D Newtonian potential: ∇²φ = 4πG·ρ
- R-field → metric perturbation h_μν in 3D
- Verify: g_μν ≈ η_μν + h_μν with |h| ≪ 1
- Test particle acceleration: a = -∇φ = -∇R

**Implementation**:
```cpp
// File: test/test_weak_field_3d.cpp
// - Initialize Gaussian R-field perturbation (3D)
// - Test particle at various radii
// - Measure acceleration a_SMFT = -∇R
// - Compare to Newtonian: a_Newton = GM/r²
// - Verify relative error < 0.1%
```

**Quality Gates**:
- ✅ |a_SMFT - a_Newton| / |a_Newton| < 0.1% (weak field)
- ✅ Energy conservation < 0.1%
- ✅ 1/r² falloff verified at multiple radii

**Deliverables**:
- `test/test_weak_field_3d.cpp`
- `config/weak_field_3d.yaml`

---

## Phase 3: Advanced EM-Gravity Coupling (G2, G3)

### 3.1 Three-Body EM Dynamics (G2)

**Current Status**: ⚠️ Not implemented in 2D - CREATE FROM SCRATCH

**Implementation Requirements**:
- 3 charged particles in 3D space
- Coulomb forces: F_ij = k·q_i·q_j·r_ij/|r_ij|³
- Verify superposition principle
- Energy and momentum conservation

**Test Configuration**:
```cpp
// File: test/test_three_body_em_3d.cpp
// Particles:
// - P1: (+q, position (x1, y1, z1))
// - P2: (+q, position (x2, y2, z2))
// - P3: (-2q, position (x3, y3, z3)) → net charge = 0
//
// Initial: Equilateral triangle in 3D space
// Evolve: Boris integrator with E-field coupling
// Measure: Total energy, momentum, trajectories
```

**Quality Gates**:
- ✅ Energy conservation: ΔE/E < 0.1%
- ✅ Momentum conservation: Δp/p < 0.1%
- ✅ Forces match analytical 3-body Coulomb within 5%

**Deliverables**:
- `test/test_three_body_em_3d.cpp`
- `config/three_body_em_3d.yaml`

---

### 3.2 EM-Gravity Coupling (G3)

**Current Status**: ⚠️ Placeholder only (`test_em_gravity_coupling.cpp`)

**Implementation Requirements**:
- Stress-energy tensor: T_μν^(EM) from Maxwell3D
- R-field evolution: ∂R/∂t ~ -∇·T^(EM) (energy coupling)
- High-energy EM pulse → R-field perturbation

**Implementation**:
```cpp
// File: test/test_em_gravity_coupling_3d.cpp
// - Initialize high-energy EM Gaussian pulse (3D)
// - Disable Kuramoto (K=0) to isolate EM→R effect
// - Evolve Maxwell3D + R-field with coupling
// - Measure ΔR/ρ_EM (coupling constant ε)
// - Verify energy-momentum conservation
```

**Quality Gates**:
- ✅ ΔR/ρ_EM ≈ ε (coupling constant) within 10%
- ✅ Energy conservation (EM + gravity) < 0.1%
- ✅ R-field perturbation follows EM energy density

**Deliverables**:
- `test/test_em_gravity_coupling_3d.cpp`
- GPU shader: `shaders/em_stress_energy.comp`
- Extension to `SMFTCore3D` for EM coupling
- `config/em_gravity_coupling_3d.yaml`

---

## Implementation Strategy

### Week 1: EM Validation (Lorentz + Stückelberg)
**Days 1-2**: Lorentz force 3D
- Implement Boris integrator for 3D
- Test cyclotron motion in Bx, By, Bz fields
- Verify helical trajectories

**Days 3-4**: Stückelberg vortex 3D
- Extend vortex to lines (not points)
- Test vortex along x, y, z axes
- Implement vortex ring (toroidal)

**Day 5**: Testing and validation
- Run all tests, compare to 2D results
- Document deviations (if any)

---

### Week 2: GR Validation (Category A)
**Days 1-2**: Geodesic 3D
- Extend `GeodesicIntegrator` to 3D
- Test curved spacetime geodesics
- Verify tidal forces

**Days 3-4**: Weak field 3D
- 3D Newtonian potential
- Test particle trajectories
- 1/r² falloff verification

**Day 5**: Integration testing
- Cross-check with 2D results
- Performance benchmarking (32³ vs 64³)

---

### Week 3: Advanced Tests (G2, G3)
**Days 1-2**: Three-body EM
- Implement 3-charge system
- Coulomb forces in 3D
- Energy/momentum conservation

**Days 3-4**: EM-gravity coupling
- Implement stress-energy tensor coupling
- GPU shader integration
- R-field back-reaction

**Day 5**: Full validation report
- All tests passing
- Performance metrics
- Documentation

---

### Week 4: Buffer & Documentation
**Days 1-2**: QA review
- Re-run all tests
- Regression testing (2D tests still pass)
- Edge case validation

**Days 3-4**: Performance optimization
- Profile bottlenecks
- GPU optimization (if needed)
- Benchmark 64³ grid performance

**Day 5**: Final documentation
- Complete migration report
- Update TODO.md with completion status
- Archive 2D tests (mark as legacy)

---

## Quality Gates (All Tests)

**Required for Test Completion**:
- ✅ 3D tests pass with same accuracy as 2D equivalents
- ✅ Conservation laws verified (energy, momentum, charge, norm)
- ✅ No regressions in existing 2D tests
- ✅ Performance acceptable (64³ grid < 60s per test)
- ✅ Documentation complete (inline + markdown)

**Required for Phase Completion**:
- ✅ All Phase tests passing
- ✅ Cross-validation with 2D results
- ✅ Git commits with descriptive messages
- ✅ Updated CMakeLists.txt

---

## Deliverables Summary

**New Test Files** (~8 files):
1. `test/test_lorentz_force_3d.cpp` ⭐ **CREATE FROM SCRATCH**
2. `test/test_stuckelberg_vortex_3d.cpp` (extend 2D)
3. `test/test_geodesic_3d.cpp` (extend 2D)
4. `test/test_weak_field_3d.cpp` (extend 2D)
5. `test/test_three_body_em_3d.cpp` ⭐ **CREATE FROM SCRATCH**
6. `test/test_em_gravity_coupling_3d.cpp` (complete placeholder)
7. `test/test_maxwell3d_wave_enhanced.cpp` (optional enhancement)
8. `test/test_vortex_ring_stability_3d.cpp` (bonus - vortex rings)

**Config Files** (~8 YAML):
- Corresponding `config/*.yaml` for each test

**Source Extensions** (if needed):
- `src/physics/StuckelbergEM3D.cpp` (if 3D extension needed)
- `src/GeodesicIntegrator3D.cpp` (if separate from 2D)
- `shaders/em_stress_energy.comp` (GPU coupling)

**Documentation**:
- This file (`2D_TO_3D_VALIDATION_MIGRATION_PLAN.md`)
- `2D_TO_3D_VALIDATION_REPORT.md` (upon completion)
- Updated `TODO.md` with completion status

---

## Success Criteria

**All 2D validation results REPRODUCED in 3D**:
- Lorentz force: ω_cyclotron within 3% ✓
- Stückelberg: B_max = 1.567 ± 5% ✓
- Geodesics: Trajectory deviation < 1% ✓
- Weak field: Newtonian limit verified ✓
- EM coupling: Energy conservation < 0.1% ✓

**Timeline**: 3-4 weeks (comprehensive porting)

**Start Date**: 2026-01-02
**Target Completion**: 2026-01-30

---

## Immediate Next Actions

**Action 1**: Create `test/test_lorentz_force_3d.cpp`
- Boris integrator implementation
- Test all 3 field orientations (Bx, By, Bz)
- Helical motion validation

**Action 2**: Extend `test_stuckelberg_vortex_bfield.cpp` to 3D
- Vortex line geometry
- Flux quantization
- Topological stability

**Action 3**: Update CMakeLists.txt
- Add new test executables
- Link Maxwell3D, Dirac3D, SMFTCore3D

**Action 4**: Create config files
- YAML configurations for all new tests
- Document test parameters

---

**READY TO EXECUTE - AWAITING USER CONFIRMATION TO BEGIN PHASE 1**
