# Multi-Defect Interaction Infrastructure Discovery Report

**Date**: 2025-12-28
**Sprint**: 2.1 - Step 1 (Discovery)
**Objective**: Research existing multi-defect infrastructure and identify implementation gaps

## 1. Existing Infrastructure Status

### ✅ What EXISTS

#### Single Vortex Infrastructure
- **InitialConditions::vortexCore()** - Creates radial R-field profile with core depression
- **SMFTTestRunner** - Single vortex initialization with:
  - Winding number W = ±1 support
  - Core radius specification (physical units)
  - Center position control
  - R-field tanh profile: R(r) = R_min + (R_max - R_min) * tanh(r/r_core)

#### EM Field Computation
- **EMFieldComputer** - Full electromagnetic field extraction from phase gradients:
  - Gauge potential: A_μ = ∂_μ θ
  - Field strengths: E = -∇φ - ∂_t A, B = ∇×A
  - Energy density: U = ∫(E² + B²)/(8π) dV
  - Poynting vector for energy flux
  - Lorentz force computation

#### Validation Infrastructure
- **ScenarioValidator** - Validates defect localization scenarios
- **computeWindingNumber()** - Topological charge W = (1/2π) ∮ ∇θ · dl
- **checkVortexStructure()** - Verifies |W - 1| < tolerance
- **checkRFieldCore()** - Confirms R_min < threshold at core

### ⚠️ What's PARTIALLY Implemented

#### Antiparticle Separation Test (Test 3.4)
- **config/antiparticle_separation_validation.yaml** exists
- Tests particle-antiparticle separation in single vortex field
- Uses "two_particle" initialization mode
- **BUT**: This is particle pairs, NOT vortex pairs

### ❌ What's MISSING

#### Multi-Vortex Creation
- **No vortex pair initialization** (W=+1 and W=-1 at different positions)
- **No separation control** for vortex-antivortex pairs
- **No multi-vortex phase field generation**
- **No addVortex() or createVortexPair() methods**

#### Multi-Defect Validation
- **No total winding number tracking** W_total = ΣW_i
- **No vortex-antivortex annihilation detection**
- **No interaction force measurement** between vortices
- **No multi-defect energy conservation checks**

#### Test Configurations
- **No Test 2.7**: Vortex pair separation scan configs
- **No Test 2.8**: Annihilation dynamics configs
- **No Test 2.9**: EM field topology during annihilation configs

## 2. Physics Requirements Analysis

### Core Multi-Defect Phenomena

#### Vortex-Antivortex Pairs
- **Creation**: Two vortices with opposite winding (W=+1, W=-1)
- **Phase field**: θ(r) = W₁·arg(z - z₁) + W₂·arg(z - z₂)
- **Total charge**: W_total = W₁ + W₂ = 0 for pair
- **R-field**: Overlapping tanh profiles at each core

#### Interaction Dynamics
- **Force law**: F ∝ W₁·W₂/d² (attractive for opposite charges)
- **Critical separation**: d_c ~ 10-20 ℓ_P (below which annihilation occurs)
- **Energy**: E_interaction ~ -W₁·W₂·log(d/a) where a is lattice spacing

#### Annihilation Process
- **Approach**: Vortices spiral toward each other
- **Merger**: When d < 2r_core, cores overlap
- **Outcome**: W_total = 0, phase becomes smooth, R → 1

#### EM Coupling Effects
- **Flux quantization**: Each vortex carries quantized magnetic flux
- **Field interaction**: Overlapping B fields during approach
- **Energy release**: ΔE_EM during annihilation event

## 3. Required Test Configurations

### Test 2.7: Vortex Pair Separation Scan
**Purpose**: Validate interaction force vs separation
```yaml
separations: [5, 10, 20, 40]  # In Planck lengths
winding_numbers: [+1, -1]      # Vortex-antivortex
measure:
  - interaction_force
  - binding_energy
  - field_overlap
```

### Test 2.8: Annihilation Dynamics
**Purpose**: Track vortex merger process
```yaml
initial_separation: 15         # Start close enough to annihilate
evolution_time: 1000          # Long enough for complete merger
track:
  - vortex_positions(t)
  - W_total(t)
  - energy_conservation
  - R_field_evolution
```

### Test 2.9: EM Field Topology
**Purpose**: Measure EM fields during topology change
```yaml
scenarios:
  - pre_annihilation
  - during_merger
  - post_annihilation
measure:
  - magnetic_flux
  - electric_field_peaks
  - gauge_continuity
```

## 4. Implementation Gap Analysis

### Critical Missing Components

| Component | Status | Priority | Effort |
|-----------|--------|----------|--------|
| **createVortexPair()** | ❌ Missing | HIGH | 2-3 hours |
| **Multi-vortex θ field** | ❌ Missing | HIGH | 2-3 hours |
| **Overlapping R-fields** | ❌ Missing | HIGH | 1-2 hours |
| **W_total tracking** | ❌ Missing | MEDIUM | 1 hour |
| **Annihilation detection** | ❌ Missing | MEDIUM | 1-2 hours |
| **Force measurement** | ❌ Missing | LOW | 2-3 hours |
| **Test configs 2.7-2.9** | ❌ Missing | HIGH | 1 hour |
| **Multi-defect validator** | ❌ Missing | MEDIUM | 2-3 hours |

### Existing Components (Reusable)

| Component | Status | Notes |
|-----------|--------|-------|
| **EMFieldComputer** | ✅ Ready | Works for any phase field |
| **ScenarioValidator base** | ✅ Ready | Can extend for multi-vortex |
| **Winding number compute** | ✅ Ready | Already handles topology |
| **ObservableComputer** | ✅ Ready | Energy, momentum work |

## 5. Recommended Implementation Approach

### Phase 1: Core Infrastructure (4-5 hours)
1. **Extend InitialConditions class**:
   ```cpp
   static std::vector<float> vortexPair(
       int Nx, int Ny,
       float x1, float y1, int W1,  // First vortex
       float x2, float y2, int W2,  // Second vortex
       float core_radius);
   ```

2. **Add phase field superposition**:
   ```cpp
   static std::vector<float> multiVortexPhase(
       int Nx, int Ny,
       const std::vector<VortexConfig>& vortices);
   ```

3. **Create test configurations**:
   - config/test_2.7_separation_scan.yaml
   - config/test_2.8_annihilation.yaml
   - config/test_2.9_em_topology.yaml

### Phase 2: Validation & Metrics (3-4 hours)
1. **Extend ScenarioValidator**:
   - validateVortexPairSeparation()
   - validateAnnihilationDynamics()
   - validateEMTopology()

2. **Add observables**:
   - Total winding number W_total(t)
   - Vortex separation distance d(t)
   - Interaction energy E_int(t)

### Phase 3: Testing & Verification (2-3 hours)
1. Run separation scan (d = 5, 10, 20, 40 ℓ_P)
2. Verify force law F ∝ 1/d²
3. Track annihilation dynamics
4. Measure EM field evolution

## Key Findings

### Quick Wins
1. **EMFieldComputer is ready** - Will work for multi-vortex fields
2. **Validation framework extensible** - Can add multi-vortex checks
3. **Single vortex code is solid** - Good foundation to build on

### Major Gaps
1. **No vortex pair creation** - Must implement from scratch
2. **No multi-vortex tests** - Need full test suite
3. **No annihilation tracking** - Critical for physics validation

### Physics Insights
- Vortex-antivortex annihilation is THE key test for EM coupling
- Force law validation requires precise separation control
- EM field topology change during annihilation tests gauge theory

## Recommendation

**Proceed with implementation** - The infrastructure gap is significant but manageable. The existing single-vortex code provides a solid foundation. Implementation should take ~10-12 hours total.

**Priority order**:
1. Create vortex pair initialization (enables all tests)
2. Add test configurations (defines success criteria)
3. Implement validation metrics (measures physics)
4. Run comprehensive test suite (validates theory)

The multi-defect interaction tests are essential for validating SMFT's electromagnetic coupling hypothesis.