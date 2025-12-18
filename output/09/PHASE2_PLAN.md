# Phase 2: MSFT Coupling - Implementation Plan

## Objective

Couple Dirac evolution to dynamic MSFT synchronization field R(x,y,t) from Kuramoto dynamics.

**Goal**: Demonstrate particle localization in synchronization defects/gradients.

---

## Physics

### Mass Field from Synchronization
```
m(x,y,t) = Δ · R(x,y,t)
```

Where:
- R(x,y,t) = |⟨exp(iθ)⟩| is local Kuramoto synchronization
- Δ = mass gap parameter (MSFT control)
- θ(x,y,t) evolves via: dθ/dt = ω + K·Σsin(θⱼ - θᵢ)

### Force on Dirac Particle
```
F(x,y,t) = -<β> · ∇m(x,y,t) = -<β> · Δ · ∇R(x,y,t)
```

For <β> > 0 (particle state):
- Repelled from HIGH synchronization (high R)
- Attracted to LOW synchronization (low R, defects)

**Prediction**: Particles localize in synchronization defects

---

## Implementation Steps

### 1. GPU Kuramoto Evolution (Already Exists)
Location: `MSFTEngine::step(dt, K, damping)`
- ✓ Computes θ(x,y,t) on GPU
- ✓ Computes R(x,y,t) = synchronization field
- ✓ Outputs R_field via `getSyncField()`

### 2. CPU Dirac Evolution (Phase 1 Complete)
Location: `MSFTEngine::stepWithDirac(dt, lambda_coupling)`
- ✓ Computes Ψ(x,y,t) evolution
- ✓ Takes mass_field = Δ·R as input
- ✓ Tracks center of mass

### 3. Coupling (NEW - Phase 2)
```cpp
// In MSFTEngine::step()
// 1. Evolve Kuramoto (GPU)
step(dt, K, damping);

// 2. Get synchronization field
auto R = getSyncField();

// 3. Compute mass field
std::vector<float> mass_field(Nx * Ny);
for (i = 0; i < Nx*Ny; i++) {
    mass_field[i] = Delta * R[i];
}

// 4. Evolve Dirac (CPU)
stepWithDirac(dt, lambda_coupling, mass_field);

// 5. Optional: Feedback Ψ density -> θ evolution
// (Quantum-classical coupling)
```

---

## Test Scenarios

### Scenario 1: Static Defect
**Setup**:
- Initialize Kuramoto with localized desynchronized region (low R)
- Place Dirac wavepacket nearby
- Evolve coupled system

**Expected**:
- Wavepacket attracted to low-R defect
- Particle density accumulates in defect core

**Measurement**:
- Track CoM trajectory
- Plot ρ(x,y,t) vs R(x,y)
- Verify correlation: high ρ where low R

### Scenario 2: Traveling Wave
**Setup**:
- Initialize Kuramoto with phase gradient → R(x,t) wave
- Place Dirac wavepacket in wave path
- Evolve

**Expected**:
- Wavepacket "surfs" on R gradient
- Moves with synchronization wave

**Measurement**:
- CoM velocity vs wave velocity
- ρ(x,y,t) overlay with R(x,y,t)

### Scenario 3: Defect Collision
**Setup**:
- Two desynchronized regions
- Wavepacket between them
- Allow defects to merge

**Expected**:
- Particle trapped in merged defect
- Energy released during merger

**Measurement**:
- Energy E = ⟨Ψ|H|Ψ⟩ vs time
- ρ distribution during collision

---

## Deliverables

### Code
1. `test/test_msft_coupling.cpp` - Standalone coupling test
2. Updated `MSFTEngine` with coupled evolution loop
3. Analysis scripts for ρ vs R correlation

### Data (output/10/)
1. `defect_localization/` - Scenario 1 data
2. `wave_surfing/` - Scenario 2 data  
3. `defect_collision/` - Scenario 3 data

### Visualizations
1. Particle trajectory overlaid on R(x,y,t) field
2. Density ρ vs synchronization R scatter plot
3. Time evolution: side-by-side ρ and R

### Documentation
1. `PHASE2_RESULTS.md` - Findings summary
2. Physics interpretation of observed localization
3. Comparison to predictions

---

## Success Criteria

✓ Particle CoM moves toward low-R regions  
✓ Density ρ correlates with 1/R (inverse sync)  
✓ Stable evolution over 10k+ coupled steps  
✓ Visualizations clearly show localization  
✓ Results match theoretical prediction F ∝ -∇R

---

## Outstanding from Phase 1

⚠️ Dispersion E(k) full curve - Non-blocking
- Can complete during Phase 2 downtime
- Not required for coupled dynamics
- Add to Phase 1 documentation when done

---

## Timeline Estimate

- Scenario 1 (Static Defect): 1-2 hours
- Scenario 2 (Traveling Wave): 2-3 hours  
- Scenario 3 (Defect Collision): 2-3 hours
- Analysis & Visualization: 2-3 hours
- Documentation: 1-2 hours

**Total**: 1-2 days for complete Phase 2

---

## Ready to Begin

All Phase 1 prerequisites satisfied. Proceeding with Scenario 1.
