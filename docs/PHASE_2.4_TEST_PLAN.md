# Phase 2.4: Relativistic Mass Breakdown Investigation

**Objective**: Understand why Scenario 2.3 showed systematic γ factor underestimation at v>0.5c and identify fundamental SMFT limitations.

**Status**: Phase 2.3 showed **76% pass rate** (38/50 tests), revealing:
- ✓ **Success**: v≤0.3c validated (100% pass, <5% error)
- ✗ **Failure**: v≥0.5c systematic error (8.7-18.7% underestimation)
- ❓ **Mystery**: N=1 shows "mass freeze" (γ≈1 regardless of velocity)

---

## Three-Part Investigation Strategy

### 2.4A: Velocity Threshold Identification (1 week)

**Question**: At what exact velocity does the formula break down?

**Method**: Fine-grained velocity sweep
- **Velocities**: v = [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]c
- **Grid**: 64×64 (validated as "good enough" in 2.3)
- **N**: 10 only (optimal from 2.3)
- **Duration**: 10,000 steps (100 time units)

**Expected Outcomes**:
1. **Interpolate v_critical**: Exact velocity where γ error exceeds 5%
2. **Error scaling**: Does error grow linearly, quadratically, or exponentially with v?
3. **Physics hypothesis**: Test if breakdown correlates with:
   - Lorentz contraction (σ/γ < grid spacing?)
   - Wavepacket dispersion (spreading faster than mass response?)
   - Relativistic corrections needed (Klein-Gordon vs Dirac?)

**Deliverables**:
- γ(v) error curve with interpolated threshold
- Comparison: measured vs theory for 8 velocities
- Recommendation: safe velocity range for SMFT

**Config**: `config/scenario_2.4A_velocity_threshold.yaml`

---

### 2.4B: R-Field Dynamics Study (1 week)

**Question**: Why does N=1 freeze the mass field at γ≈1?

**Hypothesis**: Operator splitting (N>1) is required for R-field to respond to moving spinor. At N=1, Kuramoto evolution happens on too-slow timescale compared to particle motion.

**Method**: Direct N=1 vs N=10 comparison with full spatial R(x,y,t) snapshots
- **Velocities**: v = [0.5, 0.7]c (breakdown onset + clear failure)
- **Grid**: 128×128 (high resolution for spatial analysis)
- **N**: [1, 10] (direct comparison)
- **Duration**: 5,000 steps (shorter for detailed snapshots)
- **Snapshots**: 9 spatial snapshots at t = [0, 10, 25, 50, 100, 200, 300, 400, 500]

**Measurements**:
1. **R(x,y,t) evolution**: Full 2D fields at particle location
2. **∇R along trajectory**: Does R gradient follow moving particle?
3. **R_local vs R_avg**: Compare R at particle vs background R
4. **Response timescale**: τ_R = characteristic time for R to respond
5. **Feedback strength**: β = d(ln R)/d(ln σ) coupling coefficient

**Expected Outcomes**:
- **N=1 hypothesis test**: Does R freeze at background value?
- **N=10 dynamics**: Does R follow particle and modulate mass?
- **Mechanism identified**: Understand Born-Oppenheimer requirement

**Deliverables**:
- R(x,y,t) spatial evolution movies (N=1 vs N=10)
- Trajectory-following R(t) plots
- τ_R measurements for N=1 vs N=10
- Explanation document: "Why N>1 is required for relativistic motion"

**Config**: `config/scenario_2.4B_R_field_dynamics.yaml`

---

### 2.4C: Ultra-Relativistic + Ultra-Fine Grid (2 weeks)

**Question**: Is the v=0.7c failure due to insufficient grid resolution?

**Method**: 512×512 grid at ultra-relativistic velocities
- **Velocities**: v = [0.7, 0.8, 0.9]c
- **Grid**: 512×512 (Δx = 0.195 ℓ_P, 15 points per σ!)
- **N**: [10, 100] (optimal + ultra-high accuracy)
- **Duration**: 10,000 steps (dt reduced to 0.005 for CFL stability)

**Grid Resolution**:
- **Scenario 2.3**: 64×64 → Δx = 1.56 ℓ_P (σ = 3 ℓ_P ≈ 2 grid points)
- **Scenario 2.4C**: 512×512 → Δx = 0.195 ℓ_P (σ = 3 ℓ_P ≈ 15 grid points!)

**Expected Outcomes**:

**If grid resolution was the problem**:
- v=0.7c error should drop from 8.8% (64×64) to <5% (512×512)
- v=0.8c and v=0.9c should show reasonable γ factors

**If grid resolution is NOT the problem**:
- v=0.7c error remains ~9% (fundamental SMFT limitation)
- v=0.8c, v=0.9c show even worse errors
- Conclusion: Need Klein-Gordon corrections or different formalism

**Additional Physics Tests**:
1. **Lorentz contraction**: Does wavepacket contract as σ_∥ = σ/γ?
2. **Time dilation**: Does evolution slow as dt_eff = γ·dt?
3. **E-p relation**: E² = (mc²)² + (pc)² in ultra-relativistic limit?
4. **Dispersion**: Does wavepacket spread excessively at v→c?

**Deliverables**:
- γ(v) validation at v=[0.7, 0.8, 0.9]c on ultra-fine grid
- Grid resolution study: 64×64 vs 256×256 vs 512×512
- Fundamental limits document: "SMFT validity range for relativistic motion"
- Recommendation: Implement Klein-Gordon if needed

**Config**: `config/scenario_2.4C_ultra_relativistic.yaml`

---

## Timeline & Resource Estimate

| Scenario | Duration | Grid | Tests | Compute Time | Analysis Time |
|----------|----------|------|-------|--------------|---------------|
| **2.4A** | 1 week   | 64×64 | 8 velocities | ~2 hours | ~1 day |
| **2.4B** | 1 week   | 128×128 | 2v × 2N × 9 snapshots | ~4 hours | ~2 days |
| **2.4C** | 2 weeks  | 512×512 | 3v × 2N | ~24 hours | ~3 days |
| **Total** | 4 weeks  | - | 23 tests | ~30 hours | ~6 days |

**Critical Path**: Run 2.4A first → Informs 2.4B interpretation → 2.4C confirms/refutes

---

## Success Criteria

### 2.4A Success
✓ Identified v_critical within ±0.05c
✓ γ error curve matches polynomial/exponential model
✓ Safe velocity range documented (e.g., v<0.45c for <5% error)

### 2.4B Success
✓ R-field freeze mechanism understood (N=1 vs N=10)
✓ Quantified τ_R response timescale
✓ Explained why Born-Oppenheimer approximation required

### 2.4C Success (Grid Resolution Hypothesis)
✓ If 512×512 fixes v=0.7c: Grid resolution was the issue
✓ If 512×512 doesn't help: Fundamental SMFT limitation identified

### Overall Phase 2.4 Success
✓ Understand all three failure modes from 2.3
✓ Document SMFT validity range for relativistic physics
✓ Recommend next steps (Klein-Gordon? Different coupling? Acceptable limitations?)

---

## Failure Modes & Contingencies

### If 2.4A shows gradual breakdown (no sharp threshold)
→ Suggests intrinsic SMFT approximation (non-relativistic Dirac operator)
→ Proceed to 2.4C to test grid resolution
→ If grid doesn't help, recommend Klein-Gordon formalism

### If 2.4B shows R-field responds even at N=1
→ Hypothesis wrong - mass freeze is NOT R-field freezing
→ Alternative: Check energy balance (is kinetic energy accounting wrong?)
→ Re-examine m_eff extraction method

### If 2.4C shows no improvement with 512×512
→ Conclude: Not a grid resolution issue
→ Fundamental physics limitation of non-relativistic SMFT
→ Recommend:
  1. Implement relativistic Dirac operator (boost-invariant)
  2. Test Klein-Gordon equation for comparison
  3. Accept SMFT validity range: v<0.4c for quantitative work

---

## Next Steps After Phase 2.4

**If SMFT validated to v≈0.5c**:
- Phase 2.5: Energy-momentum relation E²=(mc²)²+(pc)²
- Phase 2.6: Lorentz transformations (boosted frames)
- Phase 3: Gravitational coupling ∇R → g_field

**If fundamental limitations found**:
- Implement Klein-Gordon corrections
- Compare SMFT vs Klein-Gordon at v>0.5c
- Document SMFT as "semi-classical relativistic theory" (valid v<0.5c)
- Pivot to quantum field theory formulation if needed

---

## Configuration Files Summary

1. **`scenario_2.4A_velocity_threshold.yaml`**
   - 8 velocities: [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]c
   - 64×64 grid, N=10 only
   - Find v_critical for 5% error threshold

2. **`scenario_2.4B_R_field_dynamics.yaml`**
   - 2 velocities: [0.5, 0.7]c
   - 128×128 grid, N=[1, 10]
   - 9 spatial snapshots, full R(x,y,t) evolution

3. **`scenario_2.4C_ultra_relativistic.yaml`**
   - 3 velocities: [0.7, 0.8, 0.9]c
   - 512×512 grid, N=[10, 100]
   - Ultra-fine resolution, v→c limit

---

**Verdict on Scenario 2.3**:
76% pass rate is NOT a "qualified pass" - it's a **discovery of SMFT regime limits**.

**Phase 2.4 Goal**:
Understand the 24% failure rate and map the **validity boundary** of SMFT relativistic physics.

---

**Author**: SMFT Research Team
**Date**: 2025-12-19
**Phase 2.3 Complete**: ✓
**Phase 2.4 Status**: Configurations ready, awaiting execution
