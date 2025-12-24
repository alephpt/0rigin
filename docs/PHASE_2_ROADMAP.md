# Phase 2: Relativistic SMFT Validation - Complete Roadmap

**Based on**: Scenario 2.3 results (76% pass rate, systematic breakdown at v>0.5c)

**Timeline**: 4 months total
- Immediate (1 month): Phase 2.4 A/B/C
- Short-term (2 months): Phase 2.5 A/B/C
- Long-term (1 month): Phase 2.6 A/B/C

---

## Phase 2.3: Relativistic Mass Validation ✓ COMPLETE

**Result**: **76% pass rate** (38/50 tests passed)

**Successes**:
- ✓ v≤0.3c: 100% pass rate, <5% error
- ✓ Grid convergence validated (256×256 < 2% error)
- ✓ Operator splitting optimal at N=10
- ✓ Conservation laws excellent (norm <0.5%, energy <1%)

**Failures**:
- ✗ v=0.5c: 8.7% error (marginal)
- ✗ v=0.7c: 18.7% error (systematic underestimation)
- ✗ N=1 "mass freeze": γ≈1 regardless of velocity

**Scientific Impact**: Discovered SMFT regime limits for relativistic physics

---

## Phase 2.4: Breakdown Investigation (Immediate - 1 month)

**Goal**: Understand WHY Scenario 2.3 failed at v>0.5c

### 2.4A: Velocity Threshold Identification
**Question**: At what exact velocity does breakdown occur?

- **Config**: `scenario_2.4A_velocity_threshold.yaml`
- **Velocities**: [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]c (fine-grained)
- **Grid**: 64×64 (fast)
- **N**: 10 only (optimal)
- **Duration**: ~2 hours compute

**Deliverables**:
- v_critical interpolated to ±0.05c
- γ error scaling law (linear/quadratic/exponential?)
- Safe velocity range recommendation

---

### 2.4B: R-Field Dynamics Study
**Question**: Why does N=1 freeze the mass field?

- **Config**: `scenario_2.4B_R_field_dynamics.yaml`
- **Velocities**: [0.5, 0.7]c (breakdown onset + failure)
- **Grid**: 128×128 (high spatial resolution)
- **N**: [1, 10] (direct comparison)
- **Duration**: ~4 hours compute
- **Snapshots**: 9 full R(x,y,t) spatial fields

**Measurements**:
- R(x,y,t) evolution at particle location
- ∇R along trajectory
- Response timescale τ_R
- Feedback strength β = d(ln R)/d(ln σ)

**Hypothesis Test**: Does N=1 freeze R-field, preventing mass modulation?

**Deliverables**:
- R(x,y,t) evolution movies (N=1 vs N=10)
- τ_R measurements
- Mechanism document: "Why Born-Oppenheimer required"

---

### 2.4C: Ultra-Relativistic + Ultra-Fine Grid
**Question**: Is v=0.7c failure due to insufficient resolution?

- **Config**: `scenario_2.4C_ultra_relativistic.yaml`
- **Velocities**: [0.7, 0.8, 0.9]c
- **Grid**: 512×512 (Δx = 0.195 ℓ_P, 15 points per σ!)
- **N**: [10, 100]
- **Duration**: ~24 hours compute (large grid)

**Grid Resolution Test**:
- 64×64: σ ≈ 2 grid points (coarse)
- 512×512: σ ≈ 15 grid points (ultra-fine)

**Outcome Decision Tree**:
- **If 512×512 fixes v=0.7c error**: Grid resolution was the problem → Use finer grids
- **If error persists**: Fundamental SMFT limitation → Need Klein-Gordon or accept v<0.4c validity

**Deliverables**:
- γ(v) at [0.7, 0.8, 0.9]c on ultra-fine grid
- Grid resolution study report
- Fundamental limits document

---

## Phase 2.5: Alternative Formulations (Short-term - 2 months)

**Goal**: Test alternative field theories and reference frames

### 2.5A: Klein-Gordon Comparison
**Question**: Does relativistic scalar field perform better than non-relativistic Dirac?

- **Config**: `scenario_2.5A_klein_gordon_comparison.yaml`
- **Field Types**: Klein-Gordon vs Dirac (direct comparison)
- **Velocities**: [0.0, 0.3, 0.5, 0.7, 0.9]c (full range)
- **Grid**: 128×128
- **Duration**: ~8 hours (2 field types)

**Physics**: Klein-Gordon: (∂²_t - ∇² + m²)φ = 0 with m(x) = Δ·R(x)

**Expected Outcome**:
- If KG shows better γ agreement at v>0.5c → Dirac operator is the bottleneck
- If KG shows same errors → Mass formula m=Δ·R is the limitation

**Deliverable**: KG vs Dirac comparison report

---

### 2.5B: Boosted Frame Analysis
**Question**: Does m(v)=γ·Δ·R hold in moving reference frames?

- **Config**: `scenario_2.5B_boosted_frame_analysis.yaml`
- **Velocities**: [0.3, 0.5, 0.7]c
- **Analysis**: Lorentz transform coordinates (x',t') = Λ(x,t)
- **Grid**: 128×128
- **Duration**: ~6 hours

**Tests**:
- Measure R'(x',t') in boosted frame
- Verify m' = m (mass invariant)
- Check E'² - p'² = E² - p² (Lorentz invariant)

**Expected**: In particle rest frame (v'=0), should see γ'=1 and m'=m_rest

**Deliverable**: Lorentz covariance validation report

---

### 2.5C: Ultra-Relativistic Limit
**Question**: What happens at v→c?

- **Config**: `scenario_2.5C_ultra_relativistic_limit.yaml`
- **Velocities**: [0.90, 0.95, 0.99]c (extreme test)
- **Grid**: 256×256
- **Duration**: ~12 hours
- **γ factors**: 2.29, 3.20, 7.09 (up to 609% mass increase!)

**Tests**:
- E≈pc limit (massless regime)
- Time dilation (evolution slows as 1/γ?)
- Lorentz contraction (σ_∥ = σ/γ → σ/7 at v=0.99c!)
- Causality (v<c respected?)

**Critical Question**: Does SMFT fail completely or gracefully approach E=pc?

**Deliverable**: v→c limit analysis, SMFT validity boundary

---

## Phase 2.6: Advanced Physics (Long-term - 1 month)

**Goal**: Extend SMFT to full relativistic + EM framework

### 2.6A: Energy-Momentum Relation
**Question**: Does SMFT satisfy Einstein's E²=(mc²)²+(pc)²?

- **Config**: `scenario_2.6A_energy_momentum_relation.yaml`
- **Velocities**: [0.0 to 0.9]c in 0.1c steps (10 velocities)
- **Grids**: [64, 128, 256] (convergence)
- **Duration**: ~15 hours

**Measurements**:
- E² - p² = m₀²c⁴ invariant
- Rest mass from E-p: m₀ = √(E²-p²)/c²
- Kinetic energy: T = (γ-1)mc²
- Low-v limit: T ≈ ½mv²?
- High-v limit: E ≈ pc?

**Deliverable**: E(p) dispersion curve, comprehensive E-p validation

---

### 2.6B: Electromagnetic Coupling
**Question**: Can A_μ = ∂_μθ produce emergent electromagnetism?

- **Config**: `scenario_2.6B_electromagnetic_coupling.yaml`
- **NEW PHYSICS**: Couple ∇θ to spinor as electromagnetic potential
- **Velocities**: [0.3, 0.5]c
- **Duration**: ~4 hours
- **Status**: **EXPLORATORY** (requires new implementation)

**Tests**:
- Compute A_μ = (φ, A) from phase gradients
- Extract E = -∇φ - ∂_tA, B = ∇×A
- Test Lorentz force F = q(E + v×B)
- Validate Maxwell equations

**MAJOR DISCOVERY POTENTIAL**: If ∇θ naturally couples as EM field, this validates SMFT as unified field theory

**Deliverable**: Electromagnetic emergence report (high-risk, high-reward)

---

### 2.6C: Soliton Stability
**Question**: Does high boost destroy self-trapping?

- **Config**: `scenario_2.6C_soliton_stability.yaml`
- **Velocities**: [0.0, 0.3, 0.5, 0.7, 0.9]c
- **Grid**: 256×256
- **Duration**: **50,000 steps** (~30 hours - long evolution!)
- **Domain**: 200 ℓ_P (large for long trajectory)

**Measurements**:
- Wavepacket width σ(t) (growth = dispersion)
- Localization L = ∫|Ψ|⁴ / (∫|Ψ|²)²
- Binding energy E_bind = ∫Ψ†(Δ·R)Ψ

**Stability Criteria**:
- STABLE: σ(t) ≈ constant, L ≈ constant
- DISPERSING: σ(t) ∝ √t, L → 0
- COLLAPSE: σ(t) → 0 (unlikely)

**Hypothesis**: At v→c, particle outruns R-field response → soliton disperses

**Deliverable**: Soliton stability phase diagram, v_critical for self-trapping

---

## Test Matrix Summary

| Phase | Scenario | Velocities | Grids | N-values | Compute Time | Status |
|-------|----------|------------|-------|----------|--------------|--------|
| **2.4A** | Velocity threshold | 8 | 128² | 1 | 2 hr | Ready |
| **2.4B** | R-field dynamics | 2 | 128² | 2 | 4 hr | Ready |
| **2.4C** | Ultra-fine grid | 3 | 512² | 2 | 24 hr | Ready |
| **2.5A** | Klein-Gordon | 5 | 128² | 1×2 fields | 8 hr | ⚠️ Solver ✓, Observables blocked |
| **2.5B** | Boosted frame | 3 | 128² | 1 | 6 hr | ✓ READY (Sprint 2 complete) |
| **2.5C** | v→c limit | 3 | 256² | 2 | 12 hr | Ready |
| **2.6A** | E-p relation | 10 | 3 grids | 1 | 15 hr | Ready |
| **2.6B** | EM coupling | 2 | 128² | 1 | 4 hr | ⛔ DEFERRED to Phase 5 |
| **2.6C** | Soliton stability | 5 | 256² | 1 | 30 hr | Ready |
| **Total** | 9 scenarios | 41 configs | - | - | **105 hr** | 7/9 ready, 1 blocked, 1 deferred |

---

## Implementation Requirements

### Already Implemented ✓
- Boosted Gaussian initialization (Scenario 2.3)
- Velocity sweep infrastructure
- Grid-independent physical parameters
- Operator splitting (N=1, 10, 100)
- Observables computation (E, p, m_eff, γ)

### Sprint 2 Complete (2025-12-23) ✓
1. **Klein-Gordon Evolution** (2.5A) - ✅ COMPLETE
   - ✅ Implement (∂²_t - ∇² + m²)φ = 0
   - ✅ Couple to R-field: m(x) = Δ·R(x)
   - ⚠️ **BLOCKER**: Compute KG observables (ObservableComputer doesn't support KG fields)
   - **Fix**: 4-8 hours to implement ⟨x⟩, ⟨p⟩, γ for Klein-Gordon

2. **Lorentz Transform** (2.5B) - ✅ COMPLETE
   - ✅ Coordinate transform: (x',t') = Λ(x,t)
   - ✅ Boost R-field: R'(x',t')
   - ✅ Compute boosted-frame observables
   - ✅ Files: src/physics/LorentzTransform.{h,cpp} (645 lines)

3. **Electromagnetic Coupling** (2.6B) - ⛔ DEFERRED TO PHASE 5
   - **Decision**: DEFER (3-4 weeks, 20-30% failure risk, delays Phase 2 by 4-6 weeks)
   - **Rationale**: 19/20 scenarios = 95% validation sufficient for publication
   - **Documentation**: Feasibility report in notepad (note_id: ccde6fcd)
   - **Future**: Phase 5 exploratory research (separate paper)

### Analysis Extensions 📊
- Interpolation for v_critical (2.4A)
- R-field spatial visualization (2.4B)
- Wavepacket width tracking (2.6C)
- E(p) dispersion curves (2.6A)

---

## Timeline & Milestones

### Month 1: Phase 2.4 (Immediate)
**Week 1**: Run 2.4A (velocity threshold)
- Identify v_critical ≈ 0.45-0.55c (estimate)
- Document safe velocity range

**Week 2**: Run 2.4B (R-field dynamics)
- Generate R(x,y,t) movies
- Measure τ_R, explain N=1 freeze

**Week 3-4**: Run 2.4C (ultra-fine grid)
- 512×512 at v=[0.7, 0.8, 0.9]c
- Final verdict on grid resolution vs fundamental limit

**Milestone 1**: Complete understanding of Scenario 2.3 failures

---

### Month 2-3: Phase 2.5 (Short-term)
**Week 5-6**: Implement & run 2.5A (Klein-Gordon)
- Code KG evolution
- Compare KG vs Dirac at all velocities

**Week 7-8**: Implement & run 2.5B (boosted frames)
- Code Lorentz transforms
- Validate frame independence

**Week 9-10**: Run 2.5C (v→c limit)
- Test extreme velocities [0.9, 0.95, 0.99]c
- Map SMFT breakdown

**Milestone 2**: Alternative formulations tested, SMFT boundaries mapped

---

### Month 4: Phase 2.6 (Long-term)
**Week 11-12**: Run 2.6A (E-p relation)
- Comprehensive E(p) curve
- Validate Einstein's relation

**Week 13**: Implement & run 2.6B (EM coupling)
- HIGH RISK: May discover emergent EM or fail
- Exploratory physics

**Week 14**: Run 2.6C (soliton stability)
- Long evolution (50k steps)
- Stability phase diagram

**Milestone 3**: Complete relativistic SMFT characterization

---

## Success Criteria

### Phase 2.4 Success ✓
- v_critical identified to ±0.05c
- R-field freeze mechanism understood
- Grid resolution vs fundamental limit resolved

### Phase 2.5 Success ✓
- Klein-Gordon comparison complete
- Lorentz covariance validated (or refuted)
- SMFT validity range documented

### Phase 2.6 Success ✓
- E²=(mc²)²+(pc)² validated
- EM emergence tested (success OR null result both valuable)
- Soliton stability understood

### Overall Phase 2 Success ✓
- **Complete map** of SMFT relativistic regime
- **Documented limitations** (if v<0.4c only, that's acceptable!)
- **Next steps** defined (Phase 3: Gravity OR pivot to QFT)

---

## Risk Assessment

### Low Risk (Will Complete)
- 2.4A, 2.4C, 2.5C, 2.6A, 2.6C (all use existing code)
- Expected timeline: 1-2 months

### Medium Risk (Implementation Needed)
- 2.5A (Klein-Gordon): Straightforward PDE, 1 week implementation
- 2.5B (Lorentz transforms): Standard SR, 1 week implementation

### High Risk (Exploratory)
- 2.6B (EM coupling): Novel physics, may not work as expected
- 2.4B (R-field dynamics): May not explain N=1 freeze

---

## Recommended Execution Order

1. **2.4A** (2 hours) → Immediate insight into breakdown
2. **2.4C** (24 hours) → Resolve grid vs fundamental question
3. **2.4B** (4 hours) → Understand N=1 mechanism
4. **2.5C** (12 hours) → Map v→c limit
5. **2.6A** (15 hours) → Comprehensive E-p validation
6. **2.5A** (8 hours + 1 week impl.) → KG comparison
7. **2.6C** (30 hours) → Long soliton evolution
8. **2.5B** (6 hours + 1 week impl.) → Frame independence
9. **2.6B** (4 hours + 2 weeks impl.) → EM exploration (last - high risk)

**Critical Path**: 2.4A → 2.4C → (decision point) → rest in parallel

---

## Configuration Files

All 9 test configurations ready:

**Phase 2.4 (Immediate)**:
- `config/scenario_2.4A_velocity_threshold.yaml`
- `config/scenario_2.4B_R_field_dynamics.yaml`
- `config/scenario_2.4C_ultra_relativistic.yaml`

**Phase 2.5 (Short-term)**:
- `config/scenario_2.5A_klein_gordon_comparison.yaml`
- `config/scenario_2.5B_boosted_frame_analysis.yaml`
- `config/scenario_2.5C_ultra_relativistic_limit.yaml`

**Phase 2.6 (Long-term)**:
- `config/scenario_2.6A_energy_momentum_relation.yaml`
- `config/scenario_2.6B_electromagnetic_coupling.yaml`
- `config/scenario_2.6C_soliton_stability.yaml`

---

**Phase 2.3 Status**: ✓ COMPLETE (76% pass - regime discovery)
**Sprint 2 Status (2025-12-23)**: ✓ COMPLETE
  - ✅ Lorentz transforms implemented (production-ready)
  - ✅ Klein-Gordon solver complete (observables blocked)
  - ✅ EM coupling evaluated (deferred to Phase 5)
**Phase 2.4-2.6 Status**: ⚙ 7/9 scenarios ready, 1 blocked (2.5A observables), 1 deferred (2.6B)
**Total Estimated Time**: 3 months (101 compute hours, EM deferred saves 4 weeks)

**Next Action**: Sprint 3 - Fix Klein-Gordon observables + Execute Phase 2.4A velocity threshold

---

**Author**: SMFT Research Team
**Date**: 2025-12-19
**Document**: Complete Phase 2 roadmap based on 2.3 results
