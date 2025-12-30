# Falsification Criteria for SMFT Metric Derivation

**Date**: 2025-12-29
**Purpose**: Establish clear, testable criteria to falsify (or verify) emergent spacetime metric claims
**Commitment**: Honest documentation of failures, no hand-waving, no retreat to analogy

---

## Executive Summary

This document defines **concrete, quantitative tests** that would falsify the claim that SMFT generates an emergent spacetime metric. We commit to:

1. **No moving goalposts**: Criteria are fixed before calculation
2. **Quantitative thresholds**: Numerical agreement/disagreement defined precisely
3. **Honest failure reporting**: If criteria fail, we document it clearly
4. **No hand-waving escapes**: Cannot retreat to "effective theory" or "analogy" if rigorous derivation fails

---

## 1. Falsification Hierarchy

### Level 1: Metric Existence
**Claim**: A unique functional g_μν = f[R, ∂R, θ, ∂θ] exists

**Falsification Criteria**:
- **FAIL A1**: Multiple inequivalent metrics reproduce SMFT physics
- **FAIL A2**: No metric reproduces key SMFT features
- **FAIL A3**: Metric requires non-local terms (not geometric)

**Test**: Compare conformal, acoustic, and effective action metrics. If they give different predictions for same SMFT configuration, claim is falsified.

**Status**: **FAILED A1** - Multiple candidate metrics (conformal vs. acoustic) with different properties. No unique derivation.

---

### Level 2: Consistency with SMFT Implementation
**Claim**: Derived metric reproduces SMFT Dirac evolution

**Falsification Criteria**:
- **FAIL B1**: Metric predicts spin connection terms ∂_μR not present in SMFT code
- **FAIL B2**: Mass scaling is wrong (metric gives m ~ 1/R, SMFT implements m ~ R)
- **FAIL B3**: Dispersion relation mismatch: metric prediction vs. SMFT evolution

**Test**: Compare Dirac equation in curved spacetime (with derived metric) against flat-space SMFT implementation.

**Status**:
- **FAILED B1** - Conformal metric g_μν = R²η_μν requires spin connection Γ_μ ~ ∂_μ ln R not in code
- **FAILED B2** - Conformal metric predicts m_conf ~ 1/R, SMFT implements m = Δ·R (opposite scaling)
- **PENDING B3** - Acoustic metric dispersion needs numerical verification

---

### Level 3: General Covariance
**Claim**: SMFT exhibits general covariance (coordinate invariance)

**Falsification Criteria**:
- **FAIL C1**: Theory depends on preferred coordinate frame (Kuramoto lattice)
- **FAIL C2**: Gamma matrices are flat-space γ^μ_flat, not curved γ^μ = e^μ_a γ^a
- **FAIL C3**: Coordinate transformations change physics

**Test**: Apply coordinate transformation x^μ → x'^μ and check if physics is invariant.

**Status**: **FAILED C1, C2** - SMFT uses flat-space Dirac matrices and Kuramoto lattice provides preferred frame. Not generally covariant.

---

### Level 4: Einstein Equations
**Claim**: Derived metric satisfies G_μν = 8πG_eff T_μν

**Falsification Criteria**:
- **FAIL D1**: Different components give different G_eff (inconsistent)
- **FAIL D2**: G_eff is negative or zero (unphysical)
- **FAIL D3**: G_eff varies with position by >10% (not constant)
- **FAIL D4**: Einstein tensor diverges or is undefined

**Test**: Compute G_μν and T_μν numerically, extract G_eff from each component, check consistency.

**Quantitative Threshold**:
```
σ(G_eff) / ⟨G_eff⟩ < 0.10   (10% variation allowed)
```

**Status**: **PENDING** - Calculation not yet performed. Framework in `EINSTEIN_EQUATIONS.md`.

---

### Level 5: Physical Predictions
**Claim**: Metric reproduces observed SMFT phenomena

**Falsification Criteria**:
- **FAIL E1**: Geodesics deviate from fermion trajectories by >5% RMS
- **FAIL E2**: Lensing angles don't match fermion deflection
- **FAIL E3**: Metric predicts superluminal propagation (causality violation)
- **FAIL E4**: Time dilation factor disagrees with SMFT implementation

**Test**: Integrate geodesics in derived metric and compare with SMFT wavepacket evolution.

**Quantitative Thresholds**:
```
RMS(x_geodesic - x_fermion) / L_domain < 0.05   (5% spatial deviation)
|t_metric - t_SMFT| / t_total < 0.05            (5% time deviation)
```

**Status**: **PENDING** - Requires geodesic integrator and SMFT trajectory data.

---

## 2. Detailed Falsification Tests

### Test 2.1: Conformal Metric Falsification

**Metric**: g_μν = R²η_μν

**Prediction 1**: Mass term in Dirac equation should be m_eff = m₀/R

**SMFT Reality**: Mass term is m_eff = Δ·R

**Comparison**:
```
R = 0.5  →  Metric predicts m = 2m₀, SMFT has m = 0.5Δ
```

**Verdict**: **FALSIFIED** - Wrong mass scaling by factor of 4 in defect regions.

---

**Prediction 2**: Spin connection should appear:
```
Γ_μ = (3i/4)(γ_μγ^ν - δ_μ^ν)(∂_ν ln R)
```

**SMFT Reality**: No spin connection in `DiracEvolution.cpp::applyPotentialStep()`

**Code Check**:
```cpp
// Lines 219-229: Only mass term, no ∂_μR terms
for (uint32_t i = 0; i < _N_points; ++i) {
    float m = mass_field[i];  // m = Δ·R
    std::complex<float> phase_upper = std::exp(-_beta_sign*m*dt*1if);
    std::complex<float> phase_lower = std::exp(+_beta_sign*m*dt*1if);
    // No gradient terms!
}
```

**Verdict**: **FALSIFIED** - Required spin connection is absent.

---

### Test 2.2: Acoustic Metric Consistency

**Metric**: g_μν = R²[η_μν + h_μν(∇θ)]

**Prediction 1**: Time dilation dτ = R√(1-v²)dt

**SMFT Reality**: dτ = R·dt (time_dilation_mode)

**Comparison**: For v² << 1, acoustic predicts dτ ≈ R(1 - v²/2)dt

**Quantitative Check**:
```
For v = 0.3: Acoustic gives dτ = 0.955R·dt
            SMFT gives dτ = R·dt
            Discrepancy: 4.5%
```

**Verdict**: **MARGINAL AGREEMENT** - Within 5% for subsonic flows, fails for v → 1.

---

**Prediction 2**: Frame-dragging term g₀ᵢ = -R²v_i

**SMFT Reality**: Not explicitly implemented, but could emerge from ∇θ coupling

**Test Method**:
1. Compute fermion angular momentum in vortex
2. Compare with geodesic angular momentum in acoustic metric
3. Check for co-rotation signature

**Status**: **PENDING** - Requires numerical implementation.

---

### Test 2.3: Einstein Equation Falsification

**Setup**: Compute G_μν and T_μν for SMFT configuration at t=10, grid 64×64

**Extract G_eff from each component**:
```
G_eff,(00) = G₀₀ / (8π T₀₀)
G_eff,(01) = G₀₁ / (8π T₀₁)
... (6 independent components)
```

**Consistency Check**:
```python
G_eff_values = [G_eff from each component]
mean_G = np.mean(G_eff_values)
std_G = np.std(G_eff_values)
variation = std_G / mean_G

if variation < 0.10:
    PASS
else:
    FAIL - Report discrepancy
```

**Falsification Scenarios**:

1. **G_eff has wrong sign**:
   ```
   If mean_G < 0:  FATAL FAILURE (negative gravity)
   ```

2. **G_eff varies wildly**:
   ```
   If variation > 0.50:  FAIL (not Einstein gravity)
   If variation > 0.10:  MARGINAL (modified gravity)
   ```

3. **G_eff is unphysical scale**:
   ```
   G_Planck = 1 (in Planck units)
   If |mean_G - 1| > 10:  FAIL (wrong quantum scale)
   ```

---

### Test 2.4: Geodesic Deviation

**Setup**: Initialize fermion wavepacket at (x₀, y₀) with momentum (p_x, p_y)

**SMFT Evolution**: Run for t = 100 time units

**Metric Prediction**: Integrate geodesic equation:
```
d²x^μ/dλ² + Γ^μ_ρσ (dx^ρ/dλ)(dx^σ/dλ) = 0
```
with initial conditions matching fermion.

**Comparison**:
```
Δx(t) = |x_fermion(t) - x_geodesic(t)|
RMS_deviation = √(⟨Δx²⟩_t)
```

**Falsification Thresholds**:
```
If RMS_deviation / L_domain > 0.05:  FAIL
If RMS_deviation / L_domain > 0.20:  CATASTROPHIC FAIL
```

**Expected Outcomes**:
- RMS < 5%: SUCCESS - Metric accurately describes fermion propagation
- RMS 5-20%: PARTIAL - Metric is approximate effective theory
- RMS > 20%: FAIL - Metric does not describe SMFT physics

---

### Test 2.5: Causality Violation Check

**Setup**: Compute light cones from metric:
```
ds² = 0  →  g_μν dx^μ dx^ν = 0
```

For acoustic metric:
```
-R²(1-v²)dt² - 2R²v·dx dt + R²|dx|² = 0
```

Solving for |dx/dt|_max:

```
v_light = [v ± √(1 + v²)] / (1 - v²)
```

**Falsification Criterion**:
```
If v_light > c (=1 in Planck units):  FAIL - Superluminal propagation
```

**Critical Points**:
- When v → 1: Denominato → 0, v_light → ∞ (sonic horizon)
- When v > 1: v_light complex (unphysical)

**Physical Interpretation**:
- Acoustic horizons at |∇θ| = Δ
- No superluminal travel if v < c everywhere

**Test**: Scan SMFT phase field, check max|∇θ|/Δ < 1

---

## 3. Acceptance Criteria (Inverse of Falsification)

For the metric derivation to be considered **successful**, ALL of the following must hold:

### 3.1 Uniqueness
✅ One functional form g_μν = f[R, ∂R, θ, ∂θ] reproduces all SMFT physics
❌ Multiple inequivalent metrics work equally well

### 3.2 Consistency
✅ Metric reproduces SMFT Dirac equation (no extra terms)
✅ Mass scaling m_eff = Δ·R emerges naturally
✅ Time dilation dτ = R·dt is exact (not approximate)

### 3.3 Einstein Equations
✅ G_μν = 8πG_eff T_μν with constant G_eff ≈ G_Planck
✅ Variation of G_eff across components < 10%
✅ G_eff > 0 (attractive gravity)

### 3.4 Geodesics
✅ Fermion trajectories match metric geodesics within 5% RMS
✅ Deflection angles around defects agree with geodesic deviation
✅ No causality violations (all v < c)

### 3.5 Theoretical Consistency
✅ Derivation from SMFT action via functional methods (not analogy)
✅ All mathematical steps justified
✅ No hand-waving or "by analogy" arguments

---

## 4. Current Status Summary

| Test | Status | Result |
|------|--------|--------|
| Unique metric | ❌ FAILED | Multiple candidates (conformal vs. acoustic) |
| Conformal metric | ❌ FAILED | Wrong mass scaling, spurious spin connection |
| Acoustic metric | ⏳ PENDING | Promising but not derived from first principles |
| General covariance | ❌ FAILED | SMFT uses preferred frame, flat γ matrices |
| Einstein equations | ⏳ PENDING | Calculation framework ready, not executed |
| Geodesic agreement | ⏳ PENDING | Requires numerical implementation |
| Causality | ✅ LIKELY PASS | Acoustic horizons form at v=c, no superluminal |

**Overall Assessment**:
- **2 FAILURES** (uniqueness, general covariance)
- **4 PENDING** (need calculation/simulation)
- **1 LIKELY PASS** (causality structure)

---

## 5. Failure Modes and Interpretation

### 5.1 Failure Mode A: No Metric Exists

**Indicators**:
- All candidate metrics fail geodesic test
- Einstein equations unsatisfiable
- Coordinate-dependent physics

**Interpretation**: SMFT is a **non-geometric** effective theory. Synchronization generates emergent phenomena analogous to geometry but not describable by a spacetime metric.

**Scientific Implication**: New physics beyond general relativity.

**Action**: Document non-geometric features, publish as counterexample to emergent gravity.

---

### 5.2 Failure Mode B: Multiple Metrics Work

**Indicators**:
- Conformal and acoustic metrics both pass geodesic test
- Different G_eff but both positive
- Gauge freedom larger than coordinate choice

**Interpretation**: SMFT has **gauge redundancy** in metric description. Need additional physical principle to select unique metric.

**Scientific Implication**: Underdetermined theory requiring extra constraints.

**Action**: Search for selection principle (e.g., minimal action, maximal entropy).

---

### 5.3 Failure Mode C: Approximate Metric Only

**Indicators**:
- Metric works for small R fluctuations (R ≈ 1)
- Fails near defects (R → 0)
- Valid only in low-velocity limit

**Interpretation**: SMFT has emergent geometry in **limited regime**. Full theory is non-geometric.

**Scientific Implication**: Effective field theory with cutoff scale.

**Action**: Identify regime of validity, compute corrections.

---

### 5.4 Failure Mode D: Wrong Physics

**Indicators**:
- Metric predicts phenomena not seen in SMFT
- Einstein equations fail completely
- G_eff negative or divergent

**Interpretation**: The R-field is **not the correct geometric variable**. Need different field (e.g., R + derivatives, or different combination).

**Scientific Implication**: Misidentified emergent geometry.

**Action**: Try alternative field combinations, topology, or abandon metric approach.

---

## 6. No-Retreat Commitment

We commit to **NOT** retreating to the following escape hatches if rigorous tests fail:

### ❌ Forbidden Escape #1: "It's an Analogy"
**Retreat**: "The acoustic metric is an *analogy* to help understand SMFT, not a literal derivation."

**Why Forbidden**: Directive explicitly requires rigorous derivation, not analogy. If derivation fails, we admit failure.

### ❌ Forbidden Escape #2: "Effective Theory Excuse"
**Retreat**: "Well, it's an effective low-energy theory, so some mismatch is expected."

**Why Forbidden**: Either the metric emerges or it doesn't. Approximate agreement is failure of unique derivation.

### ❌ Forbidden Escape #3: "Numerical Uncertainty"
**Retreat**: "The 50% discrepancy could be numerical error."

**Why Forbidden**: We have clear thresholds (5%, 10%, 20%). Beyond that, it's physical disagreement.

### ❌ Forbidden Escape #4: "Higher-Order Corrections"
**Retreat**: "Adding more terms would fix the mismatch."

**Why Forbidden**: Unlimited freedom to add terms makes theory unfalsifiable. We test the proposed metric as stated.

### ❌ Forbidden Escape #5: "Modified Gravity"
**Retreat**: "It's not Einstein gravity, but a modified theory."

**Why Forbidden**: Directive asked for *Einstein* equations. Modified gravity is a different claim requiring separate justification.

---

## 7. Success Looks Like

If the metric derivation succeeds, we will have:

### 7.1 Concrete Deliverables
1. **Explicit formula**: g_μν = f[R, ∂R, θ, ∂θ] with all functions defined
2. **Derivation**: From SMFT action via functional integration (documented)
3. **Einstein equations**: G_μν = 8πG_eff T_μν verified numerically with G_eff ≈ 1
4. **Geodesic agreement**: RMS deviation < 5%
5. **Physical interpretation**: Clear explanation of emergent gravity mechanism

### 7.2 Scientific Impact
- **Emergent GR**: First concrete model of spacetime geometry from quantum synchronization
- **Quantum gravity**: Connection between Planck scale (Δ) and gravity (G_eff)
- **Testable predictions**: SMFT makes specific predictions for gravitational phenomena
- **New paradigm**: Gravity as collective synchronization, not fundamental interaction

### 7.3 Publication
- **Title**: "Emergent Spacetime from Synchronization: Rigorous Derivation of Einstein Equations in SMFT"
- **Journal**: Physical Review Letters (if successful)
- **Impact**: Paradigm shift in understanding gravity

---

## 8. Failure Looks Like

If the metric derivation fails, we will have:

### 8.1 Honest Documentation
1. **What failed**: Clear statement of which tests failed
2. **Why it failed**: Physical/mathematical reason for failure
3. **What this means**: Implications for SMFT and emergent geometry

### 8.2 Scientific Integrity
- **No excuses**: Admit failure without hand-waving
- **Lessons learned**: What the failure teaches us
- **Alternative paths**: Other approaches to explore

### 8.3 Publication
- **Title**: "Limits of Emergent Geometry: Why SMFT Does Not Generate Spacetime Metrics"
- **Journal**: Physical Review D (negative results valuable)
- **Impact**: Clarifies boundary between geometric and non-geometric emergent phenomena

---

## 9. Timeline for Verdict

**Week 1 (Complete)**: Analytical framework and falsification criteria
**Week 2**: Symbolic Einstein equation calculation
**Week 3**: Numerical geodesic and stress-energy computation
**Week 4**: Final verdict and report

**Decision Point**: By end of Week 4, declare SUCCESS or FAILURE based on quantitative criteria.

**No Extensions**: 4 weeks is sufficient. If undecided after 4 weeks, default to FAILURE (burden of proof on success claim).

---

## 10. Final Commitment

We commit to the following:

1. **Honest reporting**: Report failures as clearly as successes
2. **No hand-waving**: Every mathematical step justified
3. **Quantitative thresholds**: Use numerical criteria, not qualitative judgments
4. **No goalpost moving**: Criteria fixed now, not adjusted post-calculation
5. **Accept failure**: If tests fail, accept SMFT does not generate metric

**Signed** (metaphorically):
The SMFT Research Team
Date: 2025-12-29

---

## Appendix: Quick Reference Card

### Falsification Quick Checks

```
FAIL if:
- Conformal metric used (already falsified)
- "By analogy" appears in derivation
- Spin connection ∂_μR present but not in code
- Mass scaling is 1/R not R
- Multiple metrics work equally well
- G_eff varies by >10%
- G_eff < 0
- Geodesic RMS > 5%
- Any superluminal propagation

PASS if:
- Unique metric derived from action
- m_eff = Δ·R emerges naturally
- G_μν = 8πG_eff T_μν with constant G_eff ≈ 1
- Geodesic RMS < 5%
- All math rigorous, no hand-waving
```

**When in doubt, choose FAIL over PASS. Science rewards honesty.**
