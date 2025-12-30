# SMFT Falsification Criteria
## Clear Thresholds for Experimental Rejection

**Document Version**: 1.0
**Date**: December 29, 2025

---

## Principle: Falsifiability is Essential for Science

**Karl Popper's Criterion**: A theory is scientific if and only if it makes predictions that can be proven wrong by experiment.

**SMFT Commitment**: We provide explicit, quantitative criteria for falsifying SMFT. If experiments meet these thresholds, **the theory is wrong** and must be abandoned or fundamentally revised.

This document defines:
1. **What observations would falsify SMFT**
2. **Statistical confidence required** (typically 5σ)
3. **Alternative explanations that must be ruled out**

---

## Area 1: BEC Phonon Scattering Near Vortex

### SMFT Prediction
```
c_eff(r) / c_s = R(r) = tanh(r/ξ)

At r = ξ: c_eff/c_s = 0.76 (24% suppression)
```

### Falsification Criterion

**SMFT is FALSIFIED if**:
```
c_eff(r) / c_s = 1.00 ± 0.05  for all r > ξ
```

**Confidence**: 3σ over at least 5 independent measurements at r = 0.5ξ, 1ξ, 2ξ

**Statistical Test**:
```
Null hypothesis H₀: c_eff(r) = c_s (standard BEC)
Alternative H₁: c_eff(r) = c_s · tanh(r/ξ) (SMFT)

χ² test: Δχ² = χ²(H₀) - χ²(H₁)

If Δχ² < 7.8 (p > 0.05, not significant)
→ CANNOT reject H₀
→ SMFT FALSIFIED
```

### Alternative Explanations to Rule Out

Before declaring SMFT falsified, must verify:

1. **Vortex is present**: Dark core visible in density imaging
   - If no vortex: experiment invalid

2. **Thermal phonons negligible**: T < 0.5 T_c
   - If T too high: thermal depletion masks signal

3. **Imaging resolution sufficient**: Δr < 0.5ξ
   - If poor resolution: cannot resolve spatial structure

4. **Phonon wavelength appropriate**: λ > 3ξ
   - If λ ~ ξ: scattering theory breaks down (non-linear regime)

**Only if all controls pass AND c_eff(r) = const → SMFT falsified**

---

## Area 2: Casimir Force Modification

### SMFT Prediction
```
F_SMFT / F_standard = ⟨R²⟩

If R = 0.95: F_SMFT/F_standard = 0.90 (10% reduction)
```

### Falsification Criterion

**SMFT is FALSIFIED if**:
```
F_measured(d, material) / F_Lifshitz(d, material) = 1.000 ± 0.010

for all:
  - separations d = 50-500 nm
  - materials (Au, Si, graphene, ...)
  - geometries (flat, corrugated, ...)
```

**Confidence**: 5σ with systematic error < 1%

**Statistical Test**:
```
Fit: F_measured = F_Lifshitz · (1 + δ)

SMFT predicts: δ = ⟨R²⟩ - 1 ≈ -0.10 (for R=0.95)

If |δ_measured| < 0.01 at 5σ
→ NO R-field coupling detected
→ SMFT FALSIFIED
```

### Alternative Explanations to Rule Out

1. **Electrostatic patches**: Apply nulling voltage
   - Measure F(V_applied) → Find V_null where ∂F/∂V = 0

2. **Surface roughness**: Use atomically flat surfaces
   - AFM characterization: RMS roughness < 1 nm

3. **Temperature effects**: Stabilize T ± 0.1 K
   - Measure F(T) → Check T-dependence matches Lifshitz theory

4. **Dielectric function uncertainty**: Independently measure ε(ω)
   - Ellipsometry or reflectance spectroscopy

**Only if all systematics controlled AND F/F_Lifshitz = 1.00 → SMFT falsified**

---

## Area 3: CMB Non-Gaussianity & Cosmic Strings

### SMFT Predictions

**Observable 1: f_NL (Non-Gaussianity)**
```
If SMFT phase transition at T > T_recombination:
f_NL ~ 1 to 10 (defect network)
```

**Observable 2: Gμ (String Tension)**
```
If winding number W conserved:
Gμ/c² ~ 10⁻⁷ to 10⁻⁶
```

**Observable 3: n_s Running**
```
String contribution: dn_s/dlnk ≠ 0
```

### Falsification Criteria

**SMFT COSMOLOGY is FALSIFIED if ANY ONE of**:

#### (A) f_NL Consistent with Zero
```
|f_NL^local| < 1 at 5σ (CMB-S4)

Current: f_NL = -0.9 ± 5.1 (Planck 2018)
Future: σ(f_NL) ~ 1 (CMB-S4)

If CMB-S4 measures: f_NL = 0.0 ± 1.0
→ 5σ upper limit: |f_NL| < 5
→ MARGINAL (depends on transition details)

If |f_NL| < 1 at 5σ:
→ NO defect network at recombination
→ SMFT transition either:
    (a) did not occur, or
    (b) occurred at T < T_recombination (too late for CMB signal)
```

**Escape Route**: SMFT transition at z < 1100 → galaxy-scale tests instead

#### (B) String Tension Below Detection
```
Gμ/c² < 10⁻⁸ at 5σ (CMB-S4)

Current: Gμ < 8.6×10⁻⁷ (Planck 2018, 3σ)
Future: Gμ ~ 10⁻⁸ sensitivity (CMB-S4)

If CMB-S4 measures: Gμ < 10⁻⁸
→ NO cosmic string network
→ SMFT winding number NOT conserved, or
→ Transition too weak (R² ≪ 1 at phase transition)
```

**Escape Route**: Rapid transition (short correlation length) → fewer strings

#### (C) Perfect Scale Invariance
```
dn_s/dlnk = 0.000 ± 0.001 (CMB-S4)

Current: dn_s/dlnk = -0.0045 ± 0.0067 (Planck)
Future: σ ~ 0.001 (CMB-S4)

If CMB-S4 measures: dn_s/dlnk = 0.000 ± 0.001
→ NO string contribution to power spectrum
→ Perfectly scale-invariant (pure inflation)
→ SMFT string network absent or negligible
```

**Escape Route**: String contribution < 1% of total power

### Alternative Explanations

Before declaring SMFT cosmology falsified:

1. **Foreground contamination**: Verify CMB cleaning (dust, synchrotron)
2. **Systematics**: Check beam, calibration, polarization systematics
3. **Theoretical uncertainty**: Recompute string signals with updated simulations

**Only if all checks pass AND (A) or (B) or (C) at 5σ → SMFT cosmology falsified**

---

## Area 4: Quantum Synchronization Critical Exponent

### SMFT Prediction
```
β = 0.099 ± 0.004 (novel universality class)

vs.
β_2D-Ising = 0.125
β_2D-XY = 0.23
```

### Falsification Criterion

**SMFT UNIVERSALITY is FALSIFIED if**:
```
β_measured = 0.125 ± 0.010 (2D Ising confirmed at 5σ)
```

**Statistical Test**:
```
H₀: β = 0.099 (SMFT)
H₁: β = 0.125 (2D Ising)

Separation: |0.125 - 0.099| = 0.026

If σ(β_measured) < 0.005:
→ Can distinguish at 0.026 / 0.005 = 5.2σ

If measurement gives:
β_measured = 0.125 ± 0.005
→ SMFT ruled out at (0.125 - 0.099) / 0.005 = 5.2σ
```

**Confidence**: 5σ with N > 100 qubits

### Finite-Size Scaling Test

**Additional Falsification**: Data collapse failure

```
SMFT predicts: R(σ, L) = L^(-β/ν) f[(σ - σ_c)L^(1/ν)]

with β/ν = 0.040, 1/ν = 0.40

Test: Plot L^(0.040) R vs. (σ - σ_c) L^(0.40)

If curves do NOT collapse:
→ χ²_collapse > threshold
→ SMFT exponents wrong
→ Try 2D Ising exponents (β/ν = 0.125, 1/ν = 1.0)

If 2D Ising exponents give good collapse:
→ SMFT universality FALSIFIED
```

### Alternative Explanations

1. **Crossover behavior**: β_eff(L) → 0.125 as L → ∞
   - Test: Measure β vs. system size
   - If β(L=50) = 0.099, β(L=200) = 0.125 → SMFT is finite-size effect

2. **First-order transition**: Not continuous critical point
   - Test: Look for discontinuity in R(σ)
   - If R jumps at σ_c → Different mechanism

3. **Quantum vs. classical**: β_quantum ≠ β_classical
   - Control: Simulate classical Kuramoto model (Monte Carlo)
   - If classical also gives β = 0.099 → Not quantum-specific

**Only if crossover, first-order, and quantum effects ruled out AND β = 0.125 → SMFT falsified**

---

## Area 5: High-Energy Collider Tests

### SMFT Prediction
```
Fermion mass: m_f = Δ · R(x,t)

If R varies: m_f(environment) ≠ const
```

### Falsification Criteria

**SMFT DYNAMICAL MASS is FALSIFIED if ALL THREE**:

#### (A) Lepton Universality Perfect
```
R_e / R_μ = 1.0000 ± 0.0001 (FCC-ee)

Current: R_e/R_μ = 1.0001 ± 0.0014
Future: σ ~ 0.0001 (FCC-ee, 10¹² Z bosons)

If FCC-ee measures: R_e/R_μ = 1.0000 ± 0.0001
→ Lepton masses IDENTICAL up to 0.01%
→ NO R-field modulation
→ SMFT dynamical mass FALSIFIED
```

#### (B) Top Mass Environment-Independent
```
|m_t(hadronic) - m_t(leptonic)| < 10 MeV (FCC-ee threshold)

Current: |Δm_t| < 1 GeV (LHC)
Future: δm_t ~ 10 MeV (FCC-ee)

If FCC-ee measures: m_t(all channels) = 172.5 ± 0.01 GeV
→ Top mass STABLE across environments
→ NO dynamical dependence
→ SMFT FALSIFIED
```

#### (C) Higgs Couplings Standard Model
```
μ_f = σ(H → ff̄) / σ_SM = 1.000 ± 0.005 for all f

Current: μ_f = 1.0 ± 0.1 (LHC)
Future: σ(μ_f) ~ 0.005 (HL-LHC, FCC-ee)

If all measured: μ_b = μ_τ = μ_t = 1.000 ± 0.005
→ Higgs couplings EXACTLY Standard Model
→ NO R² modification
→ SMFT FALSIFIED
```

**All three at 5σ → Fermion mass is NOT dynamical → SMFT wrong**

### SMFT Escape Routes

If falsified at colliders, SMFT can still survive by:

1. **R = 1 at high energies**: Phase transition complete at T ≪ T_EW
   - Predict: Low-energy tests (atomic physics, condensed matter)

2. **Composite fermions only**: SMFT applies to hadrons, not leptons
   - Predict: Proton structure modifications, quark mass variations

3. **Gravitational sector**: R couples to g_μν, not fermion mass
   - Predict: Equivalence principle violations, gravitational tests

**If all three escape routes fail → SMFT fundamentally wrong**

---

## Meta-Falsification: Theory-Wide Rejection

### SMFT as a Whole is FALSIFIED if:

**Condition**: ANY TWO of the five areas show 5σ rejection

```
Example 1: BEC phonon + Casimir force
  - BEC: c_eff(r) = c_s (no R-field modulation)
  - Casimir: F/F_Lifshitz = 1.00 (no R² coupling)
  → R-field does NOT couple to observables
  → CORE SMFT MECHANISM WRONG

Example 2: Quantum sync + CMB cosmology
  - Quantum: β = 0.125 (standard 2D Ising)
  - CMB: f_NL = 0, Gμ = 0 (no defects)
  → SMFT critical behavior is conventional
  → Phase transition exists but NOT novel

Example 3: BEC + Collider
  - BEC: c_eff(r) = c_s (no coupling)
  - Collider: m_f = const (no dynamics)
  → R-field is mathematical artifact, not physical
  → SMFT FALSIFIED
```

**Two independent falsifications → overwhelming evidence against theory**

---

## Statistical Confidence Requirements

### Standard Thresholds

| Confidence | σ (std. dev.) | p-value | Physics Usage |
|------------|---------------|---------|---------------|
| 90% CL | 1.6σ | 0.10 | Preliminary |
| 95% CL | 2.0σ | 0.05 | Suggestive |
| 99% CL | 2.6σ | 0.01 | Evidence |
| **3σ** | **3.0σ** | **0.003** | **Strong evidence** |
| 99.99% CL | 4.0σ | 10⁻⁴ | Very strong |
| **5σ** | **5.0σ** | **3×10⁻⁷** | **Discovery standard** |

**SMFT Falsification Requires**:
- **Minimum**: 3σ (99.7% confidence)
- **Standard**: 5σ (discovery threshold in particle physics)

### Why 5σ?

**Historical Precedent**:
- Higgs boson: Discovered at 5σ (July 2012)
- Neutrino oscillations: Confirmed at >5σ
- Gravitational waves: First detection at 5.1σ (GW150914)

**Reason**: False positive rate ~ 1 in 3.5 million (accounting for look-elsewhere effect)

**SMFT Standard**:
- Use 5σ for positive claims (novel physics)
- Use 5σ for falsification (excluding SMFT)
- Ensures robust conclusions resistant to statistical fluctuations

---

## Systematic Error Treatment

### Principle: Systematics ≤ Statistical Error

**Good Practice**:
```
Total error: σ_total² = σ_stat² + σ_sys²

For 5σ claim: σ_total < |signal| / 5

If σ_sys > σ_stat:
→ Experiment is systematics-limited
→ Must improve systematic control OR increase statistics
```

**Example (BEC Phonon)**:
```
Signal: c_eff(ξ)/c_s = 0.76 (24% suppression)
Target: 5σ → σ_total < 24% / 5 = 4.8%

Statistical (20 runs): σ_stat ~ 2% (√20 averaging)
Systematic budget: σ_sys ~ 5% (imaging, thermal, ...)

Total: σ_total = √(2² + 5²) = 5.4%

Signal / σ_total = 24% / 5.4% = 4.4σ

Sufficient? MARGINAL (aim for σ_sys < σ_stat)
```

### Systematic Error Mitigation

**Strategy**:
1. **Identify dominant systematics** (error budget)
2. **Control experiments** (vary one parameter at a time)
3. **Redundant measurements** (cross-check with independent methods)
4. **Blind analysis** (avoid confirmation bias)

**Example (Casimir)**:
- Electrostatic: Vary applied voltage → nulling
- Roughness: Measure multiple samples → average
- Temperature: Scan T → verify T-dependence
- Dielectric: Independent ε(ω) measurement → constrain

**Only after systematics controlled → interpret as SMFT falsification**

---

## Decision Tree: Interpreting Experimental Results

```
┌─ Experiment Performed ─┐
│                        │
├─ Is signal present?    │
│  (deviation from SM)   │
│                        │
├─ NO ──────────────────┐│
│                       ││
├─ YES ─────────────────┤│
│                       ││
│  ┌──────────────────┐││
│  │ Systematics OK?  │││
│  │ (all controls)   │││
│  └──────────────────┘││
│          │            ││
│     ┌────┴────┐      ││
│     │         │      ││
│    YES       NO      ││
│     │         │      ││
│  ┌──┴───┐ ┌──┴───┐  ││
│  │ 5σ?  │ │ Redo │  ││
│  └──────┘ └──────┘  ││
│     │                ││
│  ┌──┴───┐            ││
│  │ YES  │            ││
│  └──┬───┘            ││
│     │                ││
│  ┌──┴──────────────┐ ││
│  │ Matches SMFT?   │ ││
│  │ (quantitative)  │ ││
│  └──┬─────────┬────┘ ││
│     │         │      ││
│    YES       NO      ││
│     │         │      ││
│  ┌──┴───┐ ┌──┴────┐ ││
│  │SMFT  │ │ SMFT  │ ││
│  │PASS  │ │ FAIL  │ ││
│  └──────┘ └───────┘ ││
│                      ││
└──────────────────────┴┘
       │
       └─► NO signal at 5σ
           │
           ├─ Systematics OK? ─┐
           │                   │
           ├─ YES              NO
           │   │               │
           │   └─► SMFT        └─► Redo
           │       FALSIFIED
           │
           └─► Decision: Reject SMFT
```

---

## Publication Ethics: Pre-Registration

### Open Science Commitment

**Before running experiments, publicly pre-register**:

1. **Prediction**: What does SMFT predict? (numerical values)
2. **Null hypothesis**: What is the standard expectation?
3. **Falsification threshold**: What result would rule out SMFT? (5σ criterion)
4. **Analysis plan**: How will data be analyzed? (blind if possible)

**Prevents**:
- Post-hoc rationalization ("epicycles")
- Cherry-picking favorable results
- Moving goalposts after data seen

**Example (BEC Experiment)**:
```
Pre-registration (before measurement):
"If c_eff(r)/c_s = 1.00 ± 0.05 for all r > ξ at 3σ,
 SMFT prediction c_eff ~ tanh(r/ξ) is falsified."

After measurement:
- If c_eff(ξ) = 0.75 ± 0.05 → SMFT CONFIRMED
- If c_eff(ξ) = 1.00 ± 0.05 → SMFT FALSIFIED
- NO post-hoc adjustment of threshold
```

**Platforms**:
- OSF (Open Science Framework): https://osf.io/
- arXiv preprint: Analysis plan before data
- Public GitHub: Code + protocol

---

## Conclusion: SMFT is Falsifiable

**This document provides explicit, quantitative criteria for falsifying SMFT:**

1. **BEC phonon**: c_eff(r) = c_s at 3σ → FALSIFIED
2. **Casimir force**: F/F_Lifshitz = 1.00 at 5σ → FALSIFIED
3. **CMB cosmology**: f_NL = 0, Gμ < 10⁻⁸, or dn_s = 0 at 5σ → FALSIFIED
4. **Quantum sync**: β = 0.125 at 5σ → FALSIFIED
5. **Collider**: Lepton universality + m_t stable + Higgs SM at 5σ → FALSIFIED

**Any TWO falsifications → SMFT as a whole rejected**

**Commitment**:
- If experiments meet these criteria, **we will declare SMFT wrong**
- No post-hoc modifications to save the theory
- Accept failure gracefully, revise or abandon SMFT

**This is what makes SMFT science, not pseudoscience.**

**Next step: Execute experiments. Let nature decide.**
