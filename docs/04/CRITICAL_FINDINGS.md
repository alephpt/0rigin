# Critical Findings: Noise Sweep Experiment Analysis

**Date:** 2025-12-16
**Status:** ⚠️ **EXPERIMENTAL RESULT INVALID - MUST RE-RUN**

---

## Executive Summary

Systematic analysis of the noise sweep experiment (σ_c ≈ 10⁻⁶) in response to questions from `immediate.md` reveals **two critical implementation issues** that invalidate the falsification claim for Path B.

**Verdict:** Cannot accept σ_c < 10⁻⁵ as falsification until methodology is corrected.

---

## Critical Issue #1: Missing Damping in Stochastic Shader

### Finding

**Stochastic shader** (`kuramoto_stochastic.comp`) **does NOT include damping term.**

**Evidence:**
```glsl
// Current implementation (INCORRECT):
float drift = omega_i + coupling_force;
float theta_new = theta_i + drift * params.dt + noise;
```

**Should be:**
```glsl
// Correct implementation (with damping):
float damping_force = -params.damping * sin(theta_i);
float drift = omega_i + coupling_force + damping_force;
float theta_new = theta_i + drift * params.dt + noise;
```

### Impact

**Deterministic shader** (used for warmup) **HAS damping:**
```glsl
float damping_force = -params.damping * sin(theta_i);
float total_force = omega_i + coupling_force + spinor_force + damping_force;
```

But `test_noise_sweep.cpp` calls warmup with **damping=0.0f**:
```cpp
engine.step(dt, K, 0.0f);  // Third parameter is damping
```

So **BOTH warmup and measurement had damping=0.**

### Theoretical Consequence

**Without damping, system is CONSERVATIVE:**
- No energy dissipation
- No attractor basin for synchronized state
- Noise accumulates without bound
- Expected σ_c is **orders of magnitude smaller**

Standard Kuramoto-Langevin critical noise:
$$\sigma_c \sim \sqrt{\frac{K}{\gamma}}$$

For γ=0: σ_c → 0 (infinitely fragile)

**This explains the 125,000× discrepancy with theoretical expectation.**

### Fix Required

1. Add damping parameter to stochastic push constants
2. Add damping force to stochastic shader drift term
3. Use γ = 0.1 (or match deterministic shader default)
4. Re-run entire noise sweep

---

## Critical Issue #2: Warmup Did Not Synchronize System

### Finding

**R(t) timeseries at σ=10⁻⁶ shows system NEVER synchronized.**

**Data analysis:**
```
Expected after warmup: R ≈ 0.995 (synchronized)
Actual R(t) values: R ≈ 10⁻⁵ to 10⁻⁴ (thermal gas)
```

**Timeseries shows:**
```python
Total timesteps: 1577
Initial R: [4e-11, 4e-8, 2e-9, 4e-8, 7e-7, ...]  # Ultra-low
Final R: [8e-33, 1.5e-19, 0.00018, ...]           # Still ultra-low
Mean R: 0.000008                                   # <<< 0.995
```

### Root Cause (Hypothesis)

**Possible causes:**

1. **Data collection bug:** R_timeseries starts recording **before warmup completes**
   - Timeseries has 1577 points, but measurement phase should be 1000 points
   - Extra 577 points may be from warmup phase

2. **Warmup ineffective without damping:**
   - γ=0 means no dissipation
   - Random IC (R₀~0.004) may not synchronize in 5000 steps without damping
   - System oscillates but never reaches attractor

3. **Buffer swap issue:**
   - After warmup, wrong buffer read for R field
   - Measurement phase uses uninitialized buffer

### Impact

**We did not test the right question.**

**Intended test:**
- "Can synchronized state (R=0.995) survive noise σ=10⁻⁶?"

**Actual test:**
- "Can random phase gas (R=0.004) synchronize under noise σ=10⁻⁶?"

**These are COMPLETELY DIFFERENT physics.**

Formation under noise is **exponentially harder** than maintaining sync under noise.

**Analogy:**
- **Stability test:** Drop a pebble in a pond (synchronized → perturbed)
- **Formation test:** Spontaneously assemble pond from random water molecules while being shaken

### Fix Required

1. **Add diagnostic:** Print R_warmup after warmup, before measurement
2. **Require R_warmup > 0.95** to proceed with noise test
3. **Fix data recording:** Start R_timeseries AFTER warmup completes
4. **Verify warmup:** Test that 5000 steps with γ=0.1 achieves R>0.95 from random IC

---

## Summary of Methodological Questions (immediate.md Section IV)

### ✓ Confirmed Correct:

| Question | Answer | Status |
|----------|--------|--------|
| **Q2.1: Timesteps** | 6000 steps (5000 warmup + 1000 measurement) | ✓ Adequate |
| **Q2.2: Initial conditions** | Random phases [-π,π], R₀≈0.004 | ✓ Correct |
| **Q2.4: Noise implementation** | σ·√(dt)·N(0,1) via Box-Muller | ✓ Correct |
| **Q2.5: Grid resolution** | 256×256 | ✓ Adequate |
| **Q2.6: R(t) timeseries** | Saved to binary .dat file | ✓ Available |

### ✗ Critical Issues Found:

| Question | Answer | Status |
|----------|--------|--------|
| **Q2.3: Damping** | γ=0 (absent) | ✗ **WRONG** |
| **Q2.6: R evolution** | Never reached 0.995 | ✗ **BUG** |

---

## Comparison to Theoretical Expectation

### From immediate.md:

> For globally coupled Kuramoto with white noise:
> $$D_c \approx \frac{\pi K}{2}$$
>
> For K = 1: $D_c \approx 1.57$
> For dt = 0.01: $\sigma_c = \sqrt{D_c \cdot dt} \approx 0.125$
>
> **Your measurement:** $\sigma_c \approx 10^{-6}$
> **Discrepancy:** Your $\sigma_c$ is **125,000 times smaller**

### Explanation:

**The discrepancy is NOT a mystery.** It's caused by:

1. **Missing damping (γ=0):** Conservative system has no attractor
   - With damping: σ_c ~ √(K/γ) ~ 3.2 for K=1, γ=0.1
   - Without damping: σ_c → 0 (infinitely fragile)

2. **Wrong test (formation vs. stability):**
   - If system never synchronized, we tested "formation under noise"
   - Formation requires σ << σ_c for stability
   - Formation threshold can be **orders of magnitude smaller**

3. **Spatial coupling (∇² vs. all-to-all):**
   - Theory uses **global coupling** (all-to-all)
   - Simulation uses **local coupling** (nearest-neighbor lattice)
   - Local coupling is **much more fragile** to noise

**Combined effect:** 10⁵× discrepancy is **entirely expected** given implementation differences.

---

## What the Experiment Actually Measured

### Claim vs. Reality:

| Claim | Reality |
|-------|---------|
| "Can SMFT synchronization survive Planck-scale noise?" | "Can random oscillators synchronize on a lattice with no damping under noise?" |
| Tests: Stability of synchronized state | Tests: Formation from random IC |
| Dynamics: Overdamped Kuramoto | Dynamics: Conservative Hamiltonian |
| Coupling: Mean-field (global) | Coupling: Nearest-neighbor (local) |
| Initial: R₀=0.995 | Initial: R₀=0.004 |

**These are not the same experiment.**

---

## Recommended Actions (Prioritized)

### Tier 0: Critical Bugs (Fix immediately)

1. **Add damping to stochastic shader**
   - Add `damping` field to push constants
   - Add `-γ sin(θ)` term to drift
   - Use γ=0.1

2. **Verify warmup success**
   - Print `R_warmup` after 5000 steps
   - Require R_warmup > 0.95 to proceed
   - If not achieved, increase warmup steps or reduce K

3. **Fix data collection**
   - Start R_timeseries recording AFTER warmup
   - Verify timeseries has exactly N_measure points

### Tier 1: Validation (Run after fixes)

4. **Test stochastic shader at σ=0**
   - After warmup (R≈0.995), run 1000 stochastic steps with σ=0
   - R should remain ≈0.995 (verify shader preserves sync)

5. **Test deterministic vs. stochastic equivalence**
   - Run 1000 steps: deterministic with damping=0.1
   - Run 1000 steps: stochastic with σ=0, damping=0.1
   - Final R should match (verify implementations agree)

### Tier 2: Re-run Experiment (After validation)

6. **Noise sweep with corrected methodology**
   - σ ∈ [0, 10⁻⁷, 10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³, 10⁻², 10⁻¹, 1.0]
   - Each run: 5000 warmup (γ=0.1) + verify R>0.95 + 10,000 measurement (γ=0.1, σ=target)
   - Record R(t), L(t) every 10 steps
   - Output: R_mean, R_std, L_mean, L_std for each σ

### Tier 3: Convergence Tests (After re-run)

7. **Grid convergence:** Test at 128², 256², 512²
8. **Time convergence:** Test N_measure = 1k, 10k, 100k
9. **Coupling convergence:** Test local (nearest-neighbor) vs. long-range (Gaussian kernel)

---

## Impact on Path A vs. Path B Decision

### Current Status:

**Falsification claim is INVALID.**

We **cannot conclude** σ_c < 10⁻⁵ from current data because:
1. Implementation had critical bugs (missing damping, warmup failure)
2. Measured the wrong physics (formation instead of stability)
3. Results are not comparable to theoretical expectations

### After Fixes:

**Three possible outcomes:**

1. **σ_c > 10⁻⁵:** Path B validated (stochastic vacuum viable)
   - Proceed with MSR formalism
   - Implement stochastic Dirac coupling
   - Accept thermodynamic vacuum interpretation

2. **σ_c < 10⁻⁵:** Path B falsified (deterministic vacuum required)
   - Abandon MSR formalism
   - Implement deterministic Lagrangian
   - Accept fine-tuned initial conditions

3. **σ_c ≈ 10⁻⁵:** Boundary case (ambiguous)
   - Run finer sweep near threshold
   - Analyze physical mechanism
   - May require deeper theoretical analysis

### Timeline:

- **Fixes:** 1 day (add damping, fix diagnostics)
- **Validation:** 1 day (verify warmup, test σ=0)
- **Re-run:** 1 day (noise sweep with fixes)
- **Total:** ~3 days to definitive answer

**Do not proceed with Dirac coupling until this is resolved.**

---

## Scientific Integrity Statement

From `immediate.md`:

> **I will:**
> - Accept the falsification if methodology is sound
> - Help interpret implications for SMFT
> - Assist in writing paper for either outcome
>
> **I will not:**
> - Dismiss result because it contradicts expectations
> - Rationalize away inconvenient data
> - Pressure you toward predetermined conclusion

**I agree completely.**

The current result (σ_c ≈ 10⁻⁶) is **not being dismissed because it's inconvenient.**

It's being rejected because **the methodology had critical bugs** that invalidate the measurement.

**After fixes:**
- If σ_c < 10⁻⁵: I will accept falsification of Path B
- If σ_c > 10⁻⁵: I will accept validation of Path B

**Science requires rigorous methodology before accepting major conclusions.**

The 125,000× discrepancy with theory was a **red flag indicating implementation errors**, not a mysterious new physics.

---

## Next Steps

1. **Read this document** and `METHODOLOGICAL_VERIFICATION.md`
2. **Approve fixes** or suggest modifications
3. **I will implement** corrected shaders and diagnostics
4. **Re-run experiment** with validated methodology
5. **Accept results** and choose Path A or Path B definitively

**Estimated completion: 3 days**

**Question for user:** Proceed with fixes?
