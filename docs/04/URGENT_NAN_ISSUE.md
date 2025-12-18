# URGENT: NaN Values in Corrected Experiment

**Status:** 🔴 **CRITICAL - GPU COMPUTE RETURNING NaN**

---

## Problem

After implementing all corrections (R_global computation, damping), the experiment runs to completion but **all phase values are NaN**.

```
R_global after warmup: -nan
R_global = -nan ± -nan (Kuramoto order parameter)
R_local  = 0 ± 0 (Spatial avg of local sync)
L = 0 ± 0
```

**This is NOT the same issue as before** - this is a GPU compute or memory readback failure.

---

## Evidence

### 1. Timeseries shows NaN from start
```
$ head -5 output/noise_sweep/sigma_0.00e+00/R_global_timeseries.dat
0 -nan
1 -nan
2 nan
3 -nan
4 nan
```

### 2. Happens even at σ=0 (deterministic warmup)
- Warmup uses `engine.step(dt, K, 0.0f)` (deterministic shader)
- Should work perfectly - used in all previous tests
- But R_global is NaN even after warm up

### 3. R_local is also 0
- Both R_global (from phases) and R_local (from sync_field shader) broken
- Suggests phase field itself contains NaN

---

## Hypotheses

### H1: GPU readback failure
- Phases computed correctly on GPU
- Readback `getPhaseField()` returns garbage/NaN
- **Test:** Print first 10 phase values after warmup

### H2: Initial conditions not set
- `setInitialPhases()` may not actually upload to GPU
- Phases start as 0 or NaN
- Warmup operates on garbage
- **Test:** Verify upload actually happens

### H3: Deterministic shader broken by changes
- We didn't modify deterministic shader
- But push constants structure might have changed
- Mismatch between CPU and GPU structures
- **Test:** Run simple deterministic-only test

### H4: Damping value causes instability
- γ=0.1 too strong for K=27.21?
- Causes exponential divergence?
- **Test:** Try γ=0 or smaller γ

---

## Immediate Debugging Steps

### Step 1: Verify initial conditions are set
```cpp
// In test, after setInitialPhases():
std::vector<float> theta_check = engine.getPhaseField();
std::cout << "First 10 phases after IC: ";
for (int i = 0; i < 10; i++) {
    std::cout << theta_check[i] << " ";
}
std::cout << std::endl;
```

### Step 2: Print phase values during warmup
```cpp
// After each 1000 warmup steps:
std::vector<float> theta = engine.getPhaseField();
float first_phase = theta[0];
float mean_phase = 0;
for (float t : theta) mean_phase += t;
mean_phase /= theta.size();
std::cout << "  Step " << step << ": first_phase=" << first_phase
          << ", mean_phase=" << mean_phase << std::endl;
```

### Step 3: Check for inf/nan in shaders
```glsl
// In kuramoto_step.comp and kuramoto_stochastic.comp:
if (isnan(theta_new) || isinf(theta_new)) {
    theta_new = 0.0;  // Clamp to prevent propagation
}
```

---

## What Changed That Could Break This?

### Changes we made:
1. ✓ Added `compute_global_R()` - pure CPU, can't break GPU
2. ✓ Added damping to stochastic shader - **could cause instability**
3. ✓ Changed push constants structure - **could cause mismatch**
4. ✓ Changed measurement loop - pure CPU, can't break GPU

### Most likely culprits:

**1. Push constants structure mismatch:**
```cpp
// Stochastic (NEW):
struct StochasticPushConstants {
    float dt, K, sigma, damping, omega_mean;
    uint32_t Nx, Ny, time_step;
};  // 5 floats + 3 uints = 32 bytes

// Deterministic (OLD):
struct DeterministicPushConstants {
    float dt, K, damping, Delta, chiral_angle;
    float Nx, Ny, N_total, neighborhood_radius;
};  // 9 floats = 36 bytes
```

If deterministic shader expects different layout, **reads garbage for Nx, Ny**.

**2. Damping causes runaway:**
- With γ=0.1 and K=27.21, damping might be too strong
- Phases spiral to ±∞, wrap creates NaN

---

## Quick Test: Run Without Modifications

**Revert to known-working version:**
1. Don't modify push constants
2. Don't add damping to stochastic
3. Just add `compute_global_R()` to test
4. See if deterministic warmup works

**This isolates whether bug is in our changes or pre-existing.**

---

##Status: BLOCKED

Cannot proceed with noise sweep until NaN issue resolved.

**Two paths forward:**

**Path A (Quick): Revert stochastic changes, just fix R_global measurement**
- Remove damping from stochastic shader
- Remove damping from push constants
- Keep `compute_global_R()` fix
- Run sweep - will still have "no damping" issue but at least measures right thing

**Path B (Thorough): Debug NaN issue properly**
- Add diagnostics to find root cause
- Fix whatever is breaking GPU compute
- Then re-run with all corrections

**Recommendation: Path A first** (get SOME result), then Path B.

We need to know if the R_global vs R_local fix alone changes the answer, before adding more complexity.

---

## User Decision Required

**Question:** Should I:

1. **Revert damping changes**, keep only R_global fix, re-run?
   - Pro: Quick, isolates one variable
   - Con: Still has "no damping" issue from before

2. **Debug NaN issue thoroughly** before any re-run?
   - Pro: Get fully corrected result
   - Con: Takes longer, might be deep bug

3. **Something else?**

**I'm blocked pending your direction.**
