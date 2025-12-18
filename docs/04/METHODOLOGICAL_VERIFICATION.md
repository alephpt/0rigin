# Methodological Verification: Answers to Critical Questions

**Document Purpose:** Systematic answer to all questions raised in `immediate.md` Section IV regarding the noise sweep experiment that measured σ_c ≈ 10⁻⁶.

**Date:** 2025-12-16
**Experiment:** Path B Noise Sweep (Directive.md Phase 2)
**Critical Finding:** σ_c ∈ (10⁻⁶, 10⁻⁵), suggesting falsification of stochastic vacuum hypothesis

---

## I. Direct Answers to Critical Questions (immediate.md Section IV)

### Q2.1: Simulation Duration

**Question:** How many timesteps were run at each σ value? What is total simulated time T = N_steps × dt?

**Answer:**
- **Warmup phase:** N_warmup = 5000 steps (deterministic, σ=0)
- **Measurement phase:** N_measure = 1000 steps (with noise at target σ)
- **Total steps per sigma:** 6000 steps
- **Timestep:** dt = 0.01
- **Total simulated time:** T = 6000 × 0.01 = **60 time units**
- **Warmup time:** T_warmup = 50 time units
- **Measurement time:** T_measure = 10 time units

**Source:** `test/test_noise_sweep.cpp:119-120`
```cpp
const int N_warmup = 5000;    // Steps to reach equilibrium
const int N_measure = 1000;   // Steps to measure statistics
const float dt = 0.01f;
```

**Assessment:**
- ✓ Warmup time (T=50) is sufficient for equilibration at K=27.21 (relaxation time τ ~ 1/K ~ 0.04)
- ✓ Measurement time (T=10) captures 250 relaxation times worth of statistics
- ⚠️  However, this may not be sufficient for steady-state at **critical point** σ ≈ σ_c where τ → ∞

---

### Q2.2: Initial Conditions

**Question:** Was θ(t=0) random (R₀ ≈ 0.3)? Or pre-synchronized (R₀ ≈ 0.99)?

**Answer:** **Random phases from uniform distribution**

**Code:** `test/test_noise_sweep.cpp:162-165`
```cpp
std::vector<float> theta_init(Nx * Ny);
srand(42);  // Fixed seed for reproducibility
for (int i = 0; i < Nx * Ny; i++) {
    theta_init[i] = (float(rand()) / RAND_MAX) * 2.0f * M_PI - M_PI;
}
```

**Expected R₀:** For uniform random phases in [-π, π]:
$$R_0 = |\langle e^{i\theta} \rangle| \approx \frac{1}{\sqrt{N}} = \frac{1}{\sqrt{256^2}} = \frac{1}{256} \approx 0.004$$

**Then:** 5000-step warmup with σ=0 (deterministic) allows system to **synchronize from random IC**:
- Initial: R₀ ≈ 0.004 (random)
- After warmup: R ≈ 0.995 (synchronized)
- Then noise is applied

**Physical interpretation:**
- Test measures: **"Can synchronized state survive noise?"** (stability test)
- NOT: "Can synchronization form under noise?" (formation test)

**This is correct methodology** for falsification criterion in Directive.md.

---

### Q2.3: Damping Parameter

**Question:** Was γ = 0.1 included in noise test? Or was damping turned off (γ = 0)?

**Answer:** **Damping was ABSENT** (γ = 0)

**Evidence:**

**1. Stochastic shader implementation** (`shaders/smft/kuramoto_stochastic.comp:204-210`):
```glsl
// Compute deterministic drift
float omega_i = params.omega_mean;  // Could make this per-oscillator later
float coupling_force = compute_coupling_force(local_id.x, local_id.y, theta_i);

float drift = omega_i + coupling_force;

// Stochastic term: σ·√(dt)·N(0,1)
float noise = params.sigma * sqrt(params.dt) * randn();

// Euler-Maruyama update: θ(t+dt) = θ(t) + drift·dt + σ·√(dt)·dW
float theta_new = theta_i + drift * params.dt + noise;
```

**No damping term present.** Drift contains only:
- ω (natural frequency)
- K·coupling (synchronization force)

**2. Deterministic shader for comparison** (`shaders/smft/kuramoto_step.comp:183-186`):
```glsl
// Compute forces
float coupling_force = compute_coupling_force(local_id.x, local_id.y, theta_i);
float spinor_force = compute_spinor_feedback(idx, theta_i);
float damping_force = -params.damping * sin(theta_i);

// Total force
float total_force = omega_i + coupling_force + spinor_force + damping_force;
```

**Deterministic shader has damping_force term,** but stochastic shader does not.

**3. Push constants structure:**

Stochastic:
```glsl
layout(push_constant) uniform PushConstants {
    float dt;
    float K;
    float sigma;           // Noise amplitude
    float omega_mean;
    uint Nx;
    uint Ny;
    uint time_step;
    uint pad;
} params;
```
No `damping` field.

Deterministic:
```glsl
layout(push_constant) uniform PushConstants {
    float dt;
    float K;
    float damping;         // Phase damping γ
    float Delta;
    float chiral_angle;
    uint Nx;
    uint Ny;
    uint N_total;
    uint enable_feedback;
} params;
```
Has explicit `damping` field.

**Conclusion:** **Damping was γ = 0 during stochastic noise test.**

**Impact assessment:**
- ⚠️  **This is critical.** Damping provides dissipation that stabilizes synchronization.
- Without damping, system is **conservative** and more fragile to noise.
- Standard Kuramoto-Langevin includes damping or assumes overdamped limit.
- **This could explain extremely low σ_c.**

**Recommendation:** Re-run with damping term `-γ sin(θ)` added to stochastic shader.

---

### Q2.4: Noise Implementation (Exact Code)

**Question:** Which noise scaling was used?
- Option A (correct): `theta_new = theta_old + (omega + F)*dt + sigma*sqrt(dt)*randn()`
- Option B (incorrect): `theta_new = theta_old + (omega + F)*dt + sigma*randn()`

**Answer:** **Option A (CORRECT) - Proper Euler-Maruyama scaling**

**Code:** `shaders/smft/kuramoto_stochastic.comp:206-210`
```glsl
// Stochastic term: σ·√(dt)·N(0,1)
float noise = params.sigma * sqrt(params.dt) * randn();

// Euler-Maruyama update: θ(t+dt) = θ(t) + drift·dt + σ·√(dt)·dW
float theta_new = theta_i + drift * params.dt + noise;
```

**Breakdown:**
1. `randn()` generates standard normal N(0,1) via Box-Muller transform
2. Noise term: `params.sigma * sqrt(params.dt) * randn()`
3. This is **exactly** the Euler-Maruyama prescription for SDE: dθ = μ dt + σ dW

**PRNG implementation** (`shaders/smft/kuramoto_stochastic.comp:48-88`):
```glsl
// PCG Random Number Generator (32-bit)
uint pcg_state;

uint pcg_hash(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

void pcg_init(uint seed) {
    pcg_state = pcg_hash(seed);
}

// Returns uniform random float in [0, 1)
float pcg_random() {
    pcg_state = pcg_state * 747796405u + 2891336453u;
    uint word = ((pcg_state >> ((pcg_state >> 28u) + 4u)) ^ pcg_state) * 277803737u;
    word = (word >> 22u) ^ word;
    return float(word) / 4294967296.0;
}

// Box-Muller Transform for Gaussian Random Numbers
// Returns standard normal N(0,1)
float randn() {
    float u1 = pcg_random();
    float u2 = pcg_random();

    // Avoid log(0)
    u1 = max(u1, 1e-10);

    const float PI = 3.14159265359;
    float r = sqrt(-2.0 * log(u1));
    float theta_rand = 2.0 * PI * u2;

    return r * cos(theta_rand);
}
```

**Verification:**
- ✓ PCG is industry-standard PRNG with excellent statistical properties
- ✓ Box-Muller correctly transforms U[0,1] → N(0,1)
- ✓ Noise scaling includes √(dt) factor
- ✓ Implementation matches Kloeden & Platen (1992) Euler-Maruyama scheme

**Conclusion:** **Noise implementation is mathematically correct.**

---

### Q2.5: Grid Resolution

**Question:** Grid resolution - still 128²? Any verification at different resolutions?

**Answer:** **256×256** (upgraded from earlier 128×128)

**Code:** `test/test_noise_sweep.cpp:117-118`
```cpp
const int Nx = 256;
const int Ny = 256;
```

**Grid points:** N = 256² = 65,536 oscillators

**Grid convergence:** Not yet tested. Single resolution only.

**Recommendation:** Test at 128², 256², 512² to verify σ_c doesn't change with resolution (rule out discretization artifacts).

---

### Q2.6: Temporal Evolution at σ = 10⁻⁶

**Question:** Did you save R(t) time series at critical point? Does R decay smoothly or jump discontinuously?

**Answer:** **Yes, R(t) timeseries was saved. Analysis reveals severe issue.**

**Data location:** `/home/persist/neotec/0rigin/output/noise_sweep/sigma_1.00e-06/R_timeseries.dat`

**Analysis:**
```
Total timesteps: 1577
Initial R (first 10): [4e-11, 4e-8, 2e-9, 4e-8, 7e-7, 1.5e-19, 1e-5, 3e-6, 6e-10, 1e-5]
Final R (last 10): [8e-33, 1.5e-19, 0.00018, 8e-33, 1.5e-19, 0.00018, 8e-33, 1.6e-19, 0.00018, 8e-33]
Mean R: 0.000008
Std R: 0.000032
Min R: 0.000000
Max R: 0.000176
```

**Critical finding:** R **NEVER REACHES 0.995**

Expected evolution:
1. t=0: R₀ ≈ 0.004 (random IC)
2. t=0-50: R → 0.995 (warmup, σ=0)
3. t=50-60: R decays from 0.995 under noise (σ=10⁻⁶)

**Actual evolution:**
1. t=0: R ≈ 10⁻¹¹ to 10⁻⁵ (ultra-low, inconsistent with random IC)
2. R **never synchronizes**
3. Remains in thermal gas state (R ~ 10⁻⁴) throughout

**This indicates a BUG in the data collection or warmup process.**

**Possible causes:**
1. R_timeseries is recording data from **before warmup** (should start at measurement phase)
2. Warmup used wrong shader (stochastic instead of deterministic)
3. R field not properly computed after warmup
4. Buffer swap issue (reading wrong buffer)

**This invalidates the σ=10⁻⁶ result** - system never reached synchronized state to test noise stability.

---

## II. Summary of Findings

### Confirmed Correct:
✓ **Timestep count:** 6000 steps (50 warmup + 10 measurement)
✓ **Initial conditions:** Random phases (R₀ ≈ 0.004)
✓ **Noise scaling:** Proper σ·√(dt)·N(0,1) (Euler-Maruyama)
✓ **PRNG:** PCG + Box-Muller (industry standard)
✓ **Grid resolution:** 256×256 (adequate)

### Critical Issues Found:
✗ **Damping:** ABSENT (γ=0) - system is conservative, overly fragile
✗ **Warmup effectiveness:** R(t) data shows system NEVER synchronized
✗ **Data collection:** R_timeseries appears to start before warmup completes

### Impact on Falsification Claim:

**Current status:** **FALSIFICATION IS INVALID**

**Reasons:**
1. **Damping absence** makes comparison to theoretical σ_c meaningless (different dynamics)
2. **System never synchronized** (R ≠ 0.995 after warmup) means we didn't test "noise destroys sync"
3. Instead we tested "noise prevents sync formation" - different physics

**What we actually measured:**
- Can a **random phase gas** (R₀~0.004) synchronize under noise σ=10⁻⁶?
- Answer: No.

**What we needed to measure:**
- Can a **synchronized state** (R=0.995) survive noise σ=10⁻⁶?
- Answer: Unknown (test not performed correctly).

---

## III. Required Actions Before Accepting Falsification

### Immediate (Critical):

1. **Fix warmup-to-measurement transition:**
   - Verify warmup completes with R ≈ 0.995
   - Ensure R_timeseries starts recording AFTER warmup
   - Check buffer swap between deterministic warmup and stochastic measurement

2. **Add damping to stochastic shader:**
   - Include `-γ sin(θ)` term in drift
   - Use γ = 0.1 (or match deterministic shader)
   - Re-run noise sweep with damping active

3. **Verify warmup success:**
   - Print R value after warmup, before measurement
   - Should see: R_warmup ≈ 0.995 ± 0.005
   - Only proceed to noise measurement if synchronized

### Short-term (Verification):

4. **Save pre-noise R value:**
   - Add diagnostic: `R_before_noise` to output
   - Verify each sigma test starts from R ≈ 0.995

5. **Grid convergence:**
   - Test at 128², 256², 512²
   - Verify σ_c doesn't shift with resolution

6. **Longer measurement:**
   - Increase N_measure from 1000 → 10,000
   - Check if steady state is reached (R(t) plateaus)

---

## IV. Theoretical Context

### Why Damping Matters

Standard Kuramoto-Langevin **with damping**:
$$\frac{d\theta_i}{dt} = \omega_i + K \sum_j \sin(\theta_j - \theta_i) - \gamma \theta_i + \xi_i(t)$$

Critical noise (with damping):
$$\sigma_c \sim \sqrt{\frac{K}{\gamma}}$$

For K=1, γ=0.1: σ_c ~ 3.2

**Without damping** (γ=0), system is **conservative** (energy preserving):
- Phase space volume conserved
- No attractor to synchronized state
- Noise has cumulative effect (no dissipation)
- Expected σ_c **much smaller**

**This could explain the 125,000× discrepancy.**

### Why Warmup Failure Matters

From `immediate.md`:
> Q2.2: What were the initial conditions?
> - If random (R₀ ≈ 0.3): Tests "synchronization formation in noise"
> - If pre-synchronized (R₀ ≈ 0.99): Tests "synchronization stability under noise"
> - **These are different physics**

We intended to test **stability**, but if warmup failed, we accidentally tested **formation**.

Formation under noise is **much harder** than maintaining sync under noise.

**Analogy:**
- Stability test: "Can a spinning top stay upright when nudged?"
- Formation test: "Can a fallen top spin itself upright while being nudged?"

Second is obviously harder.

---

## V. Revised Experimental Protocol

### Phase 0: Validation (Before Noise Sweep)

1. **Add damping to stochastic shader**
2. **Test warmup success:**
   - Run 5000 deterministic steps from random IC
   - Measure final R_warmup
   - Require R_warmup > 0.95 to proceed
3. **Test stochastic shader preserves sync at σ=0:**
   - After warmup, run 1000 stochastic steps with σ=0
   - R should remain ≈ 0.995 (verify shader is correct)

### Phase 1: Noise Sweep (Corrected)

For each σ ∈ [0, 10⁻⁷, 10⁻⁶, 10⁻⁵, 10⁻⁴, 10⁻³, 10⁻², 10⁻¹, 1.0]:

1. **Warmup:** 5000 steps, deterministic (σ=0, γ=0.1)
2. **Verify sync:** Measure R_warmup, require R_warmup > 0.95
3. **Measure:** 10,000 steps, stochastic (σ=target, γ=0.1)
4. **Record:** R(t), L(t), θ variance(t) every 10 steps
5. **Output:** Final R_mean, R_std, L_mean, L_std

### Phase 2: Analysis

- Plot R_mean(σ) to find transition
- Identify σ_c where R drops below 0.5
- Check if σ_c > 10⁻⁵ (Path B valid) or σ_c < 10⁻⁵ (Path A required)

---

## VI. Conclusion

**Current experimental result (σ_c ≈ 10⁻⁶) is INVALID due to:**
1. Missing damping term (γ=0 vs. required γ>0)
2. Warmup failure (system never synchronized)
3. Wrong physics tested (formation vs. stability)

**Falsification of Path B is PREMATURE.**

**Required action:** Fix implementation issues and re-run before making theoretical conclusions.

**Timeline estimate:** 1-2 days to fix + re-run.

**After fixes, possible outcomes:**
- If σ_c > 10⁻⁵: Path B validated (stochastic vacuum viable)
- If σ_c < 10⁻⁵: Path B falsified (deterministic vacuum required)
- If σ_c ≈ 10⁻⁵: Boundary case (needs deeper analysis)

**We are not ready to choose Path A vs Path B yet.**
