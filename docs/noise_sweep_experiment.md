# Noise Sweep Experiment: Path B Validation

**Experiment ID**: Phase 2 - Noise Sweep
**Reference**: Directive.md, Determinism.md
**Objective**: Measure critical noise amplitude σ_c to determine Path A vs Path B
**Date**: 2025-12-16

---

## I. Theoretical Background

### The Question

Can the Synchronization Mass Field Theory (SMFT) tolerate thermal noise while maintaining localized particle-like structures (solitons)?

**Path A (Deterministic)**: System is fundamentally deterministic. Any noise destroys synchronization.
- Prediction: σ_c ≪ 10^-5 (extremely fragile)
- Theory: Pure Kuramoto-based field theory

**Path B (Stochastic)**: System has intrinsic stochastic behavior at Planck scale.
- Prediction: σ_c > 10^-5 (robust to noise)
- Theory: Martin-Siggia-Rose (MSR) field theory with noise

### Physical Motivation

From Directive.md Section I (Q1):
> "Finite bandwidth → quantization error → noise at ℓ_P scale"

If the vacuum has finite information capacity (Holographic Principle), then noise at the Planck scale is unavoidable. The question is: **how much noise can the system survive?**

### The Critical Threshold

**Falsification Criterion** (from Directive.md Section III, Q7):

```
If σ_c < 10^-5: Abandon stochastic theory → Path A required ✗
If σ_c > 10^-5: Stochastic theory validated → Path B confirmed ✓
```

---

## II. Experimental Design

### Observable: Localization L

The primary observable is the **localization parameter**:

$$L = \int |R(x)|^4 \, d^2x$$

**Physical Interpretation**:
- L → ∞: Localized soliton (particle exists)
- L ~ O(1): Diffuse cloud (no particle)
- L → 0: Empty vacuum

**Critical Value**: L = 10^5 (from previous successful runs)

### Protocol (3 Weeks)

#### Week 1: Validation
1. **Day 1-2**: Implementation
   - Integrate PRNG (PCG) into compute shader ✓
   - Add noise term to Kuramoto integrator ✓
   - Compile, debug

2. **Day 3**: Verification Tests
   - Test 1: PRNG quality (mean ≈ 0, variance ≈ 1)
   - Test 2: Noise scaling (⟨θ²⟩ ∝ σ²T, independent of dt)
   - Test 3: Deterministic limit (σ=0 reproduces previous results)

**Gate**: All three tests MUST PASS before proceeding to noise sweep

#### Week 2: Noise Sweep
- Parameters: K = 27.21, Δ = 2.5, Grid = 256×256
- Noise values: σ ∈ [0, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1.0]
- For each σ:
  1. Warmup: 5000 steps with σ=0 (equilibrate)
  2. Measure: 1000 steps with noise amplitude σ
  3. Record: L(t), R(t), phase variance

**Outputs**:
- L(σ) curve
- R(σ) order parameter evolution
- Time series for each σ
- Spatial field snapshots

#### Week 3: Analysis & Decision
1. Fit L(σ) = L₀·exp(-σ/σ_c) to extract critical threshold
2. **Decision Point**:
   - σ_c > 10^-5: **Path B validated** → proceed with MSR formalism
   - σ_c < 10^-5: **Path A required** → abandon stochastic theory

---

## III. Implementation Details

### Stochastic Integrator (Euler-Maruyama)

**Deterministic Kuramoto**:
```
dθ/dt = ω + K·Σsin(θⱼ - θᵢ)
```

**Stochastic Kuramoto**:
```
dθ = [ω + K·Σsin(θⱼ - θᵢ)]dt + σ·dW
```

**Discrete Update**:
```
θ(t+dt) = θ(t) + drift·dt + σ·√(dt)·N(0,1)
```

### PRNG: PCG Algorithm

**Why PCG?**
- Industry standard for scientific computing
- Statistical quality superior to linear congruential generators
- Passes TestU01 BigCrush battery
- 2^64 period (sufficient for our grid: 256² = 65536 threads)

**Implementation** (kuramoto_stochastic.comp):
```glsl
uint pcg_state;

void pcg_init(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    pcg_state = (word >> 22u) ^ word;
}

float pcg_random() {
    pcg_state = pcg_state * 747796405u + 2891336453u;
    uint word = ((pcg_state >> ((pcg_state >> 28u) + 4u)) ^ pcg_state) * 277803737u;
    word = (word >> 22u) ^ word;
    return float(word) / 4294967296.0;
}
```

**Seeding Strategy**:
```
seed = x + y·65536 + timestep·16777216
```
Ensures each thread has unique, time-dependent seed.

### Gaussian Transform: Box-Muller

**Algorithm**:
```glsl
float randn() {
    float u1 = pcg_random();
    float u2 = pcg_random();
    u1 = max(u1, 1e-10);  // Avoid log(0)

    float r = sqrt(-2.0 * log(u1));
    float theta = 2.0 * π * u2;

    return r * cos(theta);  // N(0,1)
}
```

**Properties**:
- Generates standard normal N(0,1)
- Mean = 0, Variance = 1 (verified in Test 1)
- Two uniform → one Gaussian (50% efficiency)

---

## IV. Expected Results

### Hypothesis: σ_c ~ O(1)

**Reasoning** (from Directive.md Section XIII):
- Kuramoto synchronization is known to be robust to moderate noise
- Critical noise typically scales with coupling: σ_c ~ K
- K = 27.21 suggests σ_c ~ 10-50

**Predicted Phase Diagram**:

```
σ = 0:      L ~ 10^6    (perfect soliton)
σ = 10^-6:  L ~ 10^6    (no change)
σ = 10^-5:  L ~ 10^5.9  (slight weakening)
σ = 10^-4:  L ~ 10^5.5  (still stable)
σ = 10^-3:  L ~ 10^5    (threshold?)
σ = 10^-2:  L ~ 10^4    (degraded)
σ = 10^-1:  L ~ 10^2    (dissolving)
σ = 1.0:    L ~ 1       (destroyed)
```

**If prediction is correct**: Path B is validated ✓

---

## V. Verification Tests (Week 1, Day 3)

### Test 1: PRNG Quality

**Procedure**:
```cpp
// Generate N=10000 samples
float sum = 0, sum_sq = 0;
for (int i = 0; i < 10000; i++) {
    float x = randn();
    sum += x;
    sum_sq += x*x;
}
float mean = sum / 10000;
float variance = sum_sq / 10000 - mean*mean;
```

**Success Criteria**:
- |mean| < 0.1
- 0.9 < variance < 1.1

**If FAILS**: PRNG is broken → do NOT proceed

### Test 2: Noise Scaling

**Theory**: For stochastic differential equation,
```
⟨θ²(T)⟩ = σ²·T  (independent of dt!)
```

**Procedure**:
- Run simulation with σ=1, T=10 for dt ∈ [0.001, 0.01, 0.1]
- Measure ⟨θ²⟩ at t=T
- Should be constant across dt values

**Success Criteria**:
- ⟨θ²⟩ varies less than 10% across dt values
- Scales linearly with σ² and T

**If FAILS**: Wrong noise scaling → integrator is broken

### Test 3: Deterministic Limit

**Procedure**:
- Set σ = 0
- Run for 1000 steps
- Compare to previous deterministic runs (test_smft_phase0)

**Success Criteria**:
- Results match within numerical precision (~10^-6)

**If FAILS**: Stochastic shader has bug even at σ=0

---

## VI. Data Analysis

### Localization Time Series

For each σ, we measure L(t) over 1000 steps:

```python
import numpy as np
import matplotlib.pyplot as plt

# Load time series
L_series = np.loadtxt('output/noise_sweep/sigma_1.00e-03/L_timeseries.dat')

# Compute statistics
L_mean = np.mean(L_series[:, 1])
L_std = np.std(L_series[:, 1])

print(f"σ = 10^-3: L = {L_mean:.2e} ± {L_std:.2e}")
```

### Critical Threshold Extraction

Fit exponential decay model:

$$L(\sigma) = L_0 \cdot e^{-\sigma/\sigma_c}$$

```python
from scipy.optimize import curve_fit

sigma_values = np.array([0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0])
L_means = [...]  # measured values

def model(sigma, L0, sigma_c):
    return L0 * np.exp(-sigma / sigma_c)

params, cov = curve_fit(model, sigma_values, L_means)
L0_fit, sigma_c_fit = params

print(f"Critical threshold: σ_c = {sigma_c_fit:.2e}")
```

### Decision Tree

```
σ_c > 10^-5?
├─ YES → Path B validated ✓
│         ├─ Formalize MSR action
│         ├─ Calculate σ_physical = ℓ_P
│         ├─ Verify σ_physical ≪ σ_c
│         └─ Publish: "Stochastic SMFT"
│
└─ NO → Path A required ✗
          ├─ Abandon MSR formalism
          ├─ Write deterministic Lagrangian
          ├─ Noise sweep → Appendix A
          └─ Publish: "Deterministic SMFT"
```

---

## VII. File Organization

### Output Structure

```
build/output/noise_sweep/
├── results.csv                    # Summary: σ, L_mean, L_std, R_mean, R_std
├── sigma_0.00e+00/
│   ├── L_timeseries.dat
│   ├── R_timeseries.dat
│   ├── phase_var_timeseries.dat
│   ├── step_0000/
│   │   ├── R_field.dat
│   │   └── theta_field.dat
│   ├── step_0500/
│   └── step_0999/
├── sigma_1.00e-06/
│   └── ...
├── sigma_1.00e-05/  (critical threshold)
│   └── ...
└── sigma_1.00e+00/
    └── ...
```

### Key Files

**Implementation**:
- `shaders/smft/kuramoto_stochastic.comp` - Stochastic integrator with PRNG
- `src/SMFTEngine.h` - Added `stepStochastic()` method
- `src/SMFTEngine.cpp` - Implementation of stochastic stepping
- `test/test_noise_sweep.cpp` - Main experiment driver

**Documentation**:
- `docs/noise_sweep_experiment.md` - This file
- `Directive.md` - Authorization and protocol
- `Determinism.md` - Theoretical analysis

---

## VIII. Success Criteria Summary

### Technical Success
- ✓ All verification tests pass (Week 1, Day 3)
- ✓ Noise sweep completes for all σ values
- ✓ Clear σ_c measurement with <20% uncertainty
- ✓ Results reproducible across multiple runs

### Scientific Success
- ✓ Hypothesis tested rigorously
- ✓ Results interpreted honestly
- ✓ Theory matches experiment
- ✓ Falsifiability demonstrated

**Either Path A or Path B can achieve scientific success.**

### Failure Modes (AVOID)
- ✗ Verification tests fail → implementation broken
- ✗ σ_c ambiguous (large error bars) → need finer σ sampling
- ✗ Results contradict theory but theory kept anyway → intellectual dishonesty
- ✗ Claim validation without evidence → unscientific

---

## IX. Timeline & Milestones

### Week 1: Implementation + Verification (Dec 16-22)
- [x] **Dec 16**: Stochastic shader created
- [x] **Dec 16**: Test driver created
- [ ] **Dec 17**: Verification Test 1 (PRNG quality)
- [ ] **Dec 18**: Verification Test 2 (noise scaling)
- [ ] **Dec 19**: Verification Test 3 (deterministic limit)
- [ ] **Dec 20**: All tests passing → GATE PASSED

### Week 2: Noise Sweep (Dec 23-29)
- [ ] **Dec 23-24**: σ ∈ [0, 10^-6, 10^-5] runs
- [ ] **Dec 25-26**: σ ∈ [10^-4, 10^-3, 10^-2] runs
- [ ] **Dec 27**: σ ∈ [10^-1, 1.0] runs
- [ ] **Dec 28-29**: Data analysis, fitting σ_c

### Week 3: Decision (Dec 30 - Jan 5)
- [ ] **Dec 30**: Path A or B decision based on σ_c
- [ ] **Dec 31 - Jan 3**: Begin formalization of chosen path
- [ ] **Jan 4-5**: Draft initial paper sections

---

## X. References

1. **Directive.md** - Experimental protocol authorization (Claude analysis)
2. **Determinism.md** - Effective stochasticity hypothesis (Claude analysis)
3. **Dirac-Anomaly.md** - Defect evolution analysis
4. Mori, H. (1965). "Transport, Collective Motion, and Brownian Motion"
5. Zwanzig, R. (1961). "Memory Effects in Irreversible Thermodynamics"
6. O'Malley, D. et al. (2014). "PCG: A Family of Better Random Number Generators"
7. Box, G. E. P., Muller, M. E. (1958). "A Note on the Generation of Random Normal Deviates"

---

## XI. Contact & Questions

**Experiment Lead**: [Your Name]
**Computational Platform**: Nova Vulkan Engine
**Hardware**: [GPU Model]
**Date Started**: 2025-12-16

**Questions to resolve before Week 2**:
1. Defect definition discrepancy (3755 defects vs "defect-free")?
2. Mass range [0.103, 1.677] vs expected [0.025, 2.5]?
3. Grid convergence: 128² vs 256² comparison?

---

**Status**: Phase 1 (Implementation) - In Progress
**Next Milestone**: Verification Test 1 (PRNG Quality) - Due Dec 17
**Expected Completion**: Jan 5, 2026
