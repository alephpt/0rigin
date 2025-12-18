# Dirac Coupling Experiment: Particle Localization Test

**Experiment ID**: Phase 3 - Dirac Field Coupling
**Reference**: Dirac-Anomaly.md Section VI-VII
**Objective**: Test if Dirac spinors localize in vacuum defects and produce discrete energy levels
**Date**: 2025-12-16

---

## I. Theoretical Prediction

### The Core MSFT Hypothesis

**Vacuum provides continuous mass landscape**:
$$m(x) = \Delta \cdot R(x)$$

where R(x) ∈ [0, 1] is the synchronization field.

**Dirac equation in position-dependent mass**:
$$[i\hbar\partial_t - c\alpha \cdot \mathbf{p} - \beta m(x)c^2]\Psi = 0$$

**Key Prediction** (from Dirac-Anomaly.md):
> "Even though m(x) is continuous, the bound state energies E_n are discrete (quantized by boundary conditions, like particle in a box)"

### Analogy: Hydrogen Atom

- **Coulomb potential** V(r) = -e²/r is continuous
- **But electron energies** E_n = -13.6 eV/n² are discrete
- **Quantization from** wavefunction boundary conditions, NOT potential discreteness

### For MSFT

- **Vacuum potential** m(x) is continuous (proven in defect evolution tests)
- **But Dirac bound states** should have discrete energies E₁, E₂, E₃, ...
- **Physical particles** have masses M_n = E_n/c² corresponding to these levels

**If correct**: This explains electron, muon, tau mass spectrum from single theory!

---

## II. Five Success Criteria

From Dirac-Anomaly.md Section VI, we require:

### Criterion 1: Localization at Defects

**Measurement**: Overlap integral
$$O = \frac{\int (1-R) \cdot |\Psi|^2 \, dx}{\int |\Psi|^2 \, dx}$$

**Success**: O > 0.7 (70% of Ψ density at defects)
**Failure**: O < 0.3 (Ψ spread uniformly)

**Physical meaning**: Dirac wavefunctions should preferentially localize where R is low (defects), not spread uniformly across the grid.

### Criterion 2: Stabilization of Deep Defects

**Measurement**: Survival rate of defects with R_core < 0.3

From vacuum-only runs (Dirac-Anomaly.md):
- **Without Dirac**: S₀ ≈ 70% (1564/5319 annihilated)
- **With Dirac**: S_D = ?

**Success**: S_D > S₀ + 10% (stabilization effect observed)
**Failure**: S_D ≈ S₀ (no effect)

**Physical meaning**: Localized Dirac field should stabilize the defects via feedback, preventing annihilation.

### Criterion 3: Discrete Energy Levels

**Measurement**: Histogram of bound state energies
$$E_n = \frac{\int_{particle\, n} \bar{\Psi}H\Psi}{\int_{particle\, n}|\Psi|^2}$$

**Success**: 2-5 distinct peaks in histogram
**Failure**: Smooth continuous distribution

**Physical meaning**: THIS IS THE KEY TEST. If we see discrete peaks, it proves quantization emerges from Dirac dynamics despite continuous m(x).

### Criterion 4: Particle Number ≪ Defect Number

**Measurement**:
- N_particles: count localized |Ψ|² peaks
- N_defects: ~3,755 cells or ~75 regions (needs clarification from previous tests)

**Success**: N_particles ∈ [10, 200] (selective binding to deep wells)
**Failure**: N_particles > 1000 (every defect hosts Ψ) or N_particles < 5 (too selective)

**Physical meaning**: Not every defect should trap a particle. Only deep potential wells (R_core ≈ 0) should bind.

### Criterion 5: Long-Time Stability

**Measurement**: Does ⟨|Ψ|²⟩ remain constant or decay?

**Success**: Fluctuates around constant (particles are stable)
**Failure**: Decays exponentially (particles dissolve)

**Physical meaning**: Bound states should be stable for t ≫ 100 timesteps. If they decay, the coupling is too weak or mechanism is wrong.

---

## III. Implementation: Dirac RK4 Shader

### Dirac Equation (2D+1)

In 2D, the Dirac equation becomes:
$$i\hbar\frac{\partial\Psi}{\partial t} = [c\alpha_x p_x + c\alpha_y p_y + \beta m(x,y)c^2]\Psi$$

where Ψ is a 2-component spinor (reduced from 4D).

### Discretization

**Spatial derivatives** (finite difference):
$$p_x \Psi = -i\hbar\frac{\Psi(x+\Delta x) - \Psi(x-\Delta x)}{2\Delta x}$$

**Time evolution** (RK4):
```
k1 = f(Ψ_n)
k2 = f(Ψ_n + dt/2 * k1)
k3 = f(Ψ_n + dt/2 * k2)
k4 = f(Ψ_n + dt * k3)

Ψ_{n+1} = Ψ_n + dt/6 * (k1 + 2k2 + 2k3 + k4)
```

where f(Ψ) = (-i/ℏ)H Ψ.

### Pauli Matrices (2D Spinor)

$$\alpha_x = \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$$

$$\alpha_y = \sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}$$

$$\beta = \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$$

### Chiral Coupling (MSFT-Specific)

From R.Christopher papers, the mass couples via:
$$m_{eff} = \Delta \cdot R(x) \cdot e^{i\theta(x)\gamma^5}$$

This introduces phase dependence, creating "chiral mass".

---

## IV. Feedback Mechanism

### Vacuum → Matter (Forward)

1. Kuramoto phases θ(x,y) evolve
2. Synchronization R(x,y) = |⟨e^(iθ)⟩| computed
3. Mass field m(x,y) = Δ·R(x,y) provides Dirac potential
4. Dirac equation evolves Ψ(x,y,t)

### Matter → Vacuum (Backward)

5. Spinor density ρ(x,y) = |Ψ|² computed
6. Feedback force F = -λ·ρ·sin(θ) added to Kuramoto equation
7. High ρ → suppress phase synchronization → maintain low R → stabilize defect

**This is the "soliton handoff" mechanism.**

---

## V. Expected Results (λ_D = 1.0)

From Dirac-Anomaly.md Section VI Table:

| Observable | Predicted | Acceptable Range | Failure Mode |
|------------|-----------|------------------|--------------|
| Localization | O = 0.8 | [0.6, 0.95] | O < 0.5 → dispersion |
| Particle count | N = 50 | [20, 150] | N < 10 or N > 500 |
| Energy levels | n = 3 | [2, 5] | n = 1 or n > 10 |
| Stabilization | ΔS = +15% | [+5%, +30%] | ΔS < 0 → destabilizing |
| Lifetime | τ > 1000 | [500, ∞] | τ < 100 → unstable |

---

## VI. Experimental Protocol

### Phase 1: Initial Vacuum Equilibration (500 steps)

1. Run vacuum-only (no Dirac field)
2. Let R(x,y) reach equilibrium
3. Identify defects (regions with R < 0.3)
4. Save final state as initial condition for Dirac coupling

### Phase 2: Dirac Field Initialization

**Option A: Ground state seed**
- Place small Gaussian wavepackets at defect centers
- $$\Psi(x,y,0) = A \cdot e^{-(x-x_0)^2/2\sigma^2} \begin{pmatrix} 1 \\ 0 \end{pmatrix}$$
- Let system relax

**Option B: Plane wave**
- Uniform Ψ(x,y,0) = const
- Let localization emerge spontaneously
- More rigorous but slower convergence

**Recommendation**: Option A for initial test (faster), Option B for validation.

### Phase 3: Coupled Evolution (2000 steps)

1. Evolve both fields together:
   - Kuramoto: θ(t) with spinor feedback
   - Dirac: Ψ(t) with mass m(x) = Δ·R(x)

2. Measure every 10 steps:
   - Localization O(t)
   - Particle count N(t)
   - Energy histogram
   - Defect survival rate

3. Save snapshots at t = [0, 500, 1000, 1500, 2000]

### Phase 4: Analysis

**Energy Spectrum Analysis**:
```python
# Identify localized peaks in |Ψ|²
peaks = find_peaks(psi_density, threshold=0.1)

# For each peak, compute local energy
energies = []
for peak in peaks:
    region = extract_region(peak, radius=5)
    E = compute_energy(psi[region], m[region])
    energies.append(E)

# Histogram
plt.hist(energies, bins=50)
plt.xlabel('Energy (eV)')
plt.ylabel('Count')
plt.title('Bound State Energy Spectrum')
```

**Look for discrete peaks** at E₁, E₂, E₃, ...

---

## VII. Quantitative Predictions

### Energy Level Spacing

For 2D potential well with depth V₀ ≈ Δ·1 = 2.5 eV and size a ≈ 10 grid cells:

$$E_n \approx \frac{\hbar^2 \pi^2 n^2}{2m_e a^2} - V_0$$

**Estimate** (using grid spacing Δx = 1 nm):
- E₁ ≈ -2.0 eV (ground state)
- E₂ ≈ -1.5 eV (first excited)
- E₃ ≈ -0.8 eV (second excited)

**Spacing**: ΔE ≈ 0.5-0.7 eV

**Test**: Measured histogram should show peaks separated by ~0.5 eV.

### Localization Length

**Theory**: Bound state size ξ ∼ ℏ/√(2m·V₀)

For V₀ = 2.5 eV, m = m_e:
$$\xi \approx 0.5 \text{ nm} \approx 0.5 \text{ grid cells}$$

**Very tight localization expected!**

### Particle Density

From Dirac-Anomaly.md analysis:
- Total defect cells: ~3,755
- Deep defects (R < 0.3): ~1,000 cells
- Number of distinct deep wells: ~50

**Prediction**: N_particles ≈ 50 (one per deep well)

---

## VIII. Failure Modes and Diagnostics

### Failure Mode 1: No Localization (O < 0.3)

**Possible causes**:
1. Dirac coupling λ_D too weak
2. Mass contrast insufficient (Δ too small)
3. Defects too shallow (no deep wells)
4. Numerical instability (dt too large)

**Diagnostic**: Increase Δ from 2.5 → 5.0, re-run

### Failure Mode 2: Continuous Energy Spectrum (No Peaks)

**Possible causes**:
1. Ψ is delocalized (see Failure 1)
2. Potential wells too shallow (E_binding < thermal fluctuations)
3. Grid resolution too coarse (bound states unresolved)

**Diagnostic**: Increase grid resolution 256² → 512²

### Failure Mode 3: All Defects Trap Particles (N > 1000)

**Possible causes**:
1. Binding energy too high (all wells, even shallow, trap Ψ)
2. Initialization placed Ψ everywhere
3. No energetic barrier to prevent trapping

**Diagnostic**: Reduce initial Ψ amplitude, check binding energy calculation

### Failure Mode 4: Particles Decay (τ < 100)

**Possible causes**:
1. Feedback destabilizes instead of stabilizes
2. Numerical instability
3. Wells are metastable, not truly bound

**Diagnostic**: Check sign of feedback term, reduce dt, verify energy conservation

---

## IX. Validation Against Defect Evolution Data

From Dirac-Anomaly.md, the vacuum-only run gave:
- Initial defects: 5,319
- Final defects: 3,755
- Survival rate: 70.6%

**Prediction with Dirac coupling**:
- Deep defects stabilized: ~85-95% survival (ΔS = +15-25%)
- Shallow defects still annihilate: ~50% survival

**Measured**:
```
N_deep_initial = count_defects(R < 0.3, t=0)
N_deep_final = count_defects(R < 0.3, t=2000)
S_deep = N_deep_final / N_deep_initial

ΔS = S_deep - 0.706
```

**Success**: ΔS > 0.10 (10% improvement in survival)

---

## X. Timeline & Milestones

### Week 1: Implementation (Dec 16-22)
- [x] **Dec 16**: Created documentation
- [ ] **Dec 17**: Implement dirac_rk4.comp shader
- [ ] **Dec 18**: Implement spinor_feedback.comp shader
- [ ] **Dec 19**: Create test_dirac_coupling.cpp driver
- [ ] **Dec 20**: Compile and debug

### Week 2: Initial Tests (Dec 23-29)
- [ ] **Dec 23**: Run Phase 1 (vacuum equilibration)
- [ ] **Dec 24**: Run Phase 2 (Dirac initialization)
- [ ] **Dec 25**: Run Phase 3 (coupled evolution, λ_D = 0.1)
- [ ] **Dec 26**: Run Phase 3 (λ_D = 1.0)
- [ ] **Dec 27**: Run Phase 3 (λ_D = 10.0)
- [ ] **Dec 28-29**: Analyze energy spectra

### Week 3: Validation (Dec 30 - Jan 5)
- [ ] **Dec 30**: Check all 5 success criteria
- [ ] **Dec 31**: Grid convergence test (256² → 512²)
- [ ] **Jan 1-2**: Defect survival analysis
- [ ] **Jan 3-4**: Energy level identification
- [ ] **Jan 5**: Final report + decision

---

## XI. Output Structure

```
build/output/dirac_coupling/
├── results_summary.txt           # All 5 criteria: PASS/FAIL
├── lambda_0.1/
│   ├── psi_density_t0000.dat
│   ├── psi_density_t0500.dat
│   ├── psi_density_t1000.dat
│   ├── energy_spectrum.dat
│   ├── localization_timeseries.dat
│   └── defect_survival.dat
├── lambda_1.0/
│   └── ...
├── lambda_10.0/
│   └── ...
└── plots/
    ├── energy_histogram.png      # KEY: Shows discrete levels?
    ├── localization_vs_time.png
    ├── psi_density_overlay.png   # |Ψ|² overlaid on R(x,y)
    └── defect_survival.png
```

---

## XII. Success Decision Tree

```
Criterion 1 (Localization): O > 0.7?
├─ NO → FAIL: Dirac doesn't localize at defects
│         └─ Re-examine coupling strength, potential depth
│
└─ YES → Proceed to Criterion 2

Criterion 2 (Stabilization): ΔS > 10%?
├─ NO → FAIL: No feedback stabilization observed
│         └─ Check feedback implementation, sign, amplitude
│
└─ YES → Proceed to Criterion 3

Criterion 3 (Discrete Energies): 2-5 peaks in histogram?
├─ NO → FAIL: No quantization despite localization
│         └─ THIS WOULD BE SURPRISING - investigate carefully
│
└─ YES → Proceed to Criterion 4

Criterion 4 (Particle Count): 10 < N < 200?
├─ NO → FAIL: Either too promiscuous or too selective
│         └─ Adjust potential depth or binding criterion
│
└─ YES → Proceed to Criterion 5

Criterion 5 (Stability): τ > 500?
├─ NO → FAIL: Particles are transient, not stable
│         └─ Increase coupling or check energy conservation
│
└─ YES → SUCCESS ✓✓✓
          ALL 5 CRITERIA MET
          → MSFT mechanism validated
          → Particle generation from vacuum confirmed
          → Publish results
```

---

## XIII. References

1. **Dirac-Anomaly.md** - Critical predictions and success criteria
2. **Directive.md** - Experimental methodology
3. **Determinism.md** - Theoretical foundation
4. Dirac, P. A. M. (1928). "The Quantum Theory of the Electron"
5. Christopher, R. (2023). "Soliton Handoff Mechanism in MSFT"
6. Thaller, B. (1992). "The Dirac Equation" (textbook reference)

---

**Status**: Phase 1 (Implementation) - In Progress
**Next Milestone**: Dirac RK4 shader implementation - Due Dec 17
**Critical Test**: Energy spectrum discrete peaks (Criterion 3)
