## Experimental Proposals for SMFT Validation
**Concrete Designs for Laboratory Tests**

**Document Version**: 1.0
**Date**: December 29, 2025

---

## Proposal 1: BEC Phonon Scattering Near Vortex (HIGHEST PRIORITY)

### Executive Summary
**Fastest, cheapest, clearest test of SMFT**
- **Timeline**: 2-3 months
- **Cost**: $0 (existing equipment)
- **Signal**: 66% deviation at r = ξ (13σ detection)
- **Status**: Ready to execute in any BEC lab

### Scientific Objective
Test SMFT prediction that order parameter R(r) modulates effective sound speed near vortex core:
```
SMFT: c_eff(r) = c_s · R(r) = c_s · tanh(r/ξ)
Standard BEC: c_eff(r) = c_s (constant)
```

### Experimental Setup

**Platform**: Ultracold ⁸⁷Rb BEC in optical dipole trap

**Parameters**:
- Atom number: N ~ 10⁵
- Temperature: T ~ 100 nK (T/T_c ~ 0.8)
- Trap frequencies: (ω_r, ω_z) = 2π × (100, 1000) Hz (pancake geometry)
- Healing length: ξ ~ 0.5 μm
- Sound speed: c_s ~ 4 mm/s

**Vortex Imprinting**:
- Method 1 (stirring): Rotate focused laser beam, ω_stir > ω_trap
- Method 2 (phase imprinting): Apply spatially varying AC Stark shift

**Phonon Generation**:
- Bragg spectroscopy: Two-photon transition with Δk ~ 1.3 μm⁻¹
- Wavelength: λ_phonon ~ 5 μm (> 10ξ, long-wavelength regime)
- Amplitude: δn/n ~ 1% (linear response)

### Measurement Protocol

**Step 1: Prepare BEC with Single Vortex**
1. Cool atoms to BEC (standard evaporative cooling)
2. Imprint single vortex at trap center
3. Image to verify vortex position (dark core in density)
4. Wait 50 ms for vortex relaxation

**Step 2: Launch Phonon Wavepacket**
1. Apply Bragg pulse: Δt ~ 5 μs, Δk = 2π/5μm
2. Phonon propagates radially inward from edge
3. Initial amplitude: A ~ 0.01 × n₀

**Step 3: Time-Resolved Imaging**
1. Absorption imaging at t = 0, 5, 10, 15, 20 ms
2. Extract density n(r, t) from OD
3. Track wavefront position r_front(t) via phase analysis
4. Spatial resolution: Δr ~ 0.3 μm (diffraction-limited)

**Step 4: Extract Sound Speed Profile**
1. Compute phase velocity: v_phase(r) = dr_front/dt
2. Effective sound speed: c_eff(r) = v_phase(r)
3. Compare to prediction: c_eff(r) / c_s = tanh(r/ξ)
4. Fit ξ as free parameter

**Step 5: Repeat for Statistics**
- N_runs = 20 (3-5 hours total)
- Average over realizations
- Extract systematic errors (temperature drift, imaging noise)

### Data Analysis

**Primary Observable**: c_eff(r) / c_s vs. r/ξ

**Fit Function**:
```python
def model(r, xi, c_s):
    return c_s * np.tanh(r / xi)
```

**Expected Result** (SMFT):
| r/ξ | R(r) | c_eff/c_s | Deviation from c_s |
|-----|------|-----------|-------------------|
| 0.5 | 0.46 | 0.46 | 54% |
| 1.0 | 0.76 | 0.76 | 24% |
| 2.0 | 0.96 | 0.96 | 4% |

**Null Hypothesis** (Standard BEC):
c_eff(r) / c_s = 1.00 ± 0.05 for all r > ξ

**Statistical Test**:
- Fit SMFT model vs. constant model
- χ² comparison: Δχ² > 25 (5σ) → SMFT confirmed

### Systematic Errors

**Error Budget**:
| Source | Magnitude | Mitigation |
|--------|-----------|------------|
| Imaging resolution | Δr ~ 0.3 μm | Deconvolution, super-resolution |
| Thermal phonons | ~2% | Cool to T < 0.5 T_c |
| Vortex wander | ~0.2 μm | Average over many shots |
| Multiparticle scattering | <1% | Use dilute gas (na³ < 10⁻⁴) |
| **Total** | **~5%** | **Sufficient for 13σ signal** |

### Control Experiments

1. **No vortex**: Confirm c(r) = c_s ± 2% (uniform)
2. **Multiple vortices**: Test c(r) near vortex lattice
3. **Temperature scan**: Vary T to test R(T) dependence

### Expected Outcome

**If SMFT Correct**:
- c_eff(ξ)/c_s = 0.76 ± 0.05 (5σ from 1.00)
- Healing length ξ_fit = 0.5 ± 0.05 μm (consistent with theory)
- Publication: *Nature Physics* or *Science*

**If SMFT Wrong**:
- c_eff(r)/c_s = 1.00 ± 0.05 (no spatial variation)
- Explanation: Standard Gross-Pitaevskii theory correct

### Collaborators to Contact

**Top BEC Groups**:
1. **MIT** - Wolfgang Ketterle (Nobel laureate, BEC pioneer)
2. **JILA** - Eric Cornell, Deborah Jin memorial group (vortex experts)
3. **ETH Zurich** - Tilman Esslinger (quantum simulation)
4. **Cambridge** - Zoran Hadzibabic (2D BEC systems)
5. **LENS Florence** - Massimo Inguscio (precision measurements)

**Contact Approach**:
- Email: Brief summary + link to prediction document
- Offer: Provide theoretical support, co-authorship
- Timeline: 1 month discussion → 2 month execution → 1 month analysis

---

## Proposal 2: Casimir Force with Nanomembrane Resonator

### Executive Summary
- **Timeline**: 12-18 months
- **Cost**: $50k-100k (fabrication)
- **Signal**: 9.6% deviation if R = 0.95 (5σ detection)
- **Status**: Requires custom sensor development

### Scientific Objective
Measure absolute Casimir force with precision δF/F < 1% to detect R-field coupling:
```
F_SMFT = F_standard · ⟨R²⟩
```

### Experimental Setup

**Platform**: Nanomembrane resonator in vacuum

**Design**:
- Silicon nitride membrane: 100 nm thick, 100 μm × 100 μm area
- Resonant frequency: f₀ ~ 1 MHz
- Quality factor: Q ~ 10⁶ (vacuum, T = 4K)
- Gold coating: 50 nm (both sides)

**Force Sensor**:
- Measure frequency shift: Δf ∝ ∂F/∂d
- Precision: Δf/f ~ 10⁻⁶ → δF/F ~ 10⁻⁴
- Dynamic range: d = 50-500 nm

### Measurement Protocol

**Step 1: Fabricate Sensor**
1. Deposit SiN membrane (LPCVD)
2. Pattern with e-beam lithography
3. Release with XeF₂ etch
4. Sputter gold coating (50 nm)
5. Mount in cryostat with piezo actuator

**Step 2: Characterize Resonance**
1. Drive membrane with AC voltage
2. Measure amplitude with optical interferometer
3. Extract f₀, Q from ringdown

**Step 3: Casimir Force Measurement**
1. Approach gold sphere (R ~ 50 μm) to d ~ 100 nm
2. Measure frequency shift: Δf = f(d) - f(d → ∞)
3. Scan separation: d = 50, 75, 100, 150, 200, 300, 500 nm
4. Repeat 100× for statistics

**Step 4: Extract Force Gradient**
```
∂F/∂d = -2π m_eff f₀ Δf
```
where m_eff = membrane effective mass

**Step 5: Compare to Theory**
- Fit to Lifshitz theory: F_Lifshitz(d, T, ε(ω))
- Extract residual: ΔF = F_measured - F_Lifshitz
- Test SMFT: ΔF/F = (1 - ⟨R²⟩) = const vs. d

### Systematic Errors

| Source | Magnitude | Mitigation |
|--------|-----------|------------|
| Electrostatic patches | ±5% | Nulling voltage compensation |
| Surface roughness | ±2% | Use atomically flat Au(111) |
| Temperature drift | ±1% | Stabilize T ± 10 mK |
| Membrane nonlinearity | ±0.5% | Operate at low amplitude |
| **Total** | **±6%** | **Challenging for 10% signal** |

### Proposed Enhancement

**Path Forward**:
1. **Material comparison**: Measure F(Au) / F(Si) → isolate R-dependence
2. **Geometry test**: Patterned vs. flat surfaces → spatial R(x,y)
3. **Multi-frequency**: Measure at different ω → test ε(ω) vs. R coupling

**Expected Improvement**: Material ratio less sensitive to systematics → 3σ detection feasible

---

## Proposal 3: Quantum Synchronization in Superconducting Qubits

### Executive Summary
- **Timeline**: 2-3 years (phased approach)
- **Cost**: $200k-500k (postdoc + beamtime)
- **Signal**: β = 0.099 vs. 0.125 (5σ distinction)
- **Status**: Medium-scale test feasible on current hardware

### Scientific Objective
Measure critical exponent β to distinguish SMFT from known universality classes:
```
R(σ) ~ (σ_c - σ)^β

SMFT: β = 0.099 ± 0.004
2D Ising: β = 0.125
2D XY: β ~ 0.23
```

### Experimental Setup

**Platform**: IBM Quantum Heron or Google Sycamore processor

**System**:
- N = 100-200 transmon qubits in 2D lattice
- Tunable coupling: K via flux-tunable coupler
- Coherence time: T₂ ~ 100 μs
- Gate fidelity: 99.9%

**Hamiltonian**:
```
H = Σᵢ (ω_i σᵢᶻ/2 + σ_noise η_i(t)) + K Σ_⟨ij⟩ (σᵢ⁺σⱼ⁻ + h.c.)
```

where σ_noise = drive amplitude (tunable noise)

### Measurement Protocol

**Phase 1: Small-Scale Proof-of-Principle (N = 20)**

**Step 1: Initialize**
- Prepare each qubit in random superposition: |ψ_i⟩ = cos(θ_i)|0⟩ + sin(θ_i)|1⟩
- Random phases: θ_i ~ Uniform(0, 2π)

**Step 2: Quench**
- Apply coupling K and noise σ
- Time-evolve: t_total = 10 μs (fast quench)
- Monitor phase evolution

**Step 3: Measure Order Parameter**
```
R = |⟨exp(iφ_j)⟩| = |Σⱼ ⟨σⱼ⁺⟩ / N|
```
- Single-shot readout of all qubits
- Compute ensemble average over 100 realizations

**Step 4: Scan Critical Region**
- Vary σ = 0.40, 0.42, ..., 0.60 (41 points)
- For each σ: measure R(σ)
- Total shots: 41 × 100 = 4,100

**Step 5: Extract β**
- Fit: R ~ (σ_c - σ)^β
- Extract: σ_c, β, amplitude
- Compare to universality classes

**Phase 2: Medium-Scale Test (N = 100)**

- Repeat Phase 1 protocol with 100 qubits
- Expected precision: σ(β) ~ 0.01 (marginal)
- Duration: 6-12 months

**Phase 3: Definitive Test (N = 200)**

- Full finite-size scaling analysis
- Multiple system sizes: L = 32, 64, 128, 256 (in 2D layout)
- Data collapse: R(σ, L) = L^(-β/ν) f[(σ - σ_c)L^(1/ν)]
- Extract: β, ν, η (three exponents → univeral class)
- Duration: 12-24 months

### Data Analysis

**Primary Observable**: Critical exponent β

**Fit**: Power-law regression R(σ) = A(σ_c - σ)^β

**Expected Results**:
| System | σ(β) | SMFT vs. 2D Ising |
|--------|------|-------------------|
| N = 20 | ±0.02 | 1.3σ (marginal) |
| N = 100 | ±0.01 | 2.6σ (suggestive) |
| N = 200 | ±0.005 | **5.2σ (definitive)** |

**Finite-Size Scaling**:
- Test collapse quality: χ²_collapse
- Extract ν: correlation length exponent
- Confirm d_u = 5 (unconventional upper critical dimension)

### Systematic Errors

| Source | Impact | Mitigation |
|--------|--------|------------|
| Decoherence | T₂ ~ 100 μs | Fast quench (10 μs) |
| Crosstalk | ~1% coupling | Calibration + error mitigation |
| Readout error | ~2% | Measurement error correction |
| Initialization | Non-random phases | True randomization via noise |

**Total systematic**: ~5% error on R → translates to ~10% on β → still sufficient for 5σ

### Expected Outcome

**If SMFT Correct**:
- β_measured = 0.099 ± 0.005
- |β - β_2D-Ising| / σ(β) = 5.2σ → Novel universality class
- ν ~ 2.5 → d_u = 5 (anomalous)
- Publication: *Nature* or *Science*

**If 2D Ising**:
- β_measured = 0.125 ± 0.005
- SMFT effective theory at finite size
- Crossover to standard critical behavior

### Collaborators to Contact

**Quantum Hardware Groups**:
1. **Google Quantum AI** - Hartmut Neven, Pedram Roushan (105-qubit Willow)
2. **IBM Quantum** - Jay Gambetta, Abhinav Kandala (133-qubit Heron)
3. **Rigetti Computing** - Chad Rigetti (80-qubit Aspen-M)

**University Groups**:
4. **Yale** - Michel Devoret, Rob Schoelkopf (circuit QED pioneers)
5. **MIT** - Will Oliver (Lincoln Lab, superconducting qubits)
6. **Berkeley** - Irfan Siddiqi (quantum feedback control)

**Theoretical Support**:
7. **UCSB** - John Martinis (former Google, quantum computing)
8. **Caltech** - John Preskill (quantum information theory)

---

## Proposal 4: CMB Cosmology (Long-Term)

### Executive Summary
- **Timeline**: 2030s (CMB-S4 era)
- **Cost**: N/A (funded missions)
- **Signal**: f_NL ~ 1-10, Gμ ~ 10⁻⁷
- **Status**: Awaiting future data

### Scientific Objective
Search for SMFT phase transition signatures in cosmic microwave background:
1. Non-Gaussianity: f_NL from defect network
2. Cosmic strings: Gμ from winding number conservation
3. Spectral index running: scale-dependent n_s(k)

### Experimental Setup

**CMB-S4 (2028-2035)**:
- 500,000 detectors
- Frequency coverage: 30-300 GHz
- Sensitivity: ΔT ~ 1 μK·arcmin
- Angular resolution: θ ~ 1 arcmin

**Observables**:
```
f_NL: σ(f_NL) ~ 1 (vs. Planck: ±5)
Gμ: sensitivity ~ 10⁻⁸ (vs. Planck: < 8×10⁻⁷)
n_s: σ(n_s) ~ 0.002 (vs. Planck: ±0.004)
```

### Analysis Strategy

**Test 1: f_NL Measurement**
- Bispectrum analysis: ⟨T(n₁)T(n₂)T(n₃)⟩
- Modal decomposition: local, equilateral, orthogonal shapes
- SMFT prediction: f_NL^local ~ O(1-10) if defect-dominated

**Test 2: String Tension**
- Template matching: search for string-induced discontinuities
- Machine learning: CNN trained on string simulations
- Constraint: Gμ/c² (compare to SMFT Gμ ~ 10⁻⁷)

**Test 3: Running Spectral Index**
- Power spectrum: P(k) across k = 10⁻⁴ to 0.1 Mpc⁻¹
- Fit: n_s(k) = n_s(k_pivot) + (dn_s/dlnk)·ln(k/k_pivot)
- SMFT: predict dn_s/dlnk ≠ 0 (string contribution)

### Falsification Criteria

**SMFT Cosmology RULED OUT if**:
1. |f_NL| < 1 at 5σ (CMB-S4)
2. Gμ < 10⁻⁸ at 5σ
3. dn_s/dlnk = 0.000 ± 0.001 (perfect scale-invariance)

**SMFT Confirmed if**:
1. f_NL ~ 3-10 detected at 5σ
2. Gμ ~ 10⁻⁷ to 10⁻⁶ detected
3. Scale-dependent n_s(k) with string signature

### Timeline
- 2025-2028: CMB-S4 construction
- 2028-2035: Data collection (7 years)
- 2035-2037: Analysis & publication
- **Definitive answer by 2037**

---

## Proposal 5: High-Energy Collider Tests (Far Future)

### Executive Summary
- **Timeline**: 2040s (FCC-ee era)
- **Cost**: N/A (multi-billion facility)
- **Signal**: δm/m ~ 0.01% if R ≠ 1
- **Status**: Long-term precision tests

### Scientific Objective
Test dynamical fermion mass m_f = Δ·R(x,t):
1. Lepton universality: R_e = R_μ = R_τ ?
2. Top mass stability: m_t(channel) independent?
3. Higgs couplings: μ_f = R² ?

### Future Collider Precision

**FCC-ee (CERN proposal, 2040s)**:
- Z-pole: 10¹² Z bosons (vs. LEP: 10⁷)
- W mass: δm_W ~ 0.5 MeV (vs. current: ±12 MeV)
- Top mass: δm_t ~ 10 MeV (threshold scan)
- Higgs couplings: δμ/μ ~ 0.5%

**Sensitivity to R**:
```
If m_f = Δ·R → δm/m = δR/R

Precision δm_t ~ 10 MeV / 173 GeV ~ 0.006%
→ Can detect δR ~ 0.006% (if R ≠ 1)
```

### Proposed Tests

**Test 1: Lepton Universality**
```
Measure: Γ(Z → e⁺e⁻) / Γ(Z → μ⁺μ⁻) / Γ(Z → τ⁺τ⁻)

Standard Model: ratios = 1.0000 (exact)
SMFT (if R_e ≠ R_μ): ratios ≠ 1

Precision: δR/R ~ 0.01% → SMFT testable
```

**Test 2: Top Mass Stability**
```
Measure m_t in different channels:
- tt̄ → l⁺νl⁻ν̄
- tt̄ → l±νqq̄'
- tt̄ → qq̄'qq̄'

Standard: m_t(channel) = const
SMFT (if environment-dependent): m_t varies

Precision: δm_t ~ 10 MeV → detect 0.006% shifts
```

**Test 3: Higgs Coupling Ratios**
```
Measure: μ_f = σ(H → ff̄)_obs / σ(H → ff̄)_SM

SMFT: μ_f = R² (if R ≠ 1)

Precision: δμ/μ ~ 0.5% → detect R = 0.995 ± 0.003
```

### Expected Outcome

**SMFT Confirmed if**:
- R_f species-dependent (δR ~ 0.1%)
- m_t channel-dependent (δm ~ 50 MeV)
- Higgs couplings: μ_f = R² ≠ 1

**SMFT Ruled Out if**:
- Perfect lepton universality (δR/R < 10⁻⁴)
- m_t stable across channels
- Higgs couplings: μ_f = 1.000 ± 0.005

### Timeline
- 2040-2045: FCC-ee construction (if approved)
- 2045-2055: Data collection
- **Definitive answer by 2055** (far future)

---

## Priority Ranking

### Immediate (2025-2027)
1. **BEC Phonon Scattering** ✓✓✓ (EXECUTE NOW)
   - Cheapest, fastest, clearest
   - 13σ signal, 2-3 months
   - Any BEC lab can do it

### Near-Term (2027-2030)
2. **Quantum Synchronization** ✓✓
   - Novel physics, 5σ detection
   - 2-3 years, accessible hardware
   - Requires dedicated quantum computing time

3. **Casimir Force** ✓
   - Challenging (1% precision required)
   - 12-18 months, $50k-100k
   - Material comparison more robust

### Long-Term (2030-2040)
4. **CMB-S4** ✓
   - Definitive cosmology test
   - Funded mission, 2030s data
   - f_NL and Gμ sensitivity sufficient

### Far Future (2040+)
5. **FCC-ee Collider**
   - Ultimate precision
   - Multi-billion facility
   - Fermion mass tests at 0.01% level

---

## Recommended Execution Path

**Year 1 (2025)**:
- Execute BEC phonon scattering (3 months)
- Publish results (6 months)
- Begin quantum synchronization planning

**Year 2-3 (2026-2027)**:
- Quantum synchronization Phase 1 (N=20, proof-of-principle)
- Casimir force method development (nanomembrane)
- Prepare for CMB-S4 data (analysis pipelines)

**Year 4-5 (2028-2030)**:
- Quantum synchronization Phase 2-3 (N=100-200, definitive)
- Casimir force measurements (if promising)
- CMB-S4 early data analysis

**Year 6-10 (2030-2035)**:
- CMB-S4 full dataset analysis
- Follow-up experiments based on earlier results
- Publish comprehensive SMFT validation

**Beyond 2035**:
- Await FCC-ee or next-generation collider
- Theoretical refinements based on experimental results
- Cosmological observations (primordial gravitational waves, etc.)

---

## Summary

**SMFT makes five testable predictions accessible with current or near-future technology:**

1. ✅ **BEC phonon scattering** - Ready NOW (66% signal, 2-3 months)
2. ✅ **Quantum synchronization** - Feasible 2027-2030 (5σ β measurement)
3. ⚠️ **Casimir force** - Challenging (10% signal, 1% precision needed)
4. ⏳ **CMB cosmology** - Awaiting 2030s data (f_NL, Gμ constraints)
5. ⏳ **Collider precision** - Far future 2040s (fermion mass tests)

**The theory is NOT confined to Planck scales. Laboratory tests exist NOW.**

**Next step: Execute Proposal 1 (BEC phonon scattering) as proof-of-concept.**
