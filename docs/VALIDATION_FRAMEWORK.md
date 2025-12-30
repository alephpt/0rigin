# SMFT Validation Framework

## I. Global Requirements (ALL Simulations)

These are fundamental physics constraints that MUST hold for any SMFT simulation:

### 1. Probability Conservation
```
d/dt ∫|ψ|² dx dy = 0
```
- **Tolerance**: |Δ||ψ||²| < 0.5% over simulation
- **Violation indicates**: Numerical instability, non-unitary evolution, boundary issues
- **Action if violated**: Simulation INVALID - discard immediately

### 2. Energy Conservation (Isolated System)
```
dE_total/dt = 0
E_total = E_Dirac + E_Kuramoto + E_coupling
```
- **Tolerance**: |ΔE/E₀| < 1% over simulation
- **Violation indicates**: Time-step too large, operator splitting errors
- **Action if violated**: Results unreliable - discard

### 3. Order Parameter Bounds
```
0 ≤ R(x,y,t) ≤ 1  ∀ x, y, t
```
- **Physical meaning**: R is synchronization measure
- **Violation indicates**: Kuramoto update error, overflow, unphysical configuration
- **Action if violated**: CRITICAL FAILURE - immediate stop

### 4. Gauge Invariance
```
⟨ψ|O|ψ⟩ independent of θ → θ + α(t)
```
- **Violation indicates**: Coupling not gauge-invariant, measurement sees unphysical modes
- **Action if violated**: Theory wrong, not just implementation

### 5. Causality
```
|dx/dt| ≤ c (c=1 in Planck units)
```
- **Violation indicates**: Superluminal propagation, unphysical initial conditions
- **Action if violated**: Relativistic framework broken

### 6. Numerical Stability
```
All fields finite and non-divergent
```
- **Action if violated**: Reduce dt, check CFL condition

---

## II. Test-Specific Requirements

### Scenario 2.1: Defect Localization
**Purpose**: Validate vortex core structure in Kuramoto field

**Required**:
- ✓ Vortex present: W = ±1
- ✓ R-field core: R_min < 0.5
- ✓ Core position stable: |dx_core/dt| < 0.01

**NOT required**:
- ✗ Particle present (vacuum only)
- ✗ Boosted conditions
- ✗ Specific γ factor

---

### Scenario 2.2: Traveling Wave Surfing
**Purpose**: Validate particle-defect coupling (vortex + boosted particle)

**Required**:
- ✓ θ(x,y,t=0) vortex: W = ±1
- ✓ R(x,y,t=0) core: R_min < 0.5
- ✓ ψ(x,y,t=0) boosted Gaussian at offset
- ✓ ⟨p⟩(t=0) = γmv ± 5%
- ✓ Particle tracks vortex motion
- ✓ γ_measured within 5% of theory (if measuring mass)

**NOT required**:
- ✗ Multiple vortices
- ✗ Grid convergence (single config test)

---

### Scenario 2.3: Relativistic Mass Validation
**Purpose**: Validate m(v) = γ·Δ·R across velocities and grids

**Required**:
- ✓ θ(x,y,t=0) vortex: W = ±1
- ✓ R(x,y,t=0) core: R_min < 0.5
- ✓ ψ(x,y,t=0) boosted Gaussian: ⟨p⟩ = γmv ± 5%
- ✓ γ_measured within 5% of theory
- ✓ Grid convergence demonstrated
- ✓ Operator splitting convergence (N=1, 10, 100)

**NOT required**:
- ✗ Particle-vortex interaction (Gaussian can be anywhere)
- ✗ Long-time evolution (testing instantaneous mass)

---

### Scenario 2.4A-C: Breakdown Investigation
**Purpose**: Understand WHY Scenario 2.3 failed at high velocities

**Required** (same as 2.3):
- ✓ All Scenario 2.3 requirements
- ✓ Detailed diagnostics (R-field evolution, spatial snapshots)

**Additional for 2.4B**:
- ✓ R-field tracking particle position
- ✓ Mass modulation observable

---

### Future: EM Coupling Tests
**Purpose**: Validate A_μ = ∂_μθ emergent electromagnetism

**Required**:
- ✓ Magnetic flux quantization: ∮A·dl = 2πn
- ✓ Lorentz force: F = q(E + v×B)
- ✓ Gauge invariance: A → A + ∇χ

**NOT required**:
- ✗ Vortex present (could test uniform A)
- ✗ Relativistic boost
- ✗ Specific mass values

---

## III. Implementation Structure

### File Organization
```
src/validation/
├── global_checks.cpp       # Always run (6 global requirements)
├── defect_checks.cpp       # Scenario 2.1 (vortex structure)
├── traveling_wave_checks.cpp  # Scenario 2.2 (particle-defect)
├── relativistic_checks.cpp # Scenario 2.3/2.4 (mass validation)
├── em_checks.cpp          # Future EM tests
└── validation.h           # Common interface
```

### Configuration Integration
Add to YAML configs:
```yaml
validation:
  # GLOBAL (always enforced)
  norm_tolerance: 0.005        # 0.5%
  energy_tolerance: 0.01       # 1%
  R_bounds_check: true         # 0 ≤ R ≤ 1
  causality_check: true        # v ≤ c

  # TEST-SPECIFIC
  scenario: "relativistic_mass"  # Determines which checks to run

  # Relativistic mass specific
  gamma_tolerance: 0.05        # 5% for γ factor
  require_vortex: true         # W = ±1
  require_core: true           # R_min < 0.5
  require_boost: true          # ⟨p⟩ = γmv
  convergence_tolerance: 0.05  # Grid/N convergence
```

---

## IV. Validation Protocol

### Pre-Flight Check (Before Evolution)
```python
def global_preflight_check(config, initial_state):
    """Run before EVERY simulation"""

    # 1. Probability normalization
    norm = integrate(|initial_state.psi|²) * dx * dy
    assert abs(norm - 1.0) < 1e-6

    # 2. Order parameter bounds
    R = initial_state.R
    assert all(0 <= R <= 1)

    # 3. Causality (if boosted)
    if config.has_boost:
        v = compute_velocity(initial_state.psi)
        assert ||v|| < 1.0  # c = 1

    return PASS
```

### Scenario-Specific Initial Validation
```python
def validate_relativistic_mass_initial(config, state):
    """Only for Scenario 2.3/2.4"""

    # 1. Vortex structure
    W = compute_winding_number(state.theta)
    assert |W - 1| < 0.2

    # 2. Defect core
    R_min = min(state.R)
    assert R_min < 0.5

    # 3. Boosted momentum
    p0 = compute_momentum(state.psi)
    p_expected = config.gamma * config.mass * config.velocity
    assert |p0 - p_expected| < 0.05 * p_expected

    return PASS
```

### Runtime Monitoring (During Evolution)
```python
def monitor_global_constraints(state, t):
    """Check every N steps"""

    # 1. Probability conservation
    norm_drift = |integrate(|state.psi|²) - 1.0|
    if norm_drift > 0.005:
        WARNING("Norm drift: {norm_drift}")

    # 2. Energy conservation
    E_drift = |state.E_total - state.E_initial| / state.E_initial
    if E_drift > 0.01:
        WARNING("Energy drift: {E_drift}")

    # 3. R bounds
    if any(state.R < 0 or state.R > 1):
        CRITICAL_FAILURE("R out of bounds")
```

### Final Validation
```python
def validate_relativistic_mass_final(config, final_state):
    """Only for Scenario 2.3/2.4"""

    # Extract measured γ factor
    gamma_measured = extract_gamma_from_m_eff(final_state)
    gamma_theory = config.gamma_theory

    # Check within tolerance
    error = |gamma_measured - gamma_theory| / gamma_theory
    assert error < config.gamma_tolerance

    return PASS
```

---

## V. Why Scenario 2.4A/2.4B Failed

### They Failed Test-Specific Requirements:
```python
# Scenario 2.4B would have caught:
R_min = state.R.min()
assert R_min < 0.5  # ❌ FAIL: R_min = 0.9998 (no core!)

p0 = compute_momentum(state.psi)
assert p0 > 0.3  # ❌ FAIL: p0 = 0.0004 (not boosted!)
```

### But Passed Global Checks:
```python
norm = integrate(|psi|²)      # ✓ PASS: norm ≈ 1.000
energy_drift = |ΔE|/E         # ✓ PASS: drift < 1%
R_bounds = all(0 <= R <= 1)   # ✓ PASS: R ∈ [0.995, 1.000]
```

**This is exactly why we need both global AND test-specific validation.**

---

## VI. Recommended Action Plan

### 1. Implement Global Validation (HIGH PRIORITY)
- Add to ObservableComputer.cpp
- Check at initialization, every 100 steps, and finalization
- Auto-fail simulation if violated

### 2. Implement Scenario-Specific Validation
- Create validation module per scenario type
- Load based on `validation.scenario` in YAML
- Report violations clearly in output

### 3. Fix Phase 2.2 Before Proceeding
Run new validation framework on Phase 2.2:
```bash
./build/bin/smft --test config/traveling_wave_N10.yaml --validate
```

Expected failures:
- ❌ No vortex (W ≈ 0)
- ❌ No core (R_min > 0.95)
- ❌ No boost (p ≈ 0)

Action: **DEBUG until all 6 criteria pass**, then proceed to 2.3

### 4. Re-run ALL Phase 2 Tests
Only after Phase 2.2 validation passes.

---

## VII. Success Criteria

### Global (Always)
✓ ||ψ||² = 1 ± 0.5%
✓ |ΔE/E| < 1%
✓ R ∈ [0,1]
✓ v ≤ c
✓ Fields finite

### Phase 2.2 (Traveling Wave)
✓ W = ±1
✓ R_min < 0.5
✓ p(t=0) = 0.314 ± 5%
✓ γ_measured = 1.048 ± 5%

### Phase 2.3 (Relativistic Mass)
✓ All Phase 2.2 requirements
✓ Grid convergence (256×256 < 2% error)
✓ N convergence (N=10 optimal)

---

**Status**: Framework designed, needs implementation in C++ engine

**Next**: Implement global_checks.cpp and integrate with SMFTEngine
