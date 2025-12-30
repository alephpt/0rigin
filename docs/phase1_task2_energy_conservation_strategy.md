# Phase 1 Task 2: Energy Conservation Strategy
## Distinguishing Numerical Error from Physical Damping

### Executive Summary
Current SMFT simulations exhibit ~1% energy drift over long timescales. This analysis develops a comprehensive strategy to determine whether this drift represents numerical error, physical damping, or genuine conservation with measurement artifacts. The analysis reveals that **Strang splitting is already implemented** in the codebase, providing second-order symplectic integration that should preserve energy to machine precision in the absence of dissipation.

---

## 1. Current Implementation Analysis

### 1.1 Strang Splitting Implementation
The codebase already implements **Strang splitting** (second-order operator splitting) in `DiracEvolution::evolve()`:

```cpp
// src/DiracEvolution.cpp (lines 203-206)
// Strang splitting: K/2 - V - K/2
applyKineticHalfStep(dt_effective / 2.0f);
applyPotentialStep(mass_field, dt_effective);
applyKineticHalfStep(dt_effective / 2.0f);
```

This is superior to first-order Euler or even RK4 for Hamiltonian systems because:
- **Symplectic structure preservation**: Conserves phase space volume
- **Second-order accuracy**: O(dt²) global error
- **Energy conservation**: Bounded energy error (oscillates but doesn't drift)

### 1.2 Energy Computation Architecture

The energy is computed in `ObservableComputer` with three components:

1. **Kinetic Energy (T)**: From Dirac operator -iα·∇
   - Uses FFT for momentum-space computation in `DiracEvolution::getEnergy()`
   - Formula: `T = Σ_k |Ψ̃(k)|² ω(k)` where `ω(k) = |k|` (massless limit)

2. **Potential Energy (V)**: From mass coupling β·m(x)
   - Formula: `V = ∫ Ψ†(β·m(x))Ψ dx` where `m(x) = Δ·R(x)`
   - Upper spinor components: +m contribution
   - Lower spinor components: -m contribution

3. **Kuramoto Energy**: Currently **NOT included** in energy budget
   - Phase gradients: `∫|∇θ|² dx`
   - Synchronization potential: Function of R-field

### 1.3 Damping Parameter Usage

The damping parameter `γ` appears in multiple locations:

1. **Kuramoto Evolution** (`SMFTEngine::step()`):
   - Applied to phase dynamics: `dθ/dt = ω + K·Σsin(θⱼ - θᵢ) - γ·θ`
   - Direct energy dissipation mechanism

2. **Coupled Evolution** (`SMFTEngine::stepWithDirac()`):
   - Uses Strang splitting with damping in Kuramoto substeps
   - Default: `damping = 0.1` (10% dissipation rate)

3. **Test Configurations**:
   - Most tests use `damping: 0.1` (non-zero)
   - `phase_2.5A_zero_damping.yaml` specifically tests `damping: 0.0`

---

## 2. Theoretical Framework for Energy Budget

### 2.1 Sources of Energy Non-Conservation

#### A. Physical Dissipation
- **Kuramoto damping**: `-γ·θ` term directly removes energy
- **Expected signature**: Exponential decay `E(t) = E₀·exp(-γt)`
- **Test**: Run with `γ = 0` should eliminate this contribution

#### B. Numerical Truncation Error
- **Strang splitting error**: O(dt³) per step, O(dt²) global
- **Grid discretization**: O(dx²) for second-order finite differences
- **FFT aliasing**: High-frequency mode coupling
- **Expected signature**: Systematic drift proportional to dt²

#### C. Operator Splitting Error
- **Timescale separation**: Fast Kuramoto vs slow Dirac dynamics
- **Adiabatic approximation**: Assumes R-field quasi-static during Dirac evolution
- **Expected signature**: Oscillations at frequency ~1/N where N = substep ratio

#### D. Missing Energy Components
- **Kuramoto field energy**: Not currently tracked
- **EM field energy**: Computed but not always included
- **Interaction energy**: Kuramoto-Dirac coupling terms

### 2.2 Energy Conservation Equation

Complete energy budget should satisfy:

```
dE_total/dt = dE_Dirac/dt + dE_Kuramoto/dt + dE_EM/dt + dE_interaction/dt = -P_dissipated
```

Where:
- `E_Dirac`: Dirac kinetic + potential energy
- `E_Kuramoto`: Phase gradient + synchronization energy
- `E_EM`: Electromagnetic field energy (if enabled)
- `E_interaction`: Coupling terms between fields
- `P_dissipated`: Power lost to damping = γ·∫|∇θ|² dx

---

## 3. Detailed Verification Plan

### 3.1 Test Matrix

| Test ID | γ | dt | dx | N | Duration | Expected Result |
|---------|---|-----|-----|---|----------|-----------------|
| T1 | 0.0 | 0.01 | L/128 | 10 | 10⁴ steps | Energy conserved to ~10⁻⁶ |
| T2 | 0.1 | 0.01 | L/128 | 10 | 10⁴ steps | Exponential decay E ~ E₀e⁻⁰·¹ᵗ |
| T3 | 0.0 | 0.005 | L/128 | 10 | 2×10⁴ steps | 4× better conservation (dt²) |
| T4 | 0.0 | 0.01 | L/256 | 10 | 10⁴ steps | 4× better spatial accuracy |
| T5 | 0.0 | 0.01 | L/128 | 100 | 10⁴ steps | Reduced splitting error |
| T6 | 0.0 | 0.01 | L/128 | 10 | 10⁵ steps | Long-time stability test |

### 3.2 Energy Decomposition Analysis

Implement detailed energy tracking:

```cpp
struct EnergyComponents {
    double dirac_kinetic;      // T_Dirac = Σ_k |Ψ̃(k)|² ω(k)
    double dirac_potential;    // V_Dirac = ∫ Ψ†(β·m)Ψ dx
    double kuramoto_gradient;  // T_Kuramoto = ∫|∇θ|² dx
    double kuramoto_sync;      // V_Kuramoto = f(R)
    double em_field;          // E_EM = ∫(E² + B²)/(8π) dx
    double coupling_interaction; // E_int = λ·∫|Ψ|²·θ dx
    double total;             // Sum of all components
};
```

### 3.3 Code Modifications Required

1. **Add Kuramoto energy tracking** in `ObservableComputer`:
   ```cpp
   double computeKuramotoEnergy(
       const std::vector<float>& theta,
       const std::vector<float>& R_field,
       float K);
   ```

2. **Track energy components separately**:
   - Modify `Observables` struct to include all energy terms
   - Update CSV output to log individual components

3. **Implement energy flux diagnostics**:
   ```cpp
   double computeDissipationRate(
       const std::vector<float>& theta,
       float damping);
   ```

4. **Add Strang vs Euler comparison mode**:
   - Implement first-order Euler stepping as baseline
   - Compare energy conservation between methods

### 3.4 Analysis Methods

1. **Energy Conservation Metric**:
   ```
   ε(t) = |E(t) - E₀| / E₀
   ```
   Plot log(ε) vs log(t) to identify power-law behavior

2. **Dissipation Rate Analysis**:
   ```
   P(t) = -dE/dt ≈ γ·∫|∇θ|² dx
   ```
   Compare measured vs predicted dissipation

3. **Spectral Analysis**:
   - FFT of E(t) to identify oscillation frequencies
   - Should see peaks at operator splitting frequency

4. **Richardson Extrapolation**:
   - Run with dt, dt/2, dt/4
   - Extrapolate to dt→0 limit
   - True conservation if E(dt→0) = const

---

## 4. Success Criteria

### 4.1 True Energy Conservation
Evidence required:
- With γ=0: |ΔE/E₀| < 10⁻⁶ over 10⁵ timesteps
- Energy error oscillates but doesn't grow (Strang signature)
- Richardson extrapolation converges to constant
- All energy components tracked and summing correctly

### 4.2 Numerical Artifact
Evidence:
- Energy drift ∝ dt² (truncation error)
- Energy drift ∝ 1/N (operator splitting error)
- Grid refinement reduces drift systematically
- Missing energy component explains drift

### 4.3 Physical Dissipation
Evidence:
- Energy decay matches exponential E(t) = E₀·exp(-γt)
- Dissipation rate P = γ·∫|∇θ|² dx verified
- Zero damping eliminates decay completely
- Decay independent of numerical parameters (dt, dx, N)

---

## 5. Expected Results and Acceptance Criteria

### 5.1 With Strang Splitting (Current Implementation)

**Expected behavior with γ=0**:
- Energy oscillates with amplitude ~O(dt²)
- No systematic drift over long times
- Bounded error: |E(t) - E₀| < C·dt² for all t

**Acceptance criterion**:
- Energy conserved to within 0.01% over 10⁴ timesteps with γ=0
- This confirms numerical scheme is fundamentally conservative

### 5.2 Damping Contribution

**Expected behavior with γ=0.1**:
- Initial exponential decay: E(t) ≈ E₀·(1 - γt) for small t
- Asymptotic decay rate matches theory
- Decay rate independent of dt, dx, N

**Acceptance criterion**:
- Measured dissipation within 5% of theoretical prediction
- Clear separation between numerical and physical effects

### 5.3 Complete Energy Budget

**Expected with all components tracked**:
```
d(E_total)/dt = -P_damping ± ε_numerical
```
Where `ε_numerical < 10⁻⁶ E₀` for proper parameters

**Acceptance criterion**:
- Energy budget closes to within 0.1%
- All energy transfers between components accounted for

---

## 6. Implementation Priority

### Phase 1: Immediate Tests (No code changes)
1. Run existing `phase_2.5A_zero_damping.yaml` configuration
2. Compare with damping=0.1 baseline
3. Analyze existing energy data for conservation

### Phase 2: Enhanced Diagnostics (Minimal changes)
1. Add Kuramoto energy computation
2. Track energy components separately
3. Output detailed energy budget

### Phase 3: Comprehensive Validation (Full implementation)
1. Implement all test configurations T1-T6
2. Add Euler baseline for comparison
3. Complete Richardson extrapolation analysis
4. Generate publication-quality validation plots

---

## 7. Conclusion

The presence of Strang splitting in the current implementation is highly significant. This second-order symplectic integrator should provide excellent energy conservation in the absence of dissipation. The observed 1% drift likely comes from:

1. **Physical damping** (γ=0.1 in most tests)
2. **Missing Kuramoto energy** in the budget
3. **Operator splitting approximation** between Kuramoto and Dirac

The verification strategy outlined here will definitively identify the source and allow us to claim either:
- **"True energy conservation"** (with γ=0 and complete budget)
- **"Controlled dissipation"** (with quantified damping effects)

The fact that Strang splitting is already implemented gives high confidence that the numerical scheme is fundamentally energy-conserving, and any observed drift is either physical or due to incomplete energy accounting.