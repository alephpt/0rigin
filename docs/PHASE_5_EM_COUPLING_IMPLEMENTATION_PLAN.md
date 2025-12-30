# Phase 5 (Scenario 2.6B): Electromagnetic Coupling Implementation Plan

**Date**: 2024-12-24  
**Status**: Planning Phase  
**Physics Goal**: Test if ∇θ (Kuramoto phase gradients) naturally couples to Dirac field as electromagnetic potential A_μ

---

## EXECUTIVE SUMMARY

**Recommended Approach**: Option B (Perturbative EM Coupling)  
**Timeline**: **2 days** (19 hours focused work)  
**Risk**: **30-40%** (medium)  
**LOC**: ~1350 total (790 new + 560 modifications)  
**Discovery Potential**: ⭐⭐⭐⭐⭐ MAJOR (EM from SMFT phase dynamics)

**Comparison to Sprint 2 estimate**:
- Original: 3-4 weeks, 20-30% risk (full lattice gauge theory)
- Updated: **2 days, 30-40% risk** (perturbative approach)
- **8× faster** with comparable risk via simpler physics

---

## 1. PHYSICS BACKGROUND

### Minimal Coupling to EM Field

**Standard Dirac equation**:
```
(iγ^μ ∂_μ - m)ψ = 0
```

**With electromagnetic field**:
```
(iγ^μ D_μ - m)ψ = 0
D_μ = ∂_μ - iq A_μ    (covariant derivative)
```

**A_μ from Kuramoto phase θ(x,y,t)**:
```
A_0(x,y,t) = ∂_t θ    (scalar potential)
A_x(x,y,t) = ∂_x θ    (vector potential x)
A_y(x,y,t) = ∂_y θ    (vector potential y)
```

**Electromagnetic fields**:
```
E = -∇A_0 - ∂_t A     (electric field)
B = ∇ × A             (magnetic field)
```

**CRITICAL**: For smooth θ, F_μν = 0 (mixed partials commute).  
**RESOLUTION**: EM lives at topological defects (vortex cores).

---

## 2. IMPLEMENTATION OPTIONS

### Option A: Full Lattice Gauge Theory (NOT RECOMMENDED)
- Timeline: 3-4 weeks
- LOC: ~1250
- Risk: 20-30%
- Pros: Rigorous, preserves gauge invariance
- Cons: Too slow for "NOW" requirement

### Option B: Perturbative EM Coupling (RECOMMENDED ✅)
- Timeline: **2 days**
- LOC: ~1350
- Risk: 30-40%
- Pros: Fast, sufficient for exploration, non-invasive
- Cons: Weak field only (q·A << 1), approximate gauge invariance

**Recommendation: Option B** meets user requirement for immediate implementation.

---

## 3. TECHNICAL APPROACH (Option B)

### Modified Hamiltonian
```
H = α·∇ + β·m(x) + q·(A_0 + α·A)
         ↑          ↑
      existing    NEW EM term
```

### Evolution Sequence
```cpp
// Current (no EM):
applyKineticHalfStep(dt/2);
applyPotentialStep(mass, dt);
applyKineticHalfStep(dt/2);

// With EM coupling:
applyKineticHalfStep(dt/2);
applyPotentialStep(mass, dt);        // β·m term
applyEMPotentialStep(A_mu, q, dt);   // NEW: q·A_μ term
applyKineticHalfStep(dt/2);
```

### EM Potential Application
```cpp
void DiracEvolution::applyEMPotentialStep(
    const std::vector<double>& A_0,    // Scalar potential
    const std::vector<double>& A_x,    // Vector potential x
    const std::vector<double>& A_y,    // Vector potential y
    double q_charge,                   // Electric charge
    double dt) {
    
    // Scalar potential: exp(-iq·A_0·dt) → phase shift all components
    // Vector potential: exp(-iq·α·A·dt) → coupling via α matrices
}
```

---

## 4. NEW FILES REQUIRED

### 4.1 EMFieldComputer.h (~80 LOC)
```cpp
#pragma once
#include <vector>
#include <tuple>

class EMFieldComputer {
public:
    // Extract A_μ = (A_0, A_x, A_y) from phase field θ
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    computeVectorPotential(
        const std::vector<float>& theta_field,
        const std::vector<float>& theta_prev,
        int Nx, int Ny, double dt);
    
    // Compute E and B fields from A_μ
    static std::tuple<std::vector<double>, std::vector<double>, double>
    computeEMFields(
        const std::vector<double>& A_0,
        const std::vector<double>& A_x,
        const std::vector<double>& A_y,
        int Nx, int Ny, double dx, double dt);
};
```

### 4.2 EMFieldComputer.cpp (~300 LOC)
- Finite difference gradients for A_μ = ∂_μ θ
- E field: E = -∇A_0 - ∂_t A
- B field: B = ∇×A = ∂_x A_y - ∂_y A_x
- Periodic boundary conditions

### 4.3 EMObservables.h (~60 LOC)
```cpp
#pragma once
#include <vector>

struct MaxwellValidation {
    double div_E;           // ∇·E (should = ρ)
    double curl_B_minus_J;  // ∇×B - J (should = ∂_t E)
    double div_B;           // ∇·B (should = 0)
    bool gauss_law_valid;
    bool ampere_law_valid;
    bool no_monopoles;
};

class EMObservables {
public:
    // Compute Lorentz force F = q(E + v×B)
    static std::tuple<double, double>
    computeLorentzForce(...);
    
    // Validate Maxwell equations
    static MaxwellValidation
    validateMaxwellEquations(...);
    
    // EM energy density u = (E² + B²)/2
    static double
    computeEMEnergyDensity(...);
};
```

### 4.4 EMObservables.cpp (~350 LOC)
- Lorentz force from E, B, velocity
- Maxwell divergence/curl computations
- Energy density integration

**Total new files: 4, ~790 LOC**

---

## 5. MODIFIED FILES

### 5.1 DiracEvolution.h (+20 LOC)
```cpp
class DiracEvolution {
    // ... existing methods ...
    
    // NEW: Apply EM potential step
    void applyEMPotentialStep(
        const std::vector<double>& A_0,
        const std::vector<double>& A_x,
        const std::vector<double>& A_y,
        double q_charge, double dt);
    
    // NEW: Enable/disable EM coupling
    void setEMCoupling(bool enable, double q_charge);
    
private:
    bool _em_coupling_enabled;
    double _q_charge;
};
```

### 5.2 DiracEvolution.cpp (+100 LOC)
- Implement applyEMPotentialStep()
- Scalar potential: phase shift exp(-iq A_0 dt)
- Vector potential: α matrix coupling exp(-iq α·A dt)
- Modify step() to call when enabled

### 5.3 ObservableComputer.h (+50 LOC)
```cpp
struct Observables {
    // ... existing fields ...
    
    // NEW: EM observables
    double A_0_avg, A_x_avg, A_y_avg;
    double E_field_mag, B_field_mag;
    double em_energy_density;
    double lorentz_force_x, lorentz_force_y;
    MaxwellValidation maxwell;
};
```

### 5.4 ObservableComputer.cpp (+150 LOC)
- Compute EM observables each timestep
- Integration with existing compute() method
- CSV output formatting

### 5.5 SMFTTestRunner.cpp (+200 LOC)
```cpp
// In main evolution loop:
if (config.em_coupling_enabled) {
    // Get theta field from SMFTEngine
    auto theta_curr = engine.getPhaseField();
    
    // Compute A_μ = ∂_μ θ
    auto [A_0, A_x, A_y] = EMFieldComputer::computeVectorPotential(
        theta_curr, theta_prev, Nx, Ny, dt);
    
    // Apply to Dirac evolution
    dirac.applyEMPotentialStep(A_0, A_x, A_y, q_charge, dt);
    
    // Compute EM observables
    auto [E, B, u_em] = EMFieldComputer::computeEMFields(...);
    // ... collect for output
}
```

### 5.6 TestConfig.h (+10 LOC)
```cpp
struct TestConfig {
    // ... existing fields ...
    
    // NEW: EM coupling
    double em_coupling_strength;  // q (charge)
    bool test_A_mu_coupling;
};
```

### 5.7 TestConfig.cpp (+30 LOC)
- Parse em_coupling_strength from YAML
- Parse EM analysis flags

**Total modifications: 7 files, ~560 LOC**

---

## 6. IMPLEMENTATION TIMELINE

### Day 1: Core Implementation (8 hours)

**Morning (4h)**:
1. Create EMFieldComputer.{h,cpp} (1.5h)
2. Implement computeVectorPotential() (1h)
3. Implement computeEMFields() (1h)
4. Unit tests for gradients (0.5h)

**Afternoon (4h)**:
5. Add applyEMPotentialStep() to DiracEvolution (2h)
6. Integrate into step() method (1h)
7. Test with constant A_μ field (0.5h)
8. Verify phase shifts correct (0.5h)

### Day 2: Integration & Testing (8 hours)

**Morning (4h)**:
9. Create EMObservables.{h,cpp} (1.5h)
10. Implement computeLorentzForce() (1h)
11. Implement validateMaxwellEquations() (1h)
12. Extend ObservableComputer (0.5h)

**Afternoon (4h)**:
13. Integrate into SMFTTestRunner (2h)
14. Parse EM config parameters (0.5h)
15. Run scenario_2.6B test (0.5h)
16. Analyze results, debug (1h)

### Day 3: Validation (3 hours if needed)
17. Test gauge invariance (θ → θ + const) (1h)
18. Energy conservation check (1h)
19. Documentation (1h)

**Total: 19 hours over 2-3 days**

---

## 7. RISK ANALYSIS

### Technical Risks

**LOW (10-20%)**:
✅ EM field computation (straightforward gradients)
✅ Observable extensions (established pattern)
✅ Config parsing (mature system)

**MEDIUM (30-40%)**:
⚠️ Phase shift implementation in applyEMPotentialStep()
- Spinor structure with α matrices
- Sign errors possible
- **Mitigation**: Test with analytic solutions first

**HIGH (50%)**:
⚠️⚠️ Gauge invariance (approximate only in Option B)
- Results may depend on θ zero point
- **Mitigation**: Explicitly test θ → θ + const

### Physics Risks

**Discovery outcomes**:
- ⭐⭐⭐⭐⭐ **Extraordinary**: α ~ 1/137 emerges naturally
- ⭐⭐⭐⭐ **Success**: Lorentz force works, Maxwell satisfied
- ⭐⭐⭐ **Partial**: Coupling works but α wrong (needs renormalization)
- ⭐ **Failure**: No EM-like behavior (rules out hypothesis)

**Expected**: Partial success (still publishable)

---

## 8. SUCCESS CRITERIA

### Minimum Viable (Day 2 end):
✅ A_μ field extracted from θ
✅ E and B fields computed
✅ EM coupling modifies Dirac evolution
✅ Lorentz force measurement works
✅ Maxwell validation runs
✅ scenario_2.6B completes

### Full Success:
✅ All minimum criteria
✅ Flux quantization Φ = 2πn around vortices
✅ Lorentz force deflects charged particles
✅ Neutral particles (q=0) unaffected
✅ Maxwell equations approximately satisfied
✅ Gauge invariance tested
✅ EM energy contributes to total

### Falsification (Rules OUT):
❌ Flux not quantized
❌ No deflection with q ≠ 0
❌ Strong θ zero-point dependence
❌ Energy not conserved
❌ Maxwell violations >10×

---

## 9. DELIVERABLES

1. **Source code**:
   - 4 new files (~790 LOC)
   - 7 modified files (~560 LOC)

2. **Test results**:
   - output/scenario_2.6B_electromagnetic_coupling/
   - A_μ, E, B field snapshots
   - Lorentz force data
   - Maxwell validation report

3. **Documentation**:
   - This implementation plan
   - Usage guide
   - Physics interpretation

---

## 10. COMPARISON TO SPRINT 2

| Aspect | Sprint 2 Original | This Plan (Option B) |
|--------|------------------|----------------------|
| Timeline | 3-4 weeks | **2 days** |
| LOC | ~1250 | ~1350 |
| Risk | 20-30% | 30-40% |
| Approach | Full lattice gauge | Perturbative |
| Gauge invariance | Exact | Approximate |
| Speed improvement | 1× | **8×** |

**Can we do better?** ✅ **YES**
- 8× faster via simpler physics
- Comparable risk
- Sufficient for discovery
- Can upgrade later if needed

---

## 11. RECOMMENDATION

✅ **PROCEED with Option B (Perturbative EM Coupling)**

**Rationale**:
1. Meets user "NOW" requirement (2 days vs 3-4 weeks)
2. Sufficient for exploratory physics
3. Non-invasive to validated code
4. Lower complexity → fewer bugs
5. Can upgrade to Option A if results warrant

**Ready to implement upon user approval.**

---

**END OF IMPLEMENTATION PLAN**
