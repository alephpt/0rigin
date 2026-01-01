# Option 2: Dynamic K-Parameter Design Document

**Date**: 2026-01-01
**Status**: Design phase - awaiting parameter scan confirmation
**Goal**: Achieve B_max>1.0 AND R_avg>0.9 simultaneously

---

## Theoretical Foundation

### Problem Statement
Static K>0 destroys vortex gradients via exponential relaxation:
```
∇θ(t) ~ ∇θ(0)·exp(-K·R·t)
B_z ∝ ∇²θ → 0
```

### Proposed Solution
K = K(x,t) adapts spatially/temporally to preserve gradients while enabling synchronization.

---

## Design Options

### Option 2A: Gradient-Suppressed Coupling
**Mechanism**: Reduce K near strong phase gradients

```cpp
K(x,y,t) = K_max · [1 - tanh(|∇θ|²/G_crit)]
```

**Parameters**:
- K_max = 1.0 (maximum coupling in uniform regions)
- G_crit = critical gradient strength where K suppression activates

**Behavior**:
- |∇θ| small → K ≈ K_max (synchronize)
- |∇θ| large → K ≈ 0 (preserve vortex)

**Physical interpretation**: Strong gradients indicate topological structure that should be preserved.

---

### Option 2B: Order-Parameter Modulated Coupling
**Mechanism**: Reduce K in desynchronized regions

```cpp
K(x,y,t) = K_max · R(x,y,t)
```

**Rationale**: Where R is low, system already desynchronized → reducing K preserves local structure.

**Problem**: Feedback loop - low R → low K → can't synchronize → R stays low.

**Modified version**:
```cpp
K(x,y,t) = K_min + (K_max - K_min) · R²(x,y,t)
```

Ensures minimum K_min > 0 to seed synchronization.

---

### Option 2C: Vortex-Core Protection (Recommended)
**Mechanism**: Detect vortex cores and suppress K locally

```cpp
// Compute vorticity
ω_z(x,y) = ∂_xθ_y - ∂_yθ_x  (vorticity)

// Detect cores (high vorticity)
is_core(x,y) = |ω_z| > ω_threshold

// Suppress coupling
K(x,y,t) = K_max · [1 - α·is_core(x,y)]
```

**Parameters**:
- α ∈ [0,1] = suppression strength
- ω_threshold = vorticity threshold for core detection

**Advantages**:
- Directly targets vortex structures
- Physical interpretation: Protect topological defects
- Literature support: Topological defects break synchronization locally

---

### Option 2D: Hybrid Gradient + Order Parameter
**Mechanism**: Combine both effects

```cpp
K(x,y,t) = K_max · R(x,y,t) · [1 - tanh(|∇θ|²/G_crit)]
```

**Behavior**:
- Synchronized + uniform → K = K_max
- Desynchronized + structured → K ≈ 0
- Intermediate states → smooth interpolation

---

## Implementation Strategy

### Phase 1: CPU Prototype
Add to `DiracEvolution.cpp` or new `AdaptiveKuramoto.cpp`:

```cpp
class AdaptiveKuramotoField {
public:
    AdaptiveKuramotoField(int Nx, int Ny, float K_max, float G_crit);

    // Update K-field based on current θ and R
    void updateKField(
        const std::vector<float>& theta,
        const std::vector<float>& R_field,
        std::vector<float>& K_field);

    // Different adaptation strategies
    enum class AdaptationMode {
        GRADIENT_SUPPRESSED,    // Option 2A
        ORDER_MODULATED,        // Option 2B
        VORTEX_PROTECTED,       // Option 2C
        HYBRID                  // Option 2D
    };

private:
    float computeGradientSquared(const std::vector<float>& theta, int idx);
    float computeVorticity(const std::vector<float>& theta, int idx);
    bool detectVortexCore(float vorticity);
};
```

### Phase 2: YAML Configuration
Add to test configs:

```yaml
physics:
  K_mode: adaptive  # static | adaptive
  K_adaptation:
    mode: vortex_protected  # gradient_suppressed | order_modulated | vortex_protected | hybrid
    K_max: 1.0
    K_min: 0.0
    G_crit: 0.5
    vortex_threshold: 0.1
    suppression_strength: 0.9
```

### Phase 3: GPU Implementation
Extend Kuramoto shader with K-field buffer:

```glsl
layout(set = 0, binding = 5) buffer KField {
    float K[];
};

void main() {
    uint idx = gl_GlobalInvocationID.x;

    // Compute local K based on theta gradients
    float grad_sq = computeGradientSquared(theta, idx);
    K[idx] = K_max * (1.0 - tanh(grad_sq / G_crit));

    // Use adaptive K in phase evolution
    float coupling_term = K[idx] * R[idx] * sin(Theta[idx] - theta[idx]);
    theta_new[idx] = theta[idx] + dt * (omega[idx] + coupling_term);
}
```

---

## Validation Tests

### Test 1: Vortex Preservation
**Config**: Single vortex, adaptive K
**Success criteria**:
- B_max > 1.0 (EM fields preserved)
- R_avg > 0.9 (synchronization achieved in bulk)
- K-field visualization shows K≈0 at vortex core, K≈K_max elsewhere

### Test 2: Multi-Vortex Dynamics
**Config**: Two opposite-charge vortices
**Success criteria**:
- Both vortices preserved
- Synchronized background between vortices
- EM fields localized to vortex cores

### Test 3: Dynamic Evolution
**Config**: Start with K=0, gradually enable adaptation
**Success criteria**:
- Smooth transition from K=0 (B_max=1.57) to adaptive state
- No numerical instabilities
- Final state has both sync and EM

### Test 4: Comparison with Static K
**Metrics**:
- B_max(t), R_avg(t), EM_energy(t) vs. time
- Compare: K=0 (baseline), K=1.0 (current), K=adaptive (proposed)

---

## Expected Outcomes

### Optimistic Scenario
- ✅ K-field adapts to protect vortices
- ✅ B_max ≈ 1.5, R_avg ≈ 0.95
- ✅ Theory salvaged: "EM emerges from synchronized background with topological structure"

### Realistic Scenario
- Partial success: B_max ≈ 0.8, R_avg ≈ 0.85
- Improvement over static K=1.0 (B_max=0), worse than K=0 (B_max=1.57)
- Theory requires nuance: "EM emerges from partially synchronized systems"

### Pessimistic Scenario
- ❌ No improvement: K-adaptation doesn't stabilize coexistence
- ❌ Numerical instabilities from spatial K-variation
- → Proceed to Option 3 (multi-component) or accept Option 1 (abandon sync)

---

## Physical Justification (Critical)

**Question**: Why should K vary in space/time?

**Possible answers**:
1. **Phenomenological**: Effective description of multi-scale synchronization dynamics
2. **Microscopic**: Underlying oscillators have position-dependent coupling (e.g., distance-dependent)
3. **Emergent**: K adapts to minimize free energy functional F[θ,R,K]
4. **Stochastic**: Thermal/quantum fluctuations drive K-variation

**TODO**: Develop rigorous justification beyond "it works numerically."

---

## Next Steps (After Parameter Scan)

1. **If scan confirms no static K works**:
   - Implement Option 2C (vortex protection) first
   - Validate with Test 1

2. **If scan shows partial coexistence at K ≈ K_c**:
   - Use K_c as reference point
   - Design K(x,t) to interpolate around K_c spatially

3. **If scan shows unexpected behavior**:
   - Re-evaluate theoretical model
   - Consult additional literature

---

**Status**: Design complete, awaiting empirical data from K-parameter scan.
