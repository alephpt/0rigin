# Dirac Coupling Implementation Plan

## Current Status

**Infrastructure**: ✅ Complete
- Spinor buffers allocated: `_spinor_field`, `_spinor_buffer`, `_spinor_density_buffer`
- Shaders written: `dirac_rk4.comp`, `spinor_feedback.comp`  
- Descriptor sets configured
- Pipelines created: `_dirac_pipeline`

**Execution**: ❌ Missing
- No dispatch calls in `SMFTEngine::step()`
- Spinor density always = 0
- No feedback loop active

---

## Required Implementation

### Location: src/SMFTEngine.cpp lines 497-499

Current code:
```cpp
vkCmdDispatch(_compute_command_buffer, workgroupsX, workgroupsY, 1);

// 12. End command buffer recording
result = vkEndCommandBuffer(_compute_command_buffer);
```

### Step 1: Add Dirac Evolution Dispatch (After gravity_field)

Insert after line 497:
```cpp
// 13. Dirac evolution: i∂ₜΨ = [cα·p + βm(x)c²]Ψ
// Input: Ψ_current (spinor_buffer), m(x) (mass_field from sync_field × Δ)
// Output: Ψ_new (spinor_buffer updated)
vkCmdBindPipeline(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, _dirac_pipeline);
vkCmdBindDescriptorSets(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                       _dirac_pipeline_layout, 0, 1, &_dirac_descriptor_set, 0, nullptr);

// Push constants for Dirac evolution
DiracPushConstants dirac_push{};
dirac_push.dt = dt;
dirac_push.dx = 1.0f;  // Lattice spacing
dirac_push.c = 1.0f;   // Speed of light (natural units)
dirac_push.hbar = 1.0f;
vkCmdPushConstants(_compute_command_buffer, _dirac_pipeline_layout, 
                  VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(dirac_push), &dirac_push);

vkCmdDispatch(_compute_command_buffer, workgroupsX, workgroupsY, 1);
```

### Step 2: Add Spinor Feedback Dispatch (After Dirac evolution)

```cpp
// 14. Spinor feedback: compute ρ(x) = Ψ̄Ψ
// Input: Ψ_new (spinor_buffer)
// Output: ρ(x) (spinor_density_buffer)
vkCmdBindPipeline(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, _spinor_feedback_pipeline);
vkCmdBindDescriptorSets(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                       _spinor_feedback_pipeline_layout, 0, 1, &_spinor_feedback_descriptor_set, 0, nullptr);
vkCmdDispatch(_compute_command_buffer, workgroupsX, workgroupsY, 1);
```

### Step 3: Verify Shader Binding Points

Check that descriptor sets connect correctly:

**kuramoto_step.comp** (reads spinor_density):
```glsl
layout(set = 0, binding = 3) readonly buffer SpinorDensityBuffer {
    float spinor_density[];  // ρ(x) from previous step
};
```

**dirac_rk4.comp** (reads mass_field, updates spinor):
```glsl
layout(set = 0, binding = 0) readonly buffer MassField {
    float m[];  // m(x) = Δ·R(x)
};
layout(set = 0, binding = 1) buffer SpinorField {
    vec2 psi[];  // 4 complex components per site
};
```

**spinor_feedback.comp** (reads spinor, writes density):
```glsl
layout(set = 0, binding = 0) readonly buffer SpinorField {
    vec2 psi[];
};
layout(set = 0, binding = 1) writeonly buffer SpinorDensity {
    float rho[];  // ρ = Ψ̄Ψ
};
```

---

## Complete Execution Loop

After implementation, the pipeline will be:

```
1. kuramoto_step(θ, ω, ρ) → θ_new        [uses previous spinor density]
2. sync_field(θ) → R(x)                  [compute synchronization]
3. gravity_field(R) → g(x), m(x)         [compute gravity and mass]
4. dirac_evolution(Ψ, m) → Ψ_new         [NEW: evolve matter in mass field]
5. spinor_feedback(Ψ) → ρ(x)             [NEW: compute density for next step]
```

This closes the loop:
- Vacuum (θ) generates mass m(x) = Δ·R(x)
- Mass modulates Dirac equation → Ψ evolves
- Spinor density ρ = |Ψ|² feeds back to vacuum dynamics
- **Full vacuum ↔ matter coupling achieved**

---

## Implementation Tasks

### Phase 1: Code Integration (1-2 days)
- [ ] Add DiracPushConstants struct to SMFTEngine.h
- [ ] Declare _spinor_feedback_pipeline and _spinor_feedback_pipeline_layout
- [ ] Add dispatch calls to SMFTEngine::step()
- [ ] Verify compilation

### Phase 2: Pipeline Creation (2-3 days)
- [ ] Create _spinor_feedback_pipeline in initVulkan()
- [ ] Load spinor_feedback.comp shader
- [ ] Set up descriptor set layout for spinor feedback
- [ ] Verify pipeline creation

### Phase 3: Initial Conditions (1 day)
- [ ] Initialize spinor field Ψ(x,y,0) with sensible values
- [ ] Options:
  - Gaussian wavepacket
  - Plane wave
  - Random quantum state
  - Ground state of free Dirac equation

### Phase 4: Testing (3-4 days)
- [ ] Test Dirac evolution without feedback (ρ→0)
- [ ] Test spinor feedback computation
- [ ] Verify feedback to Kuramoto (check ρ ≠ 0 in shader)
- [ ] Check energy conservation
- [ ] Measure Zitterbewegung frequency

### Phase 5: Physics Validation (3-4 days)
- [ ] Verify Dirac equation numerics (RK4 accuracy)
- [ ] Check spinor normalization (Ψ̄Ψ conserved)
- [ ] Measure effective noise from Zitterbewegung
- [ ] Compare to Phase 0 predictions
- [ ] Test mass feedback loop stability

---

## Shader Verification Checklist

### dirac_rk4.comp
- [ ] Reads mass field m(x) = Δ·R(x)
- [ ] Implements Dirac evolution: i∂ₜΨ = [cα·p + βm(x)c²]Ψ
- [ ] Uses RK4 for numerical integration
- [ ] Handles periodic boundaries
- [ ] Updates spinor buffer in-place or with double-buffering

### spinor_feedback.comp
- [ ] Reads spinor field Ψ(x)
- [ ] Computes density ρ(x) = Ψ†Ψ (adjoint × spinor)
- [ ] Writes to spinor_density buffer
- [ ] Normalizes appropriately for feedback strength

---

## Expected Physics After Implementation

### 1. Zitterbewegung (Rapid Oscillation)
- Dirac spinor oscillates at ω ≈ 2mc²/ℏ
- For m ≈ 0.7Δ (mean defect mass), expect high frequency
- Should see ρ(x,t) fluctuate rapidly

### 2. Effective Noise
- Zitterbewegung generates "jitter" in ρ(x,t)
- This feeds back to θ dynamics as effective stochastic force
- May replace need for explicit Langevin noise

### 3. Mass-Matter Coupling
- Regions with high R → high m → strong Dirac binding
- Regions with defects (R→0) → low m → free Dirac propagation
- Should see "matter attracted to mass"

### 4. Potential Instabilities
- Positive feedback: ρ↑ → θ evolves → R changes → m changes → ρ changes
- May need damping or coupling strength tuning
- Monitor for runaway synchronization or desynchronization

---

## Debugging Strategy

1. **Start with disabled feedback**
   - Set `enable_feedback = 0` in kuramoto_step.comp
   - Verify Dirac evolves independently
   - Check spinor density computation

2. **Enable weak feedback**
   - Set λ = 0.01 (weak coupling)
   - Monitor R(t) for stability
   - Gradually increase λ

3. **Add diagnostics**
   - Save Ψ(x,y,t) timeseries
   - Compute Ψ̄Ψ integral (should be conserved)
   - Measure energy ⟨Ψ|H|Ψ⟩
   - Track correlation ⟨θ(x)ρ(x)⟩

---

## Timeline Estimate

**Optimistic**: 1 week (if shaders are correct and pipelines work first try)

**Realistic**: 2 weeks
- Week 1: Integration, pipeline creation, initial testing
- Week 2: Physics validation, debugging, parameter tuning

**Pessimistic**: 3-4 weeks (if shader bugs or numerical instabilities found)

---

## Success Criteria

✅ **Minimum Viable**:
- Dirac and spinor_feedback dispatches execute without errors
- spinor_density buffer gets non-zero values
- Kuramoto dynamics shows feedback effect (R changes compared to vacuum-only)

✅ **Full Success**:
- Spinor normalization preserved (Ψ̄Ψ = const)
- Zitterbewegung observable in ρ(x,t) timeseries
- Mass-matter coupling stable over 10k steps
- Energy approximately conserved
- Effective noise from Zitterbewegung matches predictions

✅ **Physics Breakthrough**:
- Defect-spinor interactions create particle-like states
- Mass quantization emerges from Dirac spectrum
- Gravitational waves from spinor-vacuum coupling
- Connection to actual physics (compare to experimental data)

---

## Next Immediate Step

**Start with**: Add pipeline dispatch calls to SMFTEngine::step()
- Low risk, high visibility
- Can test if pipelines are correctly initialized
- Will immediately show if descriptor sets are bound correctly
- Foundation for all subsequent work

**File to edit**: `src/SMFTEngine.cpp` lines 497-499
