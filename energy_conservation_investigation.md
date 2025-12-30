# EM Energy Conservation Investigation

## Problem
EM verification tests show 1.77%-561% energy drift, exceeding 0.1% scientific tolerance.

## Investigation Findings

### 1. Missing Kuramoto Field Energy (Partially addressed)
- **Issue**: Total energy was missing the Kuramoto field's own gradient and synchronization energy
- **Fix Attempted**: Added `computeKuramotoFieldEnergy()` but disabled due to NaN issues with vortex configurations
- **Status**: Needs proper phase-wrapped gradient computation for vortices

### 2. Phase Wrapping in EM Field Computation (Partially fixed)
- **Issue**: When computing EM fields from phase θ:
  - Temporal derivative: ∂θ/∂t was not handling 2π wrapping
  - Spatial derivatives: ∇θ was not handling branch cuts in vortices
- **Fix Applied**: 
  - Added phase wrapping to temporal differences in `EMFieldComputer::computeFromPhase()`
  - Added phase wrapping to spatial derivatives with `is_phase_field` flag
- **Result**: Energy drift reduced from 2×10^9% to 561%, but still far from <0.1% target

### 3. Remaining Issues

#### A. Numerical Integration Scheme
- Using operator splitting (Strang or sequential)
- Each substep may accumulate energy errors
- GPU computation may have different numerical precision

#### B. EM-Matter Coupling Terms
- Current energy budget includes:
  - Dirac kinetic energy: T
  - Dirac-Kuramoto coupling (mass term): V
  - EM field energy: (E² + B²)/(8π)
- Potentially missing:
  - EM-matter interaction energy beyond the mass term
  - Gauge field self-interaction terms
  - Poynting flux contributions

#### C. Field Energy Computation
- EM field energy uses (E² + B²)/(8π)
- But E and B are derived from phase gradients
- For vortex with winding number 1, gradients can be large near core
- Even with phase wrapping, numerical derivatives may be inaccurate

### 4. Next Steps

1. **Verify Phase Evolution**: Check if phase field itself is evolving correctly
2. **Test Without Vortex**: Try simpler initial conditions without topological defects
3. **Check GPU Shaders**: Ensure GPU evolution preserves energy at each step
4. **Add Energy Monitoring**: Track energy components at each substep
5. **Consider Symplectic Integrator**: Use energy-preserving numerical scheme

### 5. Test Results Summary

| Configuration | Initial Energy | Final Energy | Drift % |
|--------------|---------------|--------------|---------|
| Original (no fix) | 0.311 | 6.28×10^6 | 2×10^9% |
| With Kuramoto energy | 0.361 | -3437 | -952488% |
| Phase wrapping only | 0.361 | 2.38 | 561% |

The phase wrapping fix helps but is insufficient. The problem appears to be fundamental to how EM fields are extracted from the Kuramoto phase in the presence of topological defects.
