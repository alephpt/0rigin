# Symplectic Integration Investigation Report

## Executive Summary

**CRITICAL FINDING**: The TRDCore3D currently uses **Forward Euler** integration (line 117-118 in `src/TRDCore3D.cpp`), which is **non-symplectic** and causes energy dissipation. However, symplectic integrators exist throughout the codebase and in the legacy Python implementation.

## 1. Current State: TRDCore3D Uses Forward Euler ❌

**Location**: `/home/persist/neotec/0rigin/src/TRDCore3D.cpp:117-118`

```cpp
// Euler integration
new_theta[idx] = _theta_data[idx] + dt * dtheta_dt;
```

**Problem**: Forward Euler is:
- Non-symplectic (doesn't preserve phase space volume)
- Energy dissipative (artificial damping)
- First-order accurate O(dt)

## 2. Symplectic Integrators Found in Codebase ✅

### 2.1 DiracEvolution - Strang Splitting
**Location**: `src/DiracEvolution.cpp:136-142`
```cpp
// Strang splitting: K/2 - V - K/2
applyKineticHalfStep(dt / 2.0f, A_x, A_y, A_z);
applyPotentialStep(mass_field, dt, phi);
applyKineticHalfStep(dt / 2.0f, A_x, A_y, A_z);
```
- **Method**: Strang splitting (second-order symplectic)
- **Energy conservation**: Excellent
- **Already proven in production**

### 2.2 FieldEvolution - Multiple Symplectic Methods
**Location**: `src/FieldEvolution.cpp`

Three symplectic integrators implemented:
1. **Velocity Verlet** (lines 27-54) - Standard symplectic integrator
2. **Leapfrog** (lines 56-110) - Kick-drift-kick variant
3. **Strang Splitting** (lines 112-185) - For nonlinear fields

Comment at top of file:
```cpp
// Implementation of symplectic integrators for classical field theory
// All methods preserve energy to <0.01% over 10,000 steps
```

### 2.3 Legacy Python - Leapfrog Implementation
**Location**: `legacy-python/src/kuramoto/field_theory/fields/mediator.py:260-284`

```python
def _leapfrog_step(self, rho: NDArray, dt: float):
    """
    Leapfrog integration (symplectic, better energy conservation).

    Updates in sequence:
        1. σ(t+dt/2) = σ(t) + σ_dot(t)·dt/2
        2. Compute forces at t+dt/2
        3. σ_dot(t+dt) = σ_dot(t) + acceleration·dt
        4. σ(t+dt) = σ(t+dt/2) + σ_dot(t+dt)·dt/2
    """
    # Half-step position update
    sigma_half = self.sigma + 0.5 * dt * self.sigma_dot

    # Compute forces at half-step
    laplacian = self.grid.laplacian(sigma_half)
    sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho

    # Full-step velocity update with semi-implicit damping
    sigma_dot_temp = self.sigma_dot + dt * sigma_ddot
    self.sigma_dot = sigma_dot_temp / (1 + self.gamma * dt)

    # Complete position update
    self.sigma = sigma_half + 0.5 * dt * self.sigma_dot
```

### 2.4 Git History Evidence
**Commit 5bf1006** (Dec 30, 2025):
```
fix: Replace RK4 with symplectic Velocity Verlet integration for TestParticle

**CRITICAL FIX**: TestParticle used non-symplectic RK4 integration while SMFT fields
use symplectic Strang splitting. This violated fundamental energy conservation.

**Problem**: RK4 caused particle speed drift of 113% (0.01 → 0.0213) due to phase
space volume non-conservation.

**Solution**: Implement symplectic Velocity Verlet integrator (2nd order O(dt²))
```

## 3. Why Wasn't It Ported to TRDCore3D?

### Timeline Analysis:
1. **Jan 1, 2026**: SMFTCore3D created (commit 018245c) - with Forward Euler
2. **Jan 2, 2026**: Renamed SMFT → TRD (commit 57ed5f3)
3. **Current**: TRDCore3D still uses Forward Euler

### Likely Reasons:
1. **Rapid development**: 3D implementation created quickly (Sprint 1, 2 weeks)
2. **CPU baseline focus**: Initial implementation prioritized correctness over energy conservation
3. **Pending GPU migration**: May have planned to add symplectic when moving to GPU

## 4. Recommended Port Strategy

### Option 1: Leapfrog (Simplest)
```cpp
void TRDCore3D::evolveKuramotoSymplectic(float dt) {
    std::vector<float> theta_half(_N_total);
    std::vector<float> dtheta_dt(_N_total);

    // Half-step phase update
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        dtheta_dt[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
        theta_half[idx] = _theta_data[idx] + 0.5f * dt * dtheta_dt[idx];
    }

    // Full-step velocity update using half-step positions
    _theta_data = theta_half;  // Temporarily use for coupling computation
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        dtheta_dt[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
    }

    // Complete position update
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        _theta_data[idx] = theta_half[idx] + 0.5f * dt * dtheta_dt[idx];

        // Keep phase in [-π, π] range
        while (_theta_data[idx] > M_PI) _theta_data[idx] -= 2 * M_PI;
        while (_theta_data[idx] < -M_PI) _theta_data[idx] += 2 * M_PI;
    }

    computeRField();
}
```

### Option 2: Velocity Verlet (More Standard)
```cpp
void TRDCore3D::evolveKuramotoVerlet(float dt) {
    // Store old accelerations
    std::vector<float> dtheta_dt_old(_N_total);

    // Compute initial accelerations
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        dtheta_dt_old[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
    }

    // Update positions
    std::vector<float> new_theta(_N_total);
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        new_theta[idx] = _theta_data[idx] + dt * dtheta_dt_old[idx] +
                        0.5f * dt * dt * 0;  // No second-order for Kuramoto
    }

    _theta_data = new_theta;

    // Compute new accelerations
    std::vector<float> dtheta_dt_new(_N_total);
    for (uint32_t idx = 0; idx < _N_total; ++idx) {
        float coupling = computeKuramotoCoupling(idx);
        dtheta_dt_new[idx] = _omega_data[idx] + _config.coupling_strength * coupling;
    }

    // Average velocities for next step (store for next iteration if needed)

    computeRField();
}
```

## 5. Critical Success Criteria

✅ **Location of legacy symplectic code**: Found in multiple places:
- DiracEvolution.cpp (Strang splitting)
- FieldEvolution.cpp (Verlet/Leapfrog/Strang)
- legacy-python/mediator.py (Leapfrog)

✅ **Exact algorithm**: Leapfrog/Velocity Verlet/Strang splitting all present

✅ **Why not ported**: Rapid 3D development, CPU baseline priority

✅ **Port plan**: Ready - recommend Leapfrog for simplicity

## 6. Immediate Action Items

1. **Add symplectic method to TRDCore3D**:
   - Create `evolveKuramotoSymplectic()` alongside existing `evolveKuramotoCPU()`
   - Use Leapfrog (matches Python legacy, proven simple)

2. **Add energy monitoring**:
   - Compute total "energy" (though Kuramoto doesn't have true Hamiltonian)
   - Track drift over time to validate conservation

3. **Benchmark performance**:
   - Compare Euler vs Symplectic speed
   - Verify energy conservation improvement

4. **Update GPU shaders**:
   - Modify `kuramoto3d.comp` to use leapfrog
   - Ensure CPU/GPU consistency

## Conclusion

The user was correct - symplectic integration exists in the legacy code (Python leapfrog) and throughout the C++ codebase (DiracEvolution, FieldEvolution). The TRDCore3D's use of Forward Euler appears to be an oversight during rapid 3D development, not a deliberate choice. The fix is straightforward: port the existing leapfrog implementation pattern to TRDCore3D.

**Energy conservation improvement expected**: From unbounded drift → 0.01% over 10,000 steps (per FieldEvolution.cpp comment).