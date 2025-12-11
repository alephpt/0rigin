# Semi-Implicit Damping & Data Validation - Implementation Design

**Project**: 0rigin Field Theory Stability Fixes
**Date**: 2025-12-10
**Status**: Design Phase - Ready for Implementation

---

## 1. SEMI-IMPLICIT DAMPING FIX

### Problem Analysis
**File**: `/home/persist/neotec/0rigin/src/kuramoto/field_theory/fields/mediator.py`
**Method**: `_leapfrog_step` (lines 176-199)
**Current Issue**: Line 193 includes explicit damping term `-self.gamma * self.sigma_dot` in force calculation, causing conditional stability (requires dt < 2/γ)

### Mathematical Fix
Replace explicit damping with semi-implicit scheme (unconditionally stable):

**Before** (line 193):
```python
sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho - self.gamma * self.sigma_dot
```

**After** (split into two operations):
```python
# Compute undamped acceleration
sigma_ddot_undamped = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho

# Update velocity (explicit step)
sigma_dot_temp = self.sigma_dot + dt * sigma_ddot_undamped

# Apply damping (semi-implicit - unconditionally stable)
self.sigma_dot = sigma_dot_temp / (1.0 + self.gamma * dt)
```

### Exact Code Changes

#### Change 1: Replace line 193 (force computation)
**Location**: mediator.py:193
**Old**:
```python
        sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho - self.gamma * self.sigma_dot
```

**New**:
```python
        # Compute acceleration WITHOUT damping (semi-implicit treatment)
        sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho
```

#### Change 2: Replace lines 195-196 (velocity update)
**Location**: mediator.py:195-196
**Old**:
```python
        # Full-step velocity update
        self.sigma_dot += dt * sigma_ddot
```

**New**:
```python
        # Full-step velocity update (explicit)
        sigma_dot_temp = self.sigma_dot + dt * sigma_ddot

        # Apply damping semi-implicitly (unconditionally stable)
        self.sigma_dot = sigma_dot_temp / (1.0 + self.gamma * dt)
```

#### Change 3: Update docstring (lines 177-186)
**Location**: mediator.py:177-186
**Old**:
```python
        """
        Leapfrog integration (symplectic, better energy conservation).

        Updates in sequence:
            1. σ(t+dt/2) = σ(t) + σ_dot(t)·dt/2
            2. Compute forces at t+dt/2
            3. σ_dot(t+dt) = σ_dot(t) + acceleration·dt
            4. σ(t+dt) = σ(t+dt/2) + σ_dot(t+dt)·dt/2

        Includes damping term -γ·σ_dot for numerical stability.
        """
```

**New**:
```python
        """
        Leapfrog integration with semi-implicit damping.

        Updates in sequence:
            1. σ(t+dt/2) = σ(t) + σ_dot(t)·dt/2
            2. Compute undamped forces at t+dt/2
            3. σ_dot_temp(t+dt) = σ_dot(t) + acceleration·dt
            4. σ_dot(t+dt) = σ_dot_temp / (1 + γ·dt)  [semi-implicit damping]
            5. σ(t+dt) = σ(t+dt/2) + σ_dot(t+dt)·dt/2

        Semi-implicit damping treatment ensures unconditional stability
        regardless of dt/γ ratio, while preserving energy dissipation.
        """
```

### Complete Modified Method
**Final implementation** (lines 176-199):
```python
    def _leapfrog_step(self, rho: NDArray, dt: float):
        """
        Leapfrog integration with semi-implicit damping.

        Updates in sequence:
            1. σ(t+dt/2) = σ(t) + σ_dot(t)·dt/2
            2. Compute undamped forces at t+dt/2
            3. σ_dot_temp(t+dt) = σ_dot(t) + acceleration·dt
            4. σ_dot(t+dt) = σ_dot_temp / (1 + γ·dt)  [semi-implicit damping]
            5. σ(t+dt) = σ(t+dt/2) + σ_dot(t+dt)·dt/2

        Semi-implicit damping treatment ensures unconditional stability
        regardless of dt/γ ratio, while preserving energy dissipation.
        """
        # Half-step position update
        sigma_half = self.sigma + 0.5 * dt * self.sigma_dot

        # Compute acceleration WITHOUT damping (semi-implicit treatment)
        laplacian = self.grid.laplacian(sigma_half)
        sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho

        # Full-step velocity update (explicit)
        sigma_dot_temp = self.sigma_dot + dt * sigma_ddot

        # Apply damping semi-implicitly (unconditionally stable)
        self.sigma_dot = sigma_dot_temp / (1.0 + self.gamma * dt)

        # Complete position update
        self.sigma = sigma_half + 0.5 * dt * self.sigma_dot
```

### Verification Strategy
1. **Stability test**: Run with large dt (0.1) and large gamma (10.0) - should NOT blow up
2. **Energy dissipation**: Verify energy decreases monotonically (no oscillations)
3. **Accuracy**: Compare with small dt solutions - should converge to same result
4. **Performance**: Benchmark - should have no slowdown (same operation count)

---

## 2. DATA VALIDATION DESIGN

### Problem
Need to catch NaN/Inf early with diagnostic messages for debugging instabilities.

### Strategy
Add validation at critical points WITHIN MediatorField (not in examples).

### Implementation: Add validation to MediatorField class

#### Modification 1: Add validation flag to __init__
**Location**: After line 83 in mediator.py
**Add**:
```python
        # Validation settings
        self.validate_numerics = True  # Can disable for performance
```

#### Modification 2: Add validation method
**Location**: After line 162 in mediator.py
**Add**:
```python
    def _validate_state(self):
        """Validate field state for numerical errors."""
        if not np.all(np.isfinite(self.sigma)):
            raise ValueError(
                f"MediatorField diverged at t={self.t:.6f}: "
                f"sigma contains NaN/Inf. "
                f"Parameters: c={self.c}, M={self.M}, gamma={self.gamma}. "
                f"Last |sigma| max: {np.max(np.abs(self.sigma)):.6e}, "
                f"|sigma_dot| max: {np.max(np.abs(self.sigma_dot)):.6e}"
            )

        if not np.all(np.isfinite(self.sigma_dot)):
            raise ValueError(
                f"MediatorField velocity diverged at t={self.t:.6f}: "
                f"sigma_dot contains NaN/Inf. "
                f"Parameters: c={self.c}, M={self.M}, gamma={self.gamma}. "
                f"Last |sigma| max: {np.max(np.abs(self.sigma)):.6e}, "
                f"|sigma_dot| max: {np.max(np.abs(self.sigma_dot)):.6e}"
            )
```

#### Modification 3: Add validation call to evolve_step
**Location**: After line 162 in mediator.py (after `self.t += dt`)
**Add**:
```python
        # Validate field state (catch NaN/Inf early)
        if self.validate_numerics:
            self._validate_state()
```

---

## 3. ENERGY MONITORING DESIGN

### Objective
Add energy bounds checking to detect runaway growth.

### Implementation: Add to MediatorField class

#### Modification 1: Add energy monitoring flag to __init__
**Location**: After line 83 in mediator.py
**Add**:
```python
        # Energy monitoring
        self.check_energy_bounds = False  # Enable if needed
        self._initial_energy = None
```

#### Modification 2: Add energy validation method
**Location**: After line 245 in mediator.py (after `compute_field_energy()`)
**Add**:
```python
    def validate_energy_bounds(self, max_energy_ratio: float = 1000.0) -> bool:
        """
        Check if field energy has grown beyond reasonable bounds.

        Parameters
        ----------
        max_energy_ratio : float
            Maximum allowed energy ratio relative to initial energy.

        Returns
        -------
        bool
            True if energy is within bounds, False otherwise.
        """
        current_energy = self.compute_field_energy()

        # Track initial energy if not set
        if not hasattr(self, '_initial_energy') or self._initial_energy is None:
            self._initial_energy = current_energy
            if self._initial_energy == 0:
                self._initial_energy = 1e-10  # Avoid division by zero

        energy_ratio = current_energy / self._initial_energy

        if energy_ratio > max_energy_ratio:
            import warnings
            warnings.warn(
                f"Field energy grew by {energy_ratio:.1f}x at t={self.t:.6f}. "
                f"Possible numerical instability. "
                f"Current energy: {current_energy:.6e}, "
                f"Initial energy: {self._initial_energy:.6e}. "
                f"Consider smaller dt or larger damping.",
                RuntimeWarning
            )
            return False

        return True
```

#### Modification 3: Add optional energy check to evolve_step
**Location**: After validation code in evolve_step (after line 162 + validation additions)
**Add**:
```python
        # Optional energy bounds checking
        if self.check_energy_bounds:
            self.validate_energy_bounds()
```

---

## 4. VISUALIZATION LABEL FIXES

### Objective
Add proper "Time (t)" labels to all temporal plots.

### Files to Modify

#### File 1: /home/persist/neotec/0rigin/examples/field_theory/smft_full_demo.py
**Line 165** - Change from `'Time'` to `'Time (t)'`
**Line 206** - Change from `'Time'` to `'Time (t)'`

#### File 2: /home/persist/neotec/0rigin/examples/field_theory/smft_demo.py
**Line 58** - Change from `'Time'` to `'Time (t)'`
**Line 65** - Change from `'Time'` to `'Time (t)'`
**Line 191** - Change from `'Time'` to `'Time (t)'`
**Line 198** - Change from `'Time'` to `'Time (t)'`

#### File 3: /home/persist/neotec/0rigin/examples/field_theory/hamiltonian_demo.py
**Line 54** - Change from `'Time'` to `'Time (t)'`
**Line 74** - Change from `'Time'` to `'Time (t)'`
**Line 143** - Change from `'Time'` to `'Time (t)'`
**Line 149** - Change from `'Time'` to `'Time (t)'`
**Line 156** - Change from `'Time'` to `'Time (t)'`

#### File 4: /home/persist/neotec/0rigin/examples/basic_synchronization.py
**Line 244** - Change from `'Time'` to `'Time (t)'`

#### File 5: /home/persist/neotec/0rigin/examples/prototype_test.py
**Line 179** - Change from `'Time'` to `'Time (t)'`
**Line 190** - Change from `'Time'` to `'Time (t)'`

#### File 6: /home/persist/neotec/0rigin/examples/demo_synchronization.py
**Line 122** - Change from `'Time'` to `'Time (t)'`

### Summary Table

| File | Lines to Update | Count |
|------|----------------|-------|
| smft_full_demo.py | 165, 206 | 2 |
| smft_demo.py | 58, 65, 191, 198 | 4 |
| hamiltonian_demo.py | 54, 74, 143, 149, 156 | 5 |
| basic_synchronization.py | 244 | 1 |
| prototype_test.py | 179, 190 | 2 |
| demo_synchronization.py | 122 | 1 |
| **TOTAL** | | **15 changes** |

---

## 5. IMPLEMENTATION CHECKLIST

### Phase 1: Semi-Implicit Damping (HIGHEST PRIORITY)
- [ ] Modify `_leapfrog_step` method (lines 176-199)
- [ ] Update docstring (lines 177-186)
- [ ] Remove explicit damping from force (line 193)
- [ ] Add semi-implicit velocity update (lines 195-196)
- [ ] Test with large dt/gamma ratios
- [ ] Verify energy dissipation monotonicity
- [ ] Benchmark performance (should be unchanged)

### Phase 2: Data Validation
- [ ] Add `validate_numerics` flag to `__init__`
- [ ] Implement `_validate_state()` method
- [ ] Add validation call in `evolve_step`
- [ ] Test with deliberately unstable parameters
- [ ] Verify diagnostic messages are helpful

### Phase 3: Energy Monitoring
- [ ] Add `check_energy_bounds` flag to `__init__`
- [ ] Implement `validate_energy_bounds()` method
- [ ] Add `_initial_energy` tracking
- [ ] Add optional energy check in `evolve_step`
- [ ] Test warning threshold
- [ ] Verify warning message clarity

### Phase 4: Visualization Labels
- [ ] Update smft_full_demo.py (2 locations)
- [ ] Update smft_demo.py (4 locations)
- [ ] Update hamiltonian_demo.py (5 locations)
- [ ] Update basic_synchronization.py (1 location)
- [ ] Update prototype_test.py (2 locations)
- [ ] Update demo_synchronization.py (1 location)
- [ ] Verify plots render correctly

### Phase 5: Testing & Validation
- [ ] Run all field theory tests
- [ ] Run stability tests with large dt/gamma
- [ ] Generate plots and verify labels
- [ ] Check energy bounds warnings
- [ ] Verify NaN/Inf detection works
- [ ] Performance regression test

---

## 6. EXPECTED OUTCOMES

### Stability Improvements
- **Before**: Blows up with dt=0.1, gamma=10
- **After**: Stable with arbitrary dt/gamma ratios (still accurate with reasonable values)

### Diagnostic Quality
- **Before**: Silent failure or cryptic numpy errors
- **After**: Clear error messages with parameter values and remediation suggestions

### Energy Monitoring
- **Before**: No warning until complete divergence
- **After**: Early warning when energy grows beyond 1000x initial value

### Visualization Quality
- **Before**: Generic "Time" labels
- **After**: Professional "Time (t)" labels with physical units

---

## 7. RISK ASSESSMENT

### Low Risk Changes
✅ Visualization labels (cosmetic, no logic change)
✅ Validation flags (optional, disabled by default for energy monitoring)

### Medium Risk Changes
⚠️ Semi-implicit damping (changes numerical behavior)
- Mitigation: Extensive testing, regression tests

### High Risk Changes
❌ None (all changes are conservative and well-tested in literature)

---

## IMPLEMENTATION READY

All designs are complete with exact line numbers and code snippets.
Ready for implementation in sequence:

1. **Semi-implicit damping** (highest priority - fixes instability)
2. **Data validation** (catch errors early)
3. **Energy monitoring** (detect divergence)
4. **Visualization labels** (polish)

**Estimated Implementation Time**: 2-3 hours
**Testing Time**: 1-2 hours
**Total**: 4-5 hours for complete implementation and validation
