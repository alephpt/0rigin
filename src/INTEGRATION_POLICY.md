# SMFT Integration Policy

## MANDATORY: Symplectic Integration Only

**ALL time evolution in SMFT MUST use symplectic integrators.**

This is a **hard architectural requirement** with no exceptions. Energy conservation is fundamental to physical correctness.

---

## Why Symplectic Integration?

Symplectic integrators preserve phase space volume (Liouville's theorem), which guarantees:

1. **Energy Conservation**: Total energy E(t) remains constant (bounded error O(dt²))
2. **Bounded Trajectories**: Orbits stay on their energy surface, preventing secular drift
3. **Long-time Stability**: No spurious exponential growth or decay
4. **Fidelity to Physics**: Hamiltonian structure is preserved numerically

**Contrast with RK4** (non-symplectic):
- RK4 conserves energy order O(dt⁴) locally, but errors accumulate
- In long simulations, energy drifts secularly: E(t) ≈ E₀ + α·t
- Particle speeds drift: |v(t)| = |v₀| + β·t² → unstable orbits
- Completely unphysical for systems where energy should be conserved

---

## Approved Integrators

### 1. **Velocity Verlet** (2nd order, O(dt²))
**For particle dynamics and light systems**

Structure:
```
KICK: p → p + (dt/2)·F(x)
DRIFT: x → x + dt·p/m
KICK: p → p + (dt/2)·F(x')
```

Advantages:
- Simple to implement
- Very stable for short-range forces
- Excellent for cyclotron motion
- Integrates position and momentum on equal footing

Implementation: `TestParticle.cpp` - `evolveLorentzForce()`

---

### 2. **Strang Splitting** (2nd order, O(dt²))
**For field evolution with separable Hamiltonians**

Structure:
```
DRIFT: evolve free part for dt/2
KICK: evolve interaction for dt
DRIFT: evolve free part for dt/2
```

Advantages:
- Exact for quadratic kinetic energy
- Separates linear and nonlinear parts
- Used in SMFT field evolution

Implementation: `DiracEvolution.cpp` - `integrateStrangSplitting()`

---

### 3. **Leapfrog** (2nd order, O(dt²))
**For coupled systems with clear position/velocity separation**

Structure:
```
v → v + (dt/2)·a(x)
x → x + dt·v
v → v + (dt/2)·a(x')
```

Advantages:
- Nearly identical to Velocity Verlet
- Very stable numerically
- Good for gravitational N-body problems

---

## FORBIDDEN Integrators

### ❌ **RK4** (4th order, but non-symplectic)
**NEVER USE - Energy non-conservation causes unphysical drift**

Why forbidden:
- Non-symplectic: violates phase space volume conservation
- Secular energy drift in long simulations
- Particle speed drifts uncontrollably
- No physical justification for Hamiltonian systems

Evidence from this codebase:
- TestParticle with RK4: |v| drifts from 0.01 → 0.0213 (113% error!)
- Pure magnetic field should conserve speed exactly
- RK4 violates this fundamental conservation law

---

### ❌ **Euler** (1st order O(dt))
**NEVER USE - Too inaccurate and non-symplectic**

Why forbidden:
- Only 1st order: need dt << 1 for accuracy
- Non-symplectic: energy drifts
- Obsolete for modern computing

---

### ❌ **Midpoint/Crank-Nicolson** (for Hamiltonian systems)
**NEVER USE - While symplectic, not designed for our use case**

Why avoided:
- More complex than Velocity Verlet
- No advantage for our particle dynamics
- Reserved for implicit integration of stiff systems

---

## Code Review Checklist

Every code review must verify:

- [ ] **No RK4 in particle evolution** - grep for "RK4", "integrateRK4"
- [ ] **No Euler steps** - grep for "Euler", "forward_euler", "explicit_euler"
- [ ] **Velocity Verlet or Strang used** - explicit in implementation
- [ ] **Energy conservation tested** - test suite validates dE/E < 1% for test duration
- [ ] **Long-time stability verified** - simulation runs 10,000+ steps without drift
- [ ] **No TODO comments about integration** - if spotted, block merge

---

## Violation Consequences

If this policy is violated:

1. **Code Review Fails** - Will not merge
2. **Tests Fail** - Energy conservation tests will catch it
3. **Results Untrustworthy** - Any paper/publication based on violated code is invalid
4. **Physics Compromised** - Violates fundamental conservation law

---

## Adding New Integrators

If a new integrator is needed:

1. **Verify symplecticity** - Prove phase space volume preservation
2. **Test on toy problem** - Harmonic oscillator or cyclotron motion
3. **Validate energy conservation** - dE/E << 0.1% over 10,000 steps
4. **Document in this file** - Add to "Approved Integrators" section
5. **Get architecture review** - SMFT core team approval required

---

## References

- Verlet, L. (1967). Computer Experiments on Classical Fluids. *Physical Review*, 159(1), 98–103.
- Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration*. Springer.
- Candy & Rozmus (1991). Symplectic integration algorithms. *Journal of Computational Physics*, 92(1), 230–256.

---

**This policy is non-negotiable. Energy conservation is physics, not optional.**
