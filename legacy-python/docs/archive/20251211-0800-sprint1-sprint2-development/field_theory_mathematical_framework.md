# Field Theory Mathematical Framework

**Sprint 2 Discovery Output**
**Date**: December 10, 2025

---

## Core Transformation Sequence

### 1. Discrete Kuramoto → Hamiltonian Kuramoto

**Original Kuramoto**:
```
dθⱼ/dt = ωⱼ + (K/N) Σₖ sin(θₖ - θⱼ)
```

**Hamiltonian Extension**:
```
H = Σⱼ [pⱼ²/(2m) - ωⱼθⱼ] - (K/2N) Σⱼₖ cos(θⱼ - θₖ)

dθⱼ/dt = ∂H/∂pⱼ = pⱼ/m
dpⱼ/dt = -∂H/∂θⱼ - γpⱼ = ωⱼ + (K/N)Σₖ sin(θₖ - θⱼ) - γpⱼ

Overdamped limit (m/γ → 0): Recovers original Kuramoto
```

### 2. Global → Local Coupling

**Global Coupling** (all-to-all):
```
Fⱼ = (K/N) Σₖ sin(θₖ - θⱼ)
```

**Local Field Coupling** (via mediator σ):
```
Oscillator source: Jⱼ = g·exp(iθⱼ)δ(x - xⱼ)
Field equation: (∇² - M²)σ = J
Coupling force: Fⱼ = g·Im[σ(xⱼ)exp(-iθⱼ)]

Heavy field (M→∞): σ ≈ J/M² → Global coupling
Light field (M→0): σ propagates → Retarded interactions
```

### 3. Discrete → Continuum (N → ∞)

**Discrete System** (N oscillators):
```
θⱼ(t), j = 1...N
Order parameter: Z = (1/N)Σⱼ exp(iθⱼ)
```

**Continuum Limit** (Field theory):
```
Phase field: θ(x,t), x ∈ ℝᵈ
Density: ρ(x,t) = Σⱼ δ(x - xⱼ)
Order field: Z(x,t) = ∫ ρ(y,t)exp(iθ(y,t))dy

Field equations (from Ott-Antonsen):
∂²R/∂t² + γ∂R/∂t = D∇²R + μ²R - λR³
R²∂²ψ/∂t² + 2R∂R/∂t·∂ψ/∂t = D∇²ψ

where Z(x,t) = R(x,t)exp(iψ(x,t))
```

### 4. Non-relativistic → Relativistic

**Non-relativistic Action**:
```
S = ∫dt d³x [
    (1/2)(∂ₜR)² - (1/2)(∇R)² - V(R)
    + (1/2)R²(∂ₜθ)² - (1/2)R²(∇θ)²
]
```

**Lorentz-Covariant Action**:
```
S = ∫d⁴x √-g [
    (1/2)g^μν∂μR∂νR - (m²/2)R² - (λ/4)R⁴
    + (1/2)R²g^μν∂μθ∂νθ
    + ψ̄(iγ^μDμ - ΔR·exp(iθγ⁵))ψ
]

Metric: g_μν = diag(-1,+1,+1,+1)
Covariant derivative: Dμ = ∂μ + ieAμ
```

---

## Key Mathematical Structures

### Ott-Antonsen Manifold

**Ansatz**: Fourier modes satisfy fₙ(ω,t) = [α(ω,t)]ⁿ

**Reduction**: ∞-dimensional → 2-dimensional dynamics

**Evolution**:
```
∂ₜα = iωα + (K/2)(ᾱ - α|α|²)

Fixed points:
- Incoherent: α = 0 (unstable for K > Kc)
- Synchronized: |α| = √(1 - 2ω/K) (stable)
```

### Continuum vs Mean-Field Limits

| Property | Continuum Limit | Mean-Field Limit |
|----------|-----------------|------------------|
| Variable | x(t,ξ), ξ ∈ [0,1] | μₜ(x) probability measure |
| Equation | ∂ₜx = ∫sin(x(ξ) - x(z))dz | ∂ₜμ + ∇·(μv[μ]) = 0 |
| Tracks | Individual particles | Population density |
| Manifold | Unstable manifold of x₀ | Ott-Antonsen manifold |

### Phase Transitions

**Order Parameter Scaling**:
```
Near Kc: R ~ √(K - Kc)  (mean-field exponent β = 1/2)

Critical slowing: τ ~ |K - Kc|⁻¹

Correlation length: ξ ~ |K - Kc|⁻ν, ν = 1/2 (mean-field)
```

**Goldstone Mode** (broken U(1) symmetry):
```
Massless mode: δθ with ω² = Dk² as k → 0
Massive mode: δR with ω² = μ² + Dk²
```

---

## Numerical Discretization

### Spatial Grid

```python
# 2D lattice example
Nx, Ny = 100, 100
dx, dy = L/Nx, L/Ny

# Fields on grid
R[i,j] = R(i*dx, j*dy)
θ[i,j] = θ(i*dx, j*dy)

# Laplacian (5-point stencil)
∇²R[i,j] = (R[i+1,j] + R[i-1,j] + R[i,j+1] + R[i,j-1] - 4R[i,j])/dx²
```

### Time Integration

**Velocity-Verlet** (symplectic):
```
θₙ₊₁/₂ = θₙ + (dt/2)·pₙ/m
pₙ₊₁ = pₙ + dt·F(θₙ₊₁/₂) - γdt·pₙ₊₁/₂
θₙ₊₁ = θₙ₊₁/₂ + (dt/2)·pₙ₊₁/m
```

**Runge-Kutta 4** (high accuracy):
```
k₁ = f(tₙ, yₙ)
k₂ = f(tₙ + dt/2, yₙ + dt·k₁/2)
k₃ = f(tₙ + dt/2, yₙ + dt·k₂/2)
k₄ = f(tₙ + dt, yₙ + dt·k₃)
yₙ₊₁ = yₙ + (dt/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

### FFT for Convolutions

```python
# Screened Poisson: (∇² - M²)σ = J
# In Fourier space: -(k² + M²)σ̃ = J̃

J_k = fft2(J)                    # Source in k-space
sigma_k = -J_k / (k² + M²)       # Solve in k-space
sigma = ifft2(sigma_k)           # Back to real space
```

---

## Fermion Coupling (MSFT)

### Dirac Equation with Synchronization Mass

```
(iγ^μ∂μ - M(x,t))ψ = 0

Mass operator: M(x,t) = ΔR(x,t)exp(iθ(x,t)γ⁵)

Effective mass: m_eff = Δ⟨R⟩
```

### Yukawa Vertex

```
Lagrangian: ℒ_int = gψ̄(R + iθγ⁵)ψ

Vertex factor: -ig(1 + iγ⁵) in Feynman diagrams
```

### Chiral Symmetry Breaking

```
Before sync (R=0): ψ_L ↔ ψ_R symmetry
After sync (R≠0): ⟨ψ̄ψ⟩ ≠ 0, mass generated
```

---

## Validation Tests

### 1. Limit Recovery
```
γ → ∞: Hamiltonian → Original Kuramoto ✓
M → ∞: Local coupling → Global coupling ✓
D → ∞: Field theory → Mean-field ✓
```

### 2. Conservation Laws
```
Energy (undamped): dE/dt = 0
Momentum (translation invariant): dP/dt = 0
Angular momentum (rotation invariant): dL/dt = 0
```

### 3. Scaling Relations
```
R ~ (K - Kc)^β, β = 1/2
τ ~ |K - Kc|^-z, z = 1
ξ ~ |K - Kc|^-ν, ν = 1/2
```

### 4. Lorentz Invariance
```
Under boost x' = γ(x - vt), t' = γ(t - vx/c²):
R'(x',t') = R(x,t)  (scalar)
∂'μR∂'^μR = ∂μR∂^μR  (invariant)
```

---

## Implementation Checklist

### Week 1
- [ ] Hamiltonian formulation with momentum
- [ ] Verify overdamped limit
- [ ] Energy conservation test

### Week 2
- [ ] Screened Poisson for mediator field
- [ ] Local coupling implementation
- [ ] Heavy-field limit verification

### Week 3
- [ ] Continuum field equations
- [ ] Ott-Antonsen reduction
- [ ] Convergence as dx → 0

### Week 4
- [ ] Fermion coupling
- [ ] Mass generation demonstration
- [ ] Full MSFT integration

---

*This framework provides the mathematical foundation for Sprint 2 implementation*