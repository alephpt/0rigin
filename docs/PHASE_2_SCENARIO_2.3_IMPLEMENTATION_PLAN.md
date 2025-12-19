# Phase 2 Scenario 2.3: Relativistic Mass Testing - Implementation Plan

## Objective
Test the relativistic mass formula: **m(v) = γ · Δ · R** where γ = 1/√(1 - v²/c²)

## Physics Background

In Planck units (ℏ = c = G = Δ = 1):
- Rest mass: m₀ = Δ · R
- Boosted mass: m(v) = γ · m₀ = γ · Δ · R
- For v = 0.5c → γ = 1.155 → 15.5% mass increase

## Test Cases

| Velocity (v/c) | Lorentz Factor γ | Expected Mass Increase |
|----------------|------------------|------------------------|
| 0.0            | 1.000            | 0.0%                   |
| 0.3            | 1.048            | 4.8%                   |
| 0.5            | 1.155            | 15.5%                  |
| 0.7            | 1.400            | 40.0%                  |

## Implementation Tasks

### 1. Configuration (✅ DONE)
- [x] Add `boost_velocities` array to DiracInitialCondition
- [x] Parse YAML config for boost parameters
- [x] Create `config/relativistic_mass_validation.yaml`

### 2. Boosted Spinor Initialization (TODO)

Need to implement in `DiracEvolution::initializeGaussian()`:

**Boosted Dirac Spinor**:
A relativistic boost in the x-direction transforms a spinor:
```
ψ'(x) = exp(-iβ·α_x/2) ψ(x)
```
where:
- β = v/c (velocity in units of c)
- α_x = Dirac alpha matrix in x-direction
- For 2+1D Dirac: α_x = σ_x ⊗ I₂

**Implementation**:
```cpp
void DiracEvolution::initializeBoostedGaussian(
    float x0, float y0, float sigma,
    float boost_vx, float boost_vy)
{
    // 1. Initialize rest-frame Gaussian (existing code)
    initializeGaussian(x0, y0, sigma);

    // 2. Apply Lorentz boost operator
    //    ψ_boost = exp(-i β_x α_x / 2 - i β_y α_y / 2) ψ_rest

    for (int iy = 0; iy < _Ny; ++iy) {
        for (int ix = 0; ix < _Nx; ++ix) {
            int idx = iy * _Nx + ix;

            // Extract 4-component spinor
            std::complex<float> psi[4] = {
                spinor_data[4*idx + 0],
                spinor_data[4*idx + 1],
                spinor_data[4*idx + 2],
                spinor_data[4*idx + 3]
            };

            // Apply boost: exp(-i β·α/2) ≈ cos(β/2) - i(β·α)sin(β/2)
            float beta = sqrt(boost_vx*boost_vx + boost_vy*boost_vy);
            float cos_half = cos(beta/2);
            float sin_half = sin(beta/2);

            // Boost direction
            float nx = (beta > 1e-6) ? boost_vx/beta : 0;
            float ny = (beta > 1e-6) ? boost_vy/beta : 0;

            // Apply boost transformation (2+1D Dirac)
            // ... matrix multiplication ...

            // Write back
            spinor_data[4*idx + 0] = psi[0];
            // ...
        }
    }
}
```

### 3. Test Runner Modifications (TODO)

Similar to `grid_sizes`, need velocity sweep in `SMFTTestRunner`:

```cpp
bool SMFTTestRunner::run() {
    // ...

    // Check if velocity sweep is enabled
    if (!_config.dirac_initial.boost_velocities.empty()) {
        std::cout << "===== Velocity Sweep Testing =====" << std::endl;

        for (float v : _config.dirac_initial.boost_velocities) {
            std::cout << "Testing velocity v = " << v << "c" << std::endl;

            // Run test with this velocity
            runForVelocity(v);
        }

        // Validate γ factors
        validateRelativisticMass();
    }
}

bool SMFTTestRunner::runForVelocity(float v) {
    // Set boost velocity
    _current_boost_velocity = v;

    // Re-initialize Dirac field with boost
    // ... call initializeBoostedGaussian ...

    // Run dynamics
    runSingleTest(N);

    // Measure effective mass from trajectory
    float m_eff = computeEffectiveMass();

    // Store result
    _mass_measurements[v] = m_eff;
}
```

### 4. Effective Mass Measurement (TODO)

Extract mass from particle trajectory using Newton's 2nd law:

```cpp
float ObservableComputer::computeEffectiveMass(
    const std::vector<Observables>& trajectory)
{
    // Extract center-of-mass position over time
    std::vector<float> x_cm, y_cm;
    for (const auto& obs : trajectory) {
        x_cm.push_back(obs.x_center_of_mass);
        y_cm.push_back(obs.y_center_of_mass);
    }

    // Compute velocity: v = dx/dt
    std::vector<float> vx, vy;
    for (size_t i = 1; i < x_cm.size(); ++i) {
        vx.push_back((x_cm[i] - x_cm[i-1]) / dt);
        vy.push_back((y_cm[i] - y_cm[i-1]) / dt);
    }

    // Compute acceleration: a = dv/dt
    std::vector<float> ax, ay;
    for (size_t i = 1; i < vx.size(); ++i) {
        ax.push_back((vx[i] - vx[i-1]) / dt);
        ay.push_back((vy[i] - vy[i-1]) / dt);
    }

    // Force from Kuramoto field: F = -<β>Δ∇R
    std::vector<float> Fx, Fy;
    // ... extract from force_x, force_y observables ...

    // m_eff = |F| / |a|  (Newton's 2nd law)
    float m_eff_avg = 0.0f;
    int count = 0;
    for (size_t i = 0; i < Fx.size(); ++i) {
        float F_mag = sqrt(Fx[i]*Fx[i] + Fy[i]*Fy[i]);
        float a_mag = sqrt(ax[i]*ax[i] + ay[i]*ay[i]);

        if (a_mag > 1e-6 && F_mag > 1e-6) {
            m_eff_avg += F_mag / a_mag;
            count++;
        }
    }

    return (count > 0) ? m_eff_avg / count : 0.0f;
}
```

### 5. Relativistic Validation (TODO)

```cpp
bool SMFTTestRunner::validateRelativisticMass() {
    bool passed = true;

    // Extract rest mass (v = 0)
    float m_rest = _mass_measurements[0.0];

    std::cout << "\n===== Relativistic Mass Validation =====" << std::endl;
    std::cout << "Rest mass m₀ = " << m_rest << std::endl;

    for (const auto& [v, m_eff] : _mass_measurements) {
        if (v < 1e-6) continue;  // Skip v=0

        // Compute expected γ factor
        float gamma_expected = 1.0f / sqrt(1.0f - v*v);

        // Compute observed γ factor
        float gamma_observed = m_eff / m_rest;

        // Check agreement
        float error = fabs(gamma_observed - gamma_expected) / gamma_expected;
        bool test_passed = error < _config.validation.gamma_tolerance;

        std::cout << "v = " << v << "c:" << std::endl;
        std::cout << "  γ_expected = " << gamma_expected << std::endl;
        std::cout << "  γ_observed = " << gamma_observed << std::endl;
        std::cout << "  error = " << (error*100) << "%" << std::endl;
        std::cout << "  " << (test_passed ? "✓ PASS" : "✗ FAIL") << std::endl;

        passed &= test_passed;
    }

    return passed;
}
```

### 6. ObservableComputer Enhancements (TODO)

Add new observables:
- `x_center_of_mass`, `y_center_of_mass` (track particle position)
- `vx_average`, `vy_average` (average velocity)
- `effective_mass` (measured from F=ma)
- `gamma_factor` (m_eff / m_rest)
- `lorentz_invariant` (E² - p²c²)

## Expected Results

For each velocity test:
1. **Trajectory**: Particle should move with constant velocity (zero force background)
2. **Effective Mass**: Measured from trajectory should match m = γ·m₀
3. **Energy-Momentum**: E² - p² = m₀² should be invariant
4. **Validation**: γ_observed / γ_expected < 5% error

## Success Criteria

- ✓ Grid convergence: Results independent of grid resolution
- ✓ N convergence: Results independent of operator splitting ratio
- ✓ γ factor agreement: |γ_obs - γ_theory| / γ_theory < 5%
- ✓ Lorentz invariance: E² - p² = constant ± 1%

## Next Steps

1. Implement `initializeBoostedGaussian()` in DiracEvolution
2. Add center-of-mass tracking to ObservableComputer
3. Implement velocity loop in SMFTTestRunner (similar to grid_sizes)
4. Add effective mass computation
5. Add γ factor validation
6. Run test with v = [0.0, 0.3, 0.5, 0.7]c
7. Verify time dilation emerges automatically from Dirac dynamics
