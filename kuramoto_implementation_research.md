# Kuramoto Model Implementation Research Report

## Executive Summary
This report documents best practices for implementing the Kuramoto model in Python, focusing on numerical methods, validation approaches, and performance considerations for systems of 100-1000 oscillators.

## 1. Reference Implementations

### 1.1 fabridamicelli/kuramoto (GitHub)
- **Features**: Network-based implementation with adjacency matrix support
- **Integration**: Unspecified method, dt parameter exposed
- **Validation**: Matches English 2008 numerical results, validates Kc = √(8/π) × σ(ω)
- **System sizes**: Tested with 100-500 nodes
- **URL**: https://github.com/fabridamicelli/kuramoto

### 1.2 laszukdawid/Dynamical-systems (GitHub)
- **Features**: Multi-harmonic coupling, Jacobian support
- **Integration**: scipy.integrate.ode with "dopri5" (RK4/5 Dormand-Prince)
- **Timestep**: dt=0.05 typical
- **Special**: Includes Jacobian for improved stability
- **URL**: https://github.com/laszukdawid/Dynamical-systems

### 1.3 cbnfreitas/kuramoto_model_integrate_and_plot (GitHub)
- **Features**: Supports delay differential equations
- **Integration**: odeint (LSODA) for ODEs, ddeint for DDEs
- **Visualization**: Built-in plotting capabilities
- **URL**: https://github.com/cbnfreitas/kuramoto_model_integrate_and_plot

### 1.4 ritobandatta/kuramoto_model (GitHub)
- **Features**: Graph-based implementation
- **Focus**: Network topology effects
- **URL**: https://github.com/ritobandatta/kuramoto_model

### 1.5 netrd Library
- **Features**: Part of network reconstruction toolkit
- **Integration**: Professional implementation with documentation
- **URL**: https://netrd.readthedocs.io

## 2. Recommended Numerical Solvers

### 2.1 Performance Comparison (Kuramoto Oscillators)
Based on empirical benchmarks for random Kuramoto networks:
1. **odeint (LSODA)**: 7.76s - Most accurate, adaptive timestep
2. **ode with dopri5**: 6.66s - Fastest, good accuracy
3. **solve_ivp (RK45)**: 9.19s - Modern API, higher overhead

### 2.2 Solver Recommendations

#### For general use (N < 1000):
```python
from scipy.integrate import odeint

def kuramoto_deriv(theta, t, omega, K, N):
    """Kuramoto model derivative function."""
    theta_mat = theta[:, None] - theta
    return omega + (K/N) * np.sum(np.sin(theta_mat), axis=1)

# Use odeint with LSODA (adaptive timestep)
solution = odeint(kuramoto_deriv, theta0, t_span, args=(omega, K, N))
```

#### For high precision:
```python
from scipy.integrate import solve_ivp

# Use DOP853 for high accuracy requirements
sol = solve_ivp(kuramoto_rhs, [t0, tf], theta0,
                method='DOP853', rtol=1e-10, atol=1e-12)
```

#### For stiff problems:
```python
# Use implicit methods for stiff systems
sol = solve_ivp(kuramoto_rhs, [t0, tf], theta0, method='Radau')
```

### 2.3 Timestep Requirements
- **Typical values**: dt = 0.01 to 0.1
- **Adaptive solvers**: Let solver choose (recommended)
- **Fixed timestep**: dt < 0.1 × min(1/K, 1/max(ω))
- **For phase transitions**: Use smaller dt near Kc

## 3. Standard Validation Tests

### 3.1 Order Parameter Tests
```python
def order_parameter(theta):
    """Calculate Kuramoto order parameter."""
    z = np.mean(np.exp(1j * theta))
    return np.abs(z), np.angle(z)

# Test cases:
# 1. K=0: r should decay to ~1/√N (incoherent)
# 2. K>>Kc: r should approach 1 (synchronized)
# 3. K=Kc: r should show critical behavior
```

### 3.2 Critical Coupling Values

#### Lorentzian Distribution
- **Analytical**: Kc = 2γ (γ = half-width at half-maximum)
- **Test**: For γ=1, Kc=2 exactly

#### Gaussian Distribution
- **Semi-analytical**: Kc ≈ √(8/π) × σ
- **Test**: For σ=1, Kc ≈ 1.596

#### Uniform Distribution
- **Numerical**: Must determine via simulation
- **Test**: Width w, Kc ≈ 1.2 × w

### 3.3 Benchmark Tests

1. **Synchronization transition**:
   - Sweep K from 0 to 2Kc
   - Verify r(K) shows expected S-curve
   - Check hysteresis (forward vs backward sweep)

2. **Steady-state convergence**:
   - Run until dr/dt < 1e-6
   - Typical time: T > 100/K

3. **Finite-size scaling**:
   - Test N = [50, 100, 200, 500, 1000]
   - Verify r_steady scales as expected

## 4. Performance Considerations

### 4.1 Computational Complexity
- **Per timestep**: O(N²) for all-to-all coupling
- **Memory**: O(N²) for coupling matrix, O(N) for phases
- **Optimization**: Vectorize using NumPy broadcasting

### 4.2 Optimized Implementation
```python
def kuramoto_deriv_optimized(theta, t, omega, K, N):
    """Vectorized Kuramoto derivative."""
    # Use broadcasting to avoid explicit loops
    theta_diff = theta[:, np.newaxis] - theta
    coupling = np.sum(np.sin(theta_diff), axis=1)
    return omega + (K/N) * coupling
```

### 4.3 Performance Tips for N=100-1000
1. **Vectorization**: Use NumPy broadcasting, avoid Python loops
2. **Memory**: Pre-allocate arrays, reuse buffers
3. **Solver choice**: Use odeint for speed, DOP853 for accuracy
4. **Parallelization**: Consider parallel parameter sweeps
5. **JIT compilation**: Use Numba for derivative function

### 4.4 Expected Performance
- **N=100**: ~0.1s per 100 time units
- **N=500**: ~1-2s per 100 time units
- **N=1000**: ~5-10s per 100 time units
(Using vectorized NumPy on modern CPU)

## 5. Common Pitfalls and Solutions

### 5.1 Numerical Stability Issues

**Problem**: Phase wrapping causes discontinuities
**Solution**: Use modulo arithmetic or unwrap phases
```python
theta = np.mod(theta + np.pi, 2*np.pi) - np.pi
```

**Problem**: Stiff dynamics near synchronization
**Solution**: Use adaptive timestep or implicit methods

**Problem**: Initial transients affect measurements
**Solution**: Discard first T=100/K time units

### 5.2 Implementation Pitfalls

1. **Incorrect mean-field normalization**: Always divide by N
2. **Phase initialization**: Use random uniform on [-π, π]
3. **Frequency distribution**: Center at ω₀=0 for simplicity
4. **Boundary conditions**: Handle phase wrapping correctly

### 5.3 Analysis Pitfalls

1. **Finite-size effects**: Small N shows fluctuations
2. **Hysteresis**: Different results for increasing/decreasing K
3. **Metastable states**: System may get stuck in local minima
4. **Measurement timing**: Wait for steady state before measuring

## 6. Validation Checklist

- [ ] Order parameter r→0 as K→0
- [ ] Order parameter r→1 as K→∞
- [ ] Critical coupling matches theory for known distributions
- [ ] Phase coherence visualizations show expected patterns
- [ ] Computational time scales as O(N²)
- [ ] Results stable under timestep refinement
- [ ] Steady state reached before measurements
- [ ] Multiple random seeds give consistent statistics

## 7. Recommended Implementation Structure

```python
class KuramotoModel:
    def __init__(self, N, coupling, freq_dist='gaussian', freq_width=1.0):
        self.N = N
        self.K = coupling
        self.omega = self._init_frequencies(freq_dist, freq_width)
        self.theta = np.random.uniform(-np.pi, np.pi, N)

    def derivative(self, theta, t):
        theta_diff = theta[:, np.newaxis] - theta
        return self.omega + (self.K/self.N) * np.sum(np.sin(theta_diff), axis=1)

    def integrate(self, t_span, method='odeint'):
        if method == 'odeint':
            return odeint(self.derivative, self.theta, t_span)
        else:
            sol = solve_ivp(lambda t, y: self.derivative(y, t),
                          [t_span[0], t_span[-1]], self.theta,
                          t_eval=t_span, method=method)
            return sol.y.T

    def order_parameter(self, theta=None):
        if theta is None:
            theta = self.theta
        z = np.mean(np.exp(1j * theta), axis=-1)
        return np.abs(z), np.angle(z)
```

## Sources

1. Kuramoto, Y. (1984). Chemical Oscillations, Waves, and Turbulence
2. Strogatz, S. H. (2000). From Kuramoto to Crawford
3. Acebrón et al. (2005). The Kuramoto model: A simple paradigm
4. GitHub repositories: fabridamicelli, laszukdawid, cbnfreitas
5. SciPy documentation for integration methods
6. Research papers on critical coupling and synchronization transitions