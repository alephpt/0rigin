# Kuramoto Model: Python Implementation

A comprehensive Python library for simulating and analyzing the Kuramoto model of coupled phase oscillators.

## Overview

The Kuramoto model describes the synchronization dynamics of a population of coupled oscillators:

```
dθ_i/dt = ω_i + (K/N) Σ_j sin(θ_j - θ_i)
```

where:
- `θ_i` is the phase of oscillator `i`
- `ω_i` is its natural frequency
- `K` is the coupling strength
- `N` is the number of oscillators

This implementation provides a clean, extensible framework for studying synchronization phenomena.

## Features

- **Multiple frequency distributions**: Lorentzian, Gaussian, Uniform
- **Flexible coupling schemes**: All-to-all, network-based
- **Advanced solvers**: RK4, RK45 (adaptive), Euler
- **Comprehensive analysis**: Order parameter, phase coherence, convergence detection
- **Publication-ready visualizations**: Phase circles, bifurcation diagrams, time series
- **Clean architecture**: Modular design, type hints, comprehensive docstrings

## Installation

### Requirements

- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Matplotlib >= 3.3 (for visualization)

### Install dependencies

```bash
pip install -r requirements.txt
```

### Setup

Add the source directory to your Python path:

```python
import sys
sys.path.insert(0, '/path/to/0rigin/src')
import kuramoto as km
```

## Quick Start

### Basic Example

```python
import kuramoto as km
import numpy as np
import matplotlib.pyplot as plt

# Create model with Lorentzian distribution
model = km.KuramotoModel(
    N=100,                    # 100 oscillators
    coupling=3.0,             # Coupling strength K
    frequencies='lorentzian'  # Frequency distribution
)

# Simulate
solution = model.evolve(t_span=(0, 50), solver='rk45')

# Plot order parameter
plt.plot(solution['t'], solution['R'])
plt.xlabel('Time')
plt.ylabel('Order Parameter R')
plt.show()
```

### Using Custom Distributions

```python
from kuramoto.distributions import GaussianDistribution

# Create custom distribution
dist = GaussianDistribution(mean=0, std=1.0)
print(f"Critical coupling: Kc ≈ {dist.critical_coupling():.3f}")

# Create model
model = km.KuramotoModel(N=200, coupling=2.5, frequencies=dist)
solution = model.evolve((0, 100))
```

### Order Parameter Analysis

```python
from kuramoto.analysis import OrderParameter

# Analyze synchronization
op = OrderParameter(solution['phases'])

print(f"Mean R: {op.mean_amplitude():.3f}")
print(f"Steady-state R: {op.steady_state_amplitude():.3f}")
print(f"Converged: {op.is_synchronized()}")
print(f"Convergence time: {op.convergence_time(solution['t']):.2f}")
```

### Visualization

```python
from kuramoto.visualization import (
    plot_phases,
    plot_order_parameter,
    plot_bifurcation_diagram
)

# Phase distribution on unit circle
final_phases = solution['phases'][-1, :]
plot_phases(final_phases, model.frequencies)
plt.show()

# Order parameter evolution
plot_order_parameter(solution['t'], solution['R'], R_theory=0.7)
plt.show()
```

## Examples

Run the complete demonstration:

```bash
cd examples
python demo_synchronization.py
```

This will generate:
- Synchronization regime comparison (subcritical, critical, supercritical)
- Distribution comparison (Lorentzian, Gaussian, Uniform)
- Bifurcation diagram (R vs K)
- Phase visualization on unit circle

## API Documentation

### Core Classes

#### `KuramotoModel`

Main simulation class.

```python
model = km.KuramotoModel(
    N=100,                        # Number of oscillators
    coupling=2.0,                 # Coupling strength or Coupling object
    frequencies='lorentzian',     # Distribution name, object, or array
    initial_phases=None           # Initial phases (random if None)
)

# Simulate
solution = model.evolve(
    t_span=(0, 50),              # Time interval
    solver='rk45',               # Solver: 'rk45', 'rk4', 'euler'
    dt=None,                     # Time step (for fixed-step solvers)
    store_trajectory=True        # Store full trajectory
)

# Returns dict with keys: 't', 'phases', 'R', 'Psi'
```

### Distributions

All distributions inherit from `FrequencyDistribution` and implement:
- `sample(N, seed)`: Generate N frequencies
- `pdf(omega)`: Probability density function
- `critical_coupling()`: Analytical Kc (if known)
- `mean()`, `variance()`: Distribution statistics

#### `LorentzianDistribution`

```python
from kuramoto.distributions import LorentzianDistribution

dist = LorentzianDistribution(center=0, width=1.0)
Kc = dist.critical_coupling()  # Returns 2.0 (exact)
R_steady = dist.steady_state_order_parameter(K=3.0)
```

#### `GaussianDistribution`

```python
from kuramoto.distributions import GaussianDistribution

dist = GaussianDistribution(mean=0, std=1.0)
Kc = dist.critical_coupling()  # Returns ~1.596 (approximate)
```

#### `UniformDistribution`

```python
from kuramoto.distributions import UniformDistribution

dist = UniformDistribution(low=-1, high=1)
# No analytical Kc available
```

### Analysis Tools

#### `OrderParameter`

Comprehensive order parameter analysis.

```python
from kuramoto.analysis import OrderParameter

op = OrderParameter(phases)  # phases: (n_times, n_oscillators)

# Compute order parameter
R, Psi = op.time_series()

# Analysis metrics
mean_R = op.mean_amplitude()
steady_R = op.steady_state_amplitude(fraction=0.2)
is_sync = op.is_synchronized(threshold=0.5)
conv_time = op.convergence_time(t, threshold=0.01)
locked = op.phase_locked_fraction(tolerance=0.1)
```

#### Synchronization Metrics

```python
from kuramoto.analysis import (
    phase_coherence,
    phase_variance,
    metastability,
    frequency_entrainment
)

rho = phase_coherence(phases)
var = phase_variance(phases, circular=True)
M = metastability(R_timeseries)
entrainment = frequency_entrainment(frequencies)
```

### Visualization

All plotting functions return matplotlib Axes objects for customization.

```python
from kuramoto.visualization import (
    plot_phases,              # Phase circle
    plot_order_parameter,     # R(t) time series
    plot_bifurcation_diagram, # R vs K
    plot_phase_evolution      # Phase trajectories
)
```

## Theory

### Critical Coupling

The Kuramoto model exhibits a continuous phase transition from incoherent to partially synchronized states at critical coupling `Kc`.

**Lorentzian distribution**: `Kc = 2γ` (exact, from Ott-Antonsen theory)

**Gaussian distribution**: `Kc ≈ √(8/π) σ ≈ 1.596σ` (approximate)

**Uniform distribution**: No analytical formula (numerical determination required)

### Order Parameter

The complex order parameter `Z = R e^(iΨ)` measures synchronization:

```
Z = (1/N) Σ_j e^(iθ_j)
```

- `R = |Z|`: Synchronization amplitude (0 = incoherent, 1 = fully synchronized)
- `Ψ = arg(Z)`: Mean phase

### Synchronization Regimes

1. **Subcritical** (`K < Kc`): `R → 0` (incoherent)
2. **Critical** (`K ≈ Kc`): Fluctuating order parameter
3. **Supercritical** (`K > Kc`): `R → R_∞ > 0` (partial synchronization)

Near the critical point: `R ∝ √(K - Kc)`

## Architecture

```
src/kuramoto/
├── __init__.py              # Main exports
├── core/
│   ├── model.py             # KuramotoModel class
│   └── coupling.py          # Coupling schemes
├── distributions/
│   ├── base.py              # Abstract base class
│   ├── lorentzian.py        # Lorentzian distribution
│   ├── gaussian.py          # Gaussian distribution
│   └── uniform.py           # Uniform distribution
├── solvers/
│   ├── base.py              # Solver interface
│   └── runge_kutta.py       # RK4, RK45, Euler
├── analysis/
│   ├── order_parameter.py   # OrderParameter class
│   └── metrics.py           # Additional metrics
└── visualization/
    ├── phase_plot.py        # Phase visualizations
    ├── time_series.py       # Time series plots
    └── bifurcation.py       # Bifurcation diagrams
```

## Code Quality

- **PEP 8 compliant**: Clean, readable code
- **Type hints**: Full typing support throughout
- **Docstrings**: NumPy-style documentation for all public APIs
- **Modular design**: Clear separation of concerns
- **Extensible**: Easy to add new distributions, solvers, analyses

## Performance

- NumPy vectorization for efficient computation
- Adaptive RK45 solver for automatic step-size control
- Optimized order parameter calculation: O(N) complexity

## References

1. Kuramoto, Y. (1984). Chemical Oscillations, Waves, and Turbulence
2. Strogatz, S. H. (2000). From Kuramoto to Crawford
3. Ott, E., & Antonsen, T. M. (2008). Low dimensional behavior of large systems

## Contributing

This implementation follows strict development standards:
- Files < 500 lines
- Functions < 50 lines
- Comprehensive error handling
- No hardcoded values
- Full test coverage

## License

Research use - 0rigin Project

## Contact

0rigin Research Team
