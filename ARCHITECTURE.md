# Kuramoto Model Simulation Architecture

## Overview

This document describes the architecture of the implemented Python-based Kuramoto model simulation framework with statistical mean field theory (SMFT) extensions. The system supports both classical Kuramoto synchronization dynamics and field-theoretic formulations with Hamiltonian dynamics.

## Design Principles

1. **Modularity**: Clear separation between model definition, numerical solvers, analysis, and visualization
2. **Extensibility**: Easy to add new frequency distributions, coupling functions, and field theory components
3. **Performance**: Efficient numerical computation with optional parallelization
4. **Testability**: Comprehensive test coverage with validation against analytical results
5. **Scientific Rigor**: Accurate numerical methods with error control

## Project Structure

```
0rigin/
├── README.md                        # User documentation
├── ARCHITECTURE.md                  # This document
├── CONTRIBUTING.md                  # Development guidelines
├── requirements.txt                 # Dependencies
│
├── src/
│   └── kuramoto/
│       ├── __init__.py               # Package exports
│       │
│       ├── core/
│       │   ├── __init__.py
│       │   ├── model.py             # KuramotoModel class
│       │   └── coupling.py          # Coupling schemes
│       │
│       ├── distributions/
│       │   ├── __init__.py
│       │   ├── base.py              # FrequencyDistribution abstract base
│       │   ├── lorentzian.py        # LorentzianDistribution
│       │   ├── gaussian.py          # GaussianDistribution
│       │   └── uniform.py           # UniformDistribution
│       │
│       ├── solvers/
│       │   ├── __init__.py
│       │   ├── base.py              # Solver interface
│       │   └── runge_kutta.py       # RK4, RK45, Euler implementations
│       │
│       ├── analysis/
│       │   ├── __init__.py
│       │   ├── order_parameter.py   # OrderParameter class
│       │   └── metrics.py           # Synchronization metrics
│       │
│       ├── visualization/
│       │   ├── __init__.py
│       │   ├── phase_plot.py        # Phase circle plots
│       │   ├── time_series.py       # Time evolution plots
│       │   ├── bifurcation.py       # Bifurcation diagrams
│       │   └── utils.py             # Plotting utilities
│       │
│       └── field_theory/            # SMFT extensions (NEW)
│           ├── __init__.py
│           ├── hamiltonian/
│           │   ├── __init__.py
│           │   └── phase_space.py   # Hamiltonian phase space dynamics
│           ├── fields/
│           │   ├── __init__.py
│           │   ├── grid.py          # Spatial grid discretization
│           │   ├── scalar_field.py  # ScalarField base class
│           │   └── mediator.py      # MediatorField (Klein-Gordon)
│           ├── coupling/
│           │   ├── __init__.py
│           │   └── local_coupling.py # Local field-oscillator coupling
│           ├── pde_solvers/
│           │   ├── __init__.py
│           │   └── base.py          # PDE solver interface (py-pde integration)
│           └── smft_system.py       # High-level SMFT system class
│
├── tests/
│   ├── conftest.py                  # Pytest configuration
│   ├── test_model.py                # Core model tests
│   ├── test_coupling.py             # Coupling tests
│   ├── test_solvers.py              # Solver validation
│   ├── test_distributions.py        # Distribution tests
│   ├── test_order_parameter.py      # Order parameter tests
│   ├── test_metrics.py              # Analysis metrics tests
│   ├── test_phase_plot.py           # Visualization tests
│   ├── test_time_series.py
│   ├── test_bifurcation.py
│   └── field_theory/                # Field theory tests
│       ├── test_hamiltonian.py
│       ├── test_fields.py
│       ├── test_mediator.py
│       ├── test_local_coupling.py
│       └── test_smft_system.py
│
├── examples/
│   ├── demo_synchronization.py      # Classical Kuramoto demo
│   └── field_theory/
│       ├── smft_demo.py             # SMFT demonstration
│       └── hamiltonian_demo.py      # Hamiltonian dynamics demo
│
└── docs/
    └── validation/
        └── VALIDATION_REPORT.md     # Scientific validation results
```

## Module Architecture

### Core Module (`kuramoto.core`)

```python
# model.py
class KuramotoModel:
    """
    Main Kuramoto model class.

    Parameters
    ----------
    N : int
        Number of oscillators
    coupling : float or Coupling
        Coupling strength K or coupling object
    frequencies : array-like or Distribution
        Natural frequencies ω_i
    initial_phases : array-like, optional
        Initial phase configuration
    """

    def __init__(self, N, coupling, frequencies, initial_phases=None):
        pass

    def equations_of_motion(self, t, phases):
        """Return dθ/dt for current phase configuration."""
        pass

    def evolve(self, t_span, solver='rk45', **solver_kwargs):
        """Evolve system over time span."""
        pass

    def compute_order_parameter(self, phases=None):
        """Calculate Kuramoto order parameter R, Ψ."""
        pass

# oscillator.py
class Oscillator:
    """Single oscillator with phase and frequency."""

    def __init__(self, frequency, initial_phase=None):
        self.frequency = frequency
        self.phase = initial_phase or np.random.uniform(0, 2*np.pi)

    def update(self, coupling_field, dt):
        """Update phase based on coupling field."""
        pass

# coupling.py
class Coupling(ABC):
    """Abstract base for coupling functions."""

    @abstractmethod
    def compute_field(self, phases):
        """Compute coupling field from phase configuration."""
        pass

class SinusoidalCoupling(Coupling):
    """Standard sin(θ_j - θ_i) coupling."""

    def __init__(self, strength):
        self.K = strength

    def compute_field(self, phases):
        """Return coupling field for each oscillator."""
        pass

class NetworkCoupling(Coupling):
    """Coupling on arbitrary network topology."""

    def __init__(self, adjacency_matrix, strength):
        self.A = adjacency_matrix
        self.K = strength
```

### Distributions Module (`kuramoto.distributions`)

```python
# base.py
class FrequencyDistribution(ABC):
    """Abstract base for frequency distributions."""

    @abstractmethod
    def sample(self, N):
        """Sample N frequencies from distribution."""
        pass

    @abstractmethod
    def pdf(self, omega):
        """Probability density function."""
        pass

    @abstractmethod
    def critical_coupling(self):
        """Analytical Kc if known."""
        pass

# lorentzian.py
class LorentzianDistribution(FrequencyDistribution):
    """
    Lorentzian (Cauchy) distribution.

    g(ω) = γ/π / [(ω - ω₀)² + γ²]

    Parameters
    ----------
    center : float
        Center frequency ω₀
    width : float
        Half-width γ
    """

    def __init__(self, center=0, width=1):
        self.omega_0 = center
        self.gamma = width

    def critical_coupling(self):
        """Return Kc = 2γ."""
        return 2 * self.gamma
```

### Solvers Module (`kuramoto.solvers`)

```python
# base.py
class Solver(ABC):
    """Abstract base for ODE solvers."""

    @abstractmethod
    def integrate(self, func, y0, t_span, **kwargs):
        """Integrate ODE system."""
        pass

    @abstractmethod
    def step(self, func, t, y, dt):
        """Single integration step."""
        pass

# runge_kutta.py
class RK4Solver(Solver):
    """Fourth-order Runge-Kutta solver."""

    def step(self, func, t, y, dt):
        """RK4 step implementation."""
        k1 = func(t, y)
        k2 = func(t + dt/2, y + dt*k1/2)
        k3 = func(t + dt/2, y + dt*k2/2)
        k4 = func(t + dt, y + dt*k3)
        return y + dt*(k1 + 2*k2 + 2*k3 + k4)/6

class RK45Solver(Solver):
    """Adaptive Runge-Kutta-Fehlberg solver."""

    def __init__(self, rtol=1e-6, atol=1e-9):
        self.rtol = rtol
        self.atol = atol
```

### Analysis Module (`kuramoto.analysis`)

```python
# order_parameter.py
class OrderParameter:
    """
    Kuramoto order parameter calculator.

    R e^(iΨ) = (1/N) Σ e^(iθ_j)
    """

    @staticmethod
    def compute(phases):
        """
        Compute R and Ψ from phase configuration.

        Returns
        -------
        R : float
            Synchronization amplitude [0, 1]
        Psi : float
            Mean phase
        """
        z = np.mean(np.exp(1j * phases))
        return np.abs(z), np.angle(z)

    @staticmethod
    def time_series(solution):
        """Compute R(t) from solution trajectory."""
        pass

# stability.py
class LinearStability:
    """Linear stability analysis around fixed points."""

    def __init__(self, model):
        self.model = model

    def find_fixed_points(self):
        """Find synchronized and incoherent states."""
        pass

    def compute_jacobian(self, state):
        """Linearization around state."""
        pass

    def analyze_stability(self, fixed_point):
        """Eigenvalue analysis of Jacobian."""
        pass

# bifurcation.py
class BifurcationAnalysis:
    """Track synchronization transition vs coupling."""

    def scan_coupling(self, K_range, N=1000, samples=100):
        """Compute R(K) bifurcation diagram."""
        pass

    def find_critical_coupling(self, tolerance=1e-3):
        """Numerical determination of Kc."""
        pass
```

### Visualization Module (`kuramoto.visualization`)

```python
# phase_plot.py
class PhasePlot:
    """Visualize oscillators on unit circle."""

    def __init__(self, fig=None, ax=None):
        self.fig = fig or plt.figure()
        self.ax = ax or self.fig.add_subplot(111, projection='polar')

    def plot_snapshot(self, phases, colors=None):
        """Plot current phase configuration."""
        pass

    def plot_trajectory(self, solution, oscillator_idx=None):
        """Plot phase evolution over time."""
        pass

# order_plot.py
class OrderParameterPlot:
    """Visualize order parameter evolution."""

    def plot_time_series(self, t, R, Psi=None):
        """Plot R(t) and optionally Ψ(t)."""
        pass

    def plot_bifurcation(self, K_values, R_values):
        """Bifurcation diagram R vs K."""
        pass
```

## API Design

### Basic Usage Pattern

```python
import kuramoto as km

# Create model
model = km.KuramotoModel(
    N=100,
    coupling=km.SinusoidalCoupling(K=2.0),
    frequencies=km.LorentzianDistribution(center=0, width=1)
)

# Simulate
solution = model.evolve(t_span=(0, 100), solver='rk45')

# Analyze
R, Psi = km.analysis.OrderParameter.time_series(solution)

# Visualize
km.visualization.plot_synchronization(solution)
```

### Advanced Usage

```python
# Custom frequency distribution
class BimodalDistribution(km.FrequencyDistribution):
    def sample(self, N):
        # Custom implementation
        pass

# Network coupling
import networkx as nx
G = nx.watts_strogatz_graph(100, 4, 0.3)
coupling = km.NetworkCoupling(nx.adjacency_matrix(G), K=5.0)

# Parallel simulation
with km.parallel_context(n_workers=4):
    results = km.parameter_scan(
        K_range=np.linspace(0, 10, 100),
        N_values=[100, 500, 1000]
    )
```

## Testing Strategy

### Unit Tests
- Each class method tested independently
- Mock objects for dependencies
- Property-based testing for invariants

### Integration Tests
- End-to-end simulation workflows
- Solver accuracy validation
- Performance benchmarks

### Validation Tests
- Known analytical solutions (uniform distribution)
- Ott-Antonsen reduction (Lorentzian)
- Critical coupling values
- Asymptotic behavior (N → ∞)

### Test Examples

```python
# tests/test_model.py
def test_order_parameter_bounds():
    """R must be in [0, 1]."""
    model = create_test_model()
    R, _ = model.compute_order_parameter()
    assert 0 <= R <= 1

def test_synchronization_transition():
    """Verify R → 0 for K < Kc, R > 0 for K > Kc."""
    pass

def test_ott_antonsen_lorentzian():
    """Compare with exact OA solution."""
    pass

# tests/test_solvers.py
def test_rk4_accuracy():
    """Verify 4th order convergence."""
    pass

def test_energy_conservation():
    """For conservative variant, check H(t) = const."""
    pass
```

### Field Theory Module (`kuramoto.field_theory`)

**Status**: Implemented SMFT extensions with Hamiltonian dynamics and Klein-Gordon mediator fields.

```python
# field_theory/
├── hamiltonian/
│   └── phase_space.py          # Hamiltonian formulation of Kuramoto dynamics
├── fields/
│   ├── grid.py                 # Spatial grid discretization
│   ├── scalar_field.py         # ScalarField abstract base class
│   └── mediator.py             # MediatorField (Klein-Gordon PDE)
├── coupling/
│   └── local_coupling.py       # Local oscillator-field coupling
├── pde_solvers/
│   └── base.py                 # PDE solver interface (py-pde integration)
└── smft_system.py              # SMFTSystem high-level interface
```

#### Key Components

**Hamiltonian Dynamics** (`hamiltonian/phase_space.py`):
- Phase space formulation: (θ, p) coordinates
- Conjugate momentum dynamics
- Energy conservation validation
- Hamiltonian H = Σ[p_i²/2 + ω_i·θ_i] - K·coupling_term

**Scalar Fields** (`fields/`):
- Grid-based spatial discretization
- Generic ScalarField base class
- Klein-Gordon mediator field with wave propagation
- Mass generation mechanism: m_eff ∝ R (order parameter coupling)

**Local Coupling** (`coupling/local_coupling.py`):
- Spatial position-dependent coupling
- Field interpolation at oscillator positions
- Self-consistent field updates

**SMFT System** (`smft_system.py`):
- High-level interface combining oscillators + field
- Coupled ODE-PDE evolution
- Order parameter calculation with spatial dependence

#### Capabilities

✅ **Implemented**:
- Hamiltonian phase space dynamics
- Klein-Gordon mediator field evolution
- Local position-based coupling
- Self-consistent mean field dynamics
- Mass generation via order parameter (m_eff ∝ R)
- Wave propagation at speed c
- Energy conservation in Hamiltonian formulation

⏳ **Future Extensions**:
1. **Full PDE Coupling**: Integrate with py-pde for advanced PDE solvers
2. **Quantum Corrections**: Path integral formulation
3. **Relativistic Covariance**: Lorentz-covariant formulation
4. **GPU Acceleration**: CuPy/JAX backends for large-scale simulations
5. **Machine Learning**: Neural ODE solvers, parameter inference

## Performance Considerations

### Computational Complexity
- Single timestep: O(N²) for all-to-all coupling
- Network coupling: O(E) where E = edges
- Order parameter: O(N)

### Optimization Strategies
- Vectorized NumPy operations
- Just-in-time compilation (Numba)
- Parallel parameter scans
- Sparse matrix for networks
- Adaptive timestepping

### Memory Management
- Lazy evaluation for large N
- Chunked storage for long simulations
- Optional trajectory decimation

## Code Standards

### Style Guide
- PEP 8 compliance
- Type hints for all public APIs
- NumPy-style docstrings
- Black formatting

### Documentation Requirements
- All public functions documented
- Mathematical formulas in docstrings
- Usage examples in module docstrings
- References to papers

### Error Handling
- Descriptive error messages
- Input validation with clear requirements
- Graceful degradation for edge cases
- Warnings for numerical issues

## Dependencies

### Core Requirements
```txt
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.4.0
```

### Optional Dependencies
```txt
numba>=0.54.0       # JIT compilation
networkx>=2.6       # Network coupling
pandas>=1.3.0       # Data analysis
h5py>=3.0.0        # HDF5 storage
tqdm>=4.62.0       # Progress bars
pytest>=6.2.0      # Testing
black>=21.0        # Formatting
mypy>=0.9.0        # Type checking
```

## Version Management

### Semantic Versioning
- MAJOR: Breaking API changes
- MINOR: New features, backwards compatible
- PATCH: Bug fixes

### Release Checklist
1. Run full test suite
2. Update documentation
3. Verify examples run
4. Tag release in git
5. Deploy to PyPI

## Summary

This architecture provides:
- Clean separation of concerns
- Easy extensibility for future physics
- Comprehensive testing framework
- Performance optimization paths
- Clear documentation structure
- Professional Python packaging

The design supports both immediate Kuramoto research needs and future extensions to field-theoretic formulations while maintaining scientific rigor and code quality.