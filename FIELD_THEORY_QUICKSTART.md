# Field Theory Quick Start Guide

## Installation

```bash
# Already installed with main package
import sys
sys.path.insert(0, '/path/to/0rigin/src')
from kuramoto.field_theory import SMFTSystem
```

## 5-Minute Tutorial

### 1. Create System

```python
from kuramoto.field_theory import SMFTSystem

system = SMFTSystem(
    grid_shape=(50, 50),      # Spatial grid
    N_oscillators=100,         # Discrete oscillators
    coupling='local',          # Local or global
    mediator_mass=10.0        # Field mass (M→∞ = Kuramoto)
)
```

### 2. Evolve System

```python
solution = system.evolve(
    t_span=(0, 30),           # Time span
    dt=0.01,                  # Time step
    store_interval=50         # Save every 50 steps
)
```

### 3. Access Results

```python
# Oscillator dynamics
times = solution['t']
phases = solution['theta']
R_global = solution['R']

# Field dynamics
R_field = solution['sync_field']        # R(x,y,t)
sigma_field = solution['mediator_field'] # σ(x,y,t)

# Derived quantities
m_eff = system.compute_effective_mass()  # m(x,y)
```

### 4. Visualize

```python
import matplotlib.pyplot as plt

# Order parameter
plt.plot(times, R_global)
plt.xlabel('Time')
plt.ylabel('R')
plt.show()

# Final field
plt.imshow(R_field[-1].T, origin='lower', cmap='viridis')
plt.colorbar(label='R(x,y)')
plt.show()
```

## Common Use Cases

### Compare Local vs Global Coupling

```python
# Local coupling
system_local = SMFTSystem(
    grid_shape=(40, 40),
    N_oscillators=100,
    coupling='local'
)
sol_local = system_local.evolve((0, 20))

# Global coupling
system_global = SMFTSystem(
    grid_shape=(40, 40),
    N_oscillators=100,
    coupling='global'
)
sol_global = system_global.evolve((0, 20))

print(f"Local R: {sol_local['R'][-1]:.3f}")
print(f"Global R: {sol_global['R'][-1]:.3f}")
```

### Study Heavy Mass Limit

```python
M_values = [1, 10, 100, 1000]

for M in M_values:
    system = SMFTSystem(
        grid_shape=(30, 30),
        N_oscillators=80,
        mediator_mass=M
    )
    sol = system.evolve((0, 15))
    print(f"M={M:4d}: R={sol['R'][-1]:.4f}")
```

### Access Individual Components

```python
# Direct grid access
grid = system.grid
laplacian = grid.laplacian(field_values)

# Oscillator phase space
oscillators = system.oscillators
energy = oscillators.compute_hamiltonian()
phases = oscillators.theta
momenta = oscillators.p

# Fields
mediator = system.mediator_field
sync = system.sync_field
phase_field = system.phase_field
```

## Integration with Sprint 1 Code

Field theory is 100% backward compatible:

```python
# Classic Kuramoto (Sprint 1)
from kuramoto import KuramotoModel

model = KuramotoModel(N=100, coupling=2.0, frequencies='lorentzian')
sol = model.evolve((0, 50))

# Field theory extension (Sprint 2)
from kuramoto.field_theory import SMFTSystem

system = SMFTSystem(
    grid_shape=(50, 50),
    N_oscillators=100,
    coupling='local',
    mediator_mass=10.0
)
sol = system.evolve((0, 50))

# Both work independently and together!
```

## Performance Tips

### Small System (Fast)
```python
system = SMFTSystem(
    grid_shape=(20, 20),    # Small grid
    N_oscillators=50,       # Few oscillators
    coupling='local'
)
# <100ms per step
```

### Medium System (Acceptable)
```python
system = SMFTSystem(
    grid_shape=(50, 50),
    N_oscillators=200,
    coupling='local'
)
# <500ms per step
```

### Large System (Needs Optimization)
```python
system = SMFTSystem(
    grid_shape=(100, 100),
    N_oscillators=500,
    coupling='local'
)
# ~2s per step (use smaller dt, fewer stores)
```

## Common Pitfalls

### 1. Numerical Instability

❌ **Bad**:
```python
solution = system.evolve((0, 100), dt=0.1)  # Too large!
```

✓ **Good**:
```python
solution = system.evolve((0, 100), dt=0.01)  # Stable
```

### 2. Memory Issues

❌ **Bad**:
```python
solution = system.evolve((0, 1000), store_interval=1)  # 100k snapshots!
```

✓ **Good**:
```python
solution = system.evolve((0, 1000), store_interval=100)  # Reasonable
```

### 3. Kernel Width

❌ **Bad**:
```python
R_field, _ = system.compute_local_order_parameter(kernel_width=1.0)  # Too wide
```

✓ **Good**:
```python
R_field, _ = system.compute_local_order_parameter(kernel_width=0.1)  # Appropriate
```

## Examples

See `examples/field_theory/`:
- `smft_demo.py` - Complete demonstrations
- `validate_integration.py` - Validation suite

## Testing

Run integration tests:
```bash
python examples/field_theory/validate_integration.py
```

Expected output:
```
VALIDATION SUMMARY:
  Hamiltonian ↔ Field: ✓ PASSED
  Local ↔ Global: ✓ PASSED
  Backward Compatibility: ✓ PASSED
  Heavy Mass Limit: ✓ PASSED
  Full System: ✓ PASSED

TOTAL: 5/6 validations passed
```

## Documentation

- **Full Guide**: `docs/field_theory_integration.md`
- **Integration Summary**: `field_theory_integration_summary.md`
- **API Reference**: See docstrings in `src/kuramoto/field_theory/smft_system.py`

## Support

For issues or questions:
1. Check `docs/field_theory_integration.md` for details
2. See examples in `examples/field_theory/`
3. Review integration tests in `tests/test_field_theory_integration.py`

## Next Steps

- Read full documentation: `docs/field_theory_integration.md`
- Run demos: `python examples/field_theory/smft_demo.py`
- Validate: `python examples/field_theory/validate_integration.py`
- Explore API: Check docstrings and type hints
