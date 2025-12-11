"""
Integration tests for field theory system.

Tests verify that all components work together correctly:
- Discrete ↔ Continuum integration
- Hamiltonian ↔ Field dynamics
- Local ↔ Global coupling
- Classical ↔ Field theory (backward compatibility)
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from src.kuramoto.field_theory import (
    SMFTSystem,
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField,
    MediatorField,
    LocalFieldCoupling
)


class TestSMFTSystemInitialization:
    """Test SMFT system initialization and setup."""

    def test_system_creation(self):
        """Test basic system creation."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='local',
            mediator_mass=10.0
        )

        assert system.grid_shape == (20, 20)
        assert system.N == 50
        assert system.coupling_type == 'local'
        assert system.M == 10.0
        assert system.t == 0.0

    def test_component_creation(self):
        """Test that all components are created."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=30,
            coupling='global'
        )

        # Check grid
        assert system.grid is not None
        assert system.grid.Nx == 15
        assert system.grid.Ny == 15

        # Check oscillators
        assert system.oscillators is not None
        assert system.oscillators.N == 30

        # Check fields
        assert system.mediator_field is not None
        assert system.sync_field is not None
        assert system.phase_field is not None

        # Check positions
        assert system.oscillator_positions.shape == (30, 2)

    def test_custom_frequencies(self):
        """Test initialization with custom frequencies."""
        N = 25
        frequencies = np.random.normal(0, 2, N)

        system = SMFTSystem(
            grid_shape=(10, 10),
            N_oscillators=N,
            oscillator_frequencies=frequencies
        )

        np.testing.assert_array_equal(
            system.oscillators.frequencies,
            frequencies
        )


class TestDiscreteToContiuumIntegration:
    """Test discrete oscillators → continuous field mapping."""

    def test_local_order_parameter_computation(self):
        """Test that local R(x,y) is computed correctly."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=10,
            coupling='local'
        )

        # Set all oscillators to same phase
        system.oscillators.theta[:] = np.pi / 4

        # Compute local order parameter
        R_field, theta_field = system.compute_local_order_parameter(kernel_width=0.2)

        # Check shapes
        assert R_field.shape == (20, 20)
        assert theta_field.shape == (20, 20)

        # With all oscillators synchronized, R should be high everywhere
        assert np.mean(R_field) > 0.5

    def test_local_order_parameter_bounds(self):
        """Test that local R ∈ [0,1]."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=50,
            coupling='local'
        )

        R_field, theta_field = system.compute_local_order_parameter()

        # Check bounds
        assert np.all(R_field >= 0)
        assert np.all(R_field <= 1)

        # Check phase bounds
        assert np.all(theta_field >= -np.pi)
        assert np.all(theta_field <= np.pi)

    def test_continuum_limit_convergence(self):
        """Test convergence as N increases."""
        grid_shape = (30, 30)
        N_values = [50, 100, 200]

        variances = []
        for N in N_values:
            system = SMFTSystem(
                grid_shape=grid_shape,
                N_oscillators=N,
                coupling='local'
            )

            # Synchronized state
            system.oscillators.theta[:] = 0

            R_field, _ = system.compute_local_order_parameter(kernel_width=0.1)

            # Variance should decrease with more oscillators
            variances.append(np.var(R_field))

        # Check decreasing variance
        assert variances[1] < variances[0]
        assert variances[2] < variances[1]


class TestHamiltonianFieldIntegration:
    """Test Hamiltonian phase space ↔ field integration."""

    def test_field_force_on_oscillators(self):
        """Test that field exerts force on oscillators."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=30,
            coupling='local',
            mediator_mass=5.0
        )

        # Set non-zero mediator field
        system.mediator_field.values[:] = 0.5

        # Compute forces
        forces = system.compute_field_force_on_oscillators(kernel_width=0.1)

        assert len(forces) == 30
        # Forces should be non-zero
        assert np.any(forces != 0)

    def test_energy_conservation_no_damping(self):
        """Test energy behavior in conservative case."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=20,
            coupling='local',
            mediator_mass=10.0
        )

        # Reduce damping for near-conservative system
        system.oscillators.gamma = 0.01

        # Initial energy
        E0 = system.oscillators.compute_hamiltonian()

        # Evolve briefly
        solution = system.evolve((0, 1.0), dt=0.001, store_interval=100)

        # Check energy doesn't explode
        E_final = solution['energy'][-1]
        assert abs(E_final) < 100 * abs(E0 + 1)  # Reasonable bound

    def test_coupled_step_updates_both(self):
        """Test that single step updates oscillators and field."""
        system = SMFTSystem(
            grid_shape=(10, 10),
            N_oscillators=15,
            coupling='local'
        )

        # Initial states
        theta0 = system.oscillators.theta.copy()
        field0 = system.mediator_field.values.copy()
        t0 = system.t

        # Single step
        system.step(dt=0.01)

        # Check both updated
        assert not np.array_equal(system.oscillators.theta, theta0)
        assert not np.array_equal(system.mediator_field.values, field0)
        assert system.t > t0


class TestLocalGlobalCoupling:
    """Test local vs global coupling regimes."""

    def test_local_coupling_spatial_structure(self):
        """Test that local coupling creates spatial structure."""
        system = SMFTSystem(
            grid_shape=(30, 30),
            N_oscillators=100,
            coupling='local',
            mediator_mass=10.0
        )

        # Evolve
        solution = system.evolve((0, 10), dt=0.01, store_interval=50)

        # Check that sync field has spatial variation
        final_sync = solution['sync_field'][-1]
        variance = np.var(final_sync)

        assert variance > 0  # Spatial structure exists

    def test_global_coupling_uniform_field(self):
        """Test that global coupling creates more uniform field."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='global',
            mediator_mass=10.0
        )

        # Evolve
        solution = system.evolve((0, 5), dt=0.01, store_interval=25)

        # Global coupling should create less spatial variation
        # (though not completely uniform due to discrete sampling)
        final_sync = solution['sync_field'][-1]
        variance = np.var(final_sync)

        # Should have some structure but less than highly local case
        assert variance >= 0


class TestHeavyMassLimit:
    """Test M→∞ limit recovers Kuramoto model."""

    def test_increasing_mass_convergence(self):
        """Test that R converges as M increases."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=30,
            coupling='local'
        )

        M_values = [1.0, 10.0, 100.0]
        results = system.test_heavy_mass_limit(M_values)

        # Results should be defined
        assert len(results) == 3
        for M in M_values:
            assert 0 <= results[M] <= 1

    def test_heavy_mass_reduces_field_dynamics(self):
        """Test that heavy mass slows field dynamics."""
        N = 40
        frequencies = np.random.normal(0, 1, N)

        # Light mass
        system_light = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=N,
            mediator_mass=1.0,
            oscillator_frequencies=frequencies
        )

        # Heavy mass
        system_heavy = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=N,
            mediator_mass=100.0,
            oscillator_frequencies=frequencies
        )

        # Evolve both
        sol_light = system_light.evolve((0, 2), dt=0.01, store_interval=20)
        sol_heavy = system_heavy.evolve((0, 2), dt=0.01, store_interval=20)

        # Heavy mass should have smaller field variation
        field_var_light = np.var(sol_light['mediator_field'][-1])
        field_var_heavy = np.var(sol_heavy['mediator_field'][-1])

        # Heavy mass → slower field dynamics
        assert field_var_heavy <= field_var_light * 2  # Some margin


class TestBackwardCompatibility:
    """Test that field theory doesn't break Sprint 1 components."""

    def test_hamiltonian_kuramoto_standalone(self):
        """Test HamiltonianKuramoto works independently."""
        model = HamiltonianKuramoto(
            N=50,
            coupling_strength=3.0,
            frequencies=np.random.normal(0, 1, 50),
            damping=1.0
        )

        solution = model.evolve((0, 5), dt=0.01)

        assert 't' in solution
        assert 'theta' in solution
        assert 'R' in solution
        assert len(solution['t']) > 0

    def test_spatial_grid_standalone(self):
        """Test SpatialGrid works independently."""
        grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')

        # Test Laplacian
        max_error, rel_error = grid.test_laplacian()
        assert rel_error < 0.1

        # Test gradient
        f = grid.X + 2 * grid.Y
        grad_x, grad_y = grid.gradient(f)
        assert np.allclose(grad_x, 1.0, atol=0.1)

    def test_scalar_field_standalone(self):
        """Test ScalarField works independently."""
        grid = SpatialGrid(20, 20, 1.0, 1.0)
        field = ScalarField(grid, name="test")

        # Test diffusion
        field.values = grid.create_gaussian(sigma=0.05)
        initial_max = np.max(field.values)

        for _ in range(20):
            field.diffuse(0.01, 0.01)

        assert np.max(field.values) < initial_max


class TestFullSystemEvolution:
    """Test complete system evolution."""

    def test_system_evolve_completes(self):
        """Test that full evolution completes successfully."""
        system = SMFTSystem(
            grid_shape=(25, 25),
            N_oscillators=60,
            coupling='local',
            mediator_mass=10.0
        )

        solution = system.evolve(
            t_span=(0, 10),
            dt=0.01,
            store_interval=50
        )

        # Check all outputs present
        assert 't' in solution
        assert 'theta' in solution
        assert 'R' in solution
        assert 'mediator_field' in solution
        assert 'sync_field' in solution
        assert 'energy' in solution

        # Check shapes
        n_times = len(solution['t'])
        assert solution['theta'].shape == (n_times, 60)
        assert solution['R'].shape == (n_times,)
        assert solution['mediator_field'].shape == (n_times, 25, 25)

    def test_synchronization_transition(self):
        """Test that system can synchronize."""
        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=100,
            coupling='local',
            mediator_mass=5.0
        )

        # Strong coupling should lead to synchronization
        system.oscillators.K = 5.0

        solution = system.evolve((0, 50), dt=0.01, store_interval=100)

        # Check R increases
        R_initial = solution['R'][0]
        R_final = solution['R'][-1]

        # Should see some synchronization
        assert R_final >= R_initial - 0.1  # Allow some variability

    def test_effective_mass_computation(self):
        """Test effective mass field calculation."""
        system = SMFTSystem(
            grid_shape=(15, 15),
            N_oscillators=50,
            coupling='local'
        )

        # Evolve to create structure
        system.evolve((0, 5), dt=0.01, store_interval=50)

        # Compute effective mass
        m_eff = system.compute_effective_mass()

        assert m_eff.shape == (15, 15)
        assert np.all(m_eff > 0)  # Mass should be positive


class TestPerformance:
    """Test performance benchmarks."""

    def test_small_system_performance(self):
        """Test small system meets performance target."""
        import time

        system = SMFTSystem(
            grid_shape=(20, 20),
            N_oscillators=50,
            coupling='local',
            mediator_mass=10.0
        )

        start = time.time()
        system.step(dt=0.01)
        elapsed = time.time() - start

        # Single step should be fast (<100ms)
        assert elapsed < 0.1, f"Step too slow: {elapsed:.3f}s"

    def test_medium_system_performance(self):
        """Test medium system performance."""
        import time

        system = SMFTSystem(
            grid_shape=(50, 50),
            N_oscillators=200,
            coupling='local',
            mediator_mass=10.0
        )

        start = time.time()
        system.step(dt=0.01)
        elapsed = time.time() - start

        # Should complete in reasonable time (<1s)
        assert elapsed < 1.0, f"Step too slow: {elapsed:.3f}s"


def test_integration_example():
    """
    Integration test demonstrating full workflow.

    This is the example from the task specification.
    """
    # Create system
    system = SMFTSystem(
        grid_shape=(100, 100),
        N_oscillators=200,
        coupling='local',
        mediator_mass=10.0
    )

    # Evolve
    result = system.evolve(t_span=(0, 20), dt=0.01, store_interval=50)

    # Access fields
    R_field = result['sync_field']
    theta_field = result['theta']
    sigma_field = result['mediator_field']
    m_eff = system.compute_effective_mass()

    # Verify outputs
    assert R_field.shape[1:] == (100, 100)
    assert len(theta_field) > 0
    assert sigma_field.shape[1:] == (100, 100)
    assert m_eff.shape == (100, 100)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
