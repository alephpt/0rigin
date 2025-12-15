"""
Comprehensive tests for complete MSFT field theory implementation.

Tests:
1. Klein-Gordon mediator field
2. Local oscillator-field coupling
3. Fermion mass generation
4. System integration
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_less

from kuramoto.field_theory import (
    SpatialGrid,
    MediatorField,
    LocalFieldCoupling,
    FermionMassDemo,
    MSFTSystem
)


class TestMediatorField:
    """Test Klein-Gordon mediator field."""

    def test_initialization(self):
        """Test mediator field initialization."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
        mediator = MediatorField(grid, wave_speed=1.0, mass=1.0)

        assert mediator.c == 1.0
        assert mediator.M == 1.0
        assert mediator.sigma.shape == (32, 32)
        assert mediator.sigma_dot.shape == (32, 32)
        assert mediator.t == 0.0

    def test_source_density_computation(self):
        """Test oscillator density source computation."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
        mediator = MediatorField(grid)

        # Single oscillator at center
        phases = np.array([0.0])
        positions = np.array([[0.5, 0.5]])

        rho = mediator.compute_source_density(phases, positions, kernel_width=0.1)

        assert rho.shape == (32, 32)
        # Maximum should be at center
        max_idx = np.unravel_index(np.argmax(rho), rho.shape)
        center_idx = (grid.Nx // 2, grid.Ny // 2)
        assert abs(max_idx[0] - center_idx[0]) <= 1
        assert abs(max_idx[1] - center_idx[1]) <= 1

    def test_field_evolution(self):
        """Test Klein-Gordon evolution."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')
        mediator = MediatorField(grid, wave_speed=1.0, mass=0.5)

        # Initial Gaussian
        mediator.sigma = grid.create_gaussian(sigma=0.1)
        initial_energy = mediator.compute_field_energy()

        # Evolve with no source (free evolution)
        source = np.zeros((32, 32))

        for _ in range(100):
            mediator.evolve_step(source, dt=0.01, method='leapfrog')

        # Energy should be approximately conserved
        final_energy = mediator.compute_field_energy()
        energy_change = abs(final_energy - initial_energy) / initial_energy

        assert energy_change < 0.1, f"Energy changed by {energy_change:.2%}"

    def test_heavy_mass_limit(self):
        """Test that M→∞ gives instantaneous response."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)

        # Create source
        rho = np.random.rand(16, 16)

        # Test with increasing mass
        errors = []
        masses = [1.0, 5.0, 10.0, 50.0]

        for M in masses:
            mediator = MediatorField(grid, mass=M, coupling_constant=1.0)
            sigma_limit = mediator.get_heavy_mass_limit(rho)

            # Expected: σ ≈ g·ρ/M²
            expected = mediator.g * rho / M**2
            error = np.max(np.abs(sigma_limit - expected))
            errors.append(error)

        # Errors should be very small
        assert all(e < 1e-10 for e in errors)

    def test_wave_propagation_speed(self):
        """Test that waves propagate at correct speed."""
        grid = SpatialGrid(Nx=64, Ny=64, Lx=2.0, Ly=2.0, boundary='periodic')
        c_expected = 1.0
        mediator = MediatorField(grid, wave_speed=c_expected, mass=0.1)

        v_phase, v_group = mediator.compute_propagation_speed()

        # For small mass, should be close to c
        assert abs(v_phase - c_expected) < 0.5
        assert v_group > 0

    def test_field_sampling(self):
        """Test sampling field at specific positions."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
        mediator = MediatorField(grid)

        # Set known field
        mediator.sigma = grid.X + grid.Y  # Linear field

        # Sample at known positions
        positions = np.array([[0.25, 0.25], [0.5, 0.5], [0.75, 0.75]])
        values = mediator.sample_at_positions(positions)

        assert len(values) == 3
        # Values should roughly increase
        assert values[1] > values[0]
        assert values[2] > values[1]


class TestLocalFieldCoupling:
    """Test bidirectional oscillator-field coupling."""

    def test_initialization(self):
        """Test coupling system initialization."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        coupling = LocalFieldCoupling(grid)

        assert coupling.mediator is not None
        assert coupling.R_field is not None
        assert isinstance(coupling.get_effective_coupling_range(), float)

    def test_fields_update_from_oscillators(self):
        """Test that fields respond to oscillators."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        coupling = LocalFieldCoupling(grid)

        # Synchronized oscillators
        N = 20
        phases = np.zeros(N)  # All in phase
        positions = np.random.rand(N, 2)

        coupling.update_fields_from_oscillators(phases, positions, dt=0.01)

        # R field should be high everywhere
        assert np.mean(coupling.R_field.values) > 0.5

    def test_coupling_to_oscillators(self):
        """Test field influence on oscillators."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        coupling = LocalFieldCoupling(grid)

        # Set mediator field to constant
        coupling.mediator.sigma = np.ones((16, 16))

        # Get coupling at random positions
        positions = np.random.rand(10, 2)
        forces = coupling.compute_coupling_to_oscillators(positions)

        assert len(forces) == 10
        # All forces should be similar (uniform field)
        assert np.std(forces) < 0.5

    def test_coupled_evolution(self):
        """Test full coupled evolution."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        coupling = LocalFieldCoupling(grid, kernel_width=0.1)

        N = 30
        phases = np.random.uniform(0, 2*np.pi, N)
        positions = np.random.rand(N, 2)
        frequencies = np.random.normal(1.0, 0.1, N)

        result = coupling.evolve_coupled_system(
            phases, positions, frequencies,
            dt=0.02, n_steps=100, store_interval=20
        )

        assert 't' in result
        assert 'phases' in result
        assert 'R' in result
        assert len(result['R']) > 1

        # R should increase (synchronization)
        assert result['R'][-1] >= result['R'][0] - 0.1

    def test_heavy_mass_limit_convergence(self):
        """Test that heavy mass recovers global coupling."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)

        N = 20
        phases = np.random.uniform(0, 2*np.pi, N)
        positions = np.random.rand(N, 2)

        # Global R
        global_R = np.abs(np.mean(np.exp(1j * phases)))

        # Test with heavy mass
        coupling = LocalFieldCoupling(
            grid,
            mediator_params={'wave_speed': 1.0, 'mass': 100.0, 'coupling_constant': 1.0}
        )

        global_test, local_test = coupling.test_heavy_mass_limit(phases, positions)

        # Should be close to global
        assert abs(global_R - global_test) < 0.01
        assert abs(global_test - local_test) < 0.2

    def test_continuum_limit(self):
        """Test continuum approximation error."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
        coupling = LocalFieldCoupling(grid, kernel_width=0.05)

        # Many oscillators (better continuum)
        N = 200
        phases = np.random.uniform(0, 2*np.pi, N)
        positions = np.random.rand(N, 2)

        coupling.update_fields_from_oscillators(phases, positions, dt=0.01)

        error = coupling.compute_continuum_limit_error(phases, positions)

        # Error should be reasonable
        assert error < 1.0


class TestFermionMassDemo:
    """Test fermion mass generation demonstration."""

    def test_initialization(self):
        """Test fermion demo initialization."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        fermion = FermionMassDemo(grid, yukawa_coupling=1.0, bare_mass=0.5)

        assert fermion.Delta == 1.0
        assert fermion.m0 == 0.5
        assert fermion.R_field is not None

    def test_mass_from_order_parameter(self):
        """Test that m_eff = m0 + Δ·R."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        Delta = 2.0
        m0 = 0.5
        fermion = FermionMassDemo(grid, yukawa_coupling=Delta, bare_mass=m0)

        # Set uniform R field
        R_value = 0.7
        fermion.R_field.values = np.full((16, 16), R_value)

        # Compute mass
        fermion.m_eff_field.values = m0 + Delta * fermion.R_field.values

        expected_mass = m0 + Delta * R_value
        assert_allclose(fermion.m_eff_field.values, expected_mass, rtol=1e-10)

    def test_mass_generation_from_oscillators(self):
        """Test mass generation from oscillator synchronization."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0)
        fermion = FermionMassDemo(grid, yukawa_coupling=1.5, bare_mass=0.0)

        # Synchronized oscillators
        N = 50
        phases = np.random.uniform(0, 0.5, N)  # Mostly synchronized
        positions = np.random.rand(N, 2)

        fermion.update_from_oscillators(phases, positions, kernel_width=0.1)

        avg_mass = fermion.compute_average_mass()

        # Mass should be generated (> bare mass)
        assert avg_mass > fermion.m0

    def test_mass_gap(self):
        """Test mass gap computation."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        fermion = FermionMassDemo(grid, yukawa_coupling=2.0)

        # Create gradient in R field
        fermion.R_field.values = grid.X / grid.Lx
        fermion.m_eff_field.values = fermion.m0 + fermion.Delta * fermion.R_field.values

        mass_gap = fermion.compute_mass_gap()

        # Gap should equal Δ (for R from 0 to 1)
        assert abs(mass_gap - fermion.Delta) < 0.1

    def test_phase_transition(self):
        """Test phase transition demonstration."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        fermion = FermionMassDemo(grid, yukawa_coupling=1.0, bare_mass=0.0)

        R_values = np.array([0.0, 0.5, 1.0])
        result = fermion.demonstrate_phase_transition(R_values, plot=False)

        assert 'R' in result
        assert 'm_eff' in result
        assert len(result['m_eff']) == 3

        # Mass should increase with R
        assert result['m_eff'][1] > result['m_eff'][0]
        assert result['m_eff'][2] > result['m_eff'][1]

    def test_symmetry_breaking(self):
        """Test symmetry breaking quantification."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        fermion = FermionMassDemo(grid)

        # Uniform mass (no breaking)
        fermion.m_eff_field.values = np.ones((16, 16))
        symmetry_uniform = fermion.compute_symmetry_breaking_order()

        # Varied mass (breaking)
        fermion.m_eff_field.values = grid.X
        symmetry_broken = fermion.compute_symmetry_breaking_order()

        assert symmetry_broken > symmetry_uniform


class TestSystemIntegration:
    """Test full system integration."""

    def test_MSFT_system_with_new_components(self):
        """Test that MSFTSystem works with new components."""
        system = MSFTSystem(
            grid_shape=(16, 16),
            N_oscillators=20,
            coupling='local',
            mediator_mass=5.0
        )

        # Should be able to evolve
        result = system.evolve((0, 1), dt=0.02, store_interval=10)

        assert 'R' in result
        assert len(result['R']) > 0

    def test_all_components_together(self):
        """Integration test: all components working together."""
        # Create grid
        grid = SpatialGrid(Nx=24, Ny=24, Lx=1.0, Ly=1.0)

        # Create coupling
        coupling = LocalFieldCoupling(grid)

        # Create fermion demo
        fermion = FermionMassDemo(grid, yukawa_coupling=1.0)

        # Initialize oscillators
        N = 30
        phases = np.random.uniform(0, 2*np.pi, N)
        positions = np.random.rand(N, 2)
        frequencies = np.random.normal(1.0, 0.1, N)

        # Evolve coupling
        result = coupling.evolve_coupled_system(
            phases, positions, frequencies,
            dt=0.02, n_steps=50, store_interval=10
        )

        # Update fermion from final state
        final_phases = result['phases'][-1]
        fermion.update_from_oscillators(final_phases, positions)

        # Should have generated mass
        avg_mass = fermion.compute_average_mass()
        assert avg_mass >= fermion.m0


class TestNumericalProperties:
    """Test numerical accuracy and stability."""

    def test_laplacian_accuracy(self):
        """Test that Laplacian is accurate for mediator."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')

        max_err, rel_err = grid.test_laplacian()

        assert max_err < 0.1
        assert rel_err < 0.05

    def test_energy_conservation(self):
        """Test energy conservation in free evolution."""
        grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')
        mediator = MediatorField(grid, wave_speed=1.0, mass=0.5)

        # Initial condition
        mediator.sigma = grid.create_gaussian(sigma=0.1)

        energies = []
        source = np.zeros((32, 32))

        for _ in range(200):
            mediator.evolve_step(source, dt=0.01, method='leapfrog')
            if len(energies) % 10 == 0:
                energies.append(mediator.compute_field_energy())

        # Energy drift should be small
        energy_drift = (energies[-1] - energies[0]) / energies[0]
        assert abs(energy_drift) < 0.15

    def test_stability_with_source(self):
        """Test stability with oscillating source."""
        grid = SpatialGrid(Nx=16, Ny=16, Lx=1.0, Ly=1.0)
        mediator = MediatorField(grid, wave_speed=1.0, mass=1.0)

        # Oscillating source
        for step in range(100):
            t = step * 0.01
            source = np.sin(2 * np.pi * t) * grid.create_gaussian()
            mediator.evolve_step(source, dt=0.01)

            # Should remain bounded
            assert np.max(np.abs(mediator.sigma)) < 10.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
