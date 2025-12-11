"""
Unit tests for field theory prototypes.

Validates all four core components are working correctly.
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.kuramoto.field_theory import (
    HamiltonianKuramoto,
    SpatialGrid,
    ScalarField,
    PDESolver,
    PDE_AVAILABLE
)


class TestHamiltonianKuramoto:
    """Test Hamiltonian formulation of Kuramoto model."""

    def test_initialization(self):
        """Test model initialization."""
        N = 10
        K = 2.0
        frequencies = np.ones(N)

        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=K,
            frequencies=frequencies,
            damping=1.0
        )

        assert model.N == N
        assert model.K == K
        assert len(model.theta) == N
        assert len(model.p) == N
        assert model.t == 0.0

    def test_energy_conservation(self):
        """Test energy conservation in conservative case."""
        N = 20
        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=2.0,
            frequencies=np.zeros(N),
            damping=0.0  # Conservative system
        )

        # Initial energy
        E0 = model.compute_hamiltonian()

        # Evolve for short time with small timestep
        solution = model.evolve((0, 1.0), dt=0.001, method='rk4')

        # Final energy
        Ef = model.compute_hamiltonian()

        # Check conservation (relative error < 0.1%)
        relative_error = abs(Ef - E0) / abs(E0) if E0 != 0 else abs(Ef - E0)
        assert relative_error < 1e-3, f"Energy not conserved: {relative_error:.2e}"

    def test_overdamped_limit(self):
        """Test that high damping suppresses momentum."""
        N = 20
        frequencies = np.random.normal(0, 1, N)

        # High damping
        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=3.0,
            frequencies=frequencies,
            damping=100.0,
            initial_momenta=np.ones(N)  # Start with momentum
        )

        # Evolve
        solution = model.evolve((0, 1.0), dt=0.01)

        # Check that momenta decay to near zero
        final_p_max = np.max(np.abs(model.p))
        assert final_p_max < 0.1, f"Momentum not damped: {final_p_max}"

    def test_equations_of_motion(self):
        """Test that equations of motion are computed correctly."""
        N = 5
        model = HamiltonianKuramoto(
            N=N,
            coupling_strength=1.0,
            frequencies=np.ones(N),
            damping=1.0
        )

        state = np.concatenate([model.theta, model.p])
        derivatives = model.equations_of_motion(0, state)

        assert len(derivatives) == 2 * N
        # Check that dθ/dt = p
        np.testing.assert_array_equal(derivatives[:N], model.p)


class TestSpatialGrid:
    """Test spatial grid infrastructure."""

    def test_grid_creation(self):
        """Test grid initialization."""
        grid = SpatialGrid(50, 50, 2.0, 3.0, 'periodic')

        assert grid.Nx == 50
        assert grid.Ny == 50
        assert grid.Lx == 2.0
        assert grid.Ly == 3.0
        assert grid.X.shape == (50, 50)
        assert grid.Y.shape == (50, 50)

    def test_laplacian_periodic(self):
        """Test Laplacian with periodic BC."""
        grid = SpatialGrid(40, 40, 2*np.pi, 2*np.pi, 'periodic')
        max_error, rel_error = grid.test_laplacian()

        assert rel_error < 1e-2, f"Laplacian error too large: {rel_error}"

    def test_laplacian_dirichlet(self):
        """Test Laplacian with Dirichlet BC."""
        grid = SpatialGrid(40, 40, 2*np.pi, 2*np.pi, 'dirichlet')
        max_error, rel_error = grid.test_laplacian()

        # Dirichlet BC has larger errors at boundaries
        assert rel_error < 0.1, f"Laplacian error too large: {rel_error}"

    def test_gradient(self):
        """Test gradient operator."""
        grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')

        # Linear function: f = x + 2y
        f = grid.X + 2 * grid.Y

        grad_x, grad_y = grid.gradient(f)

        # Should be constant gradients
        assert np.allclose(grad_x, 1.0, atol=1e-2)
        assert np.allclose(grad_y, 2.0, atol=1e-2)

    def test_gaussian_creation(self):
        """Test Gaussian field creation."""
        grid = SpatialGrid(50, 50, 1.0, 1.0)
        gaussian = grid.create_gaussian(center=(0.5, 0.5), sigma=0.1)

        # Check peak is at center
        max_idx = np.unravel_index(np.argmax(gaussian), gaussian.shape)
        assert max_idx == (25, 25)  # Center of 50x50 grid

        # Check normalization
        assert np.max(gaussian) == pytest.approx(1.0)


class TestScalarField:
    """Test scalar field evolution."""

    def test_field_initialization(self):
        """Test field creation."""
        grid = SpatialGrid(20, 20, 1.0, 1.0)
        field = ScalarField(grid, name="test")

        assert field.name == "test"
        assert field.values.shape == (20, 20)
        assert field.t == 0.0
        assert np.all(field.values == 0)

    def test_diffusion(self):
        """Test diffusion reduces maximum value."""
        grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')
        field = ScalarField(grid)

        # Initial peak
        field.values = grid.create_gaussian(sigma=0.05)
        initial_max = np.max(field.values)

        # Apply diffusion
        for _ in range(50):
            field.diffuse(diffusion_coeff=0.01, dt=0.01)

        final_max = np.max(field.values)

        # Check diffusion occurred
        assert final_max < initial_max * 0.9

    def test_oscillator_coupling(self):
        """Test field update from oscillators."""
        grid = SpatialGrid(20, 20, 1.0, 1.0)
        field = ScalarField(grid)

        # Place oscillators
        N = 5
        phases = np.zeros(N)  # All in phase
        positions = np.array([
            [0.5, 0.5],  # Center
            [0.3, 0.3],
            [0.7, 0.3],
            [0.3, 0.7],
            [0.7, 0.7]
        ])

        field.update_from_oscillators(phases, positions, kernel_width=0.1)

        # Field should be non-zero where oscillators are
        assert np.max(field.values) > 0

        # Maximum should be near center (more oscillators nearby)
        max_idx = np.unravel_index(np.argmax(field.values), field.values.shape)
        assert 8 <= max_idx[0] <= 12  # Near center

    def test_pde_evolution(self):
        """Test PDE evolution with simple dynamics."""
        grid = SpatialGrid(20, 20, 1.0, 1.0)
        field = ScalarField(grid)

        # Initial condition
        field.values = np.random.randn(20, 20) * 0.1

        # Simple decay dynamics
        def dynamics(t, f):
            return -f  # df/dt = -f

        field.evolve_pde(dynamics, t_span=(0, 1.0), dt=0.01, store_interval=10)

        # Check decay
        assert np.max(np.abs(field.values)) < 0.5

        # Check history was stored
        assert len(field.history) > 0
        assert len(field.time_history) > 0


class TestPDESolver:
    """Test PDE solver integration."""

    @pytest.mark.skipif(not PDE_AVAILABLE, reason="py-pde not installed")
    def test_solver_creation(self):
        """Test solver initialization."""
        solver = PDESolver(
            grid_shape=(20, 20),
            grid_size=(1.0, 1.0),
            boundary='periodic'
        )

        assert solver.grid_shape == (20, 20)
        assert solver.grid_size == (1.0, 1.0)
        assert solver.boundary == 'periodic'

    @pytest.mark.skipif(not PDE_AVAILABLE, reason="py-pde not installed")
    def test_diffusion_solve(self):
        """Test diffusion equation solving."""
        solver = PDESolver(
            grid_shape=(20, 20),
            grid_size=(1.0, 1.0),
            boundary='periodic'
        )

        # Gaussian initial condition
        X, Y = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20), indexing='ij')
        initial = np.exp(-20 * ((X - 0.5)**2 + (Y - 0.5)**2))

        result = solver.solve_diffusion(
            initial,
            diffusion_coeff=0.1,
            t_range=(0, 0.1),
            dt=0.01,
            tracker=None
        )

        assert 'times' in result
        assert 'field' in result
        assert 'final_state' in result

        # Check diffusion occurred
        final_max = np.max(result['field'])
        initial_max = np.max(initial)
        assert final_max < initial_max

    def test_no_pde_fallback(self):
        """Test behavior when py-pde is not available."""
        if PDE_AVAILABLE:
            pytest.skip("py-pde is available")

        with pytest.raises(ImportError):
            solver = PDESolver(
                grid_shape=(20, 20),
                grid_size=(1.0, 1.0)
            )


def test_integration():
    """Test that all components work together."""
    # Create grid
    grid = SpatialGrid(30, 30, 1.0, 1.0, 'periodic')

    # Create field
    field = ScalarField(grid, name="order_parameter")

    # Create Hamiltonian model
    N = 50
    model = HamiltonianKuramoto(
        N=N,
        coupling_strength=3.0,
        frequencies=np.random.normal(0, 1, N),
        damping=1.0
    )

    # Evolve model
    solution = model.evolve((0, 1.0), dt=0.01)

    # Update field from oscillators (random positions for now)
    positions = np.random.rand(N, 2)
    field.update_from_oscillators(model.theta, positions)

    # Apply diffusion
    field.diffuse(0.01, 0.01)

    # All components should work together
    assert field.t > 0
    assert model.t > 0
    assert np.max(field.values) > 0


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])