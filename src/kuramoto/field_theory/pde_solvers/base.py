"""
PDE solver integration using py-pde library.

Wrapper around py-pde for field theory simulations.
"""

from typing import Optional, Callable, Dict, Any, Tuple
import numpy as np
from numpy.typing import NDArray

# Check if py-pde is available
try:
    from pde import (
        DiffusionPDE,
        PDE,
        ScalarField as PDEScalarField,
        UnitGrid,
        CartesianGrid,
        PeriodicBCValue,
        DirichletBC,
        NeumannBC
    )
    PDE_AVAILABLE = True
except ImportError:
    PDE_AVAILABLE = False
    print("Warning: py-pde not installed. Install with: pip install py-pde")


class PDESolver:
    """
    Wrapper for py-pde library to solve field evolution PDEs.

    Provides unified interface for different PDE types relevant
    to Kuramoto field theory.
    """

    def __init__(
        self,
        grid_shape: Tuple[int, int],
        grid_size: Tuple[float, float] = (1.0, 1.0),
        boundary: str = 'periodic'
    ):
        """
        Initialize PDE solver.

        Parameters
        ----------
        grid_shape : tuple
            Number of grid points (Nx, Ny).
        grid_size : tuple
            Physical domain size (Lx, Ly).
        boundary : str
            Boundary condition type.
        """
        if not PDE_AVAILABLE:
            raise ImportError(
                "py-pde library not available. "
                "Install with: pip install py-pde"
            )

        self.grid_shape = grid_shape
        self.grid_size = grid_size
        self.boundary = boundary

        # Create py-pde grid
        self.pde_grid = CartesianGrid(
            bounds=[[0, grid_size[0]], [0, grid_size[1]]],
            shape=grid_shape,
            periodic=(boundary == 'periodic')
        )

        # Set up boundary conditions
        if boundary == 'periodic':
            self.bc = 'periodic'
        elif boundary == 'dirichlet':
            self.bc = DirichletBC(self.pde_grid, 0)
        elif boundary == 'neumann':
            self.bc = NeumannBC(self.pde_grid, 0)
        else:
            raise ValueError(f"Unknown boundary: {boundary}")

    def solve_diffusion(
        self,
        initial_field: NDArray,
        diffusion_coeff: float,
        t_range: Tuple[float, float],
        dt: float = 0.01,
        tracker: Optional[str] = 'progress'
    ) -> Dict[str, Any]:
        """
        Solve diffusion equation: ∂u/∂t = D∇²u.

        Parameters
        ----------
        initial_field : NDArray
            Initial condition u(x,y,0).
        diffusion_coeff : float
            Diffusion coefficient D.
        t_range : tuple
            Time interval (t0, tf).
        dt : float
            Time step for output.
        tracker : str
            Tracker type for progress.

        Returns
        -------
        dict
            Solution with times and field evolution.
        """
        # Create PDE
        eq = DiffusionPDE(diffusion_coeff, bc=self.bc)

        # Create initial state
        state = PDEScalarField(self.pde_grid, initial_field.flatten())

        # Solve
        result = eq.solve(
            state,
            t_range=t_range,
            dt=dt,
            tracker=tracker
        )

        # Extract solution
        times = np.arange(t_range[0], t_range[1] + dt, dt)
        field_evolution = result.data.reshape(self.grid_shape)

        return {
            'times': times,
            'field': field_evolution,
            'final_state': result
        }

    def solve_custom_pde(
        self,
        pde_rhs: str,
        initial_field: NDArray,
        t_range: Tuple[float, float],
        dt: float = 0.01,
        parameters: Optional[Dict[str, float]] = None
    ) -> Dict[str, Any]:
        """
        Solve custom PDE using expression string.

        Parameters
        ----------
        pde_rhs : str
            Right-hand side expression, e.g., "D * laplace(u) - u**3 + u".
        initial_field : NDArray
            Initial condition.
        t_range : tuple
            Time interval.
        dt : float
            Time step.
        parameters : dict
            Parameters in PDE expression.

        Returns
        -------
        dict
            Solution dictionary.
        """
        # Default parameters
        if parameters is None:
            parameters = {}

        # Create PDE from expression
        eq = PDE(
            {
                "u": pde_rhs
            },
            bc=self.bc,
            user_funcs=parameters
        )

        # Initial state
        state = PDEScalarField.from_expression(
            self.pde_grid,
            initial_field.flatten()
        )

        # Solve
        result = eq.solve(
            state,
            t_range=t_range,
            dt=dt
        )

        return {
            'times': np.arange(t_range[0], t_range[1] + dt, dt),
            'field': result.data.reshape(self.grid_shape),
            'final_state': result
        }

    def benchmark_performance(
        self,
        grid_sizes: list = None
    ) -> Dict[str, float]:
        """
        Benchmark PDE solver performance.

        Parameters
        ----------
        grid_sizes : list
            Grid sizes to test. Defaults to [50, 100, 200].

        Returns
        -------
        dict
            Timing results for each grid size.
        """
        import time

        if grid_sizes is None:
            grid_sizes = [50, 100, 200]

        results = {}

        for N in grid_sizes:
            # Create solver for this grid
            solver = PDESolver(
                grid_shape=(N, N),
                grid_size=(1.0, 1.0),
                boundary='periodic'
            )

            # Random initial condition
            initial = np.random.randn(N, N)

            # Time the solve
            start_time = time.time()

            try:
                solver.solve_diffusion(
                    initial,
                    diffusion_coeff=0.1,
                    t_range=(0, 1.0),
                    dt=0.01,
                    tracker=None  # No progress bar for benchmarking
                )
                elapsed = time.time() - start_time
                results[f'{N}x{N}'] = elapsed

            except Exception as e:
                results[f'{N}x{N}'] = f"Failed: {str(e)}"

        return results

    @staticmethod
    def test_integration() -> bool:
        """
        Test py-pde integration with simple example.

        Returns
        -------
        bool
            True if test passes.
        """
        if not PDE_AVAILABLE:
            return False

        try:
            # Create small test grid
            solver = PDESolver(
                grid_shape=(20, 20),
                grid_size=(1.0, 1.0),
                boundary='periodic'
            )

            # Gaussian initial condition
            X, Y = np.meshgrid(
                np.linspace(0, 1, 20),
                np.linspace(0, 1, 20),
                indexing='ij'
            )
            initial = np.exp(-20 * ((X - 0.5)**2 + (Y - 0.5)**2))

            # Solve diffusion
            result = solver.solve_diffusion(
                initial,
                diffusion_coeff=0.01,
                t_range=(0, 0.5),
                dt=0.1,
                tracker=None
            )

            # Check that field diffused (max value decreased)
            final_max = np.max(result['field'])
            initial_max = np.max(initial)

            return final_max < initial_max * 0.9

        except Exception as e:
            print(f"Test failed: {e}")
            return False

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"PDESolver(grid={self.grid_shape}, "
            f"size={self.grid_size}, "
            f"boundary='{self.boundary}')"
        )