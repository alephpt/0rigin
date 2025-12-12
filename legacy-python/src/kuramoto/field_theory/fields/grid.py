"""
Spatial grid infrastructure for field representations.

Provides 2D grid with finite difference operators and boundary conditions.
"""

from typing import Literal, Tuple, Optional
import numpy as np
from numpy.typing import NDArray


class SpatialGrid:
    """
    2D spatial grid for field theory calculations.

    Supports finite difference operators (gradient, Laplacian) and
    various boundary conditions (periodic, Dirichlet, Neumann).
    """

    def __init__(
        self,
        Nx: int,
        Ny: int,
        Lx: float = 1.0,
        Ly: float = 1.0,
        boundary: Literal['periodic', 'dirichlet', 'neumann'] = 'periodic'
    ):
        """
        Initialize spatial grid.

        Parameters
        ----------
        Nx, Ny : int
            Number of grid points in x and y directions.
        Lx, Ly : float
            Physical domain size.
        boundary : str
            Boundary condition type.
        """
        self.Nx = Nx
        self.Ny = Ny
        self.Lx = Lx
        self.Ly = Ly
        self.boundary = boundary

        # Grid spacing
        self.dx = Lx / Nx
        self.dy = Ly / Ny

        # Create coordinate arrays
        self.x = np.linspace(0, Lx - self.dx, Nx)
        self.y = np.linspace(0, Ly - self.dy, Ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

    def laplacian(self, field: NDArray) -> NDArray:
        """
        Compute Laplacian ∇²field using finite differences.

        Parameters
        ----------
        field : NDArray
            2D field array of shape (Nx, Ny).

        Returns
        -------
        NDArray
            Laplacian of field.
        """
        assert field.shape == (self.Nx, self.Ny), "Field shape must match grid"

        lap = np.zeros_like(field)

        # Second derivative in x
        if self.boundary == 'periodic':
            # Periodic boundary conditions
            lap += (np.roll(field, 1, axis=0) - 2*field + np.roll(field, -1, axis=0)) / self.dx**2
            lap += (np.roll(field, 1, axis=1) - 2*field + np.roll(field, -1, axis=1)) / self.dy**2

        elif self.boundary == 'dirichlet':
            # Dirichlet BC (field = 0 at boundaries)
            # Interior points
            lap[1:-1, 1:-1] = (
                (field[2:, 1:-1] - 2*field[1:-1, 1:-1] + field[:-2, 1:-1]) / self.dx**2 +
                (field[1:-1, 2:] - 2*field[1:-1, 1:-1] + field[1:-1, :-2]) / self.dy**2
            )

            # Boundary points (assuming field=0 outside)
            # Left/right edges
            lap[0, 1:-1] = (field[1, 1:-1] - 2*field[0, 1:-1]) / self.dx**2
            lap[-1, 1:-1] = (field[-2, 1:-1] - 2*field[-1, 1:-1]) / self.dx**2

            # Top/bottom edges
            lap[1:-1, 0] = (field[1:-1, 1] - 2*field[1:-1, 0]) / self.dy**2
            lap[1:-1, -1] = (field[1:-1, -2] - 2*field[1:-1, -1]) / self.dy**2

            # Corners
            lap[0, 0] = (field[1, 0] - 2*field[0, 0]) / self.dx**2 + \
                        (field[0, 1] - 2*field[0, 0]) / self.dy**2
            lap[-1, 0] = (field[-2, 0] - 2*field[-1, 0]) / self.dx**2 + \
                         (field[-1, 1] - 2*field[-1, 0]) / self.dy**2
            lap[0, -1] = (field[1, -1] - 2*field[0, -1]) / self.dx**2 + \
                         (field[0, -2] - 2*field[0, -1]) / self.dy**2
            lap[-1, -1] = (field[-2, -1] - 2*field[-1, -1]) / self.dx**2 + \
                          (field[-1, -2] - 2*field[-1, -1]) / self.dy**2

        elif self.boundary == 'neumann':
            # Neumann BC (zero gradient at boundaries)
            # Interior points
            lap[1:-1, 1:-1] = (
                (field[2:, 1:-1] - 2*field[1:-1, 1:-1] + field[:-2, 1:-1]) / self.dx**2 +
                (field[1:-1, 2:] - 2*field[1:-1, 1:-1] + field[1:-1, :-2]) / self.dy**2
            )

            # Use one-sided differences at boundaries
            # Left/right edges (mirror boundary)
            lap[0, 1:-1] = (2*field[1, 1:-1] - 2*field[0, 1:-1]) / self.dx**2
            lap[-1, 1:-1] = (2*field[-2, 1:-1] - 2*field[-1, 1:-1]) / self.dx**2

            # Top/bottom edges
            lap[1:-1, 0] = (2*field[1:-1, 1] - 2*field[1:-1, 0]) / self.dy**2
            lap[1:-1, -1] = (2*field[1:-1, -2] - 2*field[1:-1, -1]) / self.dy**2

            # Corners
            lap[0, 0] = 2*(field[1, 0] - field[0, 0]) / self.dx**2 + \
                        2*(field[0, 1] - field[0, 0]) / self.dy**2
            lap[-1, 0] = 2*(field[-2, 0] - field[-1, 0]) / self.dx**2 + \
                         2*(field[-1, 1] - field[-1, 0]) / self.dy**2
            lap[0, -1] = 2*(field[1, -1] - field[0, -1]) / self.dx**2 + \
                         2*(field[0, -2] - field[0, -1]) / self.dy**2
            lap[-1, -1] = 2*(field[-2, -1] - field[-1, -1]) / self.dx**2 + \
                          2*(field[-1, -2] - field[-1, -1]) / self.dy**2

        return lap

    def gradient(self, field: NDArray) -> Tuple[NDArray, NDArray]:
        """
        Compute gradient ∇field using finite differences.

        Parameters
        ----------
        field : NDArray
            2D field array of shape (Nx, Ny).

        Returns
        -------
        grad_x, grad_y : NDArray
            x and y components of gradient.
        """
        assert field.shape == (self.Nx, self.Ny), "Field shape must match grid"

        if self.boundary == 'periodic':
            # Central differences with periodic BC
            grad_x = (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2 * self.dx)
            grad_y = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2 * self.dy)
            return grad_x, grad_y

        # Non-periodic boundaries
        grad_x = np.zeros_like(field)
        grad_y = np.zeros_like(field)

        if self.boundary in ['dirichlet', 'neumann']:
            # Central differences in interior
            grad_x[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2 * self.dx)
            grad_y[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2 * self.dy)

            if self.boundary == 'dirichlet':
                # One-sided at boundaries (field=0 outside)
                grad_x[0, :] = field[1, :] / self.dx
                grad_x[-1, :] = -field[-2, :] / self.dx
                grad_y[:, 0] = field[:, 1] / self.dy
                grad_y[:, -1] = -field[:, -2] / self.dy

            else:  # neumann
                # Zero gradient at boundaries
                grad_x[0, :] = 0
                grad_x[-1, :] = 0
                grad_y[:, 0] = 0
                grad_y[:, -1] = 0

        return grad_x, grad_y

    def divergence(self, vec_x: NDArray, vec_y: NDArray) -> NDArray:
        """
        Compute divergence ∇·vec of vector field.

        Parameters
        ----------
        vec_x, vec_y : NDArray
            x and y components of vector field.

        Returns
        -------
        NDArray
            Divergence of vector field.
        """
        grad_x, _ = self.gradient(vec_x)
        _, grad_y = self.gradient(vec_y)
        return grad_x + grad_y

    def test_laplacian(self) -> Tuple[float, float]:
        """
        Test Laplacian operator with analytical function.

        Uses f(x,y) = sin(kx*x)sin(ky*y) where ∇²f = -(kx² + ky²)f

        Returns
        -------
        max_error : float
            Maximum absolute error.
        rel_error : float
            Relative L2 error.
        """
        # Test function parameters
        kx = 2 * np.pi / self.Lx
        ky = 2 * np.pi / self.Ly

        # Analytical function
        f = np.sin(kx * self.X) * np.sin(ky * self.Y)

        # Analytical Laplacian
        lap_analytical = -(kx**2 + ky**2) * f

        # Numerical Laplacian
        lap_numerical = self.laplacian(f)

        # Errors (exclude boundaries for non-periodic BC)
        if self.boundary == 'periodic':
            max_error = np.max(np.abs(lap_numerical - lap_analytical))
            rel_error = np.linalg.norm(lap_numerical - lap_analytical) / np.linalg.norm(lap_analytical)
        else:
            # Only check interior points
            interior = slice(1, -1), slice(1, -1)
            max_error = np.max(np.abs(lap_numerical[interior] - lap_analytical[interior]))
            rel_error = (np.linalg.norm(lap_numerical[interior] - lap_analytical[interior]) /
                        np.linalg.norm(lap_analytical[interior]))

        return max_error, rel_error

    def create_gaussian(
        self,
        center: Tuple[float, float] = None,
        sigma: float = 0.1
    ) -> NDArray:
        """
        Create Gaussian field on grid.

        Parameters
        ----------
        center : tuple
            Center position (x0, y0). Defaults to domain center.
        sigma : float
            Gaussian width.

        Returns
        -------
        NDArray
            Gaussian field.
        """
        if center is None:
            center = (self.Lx / 2, self.Ly / 2)

        x0, y0 = center
        return np.exp(-((self.X - x0)**2 + (self.Y - y0)**2) / (2 * sigma**2))

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"SpatialGrid(Nx={self.Nx}, Ny={self.Ny}, "
            f"Lx={self.Lx:.2f}, Ly={self.Ly:.2f}, "
            f"boundary='{self.boundary}')"
        )