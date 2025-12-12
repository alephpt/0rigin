"""
Scalar field representations for Kuramoto order parameters.

Provides field objects that can couple to discrete oscillators.
"""

from typing import Optional, Callable, Tuple
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from .grid import SpatialGrid


class ScalarField:
    """
    Scalar field R(x,t) or θ(x,t) on spatial grid.

    Can represent order parameter amplitude, mean phase, or other
    scalar quantities that evolve according to field dynamics.
    """

    def __init__(
        self,
        grid: SpatialGrid,
        initial_values: Optional[NDArray] = None,
        name: str = "field"
    ):
        """
        Initialize scalar field.

        Parameters
        ----------
        grid : SpatialGrid
            Underlying spatial grid.
        initial_values : NDArray, optional
            Initial field values. Zeros if None.
        name : str
            Field name for identification.
        """
        self.grid = grid
        self.name = name
        self.t = 0.0

        # Initialize field values
        if initial_values is None:
            self.values = np.zeros((grid.Nx, grid.Ny))
        else:
            assert initial_values.shape == (grid.Nx, grid.Ny)
            self.values = initial_values.copy()

        # Storage for history
        self.history = []
        self.time_history = []

    def update_from_oscillators(
        self,
        oscillator_phases: NDArray,
        oscillator_positions: NDArray,
        kernel_width: float = 0.1,
        coupling_strength: float = 1.0
    ):
        """
        Update field from discrete oscillator configuration.

        Simple coupling: field responds to local oscillator density.

        Parameters
        ----------
        oscillator_phases : NDArray
            Phases of N oscillators.
        oscillator_positions : NDArray
            Positions (x,y) of oscillators, shape (N, 2).
        kernel_width : float
            Width of influence kernel.
        coupling_strength : float
            Coupling from oscillators to field.
        """
        N = len(oscillator_phases)

        # Compute local order parameter at each grid point
        new_values = np.zeros_like(self.values, dtype=complex)

        for i in range(N):
            x_osc, y_osc = oscillator_positions[i]

            # Distance from oscillator to each grid point
            dist_sq = (self.grid.X - x_osc)**2 + (self.grid.Y - y_osc)**2

            # Gaussian kernel
            kernel = np.exp(-dist_sq / (2 * kernel_width**2))

            # Add oscillator contribution
            new_values += kernel * np.exp(1j * oscillator_phases[i]) / N

        # Take magnitude for R field
        self.values = coupling_strength * np.abs(new_values)

    def diffuse(self, diffusion_coeff: float, dt: float):
        """
        Apply diffusion step: ∂field/∂t = D∇²field.

        Parameters
        ----------
        diffusion_coeff : float
            Diffusion coefficient D.
        dt : float
            Time step.
        """
        laplacian = self.grid.laplacian(self.values)
        self.values += diffusion_coeff * laplacian * dt
        self.t += dt

    def reaction_diffusion_step(
        self,
        reaction_func: Callable[[NDArray], NDArray],
        diffusion_coeff: float,
        dt: float
    ):
        """
        Single step of reaction-diffusion dynamics.

        ∂field/∂t = reaction(field) + D∇²field

        Parameters
        ----------
        reaction_func : callable
            Reaction term f(field).
        diffusion_coeff : float
            Diffusion coefficient.
        dt : float
            Time step.
        """
        # Reaction term
        reaction = reaction_func(self.values)

        # Diffusion term
        diffusion = diffusion_coeff * self.grid.laplacian(self.values)

        # Update
        self.values += (reaction + diffusion) * dt
        self.t += dt

    def evolve_pde(
        self,
        dynamics_func: Callable[[float, NDArray], NDArray],
        t_span: Tuple[float, float],
        dt: float = 0.01,
        store_interval: int = 10
    ):
        """
        Evolve field according to PDE dynamics.

        Parameters
        ----------
        dynamics_func : callable
            Function f(t, field) returning dfield/dt.
        t_span : tuple
            Time interval (t0, tf).
        dt : float
            Time step.
        store_interval : int
            Store solution every n steps.
        """
        t0, tf = t_span
        n_steps = int((tf - t0) / dt)

        # Reset history
        self.history = [self.values.copy()]
        self.time_history = [t0]

        # Time integration loop
        t = t0
        for step in range(n_steps):
            # Simple Euler for now
            dfield_dt = dynamics_func(t, self.values)
            self.values += dfield_dt * dt
            t += dt

            # Store history
            if (step + 1) % store_interval == 0:
                self.history.append(self.values.copy())
                self.time_history.append(t)

        self.t = t

    def compute_spatial_average(self) -> float:
        """Compute spatial average of field."""
        return np.mean(self.values)

    def compute_spatial_variance(self) -> float:
        """Compute spatial variance of field."""
        return np.var(self.values)

    def plot(self, ax: Optional[plt.Axes] = None, **kwargs):
        """
        Plot field as 2D heatmap.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. Creates new figure if None.
        **kwargs
            Additional arguments for imshow.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        im = ax.imshow(
            self.values.T,  # Transpose for correct orientation
            origin='lower',
            extent=[0, self.grid.Lx, 0, self.grid.Ly],
            aspect='auto',
            **kwargs
        )

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'{self.name}(x,y,t={self.t:.2f})')

        plt.colorbar(im, ax=ax)

        return ax

    def animate(
        self,
        interval: int = 50,
        save_path: Optional[str] = None
    ):
        """
        Animate field evolution from stored history.

        Parameters
        ----------
        interval : int
            Animation frame interval in ms.
        save_path : str, optional
            Path to save animation.
        """
        if not self.history:
            raise ValueError("No history stored. Run evolve_pde first.")

        fig, ax = plt.subplots(figsize=(6, 5))

        # Initial plot
        im = ax.imshow(
            self.history[0].T,
            origin='lower',
            extent=[0, self.grid.Lx, 0, self.grid.Ly],
            aspect='auto',
            cmap='viridis'
        )

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        title = ax.set_title(f'{self.name}(x,y,t=0.00)')

        plt.colorbar(im, ax=ax)

        def update(frame):
            im.set_data(self.history[frame].T)
            im.set_clim(vmin=self.history[frame].min(),
                       vmax=self.history[frame].max())
            title.set_text(f'{self.name}(x,y,t={self.time_history[frame]:.2f})')
            return [im, title]

        anim = FuncAnimation(
            fig, update,
            frames=len(self.history),
            interval=interval,
            blit=True
        )

        if save_path:
            anim.save(save_path, writer='pillow', fps=30)

        return anim

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"ScalarField(name='{self.name}', "
            f"grid={self.grid.Nx}x{self.grid.Ny}, "
            f"t={self.t:.2f})"
        )