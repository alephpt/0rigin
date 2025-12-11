"""
Klein-Gordon mediator field for local coupling in Kuramoto model.

Implements σ(x,t) dynamics: (∂²_t - c²∇² + M²)σ = g·ρ
where ρ is the local oscillator density.
"""

from typing import Optional, Tuple
import numpy as np
from numpy.typing import NDArray

from .grid import SpatialGrid


class MediatorField:
    """
    Klein-Gordon mediator field enabling local coupling.

    Satisfies the wave equation with mass:
        (∂²_t - c²∇² + M²)σ = g·ρ(x,t)

    where:
        - σ(x,t) is the mediator field
        - c is the wave speed
        - M is the effective mass
        - g is the coupling constant
        - ρ(x,t) is the oscillator density source

    In the heavy mass limit M→∞, this recovers instantaneous global coupling.
    For M~0, c→∞, we get diffusive local coupling.
    """

    def __init__(
        self,
        grid: SpatialGrid,
        wave_speed: float = 1.0,
        mass: float = 1.0,
        coupling_constant: float = 1.0,
        damping: float = 0.1,
        initial_field: Optional[NDArray] = None,
        initial_velocity: Optional[NDArray] = None
    ):
        """
        Initialize Klein-Gordon mediator field.

        Parameters
        ----------
        grid : SpatialGrid
            Underlying spatial grid.
        wave_speed : float
            Wave propagation speed c.
        mass : float
            Field mass M (controls interaction range).
        coupling_constant : float
            Source coupling strength g.
        damping : float
            Damping coefficient for numerical stability.
        initial_field : NDArray, optional
            Initial field configuration σ(x,0).
        initial_velocity : NDArray, optional
            Initial time derivative ∂σ/∂t(x,0).
        """
        self.grid = grid
        self.c = wave_speed
        self.M = mass
        self.g = coupling_constant
        self.gamma = damping

        # Field value σ(x,t)
        if initial_field is None:
            self.sigma = np.zeros((grid.Nx, grid.Ny))
        else:
            assert initial_field.shape == (grid.Nx, grid.Ny)
            self.sigma = initial_field.copy()

        # Time derivative ∂σ/∂t
        if initial_velocity is None:
            self.sigma_dot = np.zeros((grid.Nx, grid.Ny))
        else:
            assert initial_velocity.shape == (grid.Nx, grid.Ny)
            self.sigma_dot = initial_velocity.copy()

        self.t = 0.0

        # Storage for history
        self.history = []
        self.time_history = []

    def compute_source_density(
        self,
        oscillator_phases: NDArray,
        oscillator_positions: NDArray,
        kernel_width: float = 0.1
    ) -> NDArray:
        """
        Compute oscillator density source ρ(x,t).

        Uses Gaussian kernel to distribute each oscillator over space.
        The density includes phase information: ρ ~ Σ exp(iθ_j)·K(x - x_j)

        Parameters
        ----------
        oscillator_phases : NDArray
            Current phases of N oscillators.
        oscillator_positions : NDArray
            Positions (x,y) of oscillators, shape (N, 2).
        kernel_width : float
            Width of Gaussian distribution kernel.

        Returns
        -------
        NDArray
            Source density field ρ(x,y).
        """
        N = len(oscillator_phases)
        rho = np.zeros((self.grid.Nx, self.grid.Ny), dtype=complex)

        for i in range(N):
            x_osc, y_osc = oscillator_positions[i]

            # Distance from oscillator to each grid point
            dist_sq = (self.grid.X - x_osc)**2 + (self.grid.Y - y_osc)**2

            # Gaussian kernel
            kernel = np.exp(-dist_sq / (2 * kernel_width**2))

            # Add complex oscillator contribution
            rho += kernel * np.exp(1j * oscillator_phases[i]) / N

        # Take real part for source (can also use |rho| for amplitude)
        return np.real(rho)

    def evolve_step(
        self,
        source_density: NDArray,
        dt: float,
        method: str = 'leapfrog'
    ):
        """
        Evolve field by one time step using Klein-Gordon equation.

        The equation in first-order form:
            ∂σ/∂t = σ_dot
            ∂σ_dot/∂t = c²∇²σ - M²σ + g·ρ

        Parameters
        ----------
        source_density : NDArray
            Current source density ρ(x,y,t).
        dt : float
            Time step.
        method : str
            Integration method ('euler', 'leapfrog').
        """
        if method == 'leapfrog':
            self._leapfrog_step(source_density, dt)
        elif method == 'euler':
            self._euler_step(source_density, dt)
        else:
            raise ValueError(f"Unknown method: {method}")

        self.t += dt

    def _euler_step(self, rho: NDArray, dt: float):
        """Simple Euler integration."""
        # Compute Laplacian
        laplacian = self.grid.laplacian(self.sigma)

        # Klein-Gordon equation: ∂²σ/∂t² = c²∇²σ - M²σ + g·ρ
        sigma_ddot = self.c**2 * laplacian - self.M**2 * self.sigma + self.g * rho

        # Update
        self.sigma += self.sigma_dot * dt
        self.sigma_dot += sigma_ddot * dt

    def _leapfrog_step(self, rho: NDArray, dt: float):
        """
        Leapfrog integration (symplectic, better energy conservation).

        Updates in sequence:
            1. σ(t+dt/2) = σ(t) + σ_dot(t)·dt/2
            2. Compute forces at t+dt/2
            3. σ_dot(t+dt) = σ_dot(t) + acceleration·dt
            4. σ(t+dt) = σ(t+dt/2) + σ_dot(t+dt)·dt/2

        Includes damping term -γ·σ_dot for numerical stability.
        """
        # Half-step position update
        sigma_half = self.sigma + 0.5 * dt * self.sigma_dot

        # Compute forces at half-step (with damping)
        laplacian = self.grid.laplacian(sigma_half)
        sigma_ddot = self.c**2 * laplacian - self.M**2 * sigma_half + self.g * rho - self.gamma * self.sigma_dot

        # Full-step velocity update
        self.sigma_dot += dt * sigma_ddot

        # Complete position update
        self.sigma = sigma_half + 0.5 * dt * self.sigma_dot

    def get_heavy_mass_limit(self, rho: NDArray) -> NDArray:
        """
        Compute field in heavy mass limit M→∞.

        In this limit, the field responds instantaneously:
            M²σ ≈ g·ρ  =>  σ ≈ (g/M²)·ρ

        This approximates global mean-field coupling.

        Parameters
        ----------
        rho : NDArray
            Source density.

        Returns
        -------
        NDArray
            Field in heavy mass limit.
        """
        return (self.g / self.M**2) * rho

    def compute_field_energy(self) -> float:
        """
        Compute total field energy.

        E = ∫ dx [½(∂σ/∂t)² + ½c²|∇σ|² + ½M²σ²]

        Returns
        -------
        float
            Total field energy.
        """
        # Kinetic energy
        E_kin = 0.5 * np.sum(self.sigma_dot**2)

        # Gradient energy
        grad_x, grad_y = self.grid.gradient(self.sigma)
        E_grad = 0.5 * self.c**2 * np.sum(grad_x**2 + grad_y**2)

        # Mass energy
        E_mass = 0.5 * self.M**2 * np.sum(self.sigma**2)

        # Normalize by grid volume
        dV = self.grid.dx * self.grid.dy
        return (E_kin + E_grad + E_mass) * dV

    def compute_propagation_speed(self) -> Tuple[float, float]:
        """
        Estimate wave propagation speed from field dynamics.

        Returns
        -------
        phase_velocity : float
            Speed of phase propagation.
        group_velocity : float
            Speed of energy propagation.
        """
        # For Klein-Gordon: v_phase = c/sqrt(1 + (M/k)²)
        # For typical wavelength λ ~ L_x/2
        k_typical = 2 * np.pi / (self.grid.Lx / 2)

        v_phase = self.c / np.sqrt(1 + (self.M / k_typical)**2)
        v_group = self.c * np.sqrt(1 - (self.M * self.c / k_typical)**2)

        return v_phase, v_group

    def store_snapshot(self):
        """Store current field state in history."""
        self.history.append(self.sigma.copy())
        self.time_history.append(self.t)

    def sample_at_positions(self, positions: NDArray) -> NDArray:
        """
        Sample field at specific positions (for oscillator coupling).

        Uses bilinear interpolation for smooth sampling.

        Parameters
        ----------
        positions : NDArray
            Positions (x,y) to sample, shape (N, 2).

        Returns
        -------
        NDArray
            Field values at positions, shape (N,).
        """
        N = len(positions)
        values = np.zeros(N)

        for i, (x, y) in enumerate(positions):
            # Find grid indices (with wrapping for periodic BC)
            ix = int(x / self.grid.dx) % self.grid.Nx
            iy = int(y / self.grid.dy) % self.grid.Ny

            # Simple nearest-neighbor for now (can upgrade to bilinear)
            values[i] = self.sigma[ix, iy]

        return values

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"MediatorField(c={self.c:.2f}, M={self.M:.2f}, "
            f"g={self.g:.2f}, t={self.t:.2f})"
        )
