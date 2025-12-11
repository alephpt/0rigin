"""
Self-consistent Mean Field Theory (SMFT) system integration.

Couples discrete oscillators to continuous mediator fields, bridging
the gap between discrete Kuramoto dynamics and field theory.
"""

from typing import Optional, Tuple, Literal
import numpy as np
from numpy.typing import NDArray

from .fields import SpatialGrid, ScalarField, MediatorField
from .hamiltonian import HamiltonianKuramoto
from .coupling import LocalFieldCoupling, FermionMassDemo


class SMFTSystem:
    """
    Self-consistent Mean Field Theory system.

    Couples N discrete oscillators θⱼ(t) to continuous mediator field σ(x,t),
    creating a self-consistent dynamics where:
    - Oscillators influence the field through local coupling
    - Field influences oscillators through spatial averaging

    In the heavy mass limit (M→∞), recovers standard Kuramoto model.

    Parameters
    ----------
    grid_shape : tuple of int
        Spatial grid dimensions (Nx, Ny).
    N_oscillators : int
        Number of discrete oscillators.
    coupling : str
        Coupling type: 'local' or 'global'.
    mediator_mass : float
        Mass M of mediator field (M→∞ recovers Kuramoto).
    oscillator_frequencies : NDArray, optional
        Natural frequencies. Random if None.
    grid_size : tuple of float, optional
        Physical domain size (Lx, Ly).
    boundary : str, optional
        Boundary conditions for grid.
    """

    def __init__(
        self,
        grid_shape: Tuple[int, int],
        N_oscillators: int,
        coupling: Literal['local', 'global'] = 'local',
        mediator_mass: float = 10.0,
        oscillator_frequencies: Optional[NDArray] = None,
        grid_size: Tuple[float, float] = (1.0, 1.0),
        boundary: str = 'periodic'
    ):
        """Initialize SMFT system."""
        self.grid_shape = grid_shape
        self.N = N_oscillators
        self.coupling_type = coupling
        self.M = mediator_mass

        # Create spatial grid
        self.grid = SpatialGrid(
            Nx=grid_shape[0],
            Ny=grid_shape[1],
            Lx=grid_size[0],
            Ly=grid_size[1],
            boundary=boundary
        )

        # Initialize oscillator positions (random or structured)
        self.oscillator_positions = self._initialize_positions()

        # Initialize oscillator frequencies
        if oscillator_frequencies is None:
            oscillator_frequencies = np.random.normal(0, 1, N_oscillators)

        # Create Hamiltonian oscillator system
        self.oscillators = HamiltonianKuramoto(
            N=N_oscillators,
            coupling_strength=1.0,  # Coupling strength (mean-field scaling applied in step)
            frequencies=oscillator_frequencies,
            damping=5.0  # Increased damping for stability
        )

        # Create mediator field σ(x,t)
        self.mediator_field = ScalarField(
            self.grid,
            name="mediator_sigma"
        )

        # Create synchronization field R(x,t)
        self.sync_field = ScalarField(
            self.grid,
            name="sync_R"
        )

        # Create phase field θ(x,t)
        self.phase_field = ScalarField(
            self.grid,
            name="phase_theta"
        )

        # Effective mass field m_eff(x,t)
        self.mass_field = None

        # Time
        self.t = 0.0

    def _initialize_positions(self) -> NDArray:
        """Initialize oscillator positions in domain."""
        if self.coupling_type == 'local':
            # Random positions for local coupling
            positions = np.random.rand(self.N, 2)
            positions[:, 0] *= self.grid.Lx
            positions[:, 1] *= self.grid.Ly
        else:
            # Doesn't matter for global coupling
            positions = np.random.rand(self.N, 2)

        return positions

    def compute_local_order_parameter(
        self,
        kernel_width: float = 0.1
    ) -> Tuple[NDArray, NDArray]:
        """
        Compute local order parameter R(x,y) and θ(x,y).

        Projects discrete oscillator phases onto spatial grid using
        Gaussian kernel weighting.

        Parameters
        ----------
        kernel_width : float
            Spatial width of influence kernel.

        Returns
        -------
        R_field : NDArray
            Local synchronization amplitude.
        theta_field : NDArray
            Local mean phase.
        """
        # Compute complex order parameter field
        z_field = np.zeros((self.grid.Nx, self.grid.Ny), dtype=complex)

        for i in range(self.N):
            x_osc, y_osc = self.oscillator_positions[i]

            # Distance from oscillator to each grid point
            dist_sq = (self.grid.X - x_osc)**2 + (self.grid.Y - y_osc)**2

            # Gaussian kernel
            kernel = np.exp(-dist_sq / (2 * kernel_width**2))
            kernel /= np.sum(kernel)  # Normalize

            # Add oscillator contribution
            z_field += kernel * np.exp(1j * self.oscillators.theta[i])

        # Extract amplitude and phase
        R_field = np.abs(z_field)
        theta_field = np.angle(z_field)

        return R_field, theta_field

    def update_mediator_field(
        self,
        dt: float,
        diffusion_coeff: float = 0.01,
        kernel_width: float = 0.1
    ):
        """
        Update mediator field σ(x,t) from oscillator configuration.

        Field equation:
            M d²σ/dt² = -γ dσ/dt - ∂V/∂σ + D∇²σ + source(oscillators)

        Parameters
        ----------
        dt : float
            Time step.
        diffusion_coeff : float
            Spatial diffusion coefficient D.
        kernel_width : float
            Width of oscillator influence kernel.
        """
        # Compute local order parameter
        R_field, theta_field = self.compute_local_order_parameter(kernel_width)

        # Update sync and phase fields
        self.sync_field.values = R_field
        self.phase_field.values = theta_field

        # Simple mediator dynamics: relax toward local order parameter
        # In full theory, this would be Klein-Gordon equation
        damping = 5.0  # Increased damping for stability (matches oscillator damping)
        relaxation_rate = 1.0 / self.M

        # Source from oscillators
        source = relaxation_rate * R_field

        # Diffusion term
        if diffusion_coeff > 0:
            diffusion = diffusion_coeff * self.grid.laplacian(self.mediator_field.values)
        else:
            diffusion = 0

        # Semi-implicit damping for stability
        # Explicit update would be: σ_new = σ + dt * (-damping * (σ - source) + diffusion)
        # Semi-implicit: σ_new = (σ + dt * (-damping * source + diffusion)) / (1 + damping * dt)
        dσ_dt_undamped = -damping * source + diffusion
        σ_temp = self.mediator_field.values + dσ_dt_undamped * dt
        self.mediator_field.values = σ_temp / (1 + damping * dt)
        self.mediator_field.t += dt

    def compute_field_force_on_oscillators(
        self,
        kernel_width: float = 0.1
    ) -> NDArray:
        """
        Compute force on oscillators from mediator field.

        Each oscillator samples the field at its position, weighted
        by spatial kernel.

        Parameters
        ----------
        kernel_width : float
            Sampling kernel width.

        Returns
        -------
        NDArray
            Force contribution for each oscillator.
        """
        forces = np.zeros(self.N)

        for i in range(self.N):
            x_osc, y_osc = self.oscillator_positions[i]

            # Distance from oscillator to each grid point
            dist_sq = (self.grid.X - x_osc)**2 + (self.grid.Y - y_osc)**2

            # Gaussian sampling kernel
            kernel = np.exp(-dist_sq / (2 * kernel_width**2))
            kernel /= np.sum(kernel)

            # Sample field value
            field_value = np.sum(kernel * self.mediator_field.values)

            # Force is proportional to field value
            # In full theory: F_i = -∂H/∂θ_i ∝ -σ(x_i) sin(θ_i - Ψ(x_i))
            # Note: negative sign for restoring force
            phase_diff = self.oscillators.theta[i] - np.sum(kernel * self.phase_field.values)
            forces[i] = -field_value * np.sin(phase_diff)

        return forces

    def step(
        self,
        dt: float,
        diffusion_coeff: float = 0.01,
        kernel_width: float = 0.1
    ):
        """
        Single coupled time step.

        Updates both oscillators and field self-consistently:
        1. Compute field force on oscillators
        2. Advance oscillators with field coupling
        3. Update field from new oscillator configuration

        Parameters
        ----------
        dt : float
            Time step.
        diffusion_coeff : float
            Spatial diffusion in field.
        kernel_width : float
            Spatial coupling kernel width.
        """
        # Get field force on oscillators
        field_forces = self.compute_field_force_on_oscillators(kernel_width)

        # Advance oscillators (simple Euler for now)
        state = np.concatenate([self.oscillators.theta, self.oscillators.p])

        # Modified equations with field coupling
        dtheta_dt = self.oscillators.p

        # Standard Hamiltonian coupling
        theta_diff = self.oscillators.theta[:, np.newaxis] - self.oscillators.theta[np.newaxis, :]
        coupling_force = self.oscillators.K / self.N * np.sum(np.sin(theta_diff), axis=1)

        # Add field coupling
        total_force = coupling_force + field_forces

        # Semi-implicit damping for stability (same as mediator field)
        dp_dt_undamped = self.oscillators.frequencies + total_force

        # Update
        self.oscillators.theta += dtheta_dt * dt
        # Apply semi-implicit damping: p_new = (p_old + dp_dt_undamped * dt) / (1 + gamma * dt)
        p_temp = self.oscillators.p + dp_dt_undamped * dt
        self.oscillators.p = p_temp / (1 + self.oscillators.gamma * dt)
        self.oscillators.t += dt

        # Update mediator field from new oscillator state
        self.update_mediator_field(dt, diffusion_coeff, kernel_width)

        self.t += dt

    def evolve(
        self,
        t_span: Tuple[float, float],
        dt: float = 0.01,
        diffusion_coeff: float = 0.01,
        kernel_width: float = 0.1,
        store_interval: int = 10
    ) -> dict:
        """
        Evolve coupled oscillator-field system.

        Parameters
        ----------
        t_span : tuple
            Time interval (t0, tf).
        dt : float
            Time step.
        diffusion_coeff : float
            Field diffusion coefficient.
        kernel_width : float
            Spatial coupling kernel width.
        store_interval : int
            Store state every n steps.

        Returns
        -------
        dict
            Solution with time, oscillator phases, field values, energy.
        """
        t0, tf = t_span
        n_steps = int((tf - t0) / dt)

        # Storage
        n_store = n_steps // store_interval
        t_array = np.zeros(n_store)
        theta_trajectory = np.zeros((n_store, self.N))
        R_global = np.zeros(n_store)
        mediator_history = []
        sync_history = []
        energy_history = np.zeros(n_store)

        # Time evolution
        store_idx = 0
        for step in range(n_steps):
            # Coupled step
            self.step(dt, diffusion_coeff, kernel_width)

            # Store
            if (step + 1) % store_interval == 0:
                t_array[store_idx] = self.t
                theta_trajectory[store_idx] = self.oscillators.theta

                # Global order parameter
                z = np.mean(np.exp(1j * self.oscillators.theta))
                R_global[store_idx] = np.abs(z)

                # Field snapshots
                mediator_history.append(self.mediator_field.values.copy())
                sync_history.append(self.sync_field.values.copy())

                # Energy
                energy_history[store_idx] = self.oscillators.compute_hamiltonian()

                store_idx += 1

        return {
            't': t_array,
            'theta': theta_trajectory,
            'R': R_global,
            'mediator_field': np.array(mediator_history),
            'sync_field': np.array(sync_history),
            'energy': energy_history,
            'positions': self.oscillator_positions
        }

    def compute_effective_mass(self) -> NDArray:
        """
        Compute effective mass field m_eff(x,t).

        In field theory, effective mass emerges from curvature of
        potential around equilibrium.

        Returns
        -------
        NDArray
            Effective mass field.
        """
        # Simple model: m_eff ∝ 1 / R(x)
        # Where R is low (incoherent), mass is high (slow dynamics)
        epsilon = 0.01  # Regularization
        m_eff = self.M / (self.sync_field.values + epsilon)

        self.mass_field = m_eff
        return m_eff

    def test_heavy_mass_limit(self, M_values: list) -> dict:
        """
        Test that M→∞ recovers standard Kuramoto.

        Parameters
        ----------
        M_values : list
            List of mediator masses to test.

        Returns
        -------
        dict
            Results for each mass value.
        """
        results = {}

        for M in M_values:
            # Save original mass
            M_original = self.M

            # Set new mass
            self.M = M

            # Reset system
            self.oscillators.reset()
            self.mediator_field.values[:] = 0
            self.t = 0.0

            # Evolve
            solution = self.evolve((0, 20), dt=0.01, store_interval=20)

            # Store final R
            results[M] = solution['R'][-1]

            # Restore
            self.M = M_original

        return results

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"SMFTSystem(grid={self.grid_shape}, "
            f"N={self.N}, coupling='{self.coupling_type}', "
            f"M={self.M:.1f}, t={self.t:.2f})"
        )
