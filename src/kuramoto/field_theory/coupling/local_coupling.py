"""
Bidirectional coupling between oscillators and fields.

Implements the SMFT coupling mechanism:
    - Oscillators → Field: via density ρ(x,t)
    - Field → Oscillators: via local sampling σ(x_j,t)
"""

from typing import Optional, Dict, Tuple
import numpy as np
from numpy.typing import NDArray

from ..fields.grid import SpatialGrid
from ..fields.mediator import MediatorField
from ..fields.scalar_field import ScalarField


class LocalFieldCoupling:
    """
    Manages bidirectional coupling between discrete oscillators and fields.

    System dynamics:
        1. Oscillators evolve: dθ_j/dt = ω_j + λ·σ(x_j)
        2. Oscillators source field: (∂²_t - c²∇² + M²)σ = g·ρ
        3. Field feeds back to oscillators

    This creates self-consistent local dynamics.
    """

    def __init__(
        self,
        grid: SpatialGrid,
        mediator_params: Optional[Dict] = None,
        oscillator_to_field_coupling: float = 1.0,
        field_to_oscillator_coupling: float = 1.0,
        kernel_width: float = 0.1
    ):
        """
        Initialize local field coupling system.

        Parameters
        ----------
        grid : SpatialGrid
            Spatial grid for fields.
        mediator_params : dict, optional
            Parameters for mediator field (wave_speed, mass, coupling_constant).
        oscillator_to_field_coupling : float
            Strength of oscillator → field coupling (g).
        field_to_oscillator_coupling : float
            Strength of field → oscillator coupling (λ).
        kernel_width : float
            Width of spatial kernels for density.
        """
        self.grid = grid
        self.kernel_width = kernel_width
        self.g_osc_to_field = oscillator_to_field_coupling
        self.lambda_field_to_osc = field_to_oscillator_coupling

        # Initialize mediator field
        if mediator_params is None:
            mediator_params = {
                'wave_speed': 1.0,
                'mass': 1.0,
                'coupling_constant': oscillator_to_field_coupling
            }

        self.mediator = MediatorField(grid, **mediator_params)

        # Optional: track order parameter field R(x,t)
        self.R_field = ScalarField(grid, name="R_field")
        self.theta_field = ScalarField(grid, name="theta_field")

    def update_fields_from_oscillators(
        self,
        phases: NDArray,
        positions: NDArray,
        dt: float
    ):
        """
        Update all fields based on current oscillator state.

        1. Compute density ρ(x,t) from oscillators
        2. Evolve mediator field σ
        3. Update order parameter field R

        Parameters
        ----------
        phases : NDArray
            Current oscillator phases, shape (N,).
        positions : NDArray
            Oscillator positions, shape (N, 2).
        dt : float
            Time step for field evolution.
        """
        # Compute source density
        rho = self.mediator.compute_source_density(
            phases, positions, self.kernel_width
        )

        # Evolve mediator field
        self.mediator.evolve_step(rho, dt, method='leapfrog')

        # Update order parameter field
        self.R_field.update_from_oscillators(
            phases, positions,
            kernel_width=self.kernel_width,
            coupling_strength=1.0
        )

    def compute_coupling_to_oscillators(
        self,
        positions: NDArray
    ) -> NDArray:
        """
        Compute field influence on oscillators.

        For each oscillator at x_j, sample the mediator field:
            coupling_j = λ·σ(x_j)

        Parameters
        ----------
        positions : NDArray
            Oscillator positions, shape (N, 2).

        Returns
        -------
        NDArray
            Coupling contribution to each oscillator, shape (N,).
        """
        # Sample mediator at oscillator positions
        sigma_at_positions = self.mediator.sample_at_positions(positions)

        return self.lambda_field_to_osc * sigma_at_positions

    def evolve_coupled_system(
        self,
        phases: NDArray,
        positions: NDArray,
        frequencies: NDArray,
        dt: float,
        n_steps: int,
        store_interval: int = 10
    ) -> Dict:
        """
        Evolve fully coupled oscillator-field system.

        Alternates between:
            1. Update oscillators using current field
            2. Update field using current oscillators

        Parameters
        ----------
        phases : NDArray
            Initial oscillator phases, shape (N,).
        positions : NDArray
            Oscillator positions (fixed), shape (N, 2).
        frequencies : NDArray
            Natural frequencies, shape (N,).
        dt : float
            Time step.
        n_steps : int
            Number of steps to evolve.
        store_interval : int
            Store snapshots every n steps.

        Returns
        -------
        dict
            Evolution history with phases, fields, order parameter.
        """
        N = len(phases)

        # Storage
        t_array = np.arange(n_steps + 1) * dt
        phase_history = []
        R_history = []
        mediator_history = []
        energy_history = []

        # Initial state
        current_phases = phases.copy()
        phase_history.append(current_phases.copy())
        R_history.append(np.abs(np.mean(np.exp(1j * current_phases))))
        mediator_history.append(self.mediator.sigma.copy())
        energy_history.append(self.mediator.compute_field_energy())

        # Evolution loop
        for step in range(n_steps):
            # 1. Compute field coupling to oscillators
            field_coupling = self.compute_coupling_to_oscillators(positions)

            # 2. Update oscillators (simple Euler)
            # dθ/dt = ω + λ·σ(x_j)
            dphases_dt = frequencies + field_coupling
            current_phases += dphases_dt * dt
            current_phases = np.mod(current_phases, 2 * np.pi)

            # 3. Update fields from new oscillator state
            self.update_fields_from_oscillators(current_phases, positions, dt)

            # 4. Store snapshots
            if (step + 1) % store_interval == 0:
                phase_history.append(current_phases.copy())
                R_history.append(np.abs(np.mean(np.exp(1j * current_phases))))
                mediator_history.append(self.mediator.sigma.copy())
                energy_history.append(self.mediator.compute_field_energy())

        return {
            't': t_array[::store_interval],
            'phases': np.array(phase_history),
            'R': np.array(R_history),
            'mediator_field': np.array(mediator_history),
            'field_energy': np.array(energy_history),
            'R_field': self.R_field.values.copy()
        }

    def compute_continuum_limit_error(
        self,
        phases: NDArray,
        positions: NDArray
    ) -> float:
        """
        Estimate error in continuum approximation.

        Compares discrete density ρ to smooth field approximation.

        Parameters
        ----------
        phases : NDArray
            Oscillator phases.
        positions : NDArray
            Oscillator positions.

        Returns
        -------
        float
            RMS error between discrete and continuum.
        """
        # Discrete density
        rho_discrete = self.mediator.compute_source_density(
            phases, positions, self.kernel_width
        )

        # Smooth field approximation
        R_smooth = self.R_field.values

        # RMS difference
        error = np.sqrt(np.mean((rho_discrete - R_smooth)**2))

        return error

    def get_effective_coupling_range(self) -> float:
        """
        Compute effective range of coupling.

        For Klein-Gordon with mass M, the screening length is:
            ξ = c/M

        Returns
        -------
        float
            Effective coupling range.
        """
        return self.mediator.c / self.mediator.M

    def test_heavy_mass_limit(
        self,
        phases: NDArray,
        positions: NDArray
    ) -> Tuple[float, float]:
        """
        Test that M→∞ recovers global coupling.

        Parameters
        ----------
        phases : NDArray
            Oscillator phases.
        positions : NDArray
            Oscillator positions.

        Returns
        -------
        global_R : float
            Global order parameter.
        local_R_mean : float
            Mean of local order parameter field.
        """
        # Global order parameter
        global_R = np.abs(np.mean(np.exp(1j * phases)))

        # Local field average
        self.R_field.update_from_oscillators(
            phases, positions,
            kernel_width=self.kernel_width
        )
        local_R_mean = np.mean(self.R_field.values)

        return global_R, local_R_mean

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"LocalFieldCoupling(g={self.g_osc_to_field:.2f}, "
            f"λ={self.lambda_field_to_osc:.2f}, "
            f"ξ={self.get_effective_coupling_range():.2f})"
        )
