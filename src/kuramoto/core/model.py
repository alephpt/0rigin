"""
Core Kuramoto model implementation.

This module provides the main KuramotoModel class for simulating
the dynamics of coupled phase oscillators.
"""

from typing import Optional, Union, Tuple, Callable
import numpy as np
from numpy.typing import NDArray

from ..distributions.base import FrequencyDistribution
from .coupling import Coupling, SinusoidalCoupling
from ..solvers.base import Solver


class KuramotoModel:
    """
    Kuramoto model of coupled phase oscillators.

    The model simulates N coupled oscillators with dynamics:
        dθ_j/dt = ω_j + (K/N) Σ_k sin(θ_k - θ_j)

    Parameters
    ----------
    N : int
        Number of oscillators.
    coupling : float or Coupling
        Coupling strength K (if float) or coupling object.
    frequencies : array-like, FrequencyDistribution, or str
        Natural frequencies ω_i. Can be:
        - Array of N frequencies
        - FrequencyDistribution object
        - String: 'lorentzian', 'gaussian', 'uniform'
    initial_phases : array-like, optional
        Initial phase configuration. Random if None.

    Attributes
    ----------
    N : int
        Number of oscillators.
    coupling : Coupling
        Coupling object.
    frequencies : NDArray[float]
        Natural frequencies of oscillators.
    phases : NDArray[float]
        Current phase configuration.
    t : float
        Current time.

    Examples
    --------
    >>> model = KuramotoModel(N=100, coupling=2.0, frequencies='lorentzian')
    >>> solution = model.evolve((0, 50))
    >>> R, Psi = model.compute_order_parameter()
    """

    def __init__(
        self,
        N: int,
        coupling: Union[float, Coupling],
        frequencies: Union[NDArray, FrequencyDistribution, str],
        initial_phases: Optional[NDArray] = None
    ):
        """Initialize Kuramoto model."""
        self.N = N
        self.t = 0.0

        # Set up coupling
        if isinstance(coupling, (int, float)):
            self.coupling = SinusoidalCoupling(strength=coupling)
        else:
            self.coupling = coupling

        # Set up frequencies
        self.frequencies = self._initialize_frequencies(frequencies)

        # Set up initial phases
        if initial_phases is None:
            self.phases = np.random.uniform(0, 2 * np.pi, N)
        else:
            self.phases = np.asarray(initial_phases)
            assert len(self.phases) == N, "Initial phases must have length N"

        # Storage for solution history
        self._solution_t = None
        self._solution_phases = None

    def _initialize_frequencies(
        self,
        frequencies: Union[NDArray, FrequencyDistribution, str]
    ) -> NDArray:
        """Initialize frequency array from various input types."""
        if isinstance(frequencies, str):
            # Create distribution from string
            from ..distributions import (
                LorentzianDistribution,
                GaussianDistribution,
                UniformDistribution
            )
            dist_map = {
                'lorentzian': LorentzianDistribution(),
                'gaussian': GaussianDistribution(),
                'uniform': UniformDistribution(),
            }
            if frequencies not in dist_map:
                raise ValueError(
                    f"Unknown distribution: {frequencies}. "
                    f"Available: {list(dist_map.keys())}"
                )
            return dist_map[frequencies].sample(self.N)

        elif isinstance(frequencies, FrequencyDistribution):
            # Sample from distribution
            return frequencies.sample(self.N)

        else:
            # Direct array
            freq_array = np.asarray(frequencies)
            assert len(freq_array) == self.N, "Frequency array must have length N"
            return freq_array

    def equations_of_motion(self, t: float, phases: NDArray) -> NDArray:
        """
        Compute dθ/dt for current phase configuration.

        Parameters
        ----------
        t : float
            Current time (unused in autonomous system).
        phases : NDArray[float]
            Current phase configuration.

        Returns
        -------
        NDArray[float]
            Time derivatives dθ/dt for each oscillator.
        """
        # Compute coupling field
        coupling_field = self.coupling.compute_field(phases)

        # Add natural frequencies
        return self.frequencies + coupling_field

    def evolve(
        self,
        t_span: Tuple[float, float],
        solver: Union[str, Solver] = 'rk45',
        dt: Optional[float] = None,
        store_trajectory: bool = True,
        **solver_kwargs
    ) -> dict:
        """
        Evolve system over time span.

        Parameters
        ----------
        t_span : tuple of float
            Time interval (t0, tf).
        solver : str or Solver
            Solver to use. Options: 'rk45', 'rk4', 'euler'.
        dt : float, optional
            Time step (for fixed-step solvers).
        store_trajectory : bool
            Whether to store full trajectory.
        **solver_kwargs
            Additional arguments passed to solver.

        Returns
        -------
        dict
            Solution dictionary with keys:
            - 't': Time points
            - 'phases': Phase trajectories
            - 'R': Order parameter amplitude
            - 'Psi': Order parameter phase
        """
        # Initialize solver
        if isinstance(solver, str):
            solver_obj = self._get_solver(solver, **solver_kwargs)
        else:
            solver_obj = solver

        # Integrate
        t_eval, phases_trajectory = solver_obj.integrate(
            self.equations_of_motion,
            self.phases,
            t_span,
            dt=dt
        )

        # Update current state
        self.t = t_eval[-1]
        self.phases = phases_trajectory[-1]

        # Store trajectory if requested
        if store_trajectory:
            self._solution_t = t_eval
            self._solution_phases = phases_trajectory

        # Compute order parameter trajectory
        R_trajectory = []
        Psi_trajectory = []
        for phases_t in phases_trajectory:
            R, Psi = self.compute_order_parameter(phases_t)
            R_trajectory.append(R)
            Psi_trajectory.append(Psi)

        return {
            't': t_eval,
            'phases': phases_trajectory,
            'R': np.array(R_trajectory),
            'Psi': np.array(Psi_trajectory)
        }

    def _get_solver(self, name: str, **kwargs) -> Solver:
        """Get solver object by name."""
        from ..solvers import RK4Solver, RK45Solver, EulerSolver

        solvers = {
            'rk45': RK45Solver,
            'rk4': RK4Solver,
            'euler': EulerSolver
        }

        if name not in solvers:
            raise ValueError(f"Unknown solver: {name}")

        return solvers[name](**kwargs)

    def compute_order_parameter(
        self,
        phases: Optional[NDArray] = None
    ) -> Tuple[float, float]:
        """
        Calculate Kuramoto order parameter.

        The complex order parameter is defined as:
            Z = R e^(iΨ) = (1/N) Σ_j e^(iθ_j)

        Parameters
        ----------
        phases : NDArray[float], optional
            Phase configuration. Uses current phases if None.

        Returns
        -------
        R : float
            Synchronization amplitude (0 ≤ R ≤ 1).
        Psi : float
            Mean phase (-π ≤ Ψ ≤ π).
        """
        if phases is None:
            phases = self.phases

        # Complex order parameter
        z = np.mean(np.exp(1j * phases))

        return np.abs(z), np.angle(z)

    def reset(self, initial_phases: Optional[NDArray] = None):
        """
        Reset model to initial conditions.

        Parameters
        ----------
        initial_phases : NDArray[float], optional
            New initial phases. Random if None.
        """
        self.t = 0.0
        if initial_phases is None:
            self.phases = np.random.uniform(0, 2 * np.pi, self.N)
        else:
            self.phases = np.asarray(initial_phases)

        self._solution_t = None
        self._solution_phases = None

    @property
    def critical_coupling(self) -> Optional[float]:
        """
        Get theoretical critical coupling if known.

        Returns
        -------
        float or None
            Critical coupling Kc, or None if not analytically known.
        """
        if hasattr(self.coupling, 'strength'):
            # Try to get from frequency distribution
            if hasattr(self, '_freq_distribution'):
                return self._freq_distribution.critical_coupling()
        return None

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"KuramotoModel(N={self.N}, "
            f"coupling={self.coupling}, "
            f"t={self.t:.2f})"
        )