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
        coupling: Union[float, Coupling] = 1.0,
        frequencies: Union[NDArray, FrequencyDistribution, str] = 'lorentzian',
        initial_phases: Optional[NDArray] = None
    ):
        """Initialize Kuramoto model.

        Raises
        ------
        ValueError
            If N <= 0, coupling < 0, or array shapes incompatible
        TypeError
            If arguments have wrong type
        """
        # Validate N
        if not isinstance(N, (int, np.integer)):
            raise TypeError(f"N must be an integer, got {type(N).__name__}")
        if N <= 0:
            raise ValueError(f"N must be positive, got {N}")

        self.N = int(N)
        self.t = 0.0

        # Set up coupling
        if isinstance(coupling, Coupling):
            self.coupling = coupling
        elif isinstance(coupling, (int, float, np.number)):
            if coupling < 0:
                raise ValueError(f"Coupling strength must be non-negative, got {coupling}")
            self.coupling = SinusoidalCoupling(strength=float(coupling))
        else:
            raise TypeError(
                f"Coupling must be a number or Coupling object, got {type(coupling).__name__}"
            )

        # Set up frequencies
        self.frequencies = self._initialize_frequencies(frequencies)
        if len(self.frequencies) != self.N:
            raise ValueError(
                f"Frequency array length ({len(self.frequencies)}) must match N ({self.N})"
            )

        # Set up initial phases
        if initial_phases is None:
            self.phases = np.random.uniform(0, 2 * np.pi, self.N)
        else:
            # Check for complex dtype before conversion
            initial_phases_arr = np.asarray(initial_phases)
            if np.iscomplexobj(initial_phases_arr):
                raise TypeError("Initial phases must be real, not complex")

            # Convert to float
            initial_phases = np.asarray(initial_phases, dtype=float)
            if len(initial_phases) != self.N:
                raise ValueError(
                    f"Initial phases length ({len(initial_phases)}) must match N ({self.N})"
                )
            self.phases = initial_phases

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
            try:
                freq_array = np.asarray(frequencies, dtype=float)
            except (ValueError, TypeError) as e:
                raise TypeError(f"Frequencies must be numeric array, got {type(frequencies).__name__}: {e}")

            if freq_array.ndim != 1:
                raise ValueError(
                    f"Frequency array must be 1D, got shape {freq_array.shape}"
                )

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

        Raises
        ------
        ValueError
            If t_span is invalid (t_end <= t_start) or dt < 0
        TypeError
            If t_span is not a tuple
        """
        # Validate t_span
        if not isinstance(t_span, tuple) or len(t_span) != 2:
            raise TypeError("t_span must be a tuple (t_start, t_end)")

        t_start, t_end = t_span
        try:
            t_start = float(t_start)
            t_end = float(t_end)
        except (ValueError, TypeError):
            raise TypeError("t_span values must be numeric")

        if t_end <= t_start:
            raise ValueError(
                f"Invalid time span: t_end ({t_end}) must be greater than t_start ({t_start})"
            )

        # Validate dt if provided
        if dt is not None:
            try:
                dt = float(dt)
            except (ValueError, TypeError):
                raise TypeError("dt must be numeric")
            if dt < 0:
                raise ValueError(f"dt must be non-negative, got {dt}")

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