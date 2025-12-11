"""
Order parameter analysis for Kuramoto model.

Provides tools for computing and analyzing the Kuramoto order parameter
and related synchronization metrics.
"""

from typing import Tuple, Optional
import numpy as np
from numpy.typing import NDArray


class OrderParameter:
    """
    Kuramoto order parameter analyzer.

    Computes the complex order parameter:
        Z = R e^(iΨ) = (1/N) Σ_j e^(iθ_j)

    where R is the synchronization amplitude (0 ≤ R ≤ 1) and
    Ψ is the mean phase.

    Parameters
    ----------
    phases : NDArray[float]
        Phase trajectory array with shape (n_times, n_oscillators)
        or (n_oscillators,) for single snapshot.

    Attributes
    ----------
    phases : NDArray[float]
        Stored phase trajectory.
    N : int
        Number of oscillators.
    is_trajectory : bool
        True if input is time series, False if single snapshot.

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, (100, 50))  # 100 times, 50 oscillators
    >>> op = OrderParameter(phases)
    >>> R, Psi = op.time_series()
    >>> mean_R = op.mean_amplitude()
    """

    def __init__(self, phases: NDArray):
        """Initialize order parameter analyzer."""
        self.phases = np.asarray(phases)

        # Determine if trajectory or single snapshot
        if self.phases.ndim == 1:
            self.is_trajectory = False
            self.N = len(self.phases)
        elif self.phases.ndim == 2:
            self.is_trajectory = True
            self.N = self.phases.shape[1]
        else:
            raise ValueError("Phases must be 1D or 2D array")

    def compute(
        self,
        phases: Optional[NDArray] = None
    ) -> Tuple[float, float]:
        """
        Compute order parameter for given phases.

        Parameters
        ----------
        phases : NDArray[float], optional
            Phase array. Uses stored phases if None.

        Returns
        -------
        R : float
            Synchronization amplitude (0 ≤ R ≤ 1).
        Psi : float
            Mean phase (-π ≤ Ψ ≤ π).
        """
        if phases is None:
            if self.is_trajectory:
                raise ValueError(
                    "For trajectories, use time_series() method "
                    "or provide specific phase snapshot"
                )
            phases = self.phases
        else:
            phases = np.asarray(phases)

        # Complex order parameter
        z = np.mean(np.exp(1j * phases))

        return np.abs(z), np.angle(z)

    def time_series(self) -> Tuple[NDArray, NDArray]:
        """
        Compute order parameter time series.

        Returns
        -------
        R : NDArray[float]
            Synchronization amplitude at each time point.
        Psi : NDArray[float]
            Mean phase at each time point.
        """
        if not self.is_trajectory:
            # Single snapshot, return single values
            R, Psi = self.compute()
            return np.array([R]), np.array([Psi])

        # Compute for each time point
        z_trajectory = np.mean(np.exp(1j * self.phases), axis=1)

        R = np.abs(z_trajectory)
        Psi = np.angle(z_trajectory)

        return R, Psi

    def mean_amplitude(self) -> float:
        """
        Compute mean synchronization amplitude.

        Returns
        -------
        float
            Time-averaged R value.
        """
        R, _ = self.time_series()
        return np.mean(R)

    def steady_state_amplitude(
        self,
        fraction: float = 0.2
    ) -> float:
        """
        Compute steady-state order parameter amplitude.

        Averages R over the final fraction of the trajectory.

        Parameters
        ----------
        fraction : float, optional
            Fraction of trajectory to average (default: 0.2 = 20%).

        Returns
        -------
        float
            Steady-state R value.
        """
        R, _ = self.time_series()

        if len(R) == 1:
            return R[0]

        n_points = max(1, int(len(R) * fraction))
        return np.mean(R[-n_points:])

    def is_synchronized(
        self,
        threshold: float = 0.5,
        steady_state: bool = True,
        fraction: float = 0.2
    ) -> bool:
        """
        Check if system is synchronized.

        Parameters
        ----------
        threshold : float, optional
            Threshold for R to be considered synchronized (default: 0.5).
        steady_state : bool, optional
            Use steady-state R if True, mean R if False (default: True).
        fraction : float, optional
            Fraction for steady-state calculation (default: 0.2).

        Returns
        -------
        bool
            True if R > threshold.
        """
        if steady_state:
            R = self.steady_state_amplitude(fraction)
        else:
            R = self.mean_amplitude()

        return R > threshold

    def convergence_time(
        self,
        t: Optional[NDArray] = None,
        threshold: float = 0.01,
        window: int = 10
    ) -> Optional[float]:
        """
        Estimate time to reach steady state.

        Identifies when R stops changing significantly.

        Parameters
        ----------
        t : NDArray[float], optional
            Time array. Uses indices if None.
        threshold : float, optional
            Threshold for R variance to consider converged (default: 0.01).
        window : int, optional
            Window size for computing variance (default: 10).

        Returns
        -------
        float or None
            Convergence time, or None if not converged.
        """
        R, _ = self.time_series()

        if len(R) < window:
            return None

        # Compute rolling variance
        for i in range(window, len(R)):
            var = np.var(R[i-window:i])
            if var < threshold:
                return t[i] if t is not None else i

        return None

    def phase_locked_fraction(
        self,
        tolerance: float = 0.1
    ) -> float:
        """
        Compute fraction of oscillators phase-locked to mean field.

        An oscillator is considered phase-locked if its phase
        remains within tolerance of the mean phase.

        Parameters
        ----------
        tolerance : float, optional
            Phase tolerance in radians (default: 0.1).

        Returns
        -------
        float
            Fraction of phase-locked oscillators (0 to 1).
        """
        if not self.is_trajectory:
            return 0.0

        _, Psi = self.time_series()

        # Compute phase difference from mean for each oscillator
        phase_diff = self.phases - Psi[:, np.newaxis]

        # Wrap to [-π, π]
        phase_diff = np.angle(np.exp(1j * phase_diff))

        # Count oscillators within tolerance
        locked = np.mean(np.abs(phase_diff) < tolerance, axis=0)

        return np.mean(locked > 0.8)  # Locked >80% of the time

    def frequency_distribution(self) -> NDArray:
        """
        Estimate instantaneous frequency distribution.

        Computes dθ/dt by finite differences.

        Returns
        -------
        NDArray[float]
            Frequency array with shape (n_times-1, n_oscillators).
        """
        if not self.is_trajectory:
            raise ValueError("Frequency requires time series data")

        # Compute phase derivatives
        dtheta = np.diff(self.phases, axis=0)

        # Unwrap phase jumps (handle 2π wrapping)
        dtheta = np.angle(np.exp(1j * dtheta))

        return dtheta

    def __repr__(self) -> str:
        """String representation."""
        if self.is_trajectory:
            n_times = len(self.phases)
            return (
                f"OrderParameter(N={self.N} oscillators, "
                f"n_times={n_times}, "
                f"<R>={self.mean_amplitude():.3f})"
            )
        else:
            R, _ = self.compute()
            return f"OrderParameter(N={self.N} oscillators, R={R:.3f})"
