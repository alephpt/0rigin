"""
Runge-Kutta numerical integration methods.

Implements fixed-step RK4 and adaptive RK45 (Dormand-Prince) solvers.
"""

from typing import Callable, Tuple
import numpy as np
from numpy.typing import NDArray
from .base import Solver, AdaptiveSolver


class RK4Solver(Solver):
    """
    Fourth-order Runge-Kutta solver with fixed step size.

    Classic RK4 method with O(dt^4) local error.
    Suitable for smooth problems where step size is known.
    """

    def step(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> NDArray:
        """
        Perform single RK4 integration step.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        t : float
            Current time.
        y : NDArray
            Current state.
        dt : float
            Time step.

        Returns
        -------
        NDArray
            New state at time t + dt.
        """
        k1 = func(t, y)
        k2 = func(t + dt / 2, y + dt * k1 / 2)
        k3 = func(t + dt / 2, y + dt * k2 / 2)
        k4 = func(t + dt, y + dt * k3)

        return y + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


class RK45Solver(AdaptiveSolver):
    """
    Adaptive Runge-Kutta-Fehlberg 4(5) solver (Dormand-Prince).

    Uses embedded RK method to estimate error and adapt step size.
    Efficient for problems with varying time scales.

    Parameters
    ----------
    rtol : float, optional
        Relative tolerance (default: 1e-6).
    atol : float, optional
        Absolute tolerance (default: 1e-9).
    """

    # Dormand-Prince coefficients
    A = np.array([
        [0, 0, 0, 0, 0, 0, 0],
        [1/5, 0, 0, 0, 0, 0, 0],
        [3/40, 9/40, 0, 0, 0, 0, 0],
        [44/45, -56/15, 32/9, 0, 0, 0, 0],
        [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0],
        [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0],
        [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    ])

    # 5th order solution (higher accuracy)
    B = np.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0])

    # 4th order solution (for error estimation)
    B_star = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])

    # Time coefficients
    C = np.array([0, 1/5, 3/10, 4/5, 8/9, 1, 1])

    def __init__(self, rtol: float = 1e-6, atol: float = 1e-9):
        """Initialize RK45 solver with tolerances."""
        super().__init__(rtol=rtol, atol=atol)

    def step_with_error(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> Tuple[NDArray, NDArray]:
        """
        Perform RK45 step and estimate error.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        t : float
            Current time.
        y : NDArray
            Current state.
        dt : float
            Time step.

        Returns
        -------
        y_new : NDArray
            New state at time t + dt (5th order).
        error : NDArray
            Estimated error (difference between 4th and 5th order).
        """
        # Compute k values
        k = np.zeros((7, len(y)))
        k[0] = func(t, y)

        for i in range(1, 7):
            y_temp = y + dt * np.sum(self.A[i, :i, np.newaxis] * k[:i], axis=0)
            k[i] = func(t + self.C[i] * dt, y_temp)

        # 5th order solution
        y_new = y + dt * np.sum(self.B[:, np.newaxis] * k, axis=0)

        # 4th order solution
        y_star = y + dt * np.sum(self.B_star[:, np.newaxis] * k, axis=0)

        # Error estimate
        error = np.abs(y_new - y_star)

        return y_new, error


class EulerSolver(Solver):
    """
    Simple forward Euler solver.

    First-order method: y_{n+1} = y_n + dt * f(t_n, y_n)
    Low accuracy but simple and stable for small dt.
    Mainly for testing and educational purposes.
    """

    def step(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> NDArray:
        """
        Perform single Euler step.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        t : float
            Current time.
        y : NDArray
            Current state.
        dt : float
            Time step.

        Returns
        -------
        NDArray
            New state at time t + dt.
        """
        return y + dt * func(t, y)
