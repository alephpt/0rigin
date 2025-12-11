"""
Base class for ODE solvers.

This module provides the abstract interface for numerical
integration methods used to solve the Kuramoto dynamics.
"""

from abc import ABC, abstractmethod
from typing import Callable, Tuple, Optional
import numpy as np
from numpy.typing import NDArray


class Solver(ABC):
    """
    Abstract base class for ODE solvers.

    All solvers must implement methods for single-step
    integration and full trajectory integration.
    """

    @abstractmethod
    def step(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> NDArray:
        """
        Perform single integration step.

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
        pass

    def integrate(
        self,
        func: Callable[[float, NDArray], NDArray],
        y0: NDArray,
        t_span: Tuple[float, float],
        dt: Optional[float] = None,
        t_eval: Optional[NDArray] = None,
        max_steps: int = 100000
    ) -> Tuple[NDArray, NDArray]:
        """
        Integrate ODE system over time span.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        y0 : NDArray
            Initial condition.
        t_span : tuple of float
            Time interval (t0, tf).
        dt : float, optional
            Time step. If None, uses adaptive stepping.
        t_eval : NDArray, optional
            Times at which to evaluate solution.
        max_steps : int, optional
            Maximum number of integration steps.

        Returns
        -------
        t : NDArray
            Time points.
        y : NDArray
            Solution trajectory, shape (len(t), len(y0)).
        """
        t0, tf = t_span

        if t_eval is not None:
            # Evaluate at specific times
            t_out = np.asarray(t_eval)
            n_points = len(t_out)
        else:
            # Determine output times based on dt
            if dt is None:
                dt = (tf - t0) / 1000  # Default to 1000 steps
            n_points = int((tf - t0) / dt) + 1
            t_out = np.linspace(t0, tf, n_points)

        # Initialize solution arrays
        y_out = np.zeros((n_points, len(y0)))
        y_out[0] = y0

        # Integration loop
        t_current = t0
        y_current = y0.copy()
        output_idx = 1

        step_count = 0
        while t_current < tf and output_idx < n_points:
            # Check max steps
            if step_count >= max_steps:
                raise RuntimeError(
                    f"Maximum number of steps ({max_steps}) exceeded"
                )

            # Determine next output time
            t_next_output = t_out[output_idx]

            # Adaptive stepping to hit output time exactly
            dt_step = min(dt if dt else (tf - t0) / 1000, t_next_output - t_current)

            # Take integration step
            y_current = self.step(func, t_current, y_current, dt_step)
            t_current += dt_step

            # Store if we've reached an output point
            if np.abs(t_current - t_next_output) < 1e-10:
                y_out[output_idx] = y_current
                output_idx += 1

            step_count += 1

        # Trim if we didn't fill all output points
        if output_idx < n_points:
            t_out = t_out[:output_idx]
            y_out = y_out[:output_idx]

        return t_out, y_out


class AdaptiveSolver(Solver):
    """
    Base class for adaptive step-size solvers.

    Provides infrastructure for error-controlled integration
    with automatic step size adjustment.

    Parameters
    ----------
    rtol : float, optional
        Relative tolerance (default: 1e-6).
    atol : float, optional
        Absolute tolerance (default: 1e-9).
    """

    def __init__(self, rtol: float = 1e-6, atol: float = 1e-9):
        """Initialize adaptive solver."""
        self.rtol = rtol
        self.atol = atol
        self.min_dt = 1e-10
        self.max_dt = 1.0
        self.safety_factor = 0.9

    @abstractmethod
    def step_with_error(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> Tuple[NDArray, NDArray]:
        """
        Perform step and estimate error.

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
            New state at time t + dt.
        error : NDArray
            Estimated error in y_new.
        """
        pass

    def step(
        self,
        func: Callable[[float, NDArray], NDArray],
        t: float,
        y: NDArray,
        dt: float
    ) -> NDArray:
        """
        Perform adaptive step with error control.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        t : float
            Current time.
        y : NDArray
            Current state.
        dt : float
            Initial time step suggestion.

        Returns
        -------
        NDArray
            New state at time t + dt_actual.
        """
        dt_try = dt
        while True:
            # Take step and estimate error
            y_new, error = self.step_with_error(func, t, y, dt_try)

            # Compute error norm
            scale = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(y_new))
            error_norm = np.sqrt(np.mean((error / scale)**2))

            # Accept or reject step
            if error_norm <= 1:
                # Accept step and return
                return y_new

            # Reject step and reduce dt
            dt_try *= self.safety_factor * (1 / error_norm)**0.2
            dt_try = max(dt_try, self.min_dt)

            if dt_try == self.min_dt:
                # Accept step despite error
                return y_new

    def compute_initial_step(
        self,
        func: Callable[[float, NDArray], NDArray],
        t0: float,
        y0: NDArray
    ) -> float:
        """
        Estimate good initial step size.

        Parameters
        ----------
        func : callable
            Function computing dy/dt = func(t, y).
        t0 : float
            Initial time.
        y0 : NDArray
            Initial state.

        Returns
        -------
        float
            Suggested initial step size.
        """
        # Estimate based on derivative
        f0 = func(t0, y0)
        scale = self.atol + self.rtol * np.abs(y0)
        d0 = np.sqrt(np.mean((y0 / scale)**2))
        d1 = np.sqrt(np.mean((f0 / scale)**2))

        if d0 < 1e-5 or d1 < 1e-5:
            h0 = 1e-6
        else:
            h0 = 0.01 * d0 / d1

        # Take small step to estimate second derivative
        y1 = y0 + h0 * f0
        f1 = func(t0 + h0, y1)
        d2 = np.sqrt(np.mean(((f1 - f0) / scale)**2)) / h0

        if d1 <= 1e-15 and d2 <= 1e-15:
            h1 = max(1e-6, h0 * 1e-3)
        else:
            h1 = (0.01 / max(d1, d2))**(1 / 2)

        return min(100 * h0, h1)