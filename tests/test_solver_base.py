"""
Comprehensive test suite for ODE solver base classes.

Tests the Solver and AdaptiveSolver abstract base classes and their
implementations, including error handling, stability, and accuracy.
"""

import sys
import os
import numpy as np
from numpy.testing import assert_allclose
import unittest

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from kuramoto.solvers.base import Solver, AdaptiveSolver


class ConcreteSolver(Solver):
    """Concrete implementation for testing base class."""

    def step(self, func, t, y, dt):
        """Simple Euler step for testing."""
        return y + dt * func(t, y)


class ConcreteAdaptiveSolver(AdaptiveSolver):
    """Concrete adaptive solver for testing."""

    def step_with_error(self, func, t, y, dt):
        """Euler step with error estimate (half of second-order term)."""
        f = func(t, y)
        y_new = y + dt * f
        # Rough error estimate
        error = 0.5 * dt**2 * np.abs(f) * 0.1  # Simplified error
        return y_new, error


class TestSolverBase(unittest.TestCase):
    """Test base Solver class functionality."""

    def setUp(self):
        """Set up test solver."""
        self.solver = ConcreteSolver()

    def test_step_method_exists(self):
        """Test that step method is implemented."""
        self.assertTrue(hasattr(self.solver, 'step'))
        self.assertTrue(callable(self.solver.step))

    def test_simple_ode_step(self):
        """Test single step on simple ODE: dy/dt = -y."""
        def func(t, y):
            return -y

        t, y, dt = 0.0, np.array([1.0]), 0.01
        y_new = self.solver.step(func, t, y, dt)

        # Euler step: y_new = y - dt*y = y(1 - dt)
        expected = y * (1 - dt)
        assert_allclose(y_new, expected, rtol=1e-10)

    def test_integrate_constant_function(self):
        """Test integration of dy/dt = c (constant)."""
        c = 2.0

        def func(t, y):
            return np.array([c])

        y0 = np.array([0.0])
        t_span = (0, 1)
        dt = 0.01

        t, y = self.solver.integrate(func, y0, t_span, dt=dt)

        # Solution: y(t) = c*t + y0
        expected = c * t[:, np.newaxis] + y0
        assert_allclose(y, expected, rtol=0.01)  # Euler has O(dt) error

    def test_integrate_exponential_decay(self):
        """Test integration of exponential decay: dy/dt = -λy."""
        lambda_val = 2.0

        def func(t, y):
            return -lambda_val * y

        y0 = np.array([5.0])
        t_span = (0, 2)
        dt = 0.001

        t, y = self.solver.integrate(func, y0, t_span, dt=dt)

        # Analytical solution: y(t) = y0 * exp(-λt)
        expected = y0 * np.exp(-lambda_val * t[:, np.newaxis])
        assert_allclose(y, expected, rtol=0.01)  # Euler error

    def test_integrate_multidimensional(self):
        """Test integration with multiple state variables."""
        def func(t, y):
            # Simple coupled system
            return np.array([-y[0] + y[1], y[0] - y[1]])

        y0 = np.array([1.0, 0.0])
        t_span = (0, 1)
        dt = 0.01

        t, y = self.solver.integrate(func, y0, t_span, dt=dt)

        # Check shape
        self.assertEqual(y.shape[0], len(t))
        self.assertEqual(y.shape[1], len(y0))
        # Check initial condition
        assert_allclose(y[0], y0, rtol=1e-10)

    def test_integrate_with_t_eval(self):
        """Test integration with specific evaluation times."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (0, 2)
        t_eval = np.array([0, 0.5, 1.0, 1.5, 2.0])

        t, y = self.solver.integrate(func, y0, t_span, t_eval=t_eval, dt=0.001)

        # Should evaluate at requested times
        assert_allclose(t, t_eval, rtol=1e-10)
        self.assertEqual(len(y), len(t_eval))

    def test_max_steps_limit(self):
        """Test that max_steps limit is enforced."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (0, 100)
        dt = 0.0001  # Would need 1M steps

        with self.assertRaises(RuntimeError) as context:
            self.solver.integrate(func, y0, t_span, dt=dt, max_steps=1000)

        self.assertIn("Maximum number of steps", str(context.exception))

    def test_integrate_zero_span(self):
        """Test integration over zero time span."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (0, 0)

        t, y = self.solver.integrate(func, y0, t_span, dt=0.1)

        # Should return just initial condition
        self.assertEqual(len(t), 1)
        assert_allclose(y[0], y0, rtol=1e-10)

    def test_negative_time_integration(self):
        """Test backward integration (negative time)."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (2, 0)  # Backward
        dt = -0.01  # Negative step

        t, y = self.solver.integrate(func, y0, t_span, dt=dt)

        # Should integrate backward
        self.assertTrue(np.all(t[1:] < t[:-1]))  # Decreasing time


class TestAdaptiveSolverBase(unittest.TestCase):
    """Test adaptive solver base class."""

    def setUp(self):
        """Set up test adaptive solver."""
        self.solver = ConcreteAdaptiveSolver(rtol=1e-6, atol=1e-9)

    def test_initialization(self):
        """Test adaptive solver initialization."""
        self.assertEqual(self.solver.rtol, 1e-6)
        self.assertEqual(self.solver.atol, 1e-9)
        self.assertGreater(self.solver.min_dt, 0)
        self.assertLess(self.solver.max_dt, np.inf)
        self.assertLess(self.solver.safety_factor, 1.0)

    def test_step_with_error_exists(self):
        """Test that step_with_error method exists."""
        self.assertTrue(hasattr(self.solver, 'step_with_error'))
        self.assertTrue(callable(self.solver.step_with_error))

    def test_error_control_accept(self):
        """Test that step is accepted when error is small."""
        def func(t, y):
            return -y

        t, y = 0.0, np.array([1.0])
        dt = 0.001  # Small step should be accepted

        y_new = self.solver.step(func, t, y, dt)

        # Step should be taken
        self.assertIsNotNone(y_new)
        self.assertNotEqual(y_new[0], y[0])

    def test_error_control_reject_retry(self):
        """Test that large errors cause step rejection."""
        # Modify solver to have very tight tolerance
        solver = ConcreteAdaptiveSolver(rtol=1e-12, atol=1e-15)
        solver.min_dt = 1e-10

        def func(t, y):
            return -100 * y  # Stiff equation

        t, y = 0.0, np.array([1.0])
        dt = 1.0  # Way too large

        # Should still return something (after retries)
        y_new = solver.step(func, t, y, dt)
        self.assertIsNotNone(y_new)

    def test_compute_initial_step(self):
        """Test initial step size computation."""
        def func(t, y):
            return -y

        t0 = 0.0
        y0 = np.array([1.0, 2.0])

        h0 = self.solver.compute_initial_step(func, t0, y0)

        # Should return reasonable step size
        self.assertGreater(h0, 0)
        self.assertLess(h0, 1.0)  # Not too large

    def test_compute_initial_step_zero_derivative(self):
        """Test initial step when derivative is zero."""
        def func(t, y):
            return np.zeros_like(y)

        t0 = 0.0
        y0 = np.array([1.0])

        h0 = self.solver.compute_initial_step(func, t0, y0)

        # Should return small but non-zero step
        self.assertGreater(h0, 0)
        self.assertLess(h0, 0.01)

    def test_adaptive_integration(self):
        """Test full adaptive integration."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (0, 2)

        # Don't specify dt - let adaptive solver decide
        t, y = self.solver.integrate(func, y0, t_span)

        # Should complete integration
        self.assertGreater(len(t), 1)
        self.assertAlmostEqual(t[-1], t_span[1], places=5)

        # Check accuracy against analytical solution
        expected = y0 * np.exp(-t[-1])
        assert_allclose(y[-1], expected, rtol=0.01)

    def test_tolerance_effect(self):
        """Test that tighter tolerance gives more accurate results."""
        def func(t, y):
            return -y

        y0 = np.array([1.0])
        t_span = (0, 1)
        t_eval = np.array([1.0])

        # Loose tolerance
        solver_loose = ConcreteAdaptiveSolver(rtol=1e-3, atol=1e-6)
        t1, y1 = solver_loose.integrate(func, y0, t_span, t_eval=t_eval)

        # Tight tolerance
        solver_tight = ConcreteAdaptiveSolver(rtol=1e-9, atol=1e-12)
        t2, y2 = solver_tight.integrate(func, y0, t_span, t_eval=t_eval)

        # Analytical solution
        expected = y0 * np.exp(-1.0)

        # Tight tolerance should be more accurate
        error_loose = np.abs(y1[-1] - expected)
        error_tight = np.abs(y2[-1] - expected)
        self.assertLess(error_tight[0], error_loose[0])


class TestSolverStability(unittest.TestCase):
    """Test numerical stability of solvers."""

    def test_oscillatory_system(self):
        """Test solver on oscillatory system (harmonic oscillator)."""
        def harmonic_oscillator(t, y):
            # y = [position, velocity]
            # dy/dt = [velocity, -position] (ω=1)
            return np.array([y[1], -y[0]])

        solver = ConcreteSolver()
        y0 = np.array([1.0, 0.0])  # Start at x=1, v=0
        t_span = (0, 2*np.pi)
        dt = 0.01

        t, y = solver.integrate(harmonic_oscillator, y0, t_span, dt=dt)

        # Should complete one period
        # Energy should be approximately conserved
        energy = 0.5 * (y[:, 0]**2 + y[:, 1]**2)
        initial_energy = 0.5 * (y0[0]**2 + y0[1]**2)

        # Euler method has energy drift, but should be bounded
        max_drift = np.max(np.abs(energy - initial_energy))
        self.assertLess(max_drift, 0.5)  # Reasonable bound for Euler

    def test_stiff_equation(self):
        """Test solver behavior on stiff equation."""
        def stiff_func(t, y):
            return -1000 * (y - np.cos(t))

        solver = ConcreteAdaptiveSolver(rtol=1e-6, atol=1e-9)
        y0 = np.array([1.0])
        t_span = (0, 1)

        # Should handle stiff equation (adaptive stepping helps)
        t, y = solver.integrate(stiff_func, y0, t_span)

        # Solution should exist and be stable
        self.assertTrue(np.all(np.isfinite(y)))
        self.assertGreater(len(t), 10)  # Should take multiple steps

    def test_conservation_property(self):
        """Test that conservative systems maintain invariants."""
        def pendulum(t, y):
            # Nonlinear pendulum: θ'' = -sin(θ)
            # y = [θ, θ']
            return np.array([y[1], -np.sin(y[0])])

        solver = ConcreteAdaptiveSolver(rtol=1e-8, atol=1e-10)
        y0 = np.array([0.1, 0.0])  # Small angle
        t_span = (0, 10)

        t, y = solver.integrate(pendulum, y0, t_span)

        # Energy: E = 0.5*θ'^2 - cos(θ)
        energy = 0.5 * y[:, 1]**2 - np.cos(y[:, 0])
        initial_energy = 0.5 * y0[1]**2 - np.cos(y0[0])

        # Energy should be approximately conserved
        energy_variation = np.std(energy)
        self.assertLess(energy_variation, 0.01)


class TestErrorHandling(unittest.TestCase):
    """Test error handling in solvers."""

    def test_invalid_function(self):
        """Test handling of function that returns wrong shape."""
        def bad_func(t, y):
            return np.array([1.0, 2.0, 3.0])  # Wrong size

        solver = ConcreteSolver()
        y0 = np.array([1.0])
        t_span = (0, 1)

        # Should raise error or handle gracefully
        with self.assertRaises((ValueError, IndexError)):
            solver.integrate(bad_func, y0, t_span, dt=0.1)

    def test_nan_in_function(self):
        """Test handling of function that returns NaN."""
        def nan_func(t, y):
            if t > 0.5:
                return np.array([np.nan])
            return -y

        solver = ConcreteSolver()
        y0 = np.array([1.0])
        t_span = (0, 1)

        t, y = solver.integrate(nan_func, y0, t_span, dt=0.1)

        # Integration should stop or handle NaN
        # (Implementation-specific behavior)
        self.assertTrue(len(t) > 0)  # At least initial point

    def test_infinite_derivative(self):
        """Test handling of infinite derivatives."""
        def inf_func(t, y):
            if np.abs(y[0]) < 0.01:
                return np.array([np.inf])
            return -y

        solver = ConcreteAdaptiveSolver()
        y0 = np.array([0.005])  # Will trigger infinity
        t_span = (0, 1)

        # Should handle without crashing
        # (may stop early or use min_dt)
        t, y = solver.integrate(inf_func, y0, t_span)
        self.assertTrue(len(t) > 0)


class TestSolverAccuracy(unittest.TestCase):
    """Test accuracy of integration methods."""

    def test_linear_system_accuracy(self):
        """Test accuracy on linear system with known solution."""
        A = np.array([[-1, 2], [-2, -1]])

        def linear_system(t, y):
            return A @ y

        solver = ConcreteAdaptiveSolver(rtol=1e-8, atol=1e-10)
        y0 = np.array([1.0, 0.0])
        t_final = 1.0
        t_span = (0, t_final)

        t, y = solver.integrate(linear_system, y0, t_span, t_eval=[t_final])

        # Analytical solution: y(t) = exp(A*t) @ y0
        from scipy.linalg import expm
        expected = expm(A * t_final) @ y0

        # Should be accurate for linear system
        assert_allclose(y[-1], expected, rtol=0.01)

    def test_quadrature_accuracy(self):
        """Test solver as numerical integrator."""
        # Integrate f(t) by solving dy/dt = f(t), y(0) = 0
        def f(t):
            return np.sin(t)

        def ode(t, y):
            return np.array([f(t)])

        solver = ConcreteAdaptiveSolver(rtol=1e-10, atol=1e-12)
        y0 = np.array([0.0])
        t_span = (0, np.pi)

        t, y = solver.integrate(ode, y0, t_span)

        # Integral of sin(t) from 0 to π is 2
        expected = 2.0
        assert_allclose(y[-1, 0], expected, rtol=1e-6)


if __name__ == '__main__':
    unittest.main(verbosity=2)