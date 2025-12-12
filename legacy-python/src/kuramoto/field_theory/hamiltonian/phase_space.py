"""
Hamiltonian formulation of Kuramoto model.

Extends discrete Kuramoto to phase space (θ, p) with conjugate momenta.
"""

from typing import Optional, Tuple, Union
import numpy as np
from numpy.typing import NDArray


class HamiltonianKuramoto:
    """
    Hamiltonian formulation of Kuramoto model with conjugate momenta.

    The dynamics follow:
        dθ/dt = p
        dp/dt = -γp + ω + (K/N)Σsin(θ_k - θ_j)

    In the overdamped limit (γ→∞), this recovers the original Kuramoto model.
    """

    def __init__(
        self,
        N: int,
        coupling_strength: float,
        frequencies: NDArray,
        damping: float = 1.0,
        initial_phases: Optional[NDArray] = None,
        initial_momenta: Optional[NDArray] = None
    ):
        """
        Initialize Hamiltonian Kuramoto model.

        Parameters
        ----------
        N : int
            Number of oscillators.
        coupling_strength : float
            Coupling constant K.
        frequencies : NDArray
            Natural frequencies ω_i.
        damping : float
            Damping coefficient γ (γ→∞ recovers standard Kuramoto).
        initial_phases : NDArray, optional
            Initial phases θ_i(0).
        initial_momenta : NDArray, optional
            Initial momenta p_i(0).
        """
        self.N = N
        self.K = coupling_strength
        self.frequencies = np.asarray(frequencies)
        self.gamma = damping

        # Initialize phase space
        if initial_phases is None:
            self.theta = np.random.uniform(0, 2*np.pi, N)
        else:
            self.theta = np.asarray(initial_phases)

        if initial_momenta is None:
            self.p = np.zeros(N)  # Start at rest
        else:
            self.p = np.asarray(initial_momenta)

        self.t = 0.0

    def compute_hamiltonian(self) -> float:
        """
        Compute system Hamiltonian H = T + V.

        In conservative case (γ=0), this should be conserved.

        Returns
        -------
        float
            Total Hamiltonian energy.
        """
        # Kinetic energy T = Σ p_i²/2
        T = 0.5 * np.sum(self.p**2)

        # Potential from natural frequencies V_freq = -Σ ω_i θ_i
        V_freq = -np.sum(self.frequencies * self.theta)

        # Interaction potential V_int = -(K/2N) Σ_ij cos(θ_i - θ_j)
        theta_diff = self.theta[:, np.newaxis] - self.theta[np.newaxis, :]
        V_int = -0.5 * self.K / self.N * np.sum(np.cos(theta_diff))

        return T + V_freq + V_int

    def equations_of_motion(self, t: float, state: NDArray) -> NDArray:
        """
        Compute Hamilton's equations with damping.

        Parameters
        ----------
        t : float
            Current time.
        state : NDArray
            State vector [θ_1,...,θ_N, p_1,...,p_N].

        Returns
        -------
        NDArray
            Time derivatives [dθ/dt, dp/dt].
        """
        # Unpack state
        theta = state[:self.N]
        p = state[self.N:]

        # dθ/dt = ∂H/∂p = p
        dtheta_dt = p

        # Coupling forces
        theta_diff = theta[:, np.newaxis] - theta[np.newaxis, :]
        coupling_force = self.K / self.N * np.sum(np.sin(theta_diff), axis=1)

        # dp/dt = -∂H/∂θ - γp = ω + coupling - γp
        dp_dt = self.frequencies + coupling_force - self.gamma * p

        return np.concatenate([dtheta_dt, dp_dt])

    def evolve(
        self,
        t_span: Tuple[float, float],
        dt: float = 0.01,
        method: str = 'rk4'
    ) -> dict:
        """
        Integrate equations of motion.

        Parameters
        ----------
        t_span : tuple
            Time interval (t0, tf).
        dt : float
            Time step.
        method : str
            Integration method ('euler', 'rk4').

        Returns
        -------
        dict
            Solution with time, phases, momenta, energy.
        """
        t0, tf = t_span
        n_steps = int((tf - t0) / dt)

        # Storage
        t_array = np.linspace(t0, tf, n_steps + 1)
        theta_trajectory = np.zeros((n_steps + 1, self.N))
        p_trajectory = np.zeros((n_steps + 1, self.N))
        energy_trajectory = np.zeros(n_steps + 1)

        # Initial conditions
        state = np.concatenate([self.theta, self.p])
        theta_trajectory[0] = self.theta
        p_trajectory[0] = self.p
        energy_trajectory[0] = self.compute_hamiltonian()

        # Integrate
        for i in range(n_steps):
            if method == 'euler':
                state = self._euler_step(t_array[i], state, dt)
            elif method == 'rk4':
                state = self._rk4_step(t_array[i], state, dt)
            else:
                raise ValueError(f"Unknown method: {method}")

            # Store
            self.theta = state[:self.N]
            self.p = state[self.N:]
            theta_trajectory[i+1] = self.theta
            p_trajectory[i+1] = self.p
            energy_trajectory[i+1] = self.compute_hamiltonian()

        self.t = tf

        # Compute order parameter
        z = np.mean(np.exp(1j * theta_trajectory), axis=1)
        R = np.abs(z)
        Psi = np.angle(z)

        return {
            't': t_array,
            'theta': theta_trajectory,
            'p': p_trajectory,
            'energy': energy_trajectory,
            'R': R,
            'Psi': Psi
        }

    def _euler_step(self, t: float, state: NDArray, dt: float) -> NDArray:
        """Simple Euler integration step."""
        return state + dt * self.equations_of_motion(t, state)

    def _rk4_step(self, t: float, state: NDArray, dt: float) -> NDArray:
        """4th-order Runge-Kutta step."""
        k1 = self.equations_of_motion(t, state)
        k2 = self.equations_of_motion(t + dt/2, state + dt*k1/2)
        k3 = self.equations_of_motion(t + dt/2, state + dt*k2/2)
        k4 = self.equations_of_motion(t + dt, state + dt*k3)
        return state + dt * (k1 + 2*k2 + 2*k3 + k4) / 6

    def get_overdamped_limit(self) -> NDArray:
        """
        Get effective dynamics in overdamped limit (γ→∞).

        Returns phases velocities as in standard Kuramoto.
        """
        # In overdamped limit: p ≈ 0, so dp/dt ≈ 0
        # This gives: 0 = -γp + ω + coupling
        # Therefore: p ≈ (ω + coupling)/γ
        # And: dθ/dt = p ≈ (ω + coupling)/γ

        theta_diff = self.theta[:, np.newaxis] - self.theta[np.newaxis, :]
        coupling = self.K / self.N * np.sum(np.sin(theta_diff), axis=1)

        return (self.frequencies + coupling) / self.gamma