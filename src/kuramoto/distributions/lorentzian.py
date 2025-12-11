"""
Lorentzian (Cauchy) frequency distribution.

The Lorentzian distribution is special for the Kuramoto model
as it admits exact solutions via the Ott-Antonsen ansatz.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray
from .base import FrequencyDistribution


class LorentzianDistribution(FrequencyDistribution):
    """
    Lorentzian (Cauchy) frequency distribution.

    The PDF is given by:
        g(ω) = (γ/π) / [(ω - ω₀)² + γ²]

    This distribution is special because:
    1. It allows exact solution via Ott-Antonsen reduction
    2. Critical coupling is exactly Kc = 2γ
    3. The order parameter satisfies a closed equation

    Parameters
    ----------
    center : float, optional
        Center frequency ω₀ (default: 0).
    width : float, optional
        Half-width at half-maximum γ (default: 1).

    Attributes
    ----------
    omega_0 : float
        Center frequency.
    gamma : float
        Width parameter.

    Examples
    --------
    >>> dist = LorentzianDistribution(center=0, width=1)
    >>> frequencies = dist.sample(100)
    >>> Kc = dist.critical_coupling()  # Returns 2.0
    """

    def __init__(self, center: float = 0.0, width: float = 1.0):
        """Initialize Lorentzian distribution."""
        if width <= 0:
            raise ValueError("Width must be positive")
        self.omega_0 = center
        self.gamma = width
        self.center = center  # Alias
        self.width = width    # Alias

    def sample(self, N: int, seed: Optional[int] = None) -> NDArray:
        """
        Sample N frequencies from Lorentzian distribution.

        Uses the inverse transform method with the Cauchy quantile function:
            ω = ω₀ + γ * tan(π(u - 0.5))

        Parameters
        ----------
        N : int
            Number of frequencies to sample.
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        NDArray[float]
            Array of N sampled frequencies.
        """
        if seed is not None:
            np.random.seed(seed)

        # Sample uniform random numbers
        u = np.random.uniform(0, 1, N)

        # Apply inverse CDF (quantile function) of Cauchy distribution
        frequencies = self.omega_0 + self.gamma * np.tan(np.pi * (u - 0.5))

        return frequencies

    def pdf(self, omega: NDArray) -> NDArray:
        """
        Lorentzian probability density function.

        Parameters
        ----------
        omega : NDArray[float]
            Frequency values at which to evaluate PDF.

        Returns
        -------
        NDArray[float]
            PDF values g(ω).
        """
        omega = np.asarray(omega)
        return (self.gamma / np.pi) / ((omega - self.omega_0)**2 + self.gamma**2)

    def critical_coupling(self) -> float:
        """
        Return exact critical coupling for Lorentzian distribution.

        For the Lorentzian distribution, the critical coupling
        is exactly Kc = 2γ, derived from the Ott-Antonsen theory.

        Returns
        -------
        float
            Critical coupling Kc = 2γ.
        """
        return 2 * self.gamma

    def mean(self) -> float:
        """
        Mean of the Lorentzian distribution.

        Returns
        -------
        float
            Mean frequency ω₀.
        """
        return self.omega_0

    def variance(self) -> None:
        """
        Variance of the Lorentzian distribution.

        The Lorentzian distribution has infinite variance.

        Returns
        -------
        None
            Always returns None (infinite variance).
        """
        return None  # Infinite variance

    def ott_antonsen_dynamics(self, z: complex, K: float) -> complex:
        """
        Compute Ott-Antonsen dynamics for this distribution.

        The OA reduction for Lorentzian g(ω) gives:
            dz/dt = -γz + (K/2)z(1 - |z|²)

        where z is the complex order parameter.

        Parameters
        ----------
        z : complex
            Current complex order parameter.
        K : float
            Coupling strength.

        Returns
        -------
        complex
            Time derivative dz/dt.
        """
        return -self.gamma * z + (K / 2) * z * (1 - np.abs(z)**2)

    def steady_state_order_parameter(self, K: float) -> float:
        """
        Compute steady-state order parameter.

        From Ott-Antonsen theory, the steady-state value is:
            R = 0 if K ≤ Kc
            R = √(1 - (Kc/K)²) if K > Kc

        Parameters
        ----------
        K : float
            Coupling strength.

        Returns
        -------
        float
            Steady-state order parameter R.
        """
        Kc = self.critical_coupling()

        if K <= Kc:
            return 0.0
        else:
            return np.sqrt(1 - (Kc / K)**2)

    def linear_stability_eigenvalue(self, K: float) -> float:
        """
        Compute linear stability eigenvalue of incoherent state.

        The incoherent state (R=0) loses stability when
        λ = -γ + K/2 becomes positive.

        Parameters
        ----------
        K : float
            Coupling strength.

        Returns
        -------
        float
            Largest eigenvalue λ.
        """
        return -self.gamma + K / 2

    def __repr__(self) -> str:
        """String representation."""
        return f"LorentzianDistribution(center={self.omega_0}, width={self.gamma})"