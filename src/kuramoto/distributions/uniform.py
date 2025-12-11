"""
Uniform frequency distribution.

The uniform distribution represents oscillators with frequencies
uniformly distributed over an interval.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray
from .base import FrequencyDistribution


class UniformDistribution(FrequencyDistribution):
    """
    Uniform frequency distribution.

    The PDF is given by:
        g(ω) = 1/(b - a)  for a ≤ ω ≤ b
        g(ω) = 0          otherwise

    For the Kuramoto model with uniform distribution over [-Δ, Δ],
    the critical coupling must be determined numerically. There is
    no simple analytical formula like for Lorentzian distributions.

    Parameters
    ----------
    low : float, optional
        Lower bound a (default: -1).
    high : float, optional
        Upper bound b (default: 1).

    Attributes
    ----------
    a : float
        Lower bound.
    b : float
        Upper bound.

    Examples
    --------
    >>> dist = UniformDistribution(low=-1, high=1)
    >>> frequencies = dist.sample(100)
    >>> mean = dist.mean()  # Returns 0.0
    """

    def __init__(self, low: float = -1.0, high: float = 1.0):
        """Initialize uniform distribution."""
        if low >= high:
            raise ValueError("Lower bound must be less than upper bound")
        self.a = low
        self.b = high
        self.low = low    # Alias
        self.high = high  # Alias

    def sample(self, N: int, seed: Optional[int] = None) -> NDArray:
        """
        Sample N frequencies from uniform distribution.

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

        return np.random.uniform(self.a, self.b, N)

    def pdf(self, omega: NDArray) -> NDArray:
        """
        Uniform probability density function.

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
        pdf_value = 1.0 / (self.b - self.a)
        result = np.zeros_like(omega)
        mask = (omega >= self.a) & (omega <= self.b)
        result[mask] = pdf_value
        return result

    def critical_coupling(self) -> None:
        """
        Critical coupling for uniform distribution.

        The critical coupling for uniform distribution must be
        determined numerically. There is no simple analytical formula.

        Returns
        -------
        None
            Returns None (no analytical formula available).
        """
        return None

    def mean(self) -> float:
        """
        Mean of the uniform distribution.

        Returns
        -------
        float
            Mean frequency (a + b)/2.
        """
        return (self.a + self.b) / 2

    def variance(self) -> float:
        """
        Variance of the uniform distribution.

        Returns
        -------
        float
            Variance (b - a)²/12.
        """
        return ((self.b - self.a)**2) / 12

    def __repr__(self) -> str:
        """String representation."""
        return f"UniformDistribution(low={self.a}, high={self.b})"
