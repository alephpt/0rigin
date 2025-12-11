"""
Gaussian (Normal) frequency distribution.

The Gaussian distribution is commonly used in Kuramoto model studies
and has a well-defined critical coupling constant.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray
from .base import FrequencyDistribution


class GaussianDistribution(FrequencyDistribution):
    """
    Gaussian (Normal) frequency distribution.

    The PDF is given by:
        g(ω) = (1/(σ√(2π))) exp[-(ω - μ)²/(2σ²)]

    For the Kuramoto model with Gaussian distributed frequencies,
    the critical coupling is approximately:
        Kc ≈ √(8/π) σ ≈ 1.596 σ

    This approximation comes from self-consistency analysis near
    the bifurcation point.

    Parameters
    ----------
    mean : float, optional
        Mean frequency μ (default: 0).
    std : float, optional
        Standard deviation σ (default: 1).

    Attributes
    ----------
    mu : float
        Mean frequency.
    sigma : float
        Standard deviation.

    Examples
    --------
    >>> dist = GaussianDistribution(mean=0, std=1)
    >>> frequencies = dist.sample(100)
    >>> Kc = dist.critical_coupling()  # Returns ~1.596
    """

    def __init__(self, mean: float = 0.0, std: float = 1.0):
        """Initialize Gaussian distribution."""
        if std <= 0:
            raise ValueError("Standard deviation must be positive")
        self.mu = mean
        self.sigma = std
        self.mean_val = mean  # Alias
        self.std_val = std    # Alias

    def sample(self, N: int, seed: Optional[int] = None) -> NDArray:
        """
        Sample N frequencies from Gaussian distribution.

        Uses the Box-Muller transform (via numpy's normal distribution).

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

        return np.random.normal(self.mu, self.sigma, N)

    def pdf(self, omega: NDArray) -> NDArray:
        """
        Gaussian probability density function.

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
        normalization = 1.0 / (self.sigma * np.sqrt(2 * np.pi))
        exponent = -((omega - self.mu)**2) / (2 * self.sigma**2)
        return normalization * np.exp(exponent)

    def critical_coupling(self) -> float:
        """
        Return approximate critical coupling for Gaussian distribution.

        The critical coupling is approximately:
            Kc ≈ √(8/π) σ ≈ 1.596 σ

        This is derived from self-consistency analysis assuming
        the order parameter is small near the transition.

        Returns
        -------
        float
            Critical coupling Kc ≈ 1.596σ.
        """
        return np.sqrt(8 / np.pi) * self.sigma

    def mean(self) -> float:
        """
        Mean of the Gaussian distribution.

        Returns
        -------
        float
            Mean frequency μ.
        """
        return self.mu

    def variance(self) -> float:
        """
        Variance of the Gaussian distribution.

        Returns
        -------
        float
            Variance σ².
        """
        return self.sigma**2

    def __repr__(self) -> str:
        """String representation."""
        return f"GaussianDistribution(mean={self.mu}, std={self.sigma})"
