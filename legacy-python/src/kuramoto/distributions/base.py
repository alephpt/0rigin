"""
Base class for frequency distributions.

This module provides the abstract interface that all frequency
distributions must implement.
"""

from abc import ABC, abstractmethod
from typing import Optional
import numpy as np
from numpy.typing import NDArray


class FrequencyDistribution(ABC):
    """
    Abstract base class for frequency distributions.

    All frequency distributions must implement methods for
    sampling, computing the PDF, and (if known) providing
    the critical coupling value.
    """

    @abstractmethod
    def sample(self, N: int, seed: Optional[int] = None) -> NDArray:
        """
        Sample N frequencies from the distribution.

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
        pass

    @abstractmethod
    def pdf(self, omega: NDArray) -> NDArray:
        """
        Probability density function.

        Parameters
        ----------
        omega : NDArray[float]
            Frequency values at which to evaluate PDF.

        Returns
        -------
        NDArray[float]
            PDF values g(ω).
        """
        pass

    def critical_coupling(self) -> Optional[float]:
        """
        Analytical critical coupling if known.

        Returns
        -------
        float or None
            Critical coupling Kc, or None if not known analytically.
        """
        return None

    def mean(self) -> float:
        """
        Mean frequency of the distribution.

        Returns
        -------
        float
            Mean frequency ⟨ω⟩.
        """
        return 0.0  # Default to centered

    def variance(self) -> Optional[float]:
        """
        Variance of the distribution.

        Returns
        -------
        float or None
            Variance σ², or None if infinite.
        """
        return None

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}()"