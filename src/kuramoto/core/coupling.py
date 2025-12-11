"""
Coupling functions for the Kuramoto model.

This module provides various coupling schemes including all-to-all
sinusoidal coupling and network-based coupling.
"""

from abc import ABC, abstractmethod
from typing import Optional
import numpy as np
from numpy.typing import NDArray


class Coupling(ABC):
    """
    Abstract base class for coupling functions.

    All coupling classes must implement the compute_field method
    that returns the coupling contribution to each oscillator's
    frequency.
    """

    @abstractmethod
    def compute_field(self, phases: NDArray) -> NDArray:
        """
        Compute coupling field from phase configuration.

        Parameters
        ----------
        phases : NDArray[float]
            Current phases of all oscillators.

        Returns
        -------
        NDArray[float]
            Coupling contribution to dθ/dt for each oscillator.
        """
        pass


class SinusoidalCoupling(Coupling):
    """
    Standard all-to-all sinusoidal coupling.

    Implements the classic Kuramoto coupling:
        H_j = (K/N) Σ_k sin(θ_k - θ_j)

    Parameters
    ----------
    strength : float
        Coupling strength K.
    """

    def __init__(self, strength: float):
        """Initialize sinusoidal coupling."""
        self.strength = strength
        self.K = strength  # Alias for compatibility

    def compute_field(self, phases: NDArray) -> NDArray:
        """
        Compute sinusoidal coupling field.

        Parameters
        ----------
        phases : NDArray[float]
            Current phases of all oscillators.

        Returns
        -------
        NDArray[float]
            Coupling field H_j for each oscillator.
        """
        N = len(phases)

        # Vectorized computation using broadcasting
        # phase_diff[i,j] = phases[j] - phases[i]
        phase_diff = phases[np.newaxis, :] - phases[:, np.newaxis]

        # Sum over columns (axis=1) to get coupling for each oscillator
        coupling_field = (self.K / N) * np.sum(np.sin(phase_diff), axis=1)

        return coupling_field

    def __repr__(self) -> str:
        """String representation."""
        return f"SinusoidalCoupling(K={self.strength})"


class NetworkCoupling(Coupling):
    """
    Coupling on arbitrary network topology.

    Implements network-coupled Kuramoto model:
        H_j = (K/k_j) Σ_{k∈N(j)} sin(θ_k - θ_j)

    where N(j) are the neighbors of node j and k_j is its degree.

    Parameters
    ----------
    adjacency_matrix : NDArray
        Network adjacency matrix (can be weighted).
    strength : float
        Coupling strength K.
    normalize : bool, optional
        Whether to normalize by node degree (default: True).
    """

    def __init__(
        self,
        adjacency_matrix: NDArray,
        strength: float,
        normalize: bool = True
    ):
        """Initialize network coupling."""
        self.adjacency = np.asarray(adjacency_matrix)
        self.strength = strength
        self.normalize = normalize
        self.N = len(self.adjacency)

        # Precompute node degrees for normalization
        if normalize:
            self.degrees = np.sum(self.adjacency, axis=1)
            # Avoid division by zero for isolated nodes
            self.degrees = np.where(self.degrees > 0, self.degrees, 1)
        else:
            self.degrees = np.ones(self.N)

    def compute_field(self, phases: NDArray) -> NDArray:
        """
        Compute network coupling field.

        Parameters
        ----------
        phases : NDArray[float]
            Current phases of all oscillators.

        Returns
        -------
        NDArray[float]
            Coupling field H_j for each oscillator.
        """
        # Phase differences for all pairs
        phase_diff = phases[np.newaxis, :] - phases[:, np.newaxis]

        # Apply adjacency mask and compute weighted sum
        coupling_field = np.sum(
            self.adjacency * np.sin(phase_diff),
            axis=1
        )

        # Normalize by degree and apply strength
        coupling_field = self.strength * coupling_field / self.degrees

        return coupling_field

    def __repr__(self) -> str:
        """String representation."""
        return f"NetworkCoupling(N={self.N}, K={self.strength})"


class PulseCoupling(Coupling):
    """
    Pulse-coupled oscillators (Peskin model variant).

    Instead of continuous sinusoidal coupling, oscillators
    interact via discrete pulses when they fire (reach θ = 2π).

    Parameters
    ----------
    strength : float
        Pulse strength ε.
    pulse_shape : callable, optional
        Function defining pulse response.
    """

    def __init__(
        self,
        strength: float,
        pulse_shape: Optional[callable] = None
    ):
        """Initialize pulse coupling."""
        self.strength = strength
        self.pulse_shape = pulse_shape or self._default_pulse

    def _default_pulse(self, phase: float) -> float:
        """Default pulse response function."""
        # Jump by fixed amount, with saturation at 2π
        return min(phase + self.strength, 2 * np.pi)

    def compute_field(self, phases: NDArray) -> NDArray:
        """
        Compute pulse coupling field.

        Note: This is a simplified continuous approximation.
        True pulse coupling requires event-driven simulation.

        Parameters
        ----------
        phases : NDArray[float]
            Current phases of all oscillators.

        Returns
        -------
        NDArray[float]
            Coupling field (pulse approximation).
        """
        # Detect oscillators near firing (θ ≈ 2π)
        firing = np.abs(phases - 2 * np.pi) < 0.1

        # Approximate pulse influence
        N = len(phases)
        n_firing = np.sum(firing)

        if n_firing > 0:
            # All non-firing oscillators receive pulses
            field = np.where(
                ~firing,
                self.strength * n_firing / N,
                0
            )
        else:
            field = np.zeros(N)

        return field


class CustomCoupling(Coupling):
    """
    Custom coupling defined by user-provided function.

    Parameters
    ----------
    coupling_func : callable
        Function (phases) -> coupling_field.
    """

    def __init__(self, coupling_func: callable):
        """Initialize custom coupling."""
        self.coupling_func = coupling_func

    def compute_field(self, phases: NDArray) -> NDArray:
        """Apply custom coupling function."""
        return self.coupling_func(phases)


class HigherHarmonicCoupling(Coupling):
    """
    Coupling with higher harmonics.

    Generalized coupling:
        H_j = Σ_n (K_n/N) Σ_k sin(n(θ_k - θ_j) + α_n)

    Parameters
    ----------
    harmonics : dict
        Dictionary {n: (K_n, α_n)} for each harmonic.
    """

    def __init__(self, harmonics: dict):
        """Initialize higher harmonic coupling."""
        self.harmonics = harmonics

    def compute_field(self, phases: NDArray) -> NDArray:
        """
        Compute coupling with higher harmonics.

        Parameters
        ----------
        phases : NDArray[float]
            Current phases of all oscillators.

        Returns
        -------
        NDArray[float]
            Total coupling field from all harmonics.
        """
        N = len(phases)
        total_field = np.zeros(N)

        for n, (K_n, alpha_n) in self.harmonics.items():
            # Phase differences for harmonic n
            phase_diff = phases[np.newaxis, :] - phases[:, np.newaxis]
            harmonic_coupling = (K_n / N) * np.sum(
                np.sin(n * phase_diff + alpha_n),
                axis=1
            )
            total_field += harmonic_coupling

        return total_field