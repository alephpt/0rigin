"""
Additional synchronization metrics for Kuramoto model analysis.

Provides various measures of synchronization beyond the standard
order parameter.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray


def phase_coherence(phases: NDArray) -> float:
    """
    Compute phase coherence.

    Phase coherence is the magnitude of the mean complex phase:
        ρ = |⟨e^(iθ)⟩|

    This is equivalent to the order parameter amplitude R.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,) or (n_times, n_oscillators).

    Returns
    -------
    float or NDArray[float]
        Phase coherence value(s).

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> rho = phase_coherence(phases)
    """
    phases = np.asarray(phases)

    if phases.ndim == 1:
        # Single snapshot
        z = np.mean(np.exp(1j * phases))
        return np.abs(z)
    else:
        # Time series
        z = np.mean(np.exp(1j * phases), axis=1)
        return np.abs(z)


def phase_variance(phases: NDArray, circular: bool = True) -> float:
    """
    Compute phase variance.

    For circular=True, uses circular variance:
        V = 1 - |⟨e^(iθ)⟩|

    For circular=False, uses standard variance of unwrapped phases.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    circular : bool, optional
        Use circular variance if True (default: True).

    Returns
    -------
    float
        Phase variance.

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> var = phase_variance(phases, circular=True)
    """
    phases = np.asarray(phases)

    if circular:
        # Circular variance: V = 1 - R
        R = phase_coherence(phases)
        if isinstance(R, np.ndarray):
            return 1 - R
        else:
            return 1 - R
    else:
        # Standard variance
        return np.var(phases)


def local_order_parameter(
    phases: NDArray,
    neighbors: NDArray,
    weights: Optional[NDArray] = None
) -> NDArray:
    """
    Compute local order parameter for each oscillator.

    The local order parameter measures synchronization with neighbors:
        Z_i = Σ_j w_ij e^(iθ_j) / Σ_j w_ij

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    neighbors : NDArray[int]
        Neighbor indices for each oscillator, shape (n_oscillators, n_neighbors).
        Use -1 to indicate no neighbor in a slot.
    weights : NDArray[float], optional
        Weights for neighbors, shape (n_oscillators, n_neighbors).
        Uses uniform weights if None.

    Returns
    -------
    NDArray[float]
        Local order parameter for each oscillator.

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 10)
    >>> # Ring topology: each oscillator connected to 2 neighbors
    >>> neighbors = np.array([[9, 1], [0, 2], [1, 3], [2, 4], [3, 5],
    ...                       [4, 6], [5, 7], [6, 8], [7, 9], [8, 0]])
    >>> local_R = local_order_parameter(phases, neighbors)
    """
    phases = np.asarray(phases)
    neighbors = np.asarray(neighbors)
    N = len(phases)

    if weights is None:
        weights = np.ones_like(neighbors, dtype=float)
    else:
        weights = np.asarray(weights)

    local_z = np.zeros(N, dtype=complex)

    for i in range(N):
        # Get valid neighbors (not -1)
        valid = neighbors[i] >= 0
        if not np.any(valid):
            continue

        neighbor_idx = neighbors[i][valid]
        neighbor_weights = weights[i][valid]

        # Compute weighted sum
        z_sum = np.sum(neighbor_weights * np.exp(1j * phases[neighbor_idx]))
        w_sum = np.sum(neighbor_weights)

        local_z[i] = z_sum / w_sum if w_sum > 0 else 0

    return np.abs(local_z)


def frequency_entrainment(
    frequencies: NDArray,
    threshold: float = 0.01
) -> float:
    """
    Measure frequency entrainment.

    Computes the fraction of oscillators with similar instantaneous
    frequencies (within threshold of the mean).

    Parameters
    ----------
    frequencies : NDArray[float]
        Instantaneous frequency array, shape (n_oscillators,).
    threshold : float, optional
        Frequency tolerance (default: 0.01).

    Returns
    -------
    float
        Fraction of entrained oscillators (0 to 1).

    Examples
    --------
    >>> # Compute frequencies from phase time series
    >>> phases = np.random.uniform(0, 2*np.pi, (100, 50))
    >>> freqs = np.diff(phases, axis=0).mean(axis=0)
    >>> entrainment = frequency_entrainment(freqs)
    """
    frequencies = np.asarray(frequencies)
    mean_freq = np.mean(frequencies)

    entrained = np.abs(frequencies - mean_freq) < threshold
    return np.mean(entrained)


def metastability(R_timeseries: NDArray) -> float:
    """
    Compute metastability index.

    Metastability measures fluctuations in synchronization:
        M = σ(R) / ⟨R⟩

    High metastability indicates intermittent synchronization.

    Parameters
    ----------
    R_timeseries : NDArray[float]
        Order parameter time series.

    Returns
    -------
    float
        Metastability index.

    Examples
    --------
    >>> R = np.random.uniform(0.3, 0.7, 1000)
    >>> M = metastability(R)
    """
    R_timeseries = np.asarray(R_timeseries)

    mean_R = np.mean(R_timeseries)
    std_R = np.std(R_timeseries)

    if mean_R == 0:
        return 0.0

    return std_R / mean_R


def chimera_index(phases: NDArray, window_size: int = 10) -> float:
    """
    Compute chimera index.

    Measures coexistence of synchronized and desynchronized regions.
    Uses local order parameter variance as indicator.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    window_size : int, optional
        Window size for local order parameter (default: 10).

    Returns
    -------
    float
        Chimera index (0 = uniform, 1 = strong chimera).

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> chi = chimera_index(phases)
    """
    phases = np.asarray(phases)
    N = len(phases)

    # Compute local order parameter using sliding window
    local_R = np.zeros(N)

    for i in range(N):
        # Get neighbors in window
        start = max(0, i - window_size // 2)
        end = min(N, i + window_size // 2 + 1)
        local_phases = phases[start:end]

        # Compute local coherence
        z = np.mean(np.exp(1j * local_phases))
        local_R[i] = np.abs(z)

    # Chimera index: variance of local order parameters
    return np.var(local_R)


def kuramoto_shinomoto_order_parameter(phases: NDArray, kappa: float = 1.0) -> float:
    """
    Compute Kuramoto-Shinomoto order parameter.

    Generalized order parameter with tunable sensitivity:
        Z_κ = ⟨e^(iκθ)⟩

    For κ=1, reduces to standard order parameter.
    For κ>1, more sensitive to higher harmonics.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    kappa : float, optional
        Sensitivity parameter (default: 1.0).

    Returns
    -------
    float
        Generalized order parameter magnitude.

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> Z1 = kuramoto_shinomoto_order_parameter(phases, kappa=1)
    >>> Z2 = kuramoto_shinomoto_order_parameter(phases, kappa=2)
    """
    phases = np.asarray(phases)
    z = np.mean(np.exp(1j * kappa * phases))
    return np.abs(z)


def phase_frustration(
    phases: NDArray,
    natural_frequencies: NDArray,
    coupling_strength: float
) -> float:
    """
    Compute phase frustration.

    Measures the conflict between natural dynamics and coupling:
        F = ⟨|ω_i - K R sin(Ψ - θ_i)|⟩

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    natural_frequencies : NDArray[float]
        Natural frequencies ω_i.
    coupling_strength : float
        Coupling strength K.

    Returns
    -------
    float
        Average phase frustration.

    Examples
    --------
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> freqs = np.random.normal(0, 1, 100)
    >>> frustration = phase_frustration(phases, freqs, coupling_strength=2.0)
    """
    phases = np.asarray(phases)
    natural_frequencies = np.asarray(natural_frequencies)

    # Compute order parameter
    z = np.mean(np.exp(1j * phases))
    R = np.abs(z)
    Psi = np.angle(z)

    # Coupling field for each oscillator
    coupling_field = coupling_strength * R * np.sin(Psi - phases)

    # Frustration: difference from natural frequency
    frustration = np.abs(natural_frequencies - coupling_field)

    return np.mean(frustration)
