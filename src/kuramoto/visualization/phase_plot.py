"""
Phase visualization tools.

Provides functions for plotting phase distributions on the unit circle
and phase space representations.
"""

from typing import Optional, Tuple
import numpy as np
from numpy.typing import NDArray

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _check_matplotlib():
    """Check if matplotlib is available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "Matplotlib is required for visualization. "
            "Install it with: pip install matplotlib"
        )


def plot_phases(
    phases: NDArray,
    frequencies: Optional[NDArray] = None,
    ax: Optional['Axes'] = None,
    title: str = "Phase Distribution",
    cmap: str = "coolwarm"
) -> 'Axes':
    """
    Plot phase distribution as scatter on unit circle.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    frequencies : NDArray[float], optional
        Natural frequencies for color coding.
    ax : Axes, optional
        Matplotlib axes to plot on. Creates new if None.
    title : str, optional
        Plot title.
    cmap : str, optional
        Colormap name (default: 'coolwarm').

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> frequencies = np.random.normal(0, 1, 100)
    >>> ax = plot_phases(phases, frequencies)
    >>> plt.show()
    """
    _check_matplotlib()

    phases = np.asarray(phases)

    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)

    # Convert to Cartesian coordinates
    x = np.cos(phases)
    y = np.sin(phases)

    # Plot unit circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.3, linewidth=1)

    # Plot phases
    if frequencies is not None:
        scatter = ax.scatter(x, y, c=frequencies, cmap=cmap, s=50, alpha=0.7)
        plt.colorbar(scatter, ax=ax, label='Natural Frequency')
    else:
        ax.scatter(x, y, c='blue', s=50, alpha=0.7)

    # Compute and plot order parameter
    z = np.mean(np.exp(1j * phases))
    R, Psi = np.abs(z), np.angle(z)

    ax.arrow(0, 0, R*np.cos(Psi), R*np.sin(Psi),
             head_width=0.05, head_length=0.05, fc='red', ec='red',
             linewidth=2, label=f'R = {R:.3f}')

    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_aspect('equal')
    ax.set_xlabel('cos(θ)')
    ax.set_ylabel('sin(θ)')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_phase_circle(
    phases: NDArray,
    order_parameter: Optional[Tuple[float, float]] = None,
    ax: Optional['Axes'] = None,
    title: str = "Phase Circle",
    show_order_parameter: bool = True
) -> 'Axes':
    """
    Plot phases on polar circle representation.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    order_parameter : tuple of float, optional
        (R, Psi) order parameter. Computed if None.
    ax : Axes, optional
        Matplotlib axes (must be polar projection). Creates new if None.
    title : str, optional
        Plot title.
    show_order_parameter : bool, optional
        Whether to show order parameter arrow (default: True).

    Returns
    -------
    Axes
        Matplotlib axes object (polar projection).

    Examples
    --------
    >>> import numpy as np
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> ax = plot_phase_circle(phases)
    >>> plt.show()
    """
    _check_matplotlib()

    phases = np.asarray(phases)

    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='polar')

    # Plot oscillators on unit circle
    radii = np.ones_like(phases)
    ax.scatter(phases, radii, c='blue', s=50, alpha=0.6)

    if show_order_parameter:
        # Compute order parameter if not provided
        if order_parameter is None:
            z = np.mean(np.exp(1j * phases))
            R, Psi = np.abs(z), np.angle(z)
        else:
            R, Psi = order_parameter

        # Plot order parameter vector
        ax.arrow(0, 0, Psi, R,
                 head_width=0.2, head_length=0.1, fc='red', ec='red',
                 linewidth=2, length_includes_head=True,
                 label=f'R = {R:.3f}')
        ax.legend(loc='upper right')

    ax.set_ylim([0, 1.2])
    ax.set_title(title, pad=20)

    return ax


def plot_phase_histogram(
    phases: NDArray,
    bins: int = 36,
    ax: Optional['Axes'] = None,
    title: str = "Phase Histogram"
) -> 'Axes':
    """
    Plot histogram of phase distribution.

    Parameters
    ----------
    phases : NDArray[float]
        Phase array, shape (n_oscillators,).
    bins : int, optional
        Number of histogram bins (default: 36 for 10° bins).
    ax : Axes, optional
        Matplotlib axes. Creates new if None.
    title : str, optional
        Plot title.

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> phases = np.random.uniform(0, 2*np.pi, 100)
    >>> ax = plot_phase_histogram(phases)
    >>> plt.show()
    """
    _check_matplotlib()

    phases = np.asarray(phases)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    # Create histogram
    counts, edges, patches = ax.hist(
        phases, bins=bins, range=(0, 2*np.pi),
        density=True, alpha=0.7, color='blue', edgecolor='black'
    )

    # Add uniform distribution reference
    ax.axhline(1/(2*np.pi), color='red', linestyle='--',
               linewidth=2, label='Uniform distribution')

    ax.set_xlabel('Phase (radians)')
    ax.set_ylabel('Probability Density')
    ax.set_title(title)
    ax.set_xlim([0, 2*np.pi])
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_phase_space_2d(
    phases1: NDArray,
    phases2: NDArray,
    ax: Optional['Axes'] = None,
    title: str = "2D Phase Space"
) -> 'Axes':
    """
    Plot 2D phase space for two oscillators.

    Useful for analyzing synchronization between pairs.

    Parameters
    ----------
    phases1 : NDArray[float]
        Phase trajectory of first oscillator.
    phases2 : NDArray[float]
        Phase trajectory of second oscillator.
    ax : Axes, optional
        Matplotlib axes. Creates new if None.
    title : str, optional
        Plot title.

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 50, 1000)
    >>> phases1 = np.mod(1.0 * t, 2*np.pi)
    >>> phases2 = np.mod(1.1 * t, 2*np.pi)
    >>> ax = plot_phase_space_2d(phases1, phases2)
    >>> plt.show()
    """
    _check_matplotlib()

    phases1 = np.asarray(phases1)
    phases2 = np.asarray(phases2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # Plot trajectory
    ax.plot(phases1, phases2, 'b-', alpha=0.5, linewidth=0.5)
    ax.scatter(phases1[0], phases2[0], c='green', s=100,
               marker='o', label='Start', zorder=5)
    ax.scatter(phases1[-1], phases2[-1], c='red', s=100,
               marker='s', label='End', zorder=5)

    # Add diagonal (synchronization line)
    ax.plot([0, 2*np.pi], [0, 2*np.pi], 'k--', alpha=0.3, label='θ₁ = θ₂')

    ax.set_xlabel('Phase 1 (radians)')
    ax.set_ylabel('Phase 2 (radians)')
    ax.set_title(title)
    ax.set_xlim([0, 2*np.pi])
    ax.set_ylim([0, 2*np.pi])
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax
