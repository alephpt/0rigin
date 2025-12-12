"""
Time series visualization tools.

Provides functions for plotting order parameter evolution and
phase dynamics over time.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray

from .utils import _check_matplotlib

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def plot_order_parameter(
    t: NDArray,
    R: NDArray,
    Psi: Optional[NDArray] = None,
    R_theory: Optional[float] = None,
    ax: Optional['Axes'] = None,
    title: str = "Order Parameter Evolution"
) -> 'Axes':
    """
    Plot order parameter time series.

    Parameters
    ----------
    t : NDArray[float]
        Time array.
    R : NDArray[float]
        Order parameter amplitude R(t).
    Psi : NDArray[float], optional
        Order parameter phase Ψ(t). Not plotted if None.
    R_theory : float, optional
        Theoretical steady-state R value for comparison.
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
    >>> R = 0.5 * (1 - np.exp(-t/10))
    >>> ax = plot_order_parameter(t, R, R_theory=0.5)
    >>> plt.show()
    """
    _check_matplotlib()

    t = np.asarray(t)
    R = np.asarray(R)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot R(t)
    ax.plot(t, R, 'b-', linewidth=2, label='R(t)')

    # Plot theoretical value if provided
    if R_theory is not None:
        ax.axhline(R_theory, color='red', linestyle='--',
                   linewidth=2, label=f'Theory: R = {R_theory:.3f}')

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Order Parameter R', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_ylim([0, 1.05])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    return ax


def plot_phase_evolution(
    t: NDArray,
    phases: NDArray,
    frequencies: Optional[NDArray] = None,
    ax: Optional['Axes'] = None,
    title: str = "Phase Evolution",
    max_oscillators: int = 50,
    cmap: str = "viridis"
) -> 'Axes':
    """
    Plot phase evolution over time.

    Parameters
    ----------
    t : NDArray[float]
        Time array.
    phases : NDArray[float]
        Phase trajectories, shape (n_times, n_oscillators).
    frequencies : NDArray[float], optional
        Natural frequencies for color coding.
    ax : Axes, optional
        Matplotlib axes. Creates new if None.
    title : str, optional
        Plot title.
    max_oscillators : int, optional
        Maximum number of oscillators to plot (default: 50).
    cmap : str, optional
        Colormap name (default: 'viridis').

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 50, 1000)
    >>> phases = np.random.uniform(0, 2*np.pi, (1000, 20))
    >>> ax = plot_phase_evolution(t, phases)
    >>> plt.show()
    """
    _check_matplotlib()

    t = np.asarray(t)
    phases = np.asarray(phases)
    N = phases.shape[1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    # Subsample if too many oscillators
    if N > max_oscillators:
        indices = np.linspace(0, N-1, max_oscillators, dtype=int)
        phases_plot = phases[:, indices]
        if frequencies is not None:
            frequencies_plot = frequencies[indices]
    else:
        phases_plot = phases
        frequencies_plot = frequencies
        indices = np.arange(N)

    # Unwrap phases for continuous plotting
    phases_unwrapped = np.unwrap(phases_plot, axis=0)

    # Plot each oscillator
    if frequencies_plot is not None:
        colors = plt.cm.get_cmap(cmap)(
            (frequencies_plot - frequencies_plot.min()) /
            (frequencies_plot.max() - frequencies_plot.min() + 1e-10)
        )
        for i, idx in enumerate(indices):
            ax.plot(t, phases_unwrapped[:, i], color=colors[i],
                    alpha=0.6, linewidth=1)
    else:
        ax.plot(t, phases_unwrapped, alpha=0.6, linewidth=1)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Phase (unwrapped)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)

    return ax


def plot_synchronization_transition(
    t: NDArray,
    R: NDArray,
    K: float,
    Kc: Optional[float] = None,
    ax: Optional['Axes'] = None
) -> 'Axes':
    """
    Plot order parameter showing synchronization transition.

    Parameters
    ----------
    t : NDArray[float]
        Time array.
    R : NDArray[float]
        Order parameter amplitude R(t).
    K : float
        Coupling strength used in simulation.
    Kc : float, optional
        Critical coupling for reference.
    ax : Axes, optional
        Matplotlib axes. Creates new if None.

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 50, 1000)
    >>> R = 0.7 * (1 - np.exp(-t/5))
    >>> ax = plot_synchronization_transition(t, R, K=3.0, Kc=2.0)
    >>> plt.show()
    """
    _check_matplotlib()

    t = np.asarray(t)
    R = np.asarray(R)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Determine regime
    if Kc is not None:
        if K < Kc:
            regime = "Subcritical"
            color = 'blue'
        elif abs(K - Kc) / Kc < 0.1:
            regime = "Critical"
            color = 'orange'
        else:
            regime = "Supercritical"
            color = 'green'
    else:
        regime = ""
        color = 'blue'

    # Plot R(t)
    ax.plot(t, R, color=color, linewidth=2, label=f'K = {K:.2f}')

    # Add steady-state reference line
    if len(R) > 100:
        R_steady = np.mean(R[-100:])
        ax.axhline(R_steady, color='red', linestyle='--',
                   alpha=0.5, label=f'<R> = {R_steady:.3f}')

    title = "Synchronization Transition"
    if regime:
        title += f" ({regime} Regime)"
    if Kc is not None:
        title += f"\nKc = {Kc:.2f}"

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Order Parameter R', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_ylim([0, 1.05])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    return ax


def plot_complex_order_parameter(
    t: NDArray,
    R: NDArray,
    Psi: NDArray,
    ax: Optional['Axes'] = None
) -> 'Axes':
    """
    Plot complex order parameter trajectory.

    Shows Z(t) = R(t) e^(iΨ(t)) in complex plane.

    Parameters
    ----------
    t : NDArray[float]
        Time array.
    R : NDArray[float]
        Order parameter amplitude R(t).
    Psi : NDArray[float]
        Order parameter phase Ψ(t).
    ax : Axes, optional
        Matplotlib axes. Creates new if None.

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 50, 1000)
    >>> R = 0.5 * (1 - np.exp(-t/10))
    >>> Psi = 0.1 * t
    >>> ax = plot_complex_order_parameter(t, R, Psi)
    >>> plt.show()
    """
    _check_matplotlib()

    t = np.asarray(t)
    R = np.asarray(R)
    Psi = np.asarray(Psi)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # Convert to Cartesian
    x = R * np.cos(Psi)
    y = R * np.sin(Psi)

    # Plot unit circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k--', alpha=0.3, linewidth=1)

    # Plot trajectory with color gradient
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    from matplotlib.collections import LineCollection
    lc = LineCollection(segments, cmap='viridis', linewidth=2)
    lc.set_array(t)
    line = ax.add_collection(lc)

    # Mark start and end
    ax.scatter(x[0], y[0], c='green', s=100, marker='o',
               label='Start', zorder=5)
    ax.scatter(x[-1], y[-1], c='red', s=100, marker='s',
               label='End', zorder=5)

    ax.set_xlim([-1.1, 1.1])
    ax.set_ylim([-1.1, 1.1])
    ax.set_aspect('equal')
    ax.set_xlabel('Re(Z) = R cos(Ψ)', fontsize=12)
    ax.set_ylabel('Im(Z) = R sin(Ψ)', fontsize=12)
    ax.set_title('Complex Order Parameter Trajectory', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.colorbar(line, ax=ax, label='Time')

    return ax
