"""
Bifurcation diagram visualization.

Provides functions for plotting bifurcation diagrams showing
the synchronization transition as coupling strength varies.
"""

from typing import Optional, Callable
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


def plot_bifurcation_diagram(
    K_values: NDArray,
    R_values: NDArray,
    R_theory: Optional[NDArray] = None,
    Kc: Optional[float] = None,
    ax: Optional['Axes'] = None,
    title: str = "Bifurcation Diagram"
) -> 'Axes':
    """
    Plot bifurcation diagram showing R vs K.

    Parameters
    ----------
    K_values : NDArray[float]
        Coupling strength values.
    R_values : NDArray[float]
        Steady-state order parameter values.
    R_theory : NDArray[float], optional
        Theoretical R values for comparison.
    Kc : float, optional
        Critical coupling to mark on plot.
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
    >>> K_values = np.linspace(0, 5, 50)
    >>> Kc = 2.0
    >>> R_values = np.sqrt(np.maximum(0, 1 - (Kc/K_values)**2))
    >>> ax = plot_bifurcation_diagram(K_values, R_values, Kc=Kc)
    >>> plt.show()
    """
    _check_matplotlib()

    K_values = np.asarray(K_values)
    R_values = np.asarray(R_values)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot simulation results
    ax.plot(K_values, R_values, 'bo', markersize=6, label='Simulation')

    # Plot theory if provided
    if R_theory is not None:
        R_theory = np.asarray(R_theory)
        ax.plot(K_values, R_theory, 'r-', linewidth=2, label='Theory')

    # Mark critical coupling
    if Kc is not None:
        ax.axvline(Kc, color='gray', linestyle='--', linewidth=2,
                   alpha=0.7, label=f'Kc = {Kc:.3f}')

        # Add regime labels
        K_max = K_values.max()
        ax.text(Kc * 0.5, 0.5, 'Incoherent\nPhase',
                ha='center', va='center', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        ax.text((Kc + K_max) / 2, 0.5, 'Partially\nSynchronized',
                ha='center', va='center', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

    ax.set_xlabel('Coupling Strength K', fontsize=12)
    ax.set_ylabel('Order Parameter R', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim([0, K_values.max()])
    ax.set_ylim([0, 1.05])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    return ax


def plot_hysteresis_diagram(
    K_forward: NDArray,
    R_forward: NDArray,
    K_backward: NDArray,
    R_backward: NDArray,
    ax: Optional['Axes'] = None,
    title: str = "Hysteresis Diagram"
) -> 'Axes':
    """
    Plot hysteresis diagram for systems with bistability.

    Parameters
    ----------
    K_forward : NDArray[float]
        Coupling values for forward sweep.
    R_forward : NDArray[float]
        Order parameter for forward sweep.
    K_backward : NDArray[float]
        Coupling values for backward sweep.
    R_backward : NDArray[float]
        Order parameter for backward sweep.
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
    >>> K_forward = np.linspace(0, 5, 50)
    >>> R_forward = np.sqrt(np.maximum(0, 1 - (2/K_forward)**2))
    >>> K_backward = K_forward[::-1]
    >>> R_backward = R_forward[::-1]
    >>> ax = plot_hysteresis_diagram(K_forward, R_forward, K_backward, R_backward)
    >>> plt.show()
    """
    _check_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot forward and backward sweeps
    ax.plot(K_forward, R_forward, 'b-o', markersize=4,
            label='Forward (increasing K)', alpha=0.7)
    ax.plot(K_backward, R_backward, 'r-s', markersize=4,
            label='Backward (decreasing K)', alpha=0.7)

    # Add arrows to indicate direction
    mid_idx = len(K_forward) // 2
    ax.annotate('', xy=(K_forward[mid_idx+5], R_forward[mid_idx+5]),
                xytext=(K_forward[mid_idx], R_forward[mid_idx]),
                arrowprops=dict(arrowstyle='->', color='blue', lw=2))
    ax.annotate('', xy=(K_backward[mid_idx+5], R_backward[mid_idx+5]),
                xytext=(K_backward[mid_idx], R_backward[mid_idx]),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.set_xlabel('Coupling Strength K', fontsize=12)
    ax.set_ylabel('Order Parameter R', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_ylim([0, 1.05])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    return ax


def plot_phase_diagram(
    parameter1: NDArray,
    parameter2: NDArray,
    R_grid: NDArray,
    param1_name: str = "Parameter 1",
    param2_name: str = "Parameter 2",
    ax: Optional['Axes'] = None,
    title: str = "Phase Diagram",
    levels: Optional[NDArray] = None
) -> 'Axes':
    """
    Plot 2D phase diagram in parameter space.

    Parameters
    ----------
    parameter1 : NDArray[float]
        First parameter values (1D array).
    parameter2 : NDArray[float]
        Second parameter values (1D array).
    R_grid : NDArray[float]
        Order parameter grid, shape (len(parameter2), len(parameter1)).
    param1_name : str, optional
        Name of first parameter.
    param2_name : str, optional
        Name of second parameter.
    ax : Axes, optional
        Matplotlib axes. Creates new if None.
    title : str, optional
        Plot title.
    levels : NDArray[float], optional
        Contour levels. Auto-generated if None.

    Returns
    -------
    Axes
        Matplotlib axes object.

    Examples
    --------
    >>> import numpy as np
    >>> K = np.linspace(0, 5, 50)
    >>> gamma = np.linspace(0.5, 2, 40)
    >>> K_grid, gamma_grid = np.meshgrid(K, gamma)
    >>> R_grid = np.sqrt(np.maximum(0, 1 - (2*gamma_grid/K_grid)**2))
    >>> ax = plot_phase_diagram(K, gamma, R_grid, "K", "γ")
    >>> plt.show()
    """
    _check_matplotlib()

    parameter1 = np.asarray(parameter1)
    parameter2 = np.asarray(parameter2)
    R_grid = np.asarray(R_grid)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Create meshgrid
    P1, P2 = np.meshgrid(parameter1, parameter2)

    # Plot filled contours
    if levels is None:
        levels = np.linspace(0, 1, 11)

    contourf = ax.contourf(P1, P2, R_grid, levels=levels,
                           cmap='viridis', alpha=0.8)
    contour = ax.contour(P1, P2, R_grid, levels=levels,
                         colors='black', alpha=0.3, linewidths=0.5)

    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax, label='Order Parameter R')
    ax.clabel(contour, inline=True, fontsize=8)

    # Add critical line (R = 0.5) if present
    try:
        critical_contour = ax.contour(P1, P2, R_grid, levels=[0.5],
                                     colors='red', linewidths=2)
        ax.clabel(critical_contour, inline=True, fmt='R=0.5', fontsize=10)
    except:
        pass

    ax.set_xlabel(param1_name, fontsize=12)
    ax.set_ylabel(param2_name, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)

    return ax


def plot_finite_size_scaling(
    N_values: NDArray,
    Kc_values: NDArray,
    Kc_infinite: Optional[float] = None,
    ax: Optional['Axes'] = None,
    title: str = "Finite-Size Scaling"
) -> 'Axes':
    """
    Plot finite-size scaling of critical coupling.

    Shows how critical coupling approaches thermodynamic limit.

    Parameters
    ----------
    N_values : NDArray[int]
        System sizes.
    Kc_values : NDArray[float]
        Critical coupling for each size.
    Kc_infinite : float, optional
        Thermodynamic limit value.
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
    >>> N_values = np.array([10, 20, 50, 100, 200, 500, 1000])
    >>> Kc_infinite = 2.0
    >>> Kc_values = Kc_infinite + 0.5 / np.sqrt(N_values)
    >>> ax = plot_finite_size_scaling(N_values, Kc_values, Kc_infinite)
    >>> plt.show()
    """
    _check_matplotlib()

    N_values = np.asarray(N_values)
    Kc_values = np.asarray(Kc_values)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot scaling
    ax.plot(N_values, Kc_values, 'bo-', markersize=8,
            linewidth=2, label='Kc(N)')

    # Plot thermodynamic limit
    if Kc_infinite is not None:
        ax.axhline(Kc_infinite, color='red', linestyle='--',
                   linewidth=2, label=f'Kc(∞) = {Kc_infinite:.3f}')

    ax.set_xlabel('System Size N', fontsize=12)
    ax.set_ylabel('Critical Coupling Kc', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xscale('log')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    return ax
