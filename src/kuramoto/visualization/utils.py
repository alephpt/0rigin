"""
Visualization utilities.

Common utility functions shared across visualization modules.
"""

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _check_matplotlib():
    """
    Check if matplotlib is available.

    Raises
    ------
    ImportError
        If matplotlib is not installed.
    """
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "Matplotlib is required for visualization. "
            "Install it with: pip install matplotlib"
        )