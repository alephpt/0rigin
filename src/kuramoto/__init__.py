"""
Kuramoto Model Simulation Framework

A Python library for simulating and analyzing the Kuramoto model of coupled oscillators,
with support for various frequency distributions, coupling schemes, and analysis methods.

Basic Usage:
    >>> import kuramoto as km
    >>> model = km.KuramotoModel(N=100, coupling=2.0, frequencies='lorentzian')
    >>> solution = model.evolve(t_span=(0, 50))
    >>> R, Psi = model.compute_order_parameter()

Advanced Usage:
    >>> # Custom frequency distribution
    >>> freq_dist = km.distributions.LorentzianDistribution(center=0, width=1)
    >>> coupling = km.coupling.SinusoidalCoupling(K=3.0)
    >>> model = km.KuramotoModel(N=500, coupling=coupling, frequencies=freq_dist)
    >>>
    >>> # Analyze synchronization transition
    >>> analyzer = km.analysis.BifurcationAnalysis(model)
    >>> K_critical = analyzer.find_critical_coupling()
"""

__version__ = "0.1.0"
__author__ = "0rigin Research Team"

# Core classes
from .core.model import KuramotoModel
from .core.coupling import SinusoidalCoupling, NetworkCoupling

# Frequency distributions
from .distributions import (
    LorentzianDistribution,
    GaussianDistribution,
    UniformDistribution
)

# Solvers
from .solvers import RK4Solver, RK45Solver, EulerSolver

# Analysis tools (sub-modules)
from . import analysis
from . import visualization

__all__ = [
    # Core
    "KuramotoModel",
    "SinusoidalCoupling",
    "NetworkCoupling",
    # Distributions
    "LorentzianDistribution",
    "GaussianDistribution",
    "UniformDistribution",
    # Solvers
    "RK4Solver",
    "RK45Solver",
    "EulerSolver",
    # Sub-modules
    "analysis",
    "visualization",
]