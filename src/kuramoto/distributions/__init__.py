"""
Frequency distributions for Kuramoto model.
"""

from .base import FrequencyDistribution
from .lorentzian import LorentzianDistribution
from .gaussian import GaussianDistribution
from .uniform import UniformDistribution

__all__ = [
    "FrequencyDistribution",
    "LorentzianDistribution",
    "GaussianDistribution",
    "UniformDistribution",
]
