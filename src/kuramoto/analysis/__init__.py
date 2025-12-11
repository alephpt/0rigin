"""
Analysis tools for Kuramoto model simulations.

This module provides tools for analyzing synchronization phenomena,
computing order parameters, and extracting metrics from simulations.
"""

from .order_parameter import OrderParameter
from .metrics import (
    phase_coherence,
    phase_variance,
    local_order_parameter,
    frequency_entrainment,
    metastability,
    chimera_index,
    kuramoto_shinomoto_order_parameter,
    phase_frustration
)

__all__ = [
    "OrderParameter",
    "phase_coherence",
    "phase_variance",
    "local_order_parameter",
    "frequency_entrainment",
    "metastability",
    "chimera_index",
    "kuramoto_shinomoto_order_parameter",
    "phase_frustration",
]
