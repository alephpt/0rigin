"""
Dirac equation solver for MSFT.

Implements full Dirac spinor dynamics coupled to Kuramoto synchronization fields.
"""

from .gamma_matrices import (
    get_gamma_matrices_3plus1,
    get_gamma_matrices_1plus1,
    mass_operator_MSFT,
    chiral_projectors,
)

__all__ = [
    'get_gamma_matrices_3plus1',
    'get_gamma_matrices_1plus1',
    'mass_operator_MSFT',
    'chiral_projectors',
]
