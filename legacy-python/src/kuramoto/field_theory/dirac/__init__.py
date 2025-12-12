"""
Dirac equation solver for SMFT.

Implements full Dirac spinor dynamics coupled to Kuramoto synchronization fields.
"""

from .gamma_matrices import (
    get_gamma_matrices_3plus1,
    get_gamma_matrices_1plus1,
    mass_operator_smft,
    chiral_projectors,
)

__all__ = [
    'get_gamma_matrices_3plus1',
    'get_gamma_matrices_1plus1',
    'mass_operator_smft',
    'chiral_projectors',
]
