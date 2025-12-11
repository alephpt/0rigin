"""
Simplified demonstration of fermion mass generation from synchronization.

Qualitative model showing how effective mass emerges from order parameter.
"""

from typing import Tuple, Dict, Optional
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt

from ..fields.grid import SpatialGrid
from ..fields.scalar_field import ScalarField


class FermionMassDemo:
    """
    Demonstrate effective mass generation via Yukawa coupling.

    Simplified model:
        m_eff(x,t) = Δ · R(x,t)

    where:
        - m_eff is the effective fermion mass
        - Δ is the Yukawa coupling strength
        - R(x,t) is the local order parameter (like Higgs VEV)

    This shows how mass emerges from spontaneous symmetry breaking
    in the Kuramoto field theory.
    """

    def __init__(
        self,
        grid: SpatialGrid,
        yukawa_coupling: float = 1.0,
        bare_mass: float = 0.0
    ):
        """
        Initialize fermion mass demonstration.

        Parameters
        ----------
        grid : SpatialGrid
            Spatial grid for fields.
        yukawa_coupling : float
            Coupling Δ between fermions and order parameter.
        bare_mass : float
            Bare fermion mass (usually 0 for massless fermions).
        """
        self.grid = grid
        self.Delta = yukawa_coupling
        self.m0 = bare_mass

        # Order parameter field R(x,t)
        self.R_field = ScalarField(grid, name="R_order_parameter")

        # Effective mass field
        self.m_eff_field = ScalarField(grid, name="m_effective")

    def update_from_oscillators(
        self,
        phases: NDArray,
        positions: NDArray,
        kernel_width: float = 0.1
    ):
        """
        Update order parameter and compute effective mass.

        Parameters
        ----------
        phases : NDArray
            Oscillator phases, shape (N,).
        positions : NDArray
            Oscillator positions, shape (N, 2).
        kernel_width : float
            Kernel width for density estimation.
        """
        # Update R field from oscillators
        self.R_field.update_from_oscillators(
            phases, positions,
            kernel_width=kernel_width,
            coupling_strength=1.0
        )

        # Compute effective mass: m_eff = m0 + Δ·R
        self.m_eff_field.values = self.m0 + self.Delta * self.R_field.values

    def compute_mass_gap(self) -> float:
        """
        Compute mass gap: difference between max and min effective mass.

        Returns
        -------
        float
            Mass gap Δm = m_max - m_min.
        """
        return np.max(self.m_eff_field.values) - np.min(self.m_eff_field.values)

    def compute_average_mass(self) -> float:
        """
        Compute spatial average of effective mass.

        Returns
        -------
        float
            <m_eff>_space
        """
        return np.mean(self.m_eff_field.values)

    def plot_mass_vs_order_parameter(
        self,
        ax: Optional[plt.Axes] = None
    ) -> plt.Axes:
        """
        Plot scatter of m_eff vs R to show linear relationship.

        Parameters
        ----------
        ax : plt.Axes, optional
            Axes to plot on.

        Returns
        -------
        plt.Axes
            Axes with plot.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        # Scatter plot of all grid points
        R_flat = self.R_field.values.flatten()
        m_flat = self.m_eff_field.values.flatten()

        ax.scatter(R_flat, m_flat, alpha=0.3, s=10)

        # Theoretical line
        R_theory = np.linspace(0, np.max(R_flat), 100)
        m_theory = self.m0 + self.Delta * R_theory
        ax.plot(R_theory, m_theory, 'r--', linewidth=2,
                label=f'$m_{{eff}} = {self.m0:.2f} + {self.Delta:.2f}R$')

        ax.set_xlabel('Order Parameter R(x)')
        ax.set_ylabel('Effective Mass $m_{eff}(x)$')
        ax.set_title('Mass Generation from Synchronization')
        ax.legend()
        ax.grid(True, alpha=0.3)

        return ax

    def plot_fields_side_by_side(
        self,
        fig: Optional[plt.Figure] = None
    ) -> plt.Figure:
        """
        Plot R field and m_eff field side by side.

        Parameters
        ----------
        fig : plt.Figure, optional
            Figure to plot on.

        Returns
        -------
        plt.Figure
            Figure with plots.
        """
        if fig is None:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        else:
            axes = fig.subplots(1, 2)

        # Plot R field
        im1 = axes[0].imshow(
            self.R_field.values.T,
            origin='lower',
            extent=[0, self.grid.Lx, 0, self.grid.Ly],
            aspect='auto',
            cmap='viridis'
        )
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('y')
        axes[0].set_title('Order Parameter R(x,y)')
        plt.colorbar(im1, ax=axes[0])

        # Plot m_eff field
        im2 = axes[1].imshow(
            self.m_eff_field.values.T,
            origin='lower',
            extent=[0, self.grid.Lx, 0, self.grid.Ly],
            aspect='auto',
            cmap='plasma'
        )
        axes[1].set_xlabel('x')
        axes[1].set_ylabel('y')
        axes[1].set_title('Effective Mass $m_{eff}(x,y)$')
        plt.colorbar(im2, ax=axes[1])

        plt.tight_layout()

        return fig

    def compute_symmetry_breaking_order(self) -> float:
        """
        Quantify degree of symmetry breaking.

        Measures variance in effective mass (broken symmetry → variation).

        Returns
        -------
        float
            Variance of effective mass field.
        """
        return np.var(self.m_eff_field.values)

    def demonstrate_phase_transition(
        self,
        R_values: NDArray,
        plot: bool = True
    ) -> Dict:
        """
        Demonstrate mass generation as function of order parameter.

        Simulates how mass varies as system goes through sync transition.

        Parameters
        ----------
        R_values : NDArray
            Array of R values to test (e.g., from 0 to 1).
        plot : bool
            Whether to create plot.

        Returns
        -------
        dict
            Results with R, m_eff, mass_gap.
        """
        m_eff_avg = np.zeros_like(R_values)

        for i, R in enumerate(R_values):
            # Set uniform R field
            self.R_field.values = np.full((self.grid.Nx, self.grid.Ny), R)

            # Compute mass
            self.m_eff_field.values = self.m0 + self.Delta * self.R_field.values

            m_eff_avg[i] = self.compute_average_mass()

        results = {
            'R': R_values,
            'm_eff': m_eff_avg,
            'mass_gap': self.Delta * R_values  # For uniform R
        }

        if plot:
            fig, ax = plt.subplots(figsize=(8, 6))

            ax.plot(R_values, m_eff_avg, 'b-', linewidth=2,
                   label='Average Mass')
            ax.axhline(self.m0, color='gray', linestyle='--',
                      label=f'Bare Mass = {self.m0:.2f}')

            ax.set_xlabel('Order Parameter R', fontsize=12)
            ax.set_ylabel('Effective Mass $\\langle m_{eff} \\rangle$', fontsize=12)
            ax.set_title('Phase Transition: Mass Generation', fontsize=14)
            ax.legend()
            ax.grid(True, alpha=0.3)

            # Annotate critical regions
            ax.axvspan(0, 0.3, alpha=0.1, color='red',
                      label='Disordered (massless)')
            ax.axvspan(0.7, 1.0, alpha=0.1, color='blue',
                      label='Ordered (massive)')

            plt.tight_layout()
            results['figure'] = fig

        return results

    def estimate_coherence_length(self) -> float:
        """
        Estimate coherence length from R field variations.

        Coherence length ~ 1/|∇R|

        Returns
        -------
        float
            Typical coherence length.
        """
        grad_x, grad_y = self.grid.gradient(self.R_field.values)
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)

        # Average inverse gradient (avoid division by zero)
        grad_mag_safe = np.where(grad_mag > 1e-6, grad_mag, 1e-6)
        coherence_length = np.mean(1.0 / grad_mag_safe)

        return coherence_length

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"FermionMassDemo(Δ={self.Delta:.2f}, "
            f"m0={self.m0:.2f}, "
            f"<m_eff>={self.compute_average_mass():.3f})"
        )
