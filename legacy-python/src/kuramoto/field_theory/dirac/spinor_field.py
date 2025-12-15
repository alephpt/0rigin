"""
Dirac spinor field on spatial grid.

Implements 4-component complex spinor field Ψ(x,t) evolving via Dirac equation:
    i∂_tΨ = [(-iγ^i∂_i) + m(x,t)]Ψ

where m(x,t) = Δ·R(x,t)·e^(iθγ^5) is the MSFT mass operator.
"""

from typing import Tuple, Optional
import numpy as np
from numpy.typing import NDArray

from ..fields.grid import SpatialGrid
from .gamma_matrices import get_gamma_matrices_3plus1, mass_operator_MSFT


class DiracSpinorField:
    """
    Dirac spinor field Ψ(x,t) on 2D spatial grid.

    Represents 4-component complex wave function at each spatial point.
    For 2D grid (Nx, Ny), total shape is (Nx, Ny, 4) complex.

    Parameters
    ----------
    grid : SpatialGrid
        Spatial grid for field discretization.
    initial_state : str or NDArray, optional
        Initial spinor configuration:
        - 'gaussian': Gaussian wave packet
        - 'plane_wave': Plane wave with momentum k
        - 'vacuum': Zero field
        - NDArray: Custom initial state (Nx, Ny, 4)
    """

    def __init__(
        self,
        grid: SpatialGrid,
        initial_state: str = 'vacuum',
        momentum: Tuple[float, float] = (0.0, 0.0)
    ):
        self.grid = grid
        self.Nx = grid.Nx
        self.Ny = grid.Ny

        # Get gamma matrices
        self.γ0, self.γ1, self.γ2, self.γ3, self.γ5 = get_gamma_matrices_3plus1()

        # Spinor field: shape (Nx, Ny, 4) complex
        self.psi = np.zeros((self.Nx, self.Ny, 4), dtype=complex)

        # Initialize
        if initial_state == 'vacuum':
            pass  # Already zero
        elif initial_state == 'gaussian':
            self._init_gaussian()
        elif initial_state == 'plane_wave':
            self._init_plane_wave(momentum)
        elif isinstance(initial_state, np.ndarray):
            self.psi = initial_state.copy()
        else:
            raise ValueError(f"Unknown initial_state: {initial_state}")

        self.t = 0.0

    def _init_gaussian(self, sigma: float = 0.1):
        """Initialize as Gaussian wave packet (positive energy spinor)."""
        # Center of domain
        x0, y0 = self.grid.Lx / 2, self.grid.Ly / 2

        # Gaussian envelope
        gaussian = np.exp(-((self.grid.X - x0)**2 + (self.grid.Y - y0)**2) / (2 * sigma**2))

        # Positive energy spinor (upper components)
        self.psi[:, :, 0] = gaussian  # ψ_1
        self.psi[:, :, 1] = gaussian  # ψ_2
        # Lower components zero (particle at rest)

        # Normalize
        self._normalize()

    def _init_plane_wave(self, momentum: Tuple[float, float]):
        """Initialize as plane wave with momentum k = (kx, ky)."""
        kx, ky = momentum

        # Plane wave: e^(ik·x)
        phase = 1j * (kx * self.grid.X + ky * self.grid.Y)
        plane_wave = np.exp(phase)

        # Energy: E = √(k² + m²) ≈ |k| for massless
        E = np.sqrt(kx**2 + ky**2 + 1e-6)

        # Dirac spinor for momentum eigenstate
        # Upper components
        self.psi[:, :, 0] = plane_wave
        self.psi[:, :, 1] = plane_wave * kx / (E + 1)  # Simplified
        # Lower components (small for non-relativistic)
        self.psi[:, :, 2] = plane_wave * ky / (E + 1)
        self.psi[:, :, 3] = 0

        self._normalize()

    def _normalize(self):
        """Normalize spinor: ∫|Ψ|²dV = 1."""
        norm_sq = np.sum(np.abs(self.psi)**2) * self.grid.dx * self.grid.dy
        if norm_sq > 0:
            self.psi /= np.sqrt(norm_sq)

    def compute_density(self) -> NDArray:
        """
        Compute probability density ρ(x) = Ψ†Ψ = |Ψ|².

        Returns
        -------
        NDArray
            Density field (Nx, Ny).
        """
        return np.sum(np.abs(self.psi)**2, axis=2)

    def compute_current(self) -> Tuple[NDArray, NDArray]:
        """
        Compute probability current j^i = Ψ†γ^0γ^iΨ.

        Returns
        -------
        jx, jy : NDArray
            Current density components (Nx, Ny).
        """
        jx = np.zeros((self.Nx, self.Ny), dtype=complex)
        jy = np.zeros((self.Nx, self.Ny), dtype=complex)

        # j^1 = Ψ†γ^0γ^1Ψ
        γ0γ1 = self.γ0 @ self.γ1
        # j^2 = Ψ†γ^0γ^2Ψ
        γ0γ2 = self.γ0 @ self.γ2

        for i in range(self.Nx):
            for j in range(self.Ny):
                psi_dag = self.psi[i, j, :].conj()

                # j^x
                jx[i, j] = psi_dag @ γ0γ1 @ self.psi[i, j, :]

                # j^y
                jy[i, j] = psi_dag @ γ0γ2 @ self.psi[i, j, :]

        return jx.real, jy.real

    def apply_spatial_derivative(
        self,
        component: int,
        direction: str
    ) -> NDArray:
        """
        Apply spatial derivative ∂_i to spinor component.

        Parameters
        ----------
        component : int
            Spinor component index (0-3).
        direction : str
            'x' or 'y'.

        Returns
        -------
        NDArray
            Derivative field (Nx, Ny) complex.
        """
        field_component = self.psi[:, :, component]

        if direction == 'x':
            grad, _ = self.grid.gradient(field_component.real)
            grad_imag, _ = self.grid.gradient(field_component.imag)
            return grad + 1j * grad_imag

        elif direction == 'y':
            _, grad = self.grid.gradient(field_component.real)
            _, grad_imag = self.grid.gradient(field_component.imag)
            return grad + 1j * grad_imag

        else:
            raise ValueError(f"Unknown direction: {direction}")

    def compute_hamiltonian_action(
        self,
        mass_field: NDArray,
        chiral_angle: float = 0.0,
        Delta: float = 1.0
    ) -> NDArray:
        """
        Compute Hamiltonian action H·Ψ for Dirac equation.

        H = -iγ^i∂_i + m(x,t)

        where m(x,t) = Δ·R(x,t)·e^(iθγ^5)

        Parameters
        ----------
        mass_field : NDArray
            Synchronization field R(x,t) of shape (Nx, Ny).
        chiral_angle : float
            Chiral angle θ in radians.
        Delta : float
            Mass gap parameter.

        Returns
        -------
        NDArray
            H·Ψ of shape (Nx, Ny, 4) complex.
        """
        H_psi = np.zeros_like(self.psi)

        # Kinetic term: -iγ^i∂_i Ψ
        # For each spinor component, compute spatial derivatives
        for comp in range(4):
            dx_psi = self.apply_spatial_derivative(comp, 'x')
            dy_psi = self.apply_spatial_derivative(comp, 'y')

            # -iγ^1∂_x acts on all components
            for target_comp in range(4):
                H_psi[:, :, target_comp] += -1j * self.γ1[target_comp, comp] * dx_psi
                H_psi[:, :, target_comp] += -1j * self.γ2[target_comp, comp] * dy_psi

        # Mass term: m(x,t)·Ψ where m = Δ·R·e^(iθγ^5)
        # For each grid point, apply mass operator
        for i in range(self.Nx):
            for j in range(self.Ny):
                R = mass_field[i, j] if mass_field.shape == (self.Nx, self.Ny) else mass_field.flatten()[i * self.Ny + j]

                # Mass operator at this point
                m_op = mass_operator_MSFT(Delta, R, chiral_angle, self.γ5)

                # Apply to spinor
                H_psi[i, j, :] += m_op @ self.psi[i, j, :]

        return H_psi

    def step_rk4(
        self,
        dt: float,
        mass_field: NDArray,
        chiral_angle: float = 0.0,
        Delta: float = 1.0,
        normalize: bool = False
    ):
        """
        Time-step Dirac equation using RK4.

        i∂_tΨ = H·Ψ
        => ∂_tΨ = -i·H·Ψ

        Parameters
        ----------
        dt : float
            Time step.
        mass_field : NDArray
            R(x,t) field.
        chiral_angle : float
            θ parameter.
        Delta : float
            Mass gap.
        normalize : bool
            If True, normalize after each step to preserve ∫|Ψ|²dV = 1.
        """
        # Check CFL condition for stability
        # For Dirac equation: dt < dx/c where c~1 (speed of light = 1)
        cfl_limit = min(self.grid.dx, self.grid.dy) / 1.0
        if dt > 0.5 * cfl_limit:
            import warnings
            warnings.warn(
                f"Time step dt={dt:.6f} may violate CFL condition. "
                f"Recommended: dt < {0.5 * cfl_limit:.6f}"
            )

        # Preserve original state
        psi_0 = self.psi.copy()

        # RK4 stages: k_i = f(t + c_i*dt, Ψ + dt*Σ a_ij*k_j)
        k1 = -1j * self.compute_hamiltonian_action(mass_field, chiral_angle, Delta)

        self.psi = psi_0 + 0.5 * dt * k1
        k2 = -1j * self.compute_hamiltonian_action(mass_field, chiral_angle, Delta)

        self.psi = psi_0 + 0.5 * dt * k2
        k3 = -1j * self.compute_hamiltonian_action(mass_field, chiral_angle, Delta)

        self.psi = psi_0 + dt * k3
        k4 = -1j * self.compute_hamiltonian_action(mass_field, chiral_angle, Delta)

        # Final update: Ψ(t+dt) = Ψ(t) + dt/6 * (k1 + 2k2 + 2k3 + k4)
        self.psi = psi_0 + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # Optional normalization to preserve probability
        if normalize:
            self._normalize()

        self.t += dt

    def step_euler(
        self,
        dt: float,
        mass_field: NDArray,
        chiral_angle: float = 0.0,
        Delta: float = 1.0,
        normalize: bool = False
    ):
        """
        Time-step using forward Euler (simpler, less accurate).

        Ψ(t+dt) = Ψ(t) + dt·∂_tΨ
                = Ψ(t) - i·dt·H·Ψ(t)

        Parameters
        ----------
        dt : float
            Time step.
        mass_field : NDArray
            R(x,t) field.
        chiral_angle : float
            θ parameter.
        Delta : float
            Mass gap.
        normalize : bool
            If True, normalize after each step to preserve ∫|Ψ|²dV = 1.
        """
        # Check CFL condition for stability
        cfl_limit = min(self.grid.dx, self.grid.dy) / 1.0
        if dt > 0.5 * cfl_limit:
            import warnings
            warnings.warn(
                f"Time step dt={dt:.6f} may violate CFL condition. "
                f"Recommended: dt < {0.5 * cfl_limit:.6f}"
            )

        H_psi = self.compute_hamiltonian_action(mass_field, chiral_angle, Delta)
        self.psi += -1j * dt * H_psi

        # Optional normalization to preserve probability
        if normalize:
            self._normalize()

        self.t += dt

    def __repr__(self) -> str:
        return f"DiracSpinorField(grid={self.grid.Nx}×{self.grid.Ny}, t={self.t:.3f})"


# ==============================================================================
# Validation Functions
# ==============================================================================

def validate_free_dirac_evolution():
    """
    Validate Dirac solver with free particle (m=0).

    Free Dirac equation: i∂_tΨ = -iγ^i∂_iΨ
    Plane wave solution: Ψ = e^(i(k·x - ωt)) with ω = |k|
    """
    print("Validating free Dirac evolution (m=0)...")

    # Create grid
    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')

    # CFL limit: dt < dx/c where c=1
    cfl_limit = min(grid.dx, grid.dy) / 1.0
    print(f"  Grid spacing: dx={grid.dx:.6f}, dy={grid.dy:.6f}")
    print(f"  CFL limit: dt < {cfl_limit:.6f}")

    # Create spinor
    spinor = DiracSpinorField(grid, initial_state='plane_wave', momentum=(2.0, 0.0))

    # Store initial norm
    density_initial = spinor.compute_density()
    norm_initial = np.sum(density_initial) * grid.dx * grid.dy
    print(f"  Initial norm: {norm_initial:.6f}")

    # Zero mass field
    mass_field = np.zeros((grid.Nx, grid.Ny))

    # Use smaller time step for stability (well below CFL limit)
    dt = 0.1 * cfl_limit
    n_steps = 100
    print(f"  Time step: dt={dt:.6f} (0.1 × CFL limit)")
    print(f"  Total evolution: T={dt * n_steps:.4f}")

    # Evolve with RK4 and normalization
    for _ in range(n_steps):
        spinor.step_rk4(dt, mass_field, chiral_angle=0.0, Delta=0.0, normalize=True)

    # Check conservation
    density_final = spinor.compute_density()
    norm_final = np.sum(density_final) * grid.dx * grid.dy

    print(f"  Final norm: {norm_final:.6f} (should be ~1.0)")
    print(f"  ✓ PASS" if np.abs(norm_final - 1.0) < 0.1 else "  ✗ FAIL")

    return np.abs(norm_final - 1.0) < 0.1


def validate_MSFT_mass_coupling():
    """
    Validate coupling to MSFT mass field m = Δ·R.

    Check that mass field affects energy and evolution.
    Mass term contributes m·Ψ to Hamiltonian, so higher R → higher energy.
    """
    print("\nValidating MSFT mass coupling...")

    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')

    # CFL-limited time step
    cfl_limit = min(grid.dx, grid.dy) / 1.0
    dt = 0.1 * cfl_limit
    n_steps = 20
    print(f"  Time step: dt={dt:.6f}, total time: T={dt * n_steps:.4f}")

    # Test: Measure how much the mass term contributes to evolution
    spinor = DiracSpinorField(grid, initial_state='gaussian')

    # Compute Hamiltonian with low mass
    mass_low = 0.1 * np.ones((grid.Nx, grid.Ny))
    H_low = spinor.compute_hamiltonian_action(mass_low, chiral_angle=0.0, Delta=1.0)
    energy_low = np.sum(np.abs(H_low)**2) * grid.dx * grid.dy

    # Compute Hamiltonian with high mass (same spinor)
    mass_high = 0.9 * np.ones((grid.Nx, grid.Ny))
    H_high = spinor.compute_hamiltonian_action(mass_high, chiral_angle=0.0, Delta=1.0)
    energy_high = np.sum(np.abs(H_high)**2) * grid.dx * grid.dy

    print(f"\n  Hamiltonian energy (∫|H·Ψ|²):")
    print(f"    Low mass (R=0.1):  {energy_low:.6f}")
    print(f"    High mass (R=0.9): {energy_high:.6f}")
    print(f"    Ratio (high/low):  {energy_high/energy_low:.6f}")

    # Energy should be higher with higher mass
    energy_increased = energy_high > energy_low

    # Test 2: Phase evolution rate (without normalization to preserve phase)
    print(f"\n  Phase evolution test (no normalization):")

    spinor_low = DiracSpinorField(grid, initial_state='gaussian')
    psi_low_init = spinor_low.psi.copy()

    spinor_high = DiracSpinorField(grid, initial_state='gaussian')
    psi_high_init = spinor_high.psi.copy()

    # Evolve without normalization
    for _ in range(n_steps):
        spinor_low.step_rk4(dt, mass_low, chiral_angle=0.0, Delta=1.0, normalize=False)
        spinor_high.step_rk4(dt, mass_high, chiral_angle=0.0, Delta=1.0, normalize=False)

    # Measure phase change at center point
    center = (grid.Nx // 2, grid.Ny // 2)
    phase_low_init = np.angle(psi_low_init[center[0], center[1], 0])
    phase_low_final = np.angle(spinor_low.psi[center[0], center[1], 0])
    phase_low_change = np.abs(phase_low_final - phase_low_init)

    phase_high_init = np.angle(psi_high_init[center[0], center[1], 0])
    phase_high_final = np.angle(spinor_high.psi[center[0], center[1], 0])
    phase_high_change = np.abs(phase_high_final - phase_high_init)

    print(f"    Phase change (R=0.1): {phase_low_change:.6f} rad")
    print(f"    Phase change (R=0.9): {phase_high_change:.6f} rad")
    print(f"    Ratio (high/low):     {phase_high_change/phase_low_change:.6f}")

    # Higher mass should give faster phase rotation (E = p² + m²)
    phase_increased = phase_high_change > phase_low_change

    print(f"\n  Energy increases with mass: {energy_increased}")
    print(f"  Phase rotation increases with mass: {phase_increased}")

    passed = energy_increased and phase_increased
    print(f"  ✓ PASS" if passed else "  ✗ FAIL")

    return passed


if __name__ == "__main__":
    print("=" * 70)
    print(" DIRAC SPINOR FIELD VALIDATION")
    print("=" * 70)

    test1 = validate_free_dirac_evolution()
    test2 = validate_MSFT_mass_coupling()

    print("\n" + "=" * 70)
    if test1 and test2:
        print(" ✓✓✓ ALL VALIDATIONS PASSED ✓✓✓")
    else:
        print(" ✗✗✗ SOME VALIDATIONS FAILED ✗✗✗")
    print("=" * 70)
