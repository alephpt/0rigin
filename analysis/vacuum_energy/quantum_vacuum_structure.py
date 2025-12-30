#!/usr/bin/env python3
"""
Quantum Vacuum Structure Analysis for SMFT
Direction 1: Compute quantum corrections to ⟨R²⟩ from zero-point fluctuations and defect condensation

Physical Model:
- R(x) = R₀ + δR(x)  (quantum fluctuations around classical vacuum)
- ⟨R²⟩_quantum = ⟨R²⟩_classical + δ⟨R²⟩_ZPF + δ⟨R²⟩_defects
- Goal: Achieve ⟨R²⟩_quantum ~ 10^-123 from first principles
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path
from scipy.integrate import simpson
from scipy.optimize import curve_fit

# Physical constants (Planck units: ℏ = c = G = 1)
M_PLANCK = 1.22e19  # GeV
DELTA = M_PLANCK  # Mass gap = Planck mass


def compute_kuramoto_dispersion(k, K, Delta):
    """
    Compute dispersion relation for Kuramoto synchronization field

    Linearized around synchronized state:
    ω_k² = K²k² + Δ²  (phonon-like)

    Args:
        k: Momentum magnitude
        K: Coupling strength
        Delta: Mass gap

    Returns:
        omega_k: Mode frequency
    """
    omega_squared = K**2 * k**2 + Delta**2
    return np.sqrt(omega_squared)


def compute_zero_point_energy(grid_size, Delta, K=None):
    """
    Compute zero-point fluctuation contribution to ⟨R²⟩

    δ⟨R²⟩_ZPF = ∫d³k/(2π)³ × (ℏω_k)/(2Δ)

    In 2D (simulation):
    δ⟨R²⟩_ZPF = ∫dk k/(2π) × ω_k/(2Δ)

    Args:
        grid_size: Number of grid points (sets UV cutoff)
        Delta: Mass gap parameter
        K: Coupling strength (default: Delta)

    Returns:
        ZPF_contribution: Zero-point fluctuation to ⟨R²⟩
        mode_energies: Energy per mode [array]
    """
    if K is None:
        K = Delta

    # UV cutoff: k_max = π/a where a = L/N is lattice spacing
    # For domain size L = 100 ℓ_P and grid N:
    L_domain = 100.0
    a = L_domain / grid_size
    k_max = np.pi / a

    # Momentum grid (2D)
    n_k_points = 100
    k_values = np.linspace(0, k_max, n_k_points)

    # Dispersion relation
    omega_k = compute_kuramoto_dispersion(k_values, K, Delta)

    # Zero-point energy per mode (in 2D)
    # E_ZP = (1/2) ∫dk k ω_k / (2π)
    integrand = k_values * omega_k / (2 * np.pi)
    ZPF_energy = simpson(integrand, k_values)

    # Convert to R² contribution
    # ⟨R²⟩_ZPF = E_ZP / Δ
    ZPF_R_squared = ZPF_energy / Delta

    mode_energies = 0.5 * omega_k  # Energy per mode

    return ZPF_R_squared, mode_energies, k_values


def compute_defect_density(Delta, T_quantum=None):
    """
    Compute equilibrium density of virtual topological defects

    Vortex creation probability:
    P_vortex ~ exp(-E_vortex / T_quantum)

    where:
    E_vortex = πΔ  (2D vortex core energy)
    T_quantum = Δ  (quantum temperature scale)

    Args:
        Delta: Mass gap parameter
        T_quantum: Quantum temperature (default: Delta)

    Returns:
        n_defects: Defect density (per Planck area)
        P_creation: Creation probability
    """
    if T_quantum is None:
        T_quantum = Delta

    # Vortex core energy (2D)
    E_vortex = np.pi * Delta

    # Boltzmann factor
    P_creation = np.exp(-E_vortex / T_quantum)

    # Defect density (dimensional analysis)
    # n ~ (Δ/2π) × P_creation
    n_defects = (Delta / (2 * np.pi)) * P_creation

    return n_defects, P_creation


def compute_defect_suppression(n_defects, Delta):
    """
    Compute suppression of ⟨R²⟩ from defect condensate

    Each vortex core has R_core ≈ 0 over area A_core ~ (2π/Δ)²

    Effective suppression:
    ⟨R²⟩_eff = ⟨R²⟩_bulk × (1 - n_defects × A_core)

    Args:
        n_defects: Defect density
        Delta: Mass gap

    Returns:
        suppression_factor: Multiplicative suppression
        filling_fraction: Fraction of space occupied by defects
    """
    # Vortex core radius
    r_core = 2.0 / Delta  # ~2 Planck lengths
    A_core = np.pi * r_core**2

    # Filling fraction
    f_defects = n_defects * A_core

    # Cannot exceed 1 (complete coverage)
    f_defects = min(f_defects, 1.0)

    # Suppression factor
    suppression = 1.0 - f_defects

    return suppression, f_defects


def compute_quantum_R_squared(grid_size, Delta=M_PLANCK, K=None, include_defects=True):
    """
    Master function: Compute full quantum-corrected ⟨R²⟩

    ⟨R²⟩_quantum = ⟨R²⟩_classical + δ⟨R²⟩_ZPF × f_defects

    Args:
        grid_size: Simulation grid size
        Delta: Mass gap parameter
        K: Coupling strength (default: Delta)
        include_defects: Include topological defect suppression

    Returns:
        R_squared_quantum: Quantum-corrected ⟨R²⟩
        components: Dict with breakdown
    """
    # Classical expectation (fully synchronized vacuum)
    R_squared_classical = 1.0

    # Zero-point fluctuations
    ZPF_contribution, mode_energies, k_values = compute_zero_point_energy(grid_size, Delta, K)

    # Topological defects
    if include_defects:
        n_defects, P_creation = compute_defect_density(Delta)
        defect_suppression, f_defects = compute_defect_suppression(n_defects, Delta)
    else:
        n_defects = 0
        P_creation = 0
        defect_suppression = 1.0
        f_defects = 0

    # Total quantum correction
    R_squared_quantum = (R_squared_classical + ZPF_contribution) * defect_suppression

    # Breakdown
    components = {
        'classical': R_squared_classical,
        'ZPF': ZPF_contribution,
        'defect_density': n_defects,
        'defect_suppression': defect_suppression,
        'filling_fraction': f_defects,
        'total_quantum': R_squared_quantum,
        'mode_energies': mode_energies,
        'k_values': k_values
    }

    return R_squared_quantum, components


def analyze_simulation_data(output_dir, N_value=1):
    """
    Extract ⟨R²⟩ from SMFT simulation data and compare to quantum prediction

    Args:
        output_dir: Path to simulation output
        N_value: Operator splitting ratio (1, 10, or 100)

    Returns:
        simulation_R_squared: Measured ⟨R²⟩ from simulation
        grid_size: Grid resolution
    """
    # Find R field snapshots
    pattern = f"{output_dir}/N_{N_value}/R_field_*.csv"
    data_files = glob.glob(pattern)

    if not data_files:
        return None, None

    # Load all snapshots
    R_squared_values = []

    for file in sorted(data_files):
        try:
            R_field = pd.read_csv(file, header=None).values
            R_squared = np.mean(R_field**2)
            R_squared_values.append(R_squared)
        except Exception as e:
            print(f"Could not load {file}: {e}")

    if not R_squared_values:
        return None, None

    # Time average (last 50% of data)
    n_points = len(R_squared_values)
    late_time = R_squared_values[int(0.5*n_points):]
    simulation_R_squared = np.mean(late_time)

    # Get grid size
    grid_size = R_field.shape[0]

    return simulation_R_squared, grid_size


def plot_quantum_analysis(grid_sizes, results_dict):
    """
    Create comprehensive quantum vacuum analysis plots

    Args:
        grid_sizes: Array of grid resolutions tested
        results_dict: Dict mapping grid_size → components dict
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Extract data
    classical_vals = [results_dict[N]['classical'] for N in grid_sizes]
    ZPF_vals = [results_dict[N]['ZPF'] for N in grid_sizes]
    defect_supp_vals = [results_dict[N]['defect_suppression'] for N in grid_sizes]
    quantum_vals = [results_dict[N]['total_quantum'] for N in grid_sizes]

    # Plot 1: ⟨R²⟩ components vs grid size
    ax = axes[0, 0]
    ax.semilogy(grid_sizes, classical_vals, 'k--', label='Classical (⟨R²⟩ = 1)', linewidth=2)
    ax.semilogy(grid_sizes, np.abs(ZPF_vals), 'b-o', label='|ZPF contribution|')
    ax.semilogy(grid_sizes, quantum_vals, 'r-s', label='Total quantum', linewidth=2)
    ax.axhline(y=1e-123, color='g', linestyle=':', linewidth=2, label='Required (10⁻¹²³)')
    ax.set_xlabel('Grid Size N')
    ax.set_ylabel('⟨R²⟩')
    ax.set_title('Quantum Corrections to Vacuum Expectation Value')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Defect suppression
    ax = axes[0, 1]
    filling_fractions = [results_dict[N]['filling_fraction'] for N in grid_sizes]
    ax.plot(grid_sizes, defect_supp_vals, 'mo-', label='Suppression factor')
    ax.plot(grid_sizes, filling_fractions, 'co-', label='Defect filling fraction')
    ax.set_xlabel('Grid Size N')
    ax.set_ylabel('Suppression / Filling')
    ax.set_title('Topological Defect Contributions')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Mode energy spectrum (for largest grid)
    ax = axes[0, 2]
    largest_grid = max(grid_sizes)
    k_vals = results_dict[largest_grid]['k_values']
    mode_E = results_dict[largest_grid]['mode_energies']
    ax.loglog(k_vals[k_vals > 0], mode_E[k_vals > 0], 'b-', linewidth=2)
    ax.set_xlabel('Momentum k (1/ℓ_P)')
    ax.set_ylabel('Mode Energy (M_P)')
    ax.set_title(f'Zero-Point Mode Spectrum (N={largest_grid})')
    ax.grid(True, alpha=0.3)

    # Plot 4: Suppression factor vs required
    ax = axes[1, 0]
    suppression_achieved = [q / classical_vals[0] for q in quantum_vals]
    ax.semilogy(grid_sizes, suppression_achieved, 'ro-', linewidth=2, markersize=8, label='Achieved')
    ax.axhline(y=1e-123, color='g', linestyle='--', linewidth=2, label='Required')
    ax.fill_between(grid_sizes, 1e-130, 1e-116, alpha=0.2, color='green', label='Success zone')
    ax.set_xlabel('Grid Size N')
    ax.set_ylabel('Suppression Factor')
    ax.set_title('Vacuum Energy Suppression Achievement')
    ax.set_ylim([1e-140, 1e10])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 5: Convergence test
    ax = axes[1, 1]
    if len(grid_sizes) > 1:
        # Compute convergence: |⟨R²⟩_{N+1} - ⟨R²⟩_N| / ⟨R²⟩_N
        convergence = []
        grid_pairs = []
        for i in range(len(grid_sizes)-1):
            N1, N2 = grid_sizes[i], grid_sizes[i+1]
            conv = abs(quantum_vals[i+1] - quantum_vals[i]) / quantum_vals[i]
            convergence.append(conv)
            grid_pairs.append(N2)

        ax.semilogy(grid_pairs, convergence, 'bs-', linewidth=2, markersize=8)
        ax.axhline(y=0.01, color='r', linestyle='--', label='1% convergence')
        ax.set_xlabel('Grid Size N')
        ax.set_ylabel('Relative Change')
        ax.set_title('Grid Convergence of Quantum ⟨R²⟩')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 6: Summary table
    ax = axes[1, 2]
    ax.axis('off')

    summary_text = "Quantum Vacuum Analysis Summary\n"
    summary_text += "=" * 40 + "\n\n"

    # Best (largest grid) results
    best_idx = -1
    summary_text += f"Grid Size: {grid_sizes[best_idx]}×{grid_sizes[best_idx]}\n\n"
    summary_text += "Vacuum Expectation Values:\n"
    summary_text += f"  ⟨R²⟩_classical = {classical_vals[best_idx]:.3e}\n"
    summary_text += f"  δ⟨R²⟩_ZPF     = {ZPF_vals[best_idx]:.3e}\n"
    summary_text += f"  ⟨R²⟩_quantum   = {quantum_vals[best_idx]:.3e}\n\n"

    summary_text += "Defect Condensate:\n"
    summary_text += f"  n_defects = {results_dict[grid_sizes[best_idx]]['defect_density']:.3e} / ℓ_P²\n"
    summary_text += f"  Filling   = {filling_fractions[best_idx]:.1%}\n"
    summary_text += f"  Suppression = {defect_supp_vals[best_idx]:.3e}\n\n"

    summary_text += "Vacuum Energy:\n"
    summary_text += f"  Achieved: {suppression_achieved[best_idx]:.3e}\n"
    summary_text += f"  Required: 1.0e-123\n"
    summary_text += f"  Missing:  {suppression_achieved[best_idx]/1e-123:.3e}\n\n"

    # Verdict
    if suppression_achieved[best_idx] < 1e-120:
        verdict = "✓ SUCCESS"
        color = 'green'
    elif suppression_achieved[best_idx] < 1e-100:
        verdict = "⚠ PARTIAL"
        color = 'orange'
    else:
        verdict = "✗ FAILED"
        color = 'red'

    summary_text += f"Status: {verdict}"

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    return fig


def main():
    """Main analysis routine"""
    print("=" * 60)
    print("SMFT Quantum Vacuum Structure Analysis")
    print("Direction 1: Zero-Point Fluctuations + Defect Condensation")
    print("=" * 60)

    # Grid sizes to test
    grid_sizes = [16, 32, 64, 128, 256]

    print("\nComputing quantum corrections...")
    print("-" * 60)

    results_dict = {}

    for N in grid_sizes:
        R_sq_quantum, components = compute_quantum_R_squared(
            grid_size=N,
            Delta=DELTA,
            K=DELTA,
            include_defects=True
        )

        results_dict[N] = components

        print(f"\nGrid {N}×{N}:")
        print(f"  ⟨R²⟩_classical = {components['classical']:.3e}")
        print(f"  δ⟨R²⟩_ZPF     = {components['ZPF']:.3e}")
        print(f"  Defect supp.  = {components['defect_suppression']:.3e}")
        print(f"  ⟨R²⟩_quantum   = {components['total_quantum']:.3e}")
        print(f"  Suppression   = {components['total_quantum']/components['classical']:.3e}")

    # Create plots
    fig = plot_quantum_analysis(grid_sizes, results_dict)

    # Save
    output_path = "/home/persist/neotec/0rigin/output/quantum_vacuum_analysis.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")

    # Final assessment
    print("\n" + "=" * 60)
    print("ASSESSMENT")
    print("=" * 60)

    best_suppression = min([results_dict[N]['total_quantum'] for N in grid_sizes])
    required_suppression = 1e-123
    missing_factor = best_suppression / required_suppression

    print(f"Best achieved suppression: {best_suppression:.3e}")
    print(f"Required suppression:      {required_suppression:.3e}")
    print(f"Missing factor:            {missing_factor:.3e}")

    if missing_factor < 1e3:
        print("\n✓ SUCCESS: Quantum corrections achieve required suppression!")
    elif missing_factor < 1e50:
        print("\n⚠ PARTIAL: Quantum corrections help, but insufficient alone")
        print("  Recommendation: Combine with RG running (Direction 2)")
    else:
        print("\n✗ FAILED: Quantum corrections cannot resolve CC problem")
        print("  Recommendation: Must pursue alternative directions (2, 3, or 4)")

    print("\nKey Findings:")
    print(f"1. Zero-point fluctuations: {'ENHANCE' if results_dict[64]['ZPF'] > 0 else 'SUPPRESS'} vacuum energy")
    print(f"2. Defect condensation: {results_dict[64]['filling_fraction']:.1%} filling fraction")
    print(f"3. Net quantum effect: {best_suppression:.3e} suppression")
    print(f"4. Grid convergence: {'CONVERGED' if abs(results_dict[128]['total_quantum'] - results_dict[256]['total_quantum']) / results_dict[128]['total_quantum'] < 0.01 else 'NOT CONVERGED'}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
