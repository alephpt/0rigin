#!/usr/bin/env python3
"""
Visualize EM Field Evolution from Test Outputs

Plots key observables over time:
  - E field magnitude evolution
  - B field magnitude evolution
  - Poynting flux (energy flow)
  - Norm and energy conservation

Usage:
  python3 visualize_em_fields.py <output_directory>

Example:
  python3 visualize_em_fields.py output/20251228_044054_em_coupling_uniform_field_64x64_v0.3/
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path

def plot_em_evolution(obs_df, output_path, test_name="EM Field Evolution"):
    """
    Create comprehensive visualization of EM field evolution.

    Plots:
    1. EM field energy over time
    2. E and B field magnitudes
    3. Poynting vector (energy flux)
    4. Norm and total energy conservation
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'{test_name}\nObservables Over Time', fontsize=14, weight='bold')

    has_em = 'field_energy' in obs_df.columns or 'max_E' in obs_df.columns

    # Plot 1: EM field energy (if available)
    ax = axes[0, 0]
    if 'field_energy' in obs_df.columns and has_em:
        ax.plot(obs_df['time'], obs_df['field_energy'], 'b-', linewidth=2.5, label='EM Field Energy')
        ax.fill_between(obs_df['time'], obs_df['field_energy'], alpha=0.3)
        ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
        ax.set_ylabel('Energy (m_P c²)', fontsize=11)
        ax.set_title('Electromagnetic Field Energy', fontsize=12, weight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')

        # Add statistics
        mean_energy = obs_df['field_energy'].mean()
        max_energy = obs_df['field_energy'].max()
        ax.text(0.98, 0.97, f'Mean: {mean_energy:.4e}\nMax: {max_energy:.4e}',
               transform=ax.transAxes, fontsize=9, verticalalignment='top',
               horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax.text(0.5, 0.5, 'EM field energy\nnot recorded\n(em_coupling disabled?)',
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    # Plot 2: E and B field magnitudes
    ax = axes[0, 1]
    has_e_or_b = 'max_E' in obs_df.columns or 'max_B' in obs_df.columns
    if has_e_or_b:
        if 'max_E' in obs_df.columns:
            ax.semilogy(obs_df['time'], obs_df['max_E'], 'r-', linewidth=2.5, label='Max |E|', marker='o', markersize=2)
        if 'max_B' in obs_df.columns:
            ax.semilogy(obs_df['time'], obs_df['max_B'], 'b-', linewidth=2.5, label='Max |B|', marker='s', markersize=2)

        ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
        ax.set_ylabel('Field Strength (log scale)', fontsize=11)
        ax.set_title('Electric and Magnetic Field Magnitudes', fontsize=12, weight='bold')
        ax.grid(True, alpha=0.3, which='both')
        ax.legend(fontsize=10, loc='best')
    else:
        ax.text(0.5, 0.5, 'E/B field data\nnot recorded',
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    # Plot 3: Poynting flux (energy flow)
    ax = axes[1, 0]
    if 'total_flux' in obs_df.columns:
        ax.plot(obs_df['time'], obs_df['total_flux'], 'g-', linewidth=2.5, label='Poynting Flux', marker='o', markersize=2)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=1)
        ax.fill_between(obs_df['time'], obs_df['total_flux'], alpha=0.3, color='green')

        ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
        ax.set_ylabel('Energy Flux (m_P c³)', fontsize=11)
        ax.set_title('Poynting Vector - Energy Flow', fontsize=12, weight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')

        # Add statistics
        mean_flux = obs_df['total_flux'].mean()
        max_flux = obs_df['total_flux'].max()
        ax.text(0.98, 0.97, f'Mean: {mean_flux:.4e}\nMax: {max_flux:.4e}',
               transform=ax.transAxes, fontsize=9, verticalalignment='top',
               horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    else:
        ax.text(0.5, 0.5, 'Poynting flux\nnot recorded',
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    # Plot 4: Norm and Energy Conservation (always present)
    ax = axes[1, 1]

    # Norm on primary axis
    color_norm = 'tab:blue'
    ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
    ax.set_ylabel('Norm ||ψ||²', color=color_norm, fontsize=11, weight='bold')
    ax.plot(obs_df['time'], obs_df['norm'], color=color_norm, linewidth=2.5, label='Norm')
    ax.tick_params(axis='y', labelcolor=color_norm)
    ax.grid(True, alpha=0.3)

    # Energy on secondary axis
    ax2 = ax.twinx()
    color_energy = 'tab:red'
    ax2.set_ylabel('Total Energy E (m_P c²)', color=color_energy, fontsize=11, weight='bold')
    ax2.plot(obs_df['time'], obs_df['E_total'], color=color_energy, linewidth=2.5, label='Energy', linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color_energy)

    ax.set_title('Conservation Laws: Norm & Energy', fontsize=12, weight='bold')

    # Add legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)

    # Calculate and display conservation metrics
    norm_initial = obs_df['norm'].iloc[0]
    norm_final = obs_df['norm'].iloc[-1]
    norm_drift = abs(norm_final - norm_initial) / norm_initial * 100

    energy_initial = obs_df['E_total'].iloc[0]
    energy_final = obs_df['E_total'].iloc[-1]
    energy_drift = abs(energy_final - energy_initial) / energy_initial * 100

    conservation_text = f"""Norm drift: {norm_drift:.4f}%
Energy drift: {energy_drift:.4f}%"""

    ax2.text(0.02, 0.97, conservation_text, transform=ax2.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7),
            family='monospace', weight='bold')

    plt.tight_layout()
    plt.savefig(str(output_path), dpi=150, bbox_inches='tight')
    return str(output_path)

def plot_order_parameter_evolution(obs_df, output_path):
    """Plot order parameter (R field) statistics over time."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Order Parameter Evolution', fontsize=14, weight='bold')

    # R field statistics
    ax = axes[0]
    ax.plot(obs_df['time'], obs_df['R_avg'], 'b-', linewidth=2.5, label='⟨R⟩ (average)', marker='o', markersize=3)
    ax.fill_between(obs_df['time'], obs_df['R_min'], obs_df['R_max'], alpha=0.2, color='blue', label='R_min to R_max')
    ax.plot(obs_df['time'], obs_df['R_max'], 'g--', linewidth=1.5, label='R_max', alpha=0.7)
    ax.plot(obs_df['time'], obs_df['R_min'], 'r--', linewidth=1.5, label='R_min', alpha=0.7)

    ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
    ax.set_ylabel('Order Parameter R(x,y,t)', fontsize=11)
    ax.set_title('R Field Statistics', fontsize=12, weight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    ax.set_ylim([0, 1.1])

    # R field variance
    ax = axes[1]
    ax.plot(obs_df['time'], obs_df['R_var'], 'purple', linewidth=2.5, marker='s', markersize=3, label='R variance')
    ax.fill_between(obs_df['time'], obs_df['R_var'], alpha=0.3, color='purple')

    ax.set_xlabel('Time (ℓ_P/c)', fontsize=11)
    ax.set_ylabel('Variance σ²(R)', fontsize=11)
    ax.set_title('Order Parameter Variance', fontsize=12, weight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    # Statistics box
    r_avg_mean = obs_df['R_avg'].mean()
    r_var_mean = obs_df['R_var'].mean()
    stats_text = f"""Statistics:
⟨R⟩_avg: {r_avg_mean:.6f}
σ²(R)_avg: {r_var_mean:.6f}"""

    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.7),
           family='monospace')

    plt.tight_layout()
    plt.savefig(str(output_path), dpi=150, bbox_inches='tight')
    return str(output_path)

def analyze_output_directory(output_dir):
    """
    Analyze all N_X subdirectories in an output directory.
    """
    output_dir = Path(output_dir)

    if not output_dir.exists():
        print(f"ERROR: Directory not found: {output_dir}")
        sys.exit(1)

    # Find all N_X subdirectories
    n_dirs = sorted([d for d in output_dir.iterdir() if d.is_dir() and d.name.startswith('N_')])

    if not n_dirs:
        print(f"ERROR: No N_X subdirectories found in {output_dir}")
        sys.exit(1)

    print(f"Analyzing: {output_dir.name}")
    print(f"Found {len(n_dirs)} test runs: {[d.name for d in n_dirs]}\n")

    # Process each N value
    for n_dir in n_dirs:
        obs_path = n_dir / "observables.csv"

        if not obs_path.exists():
            print(f"SKIP {n_dir.name}: observables.csv not found")
            continue

        print(f"Processing {n_dir.name}...")

        try:
            obs_df = pd.read_csv(obs_path)

            # Generate visualizations
            em_plot = n_dir / "em_fields_evolution.png"
            plot_em_evolution(obs_df, em_plot, f"{output_dir.name} - {n_dir.name}")
            print(f"  ✓ {em_plot.name}")

            if 'R_avg' in obs_df.columns:
                r_plot = n_dir / "order_parameter_evolution.png"
                plot_order_parameter_evolution(obs_df, r_plot)
                print(f"  ✓ {r_plot.name}")

        except Exception as e:
            print(f"ERROR processing {n_dir.name}: {e}")
            continue

    print(f"\n✓ Analysis complete for {output_dir.name}")

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nExample usage:")
        print("  python3 visualize_em_fields.py output/20251228_044054_em_coupling_uniform_field_64x64_v0.3/")
        sys.exit(0)

    output_dir = sys.argv[1]
    analyze_output_directory(output_dir)

if __name__ == "__main__":
    main()
