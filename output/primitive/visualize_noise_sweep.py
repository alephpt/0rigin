#!/usr/bin/env python3
"""
MSFT Noise Sweep Visualization Tool

Visualizes noise sweep experiment results:
- Summary statistics vs sigma (phase transition detection)
- Timeseries evolution for different noise levels
- Spatial field snapshots at different noise levels
- Critical noise analysis (σ_c determination)

Data structure:
  noise_sweep/
    summary.dat         - Final statistics for each σ
    results.csv         - CSV format summary
    timeseries_sigma_*.dat - Time evolution at specific σ values
    sigma_X.XXeXX/      - Per-sigma directories with:
      R_global_timeseries.dat
      R_local_timeseries.dat
      L_timeseries.dat
      phase_var_timeseries.dat
      step_XXXX/        - Spatial snapshots
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
import os
from pathlib import Path
import glob


def load_summary_data(filepath):
    """
    Load summary.dat file with noise sweep results.

    Returns:
        dict with keys: sigma, R_initial, R_final, R_mean, R_std, L_final, L_mean
    """
    try:
        # Skip comment lines starting with #
        data = np.loadtxt(filepath, comments='#')

        if data.ndim == 1:
            data = data.reshape(1, -1)

        return {
            'sigma': data[:, 0],
            'R_initial': data[:, 1],
            'R_final': data[:, 2],
            'R_mean': data[:, 3],
            'R_std': data[:, 4],
            'L_final': data[:, 5],
            'L_mean': data[:, 6]
        }
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_timeseries_file(filepath):
    """
    Load timeseries file (4 columns: step, t, R_global, L).

    Returns:
        dict with keys: step, t, R_global, L
    """
    try:
        data = np.loadtxt(filepath, comments='#')

        if data.ndim == 1:
            data = data.reshape(1, -1)

        return {
            'step': data[:, 0],
            't': data[:, 1],
            'R_global': data[:, 2],
            'L': data[:, 3]
        }
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_sigma_directory_timeseries(sigma_dir):
    """
    Load all timeseries from a sigma_X.XXeXX directory.

    Returns:
        dict with R_global, R_local, L, phase_var timeseries
    """
    result = {}

    files = {
        'R_global': 'R_global_timeseries.dat',
        'R_local': 'R_local_timeseries.dat',
        'L': 'L_timeseries.dat',
        'phase_var': 'phase_var_timeseries.dat'
    }

    for key, fname in files.items():
        fpath = os.path.join(sigma_dir, fname)
        if os.path.exists(fpath):
            data = np.loadtxt(fpath)
            result[key] = data

    return result


def plot_phase_transition(summary, output_file=None, dpi=150):
    """
    Plot summary statistics vs sigma to identify phase transition.

    Main plot: R_final vs sigma (should show transition at σ_c)
    """
    sigma = summary['sigma']
    R_final = summary['R_final']
    R_std = summary['R_std']
    L_final = summary['L_final']
    L_mean = summary['L_mean']

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    fig.suptitle('MSFT Noise Sweep: Phase Transition Analysis',
                 fontsize=14, fontweight='bold')

    # Plot 1: R_final vs sigma (main diagnostic)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.semilogx(sigma, R_final, 'o-', linewidth=2, markersize=8, color='blue')
    ax1.set_xlabel('Noise σ', fontsize=11)
    ax1.set_ylabel('R_final (Sync Order)', fontsize=11)
    ax1.set_title('Order Parameter vs Noise\n(Phase Transition)', fontweight='bold')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='50% sync')
    ax1.legend()

    # Find approximate σ_c (where R drops below threshold)
    threshold = 0.9
    critical_idx = np.where(R_final < threshold)[0]
    if len(critical_idx) > 0:
        sigma_c_approx = sigma[critical_idx[0]]
        ax1.axvline(sigma_c_approx, color='orange', linestyle='--',
                   label=f'σ_c ≈ {sigma_c_approx:.2e}')
        ax1.legend()

    # Plot 2: R_std vs sigma (fluctuations peak at transition)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.semilogx(sigma, R_std, 'o-', linewidth=2, markersize=8, color='green')
    ax2.set_xlabel('Noise σ', fontsize=11)
    ax2.set_ylabel('R_std (Fluctuations)', fontsize=11)
    ax2.set_title('Order Fluctuations\n(Peak at σ_c)', fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')

    # Plot 3: L_final vs sigma (correlation length)
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.loglog(sigma, L_final, 'o-', linewidth=2, markersize=8, color='purple')
    ax3.set_xlabel('Noise σ', fontsize=11)
    ax3.set_ylabel('L_final (Cluster Size)', fontsize=11)
    ax3.set_title('Correlation Length\n(Diverges at σ_c)', fontweight='bold')
    ax3.grid(True, alpha=0.3, which='both')

    # Plot 4: Phase diagram (R vs sigma, colored by L)
    ax4 = fig.add_subplot(gs[1, 0])
    scatter = ax4.scatter(sigma, R_final, c=L_final, s=100, cmap='viridis',
                         edgecolors='black', linewidth=1.5)
    ax4.set_xscale('log')
    ax4.set_xlabel('Noise σ', fontsize=11)
    ax4.set_ylabel('R_final', fontsize=11)
    ax4.set_title('Phase Diagram\n(Color = Cluster Size)', fontweight='bold')
    ax4.grid(True, alpha=0.3, which='both')
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('L_final')

    # Plot 5: Susceptibility (dR/dσ)
    ax5 = fig.add_subplot(gs[1, 1])
    if len(sigma) > 1:
        # Numerical derivative
        dR_dsigma = np.gradient(R_final, sigma)
        ax5.semilogx(sigma, np.abs(dR_dsigma), 'o-', linewidth=2,
                    markersize=8, color='red')
        ax5.set_xlabel('Noise σ', fontsize=11)
        ax5.set_ylabel('|dR/dσ| (Susceptibility)', fontsize=11)
        ax5.set_title('Susceptibility\n(Peaks at σ_c)', fontweight='bold')
        ax5.grid(True, alpha=0.3, which='both')

    # Plot 6: Statistics table
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    # Create table data
    table_data = []
    for i in range(len(sigma)):
        row = [
            f"{sigma[i]:.1e}",
            f"{R_final[i]:.3f}",
            f"{R_std[i]:.3f}",
            f"{int(L_final[i])}"
        ]
        table_data.append(row)

    table = ax6.table(cellText=table_data,
                     colLabels=['σ', 'R_final', 'R_std', 'L_final'],
                     cellLoc='center',
                     loc='center',
                     bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)

    # Style header
    for i in range(4):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax6.set_title('Summary Statistics', fontweight='bold', pad=20)

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def plot_timeseries_comparison(noise_sweep_dir, sigma_values=None,
                               output_file=None, dpi=150):
    """
    Plot timeseries evolution for multiple sigma values.

    Args:
        noise_sweep_dir: Path to noise_sweep directory
        sigma_values: List of sigma values to plot (e.g., [0, 1e-6, 1e-3, 1])
                     If None, plots available timeseries_sigma_*.dat files
    """
    # Find available timeseries files
    pattern = os.path.join(noise_sweep_dir, 'timeseries_sigma_*.dat')
    timeseries_files = sorted(glob.glob(pattern))

    if not timeseries_files:
        print("No timeseries files found!")
        return

    # Load data
    datasets = []
    labels = []

    for fpath in timeseries_files:
        fname = os.path.basename(fpath)
        # Extract sigma index from filename
        sigma_idx = fname.replace('timeseries_sigma_', '').replace('.dat', '')

        data = load_timeseries_file(fpath)
        datasets.append(data)
        labels.append(f'σ index {sigma_idx}')

    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('MSFT Noise Sweep: Timeseries Evolution',
                 fontsize=14, fontweight='bold')

    colors = plt.cm.viridis(np.linspace(0, 1, len(datasets)))

    # Plot 1: R_global vs time
    ax1 = axes[0, 0]
    for i, (data, label) in enumerate(zip(datasets, labels)):
        ax1.plot(data['t'], data['R_global'], label=label,
                color=colors[i], linewidth=1.5, alpha=0.8)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('R_global (Sync Order)')
    ax1.set_title('Global Synchronization Evolution')
    ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Plot 2: L vs time
    ax2 = axes[0, 1]
    for i, (data, label) in enumerate(zip(datasets, labels)):
        ax2.plot(data['t'], data['L'], label=label,
                color=colors[i], linewidth=1.5, alpha=0.8)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('L (Cluster Size)')
    ax2.set_title('Correlation Length Evolution')
    ax2.legend(loc='best', fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Plot 3: R vs L phase space
    ax3 = axes[1, 0]
    for i, (data, label) in enumerate(zip(datasets, labels)):
        ax3.plot(data['R_global'], data['L'], label=label,
                color=colors[i], linewidth=1.5, alpha=0.7)
    ax3.set_xlabel('R_global')
    ax3.set_ylabel('L (Cluster Size)')
    ax3.set_title('Phase Space (R vs L)')
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Final values comparison
    ax4 = axes[1, 1]
    final_R = [data['R_global'][-1] for data in datasets]
    final_L = [data['L'][-1] for data in datasets]
    indices = range(len(datasets))

    ax4_twin = ax4.twinx()
    bars1 = ax4.bar([i - 0.2 for i in indices], final_R, width=0.4,
                    label='R_final', color='blue', alpha=0.7)
    bars2 = ax4_twin.bar([i + 0.2 for i in indices], final_L, width=0.4,
                         label='L_final', color='orange', alpha=0.7)

    ax4.set_xlabel('Dataset Index')
    ax4.set_ylabel('R_final', color='blue')
    ax4_twin.set_ylabel('L_final', color='orange')
    ax4.set_title('Final State Comparison')
    ax4.tick_params(axis='y', labelcolor='blue')
    ax4_twin.tick_params(axis='y', labelcolor='orange')
    ax4.set_xticks(indices)
    ax4.grid(True, alpha=0.3, axis='x')

    # Combined legend
    lines1, labels1 = ax4.get_legend_handles_labels()
    lines2, labels2 = ax4_twin.get_legend_handles_labels()
    ax4.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.tight_layout()

    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def analyze_critical_noise(summary):
    """
    Analyze data to estimate critical noise σ_c.

    Uses multiple methods:
    1. R threshold crossing
    2. Maximum susceptibility |dR/dσ|
    3. Maximum fluctuations R_std
    """
    sigma = summary['sigma']
    R_final = summary['R_final']
    R_std = summary['R_std']

    print("\n" + "="*60)
    print("CRITICAL NOISE ANALYSIS")
    print("="*60)

    # Method 1: Threshold crossing (R < 0.9)
    threshold = 0.9
    critical_idx = np.where(R_final < threshold)[0]
    if len(critical_idx) > 0:
        sigma_c_threshold = sigma[critical_idx[0]]
        print(f"\nMethod 1 (R < {threshold}):")
        print(f"  σ_c ≈ {sigma_c_threshold:.3e}")

    # Method 2: Maximum susceptibility
    if len(sigma) > 2:
        dR_dsigma = np.abs(np.gradient(R_final, sigma))
        max_idx = np.argmax(dR_dsigma)
        sigma_c_suscept = sigma[max_idx]
        print(f"\nMethod 2 (max |dR/dσ|):")
        print(f"  σ_c ≈ {sigma_c_suscept:.3e}")
        print(f"  |dR/dσ|_max = {dR_dsigma[max_idx]:.3e}")

    # Method 3: Maximum fluctuations
    if len(R_std) > 0:
        max_std_idx = np.argmax(R_std)
        sigma_c_fluct = sigma[max_std_idx]
        print(f"\nMethod 3 (max R_std):")
        print(f"  σ_c ≈ {sigma_c_fluct:.3e}")
        print(f"  R_std_max = {R_std[max_std_idx]:.3e}")

    print("\n" + "="*60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize MSFT noise sweep experiment results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Phase transition plot
  %(prog)s --phase-transition

  # Timeseries comparison
  %(prog)s --timeseries

  # Critical noise analysis
  %(prog)s --analyze

  # All analyses
  %(prog)s --all

  # Save outputs
  %(prog)s --phase-transition --output phase_diagram.png --dpi 300
        """
    )

    parser.add_argument('-d', '--data-dir', default='noise_sweep',
                       help='Noise sweep directory (default: noise_sweep)')
    parser.add_argument('-p', '--phase-transition', action='store_true',
                       help='Plot phase transition analysis')
    parser.add_argument('-t', '--timeseries', action='store_true',
                       help='Plot timeseries comparison')
    parser.add_argument('-a', '--analyze', action='store_true',
                       help='Analyze critical noise σ_c')
    parser.add_argument('--all', action='store_true',
                       help='Run all analyses')
    parser.add_argument('-o', '--output', help='Output filename (PNG/PDF/SVG)')
    parser.add_argument('--dpi', type=int, default=150,
                       help='Output resolution (default: 150)')

    args = parser.parse_args()

    # Default to all if nothing specified
    if not (args.phase_transition or args.timeseries or args.analyze):
        args.all = True

    if args.all:
        args.phase_transition = True
        args.timeseries = True
        args.analyze = True

    # Check data directory
    if not os.path.isdir(args.data_dir):
        print(f"Error: Directory not found: {args.data_dir}")
        return 1

    # Load summary data
    summary_file = os.path.join(args.data_dir, 'summary.dat')
    if not os.path.exists(summary_file):
        print(f"Error: summary.dat not found in {args.data_dir}")
        return 1

    summary = load_summary_data(summary_file)

    # Analysis
    if args.analyze:
        analyze_critical_noise(summary)

    # Visualizations
    try:
        if args.phase_transition:
            output = args.output if not args.all else 'noise_phase_transition.png'
            plot_phase_transition(summary, output, args.dpi)

        if args.timeseries:
            output = args.output if not args.all else 'noise_timeseries.png'
            plot_timeseries_comparison(args.data_dir, output_file=output, dpi=args.dpi)

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
