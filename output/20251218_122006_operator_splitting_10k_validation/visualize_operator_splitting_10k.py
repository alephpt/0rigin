#!/usr/bin/env python3
"""
Visualization script for operator splitting 10k timestep validation
Compares N=1, N=10, N=100 operator splitting approaches
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def load_timeseries(directory):
    """Load timeseries data from CSV file"""
    csv_path = Path(directory) / "timeseries.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    data = pd.read_csv(csv_path)
    return data

def plot_comparison(N1_data, N10_data, N100_data, output_dir):
    """Create comprehensive comparison plots"""

    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Operator Splitting 10k Timestep Validation\nN=1 (blue) vs N=10 (orange) vs N=100 (green)',
                 fontsize=14, fontweight='bold')

    # Plot 1: Synchronization field R_avg
    ax = axes[0, 0]
    ax.plot(N1_data['time'], N1_data['R_avg'], 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.plot(N10_data['time'], N10_data['R_avg'], 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.plot(N100_data['time'], N100_data['R_avg'], 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Average Sync Field R', fontsize=11)
    ax.set_title('Synchronization Field Evolution', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 2: Wavefunction norm
    ax = axes[0, 1]
    ax.plot(N1_data['time'], N1_data['norm'], 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.plot(N10_data['time'], N10_data['norm'], 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.plot(N100_data['time'], N100_data['norm'], 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Expected (1.0)')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Wavefunction Norm', fontsize=11)
    ax.set_title('Norm Conservation (Should be ~1.0)', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 3: Total energy
    ax = axes[1, 0]
    ax.plot(N1_data['time'], N1_data['energy'], 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.plot(N10_data['time'], N10_data['energy'], 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.plot(N100_data['time'], N100_data['energy'], 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Total Energy E', fontsize=11)
    ax.set_title('Energy Evolution', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 4: Maximum density
    ax = axes[1, 1]
    ax.plot(N1_data['time'], N1_data['max_density'], 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.plot(N10_data['time'], N10_data['max_density'], 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.plot(N100_data['time'], N100_data['max_density'], 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('Max Density', fontsize=11)
    ax.set_title('Wavepacket Dispersion', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/operator_splitting_comparison.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_dir}/operator_splitting_comparison.png")

    plt.close()

def plot_convergence_analysis(N1_data, N10_data, N100_data, output_dir):
    """Analyze convergence properties"""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Operator Splitting Convergence Analysis', fontsize=14, fontweight='bold')

    # Norm deviation from 1.0
    ax = axes[0]
    norm_dev_1 = np.abs(N1_data['norm'] - 1.0)
    norm_dev_10 = np.abs(N10_data['norm'] - 1.0)
    norm_dev_100 = np.abs(N100_data['norm'] - 1.0)

    ax.semilogy(N1_data['time'], norm_dev_1, 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.semilogy(N10_data['time'], norm_dev_10, 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.semilogy(N100_data['time'], norm_dev_100, 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel('|Norm - 1.0| (log scale)', fontsize=11)
    ax.set_title('Norm Conservation Error', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Energy stability (variance)
    ax = axes[1]
    window = 100  # Rolling window for variance

    energy_std_1 = N1_data['energy'].rolling(window=window).std()
    energy_std_10 = N10_data['energy'].rolling(window=window).std()
    energy_std_100 = N100_data['energy'].rolling(window=window).std()

    ax.plot(N1_data['time'], energy_std_1, 'b-', alpha=0.7, linewidth=1.5, label='N=1')
    ax.plot(N10_data['time'], energy_std_10, 'orange', alpha=0.7, linewidth=1.5, label='N=10')
    ax.plot(N100_data['time'], energy_std_100, 'g-', alpha=0.7, linewidth=1.5, label='N=100')
    ax.set_xlabel('Time', fontsize=11)
    ax.set_ylabel(f'Energy Std Dev (window={window})', fontsize=11)
    ax.set_title('Energy Stability', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/operator_splitting_convergence.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_dir}/operator_splitting_convergence.png")

    plt.close()

def compute_statistics(N1_data, N10_data, N100_data):
    """Compute and print summary statistics"""

    print("\n" + "="*70)
    print("OPERATOR SPLITTING VALIDATION STATISTICS")
    print("="*70)

    datasets = [
        ("N=1", N1_data),
        ("N=10", N10_data),
        ("N=100", N100_data)
    ]

    for name, data in datasets:
        print(f"\n{name} Statistics:")
        print(f"  Final time: {data['time'].iloc[-1]:.2f}")
        print(f"  R_avg:  mean={data['R_avg'].mean():.6f}, std={data['R_avg'].std():.6e}")
        print(f"  Norm:   final={data['norm'].iloc[-1]:.6f}, drift={(1.0 - data['norm'].iloc[-1]):.6e}")
        print(f"  Energy: mean={data['energy'].mean():.6f}, std={data['energy'].std():.6e}")
        print(f"  Max density: final={data['max_density'].iloc[-1]:.6e}")

    print("\n" + "="*70)
    print("CONVERGENCE ANALYSIS")
    print("="*70)

    # Compare N=10 vs N=1
    print("\nN=10 vs N=1 differences (final values):")
    print(f"  ΔR_avg:  {abs(N10_data['R_avg'].iloc[-1] - N1_data['R_avg'].iloc[-1]):.6e}")
    print(f"  ΔNorm:   {abs(N10_data['norm'].iloc[-1] - N1_data['norm'].iloc[-1]):.6e}")
    print(f"  ΔEnergy: {abs(N10_data['energy'].iloc[-1] - N1_data['energy'].iloc[-1]):.6e}")

    # Compare N=100 vs N=10
    print("\nN=100 vs N=10 differences (final values):")
    print(f"  ΔR_avg:  {abs(N100_data['R_avg'].iloc[-1] - N10_data['R_avg'].iloc[-1]):.6e}")
    print(f"  ΔNorm:   {abs(N100_data['norm'].iloc[-1] - N10_data['norm'].iloc[-1]):.6e}")
    print(f"  ΔEnergy: {abs(N100_data['energy'].iloc[-1] - N10_data['energy'].iloc[-1]):.6e}")

    print("\n" + "="*70)

def main():
    """Main visualization routine"""

    print("\n" + "="*70)
    print("OPERATOR SPLITTING 10K TIMESTEP VISUALIZATION")
    print("="*70 + "\n")

    # Load data
    print("Loading data...")
    try:
        N1_data = load_timeseries("N1_10k")
        print(f"  ✓ N=1 data loaded: {len(N1_data)} timesteps")

        N10_data = load_timeseries("N10_10k")
        print(f"  ✓ N=10 data loaded: {len(N10_data)} timesteps")

        N100_data = load_timeseries("N100_10k")
        print(f"  ✓ N=100 data loaded: {len(N100_data)} timesteps")
    except FileNotFoundError as e:
        print(f"\n✗ ERROR: {e}")
        print("\nPlease ensure the simulation has completed and data files exist:")
        print("  - N1_10k/timeseries.csv")
        print("  - N10_10k/timeseries.csv")
        print("  - N100_10k/timeseries.csv")
        return

    # Create output directory for plots
    output_dir = Path(".")
    output_dir.mkdir(exist_ok=True)

    # Generate plots
    print("\nGenerating comparison plots...")
    plot_comparison(N1_data, N10_data, N100_data, output_dir)

    print("\nGenerating convergence analysis...")
    plot_convergence_analysis(N1_data, N10_data, N100_data, output_dir)

    # Compute and display statistics
    compute_statistics(N1_data, N10_data, N100_data)

    print("\n✓ Visualization complete!\n")

if __name__ == "__main__":
    main()
