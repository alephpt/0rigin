#!/usr/bin/env python3
"""
Visualize Phase 1 Full Validation Results

Generates publication-quality plots for:
- Long-time norm conservation
- Energy conservation
- Operator splitting convergence (N=1, 10, 100)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Set publication quality
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 1.5

def load_observables(N):
    """Load observables CSV for given N ratio"""
    path = Path(f"phase1_full_validation/N_{N}/observables.csv")
    if not path.exists():
        print(f"Warning: {path} not found")
        return None
    return pd.read_csv(path)

def plot_long_time_stability():
    """Plot norm and energy evolution over 50k steps"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    for N in [1, 10, 100]:
        df = load_observables(N)
        if df is None:
            continue

        # Norm conservation
        ax1.plot(df['time'], df['norm'], label=f'N={N}', alpha=0.8)

        # Energy conservation
        E0 = df['E_total'].iloc[0]
        energy_drift = (df['E_total'] - E0) / E0
        ax2.plot(df['time'], energy_drift, label=f'N={N}', alpha=0.8)

    # Norm plot
    ax1.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label='Target')
    ax1.set_ylabel('Norm ||Ψ||²')
    ax1.set_title('Phase 1: Long-Time Norm Conservation (50k steps)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Energy plot
    ax2.axhline(y=0.0, color='k', linestyle='--', alpha=0.3)
    ax2.set_xlabel('Time (natural units)')
    ax2.set_ylabel('Relative Energy Drift ΔE/E₀')
    ax2.set_title('Energy Conservation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('phase1_long_time_stability.png')
    print("✓ Saved: phase1_long_time_stability.png")
    plt.close()

def plot_convergence():
    """Compare final observables across N ratios"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Load final observables for each N
    final_obs = {}
    for N in [1, 10, 100]:
        df = load_observables(N)
        if df is not None:
            final_obs[N] = df.iloc[-1]

    if len(final_obs) < 2:
        print("Not enough data for convergence plot")
        return

    N_values = list(final_obs.keys())

    # Plot 1: Norm convergence
    ax = axes[0, 0]
    norms = [final_obs[N]['norm'] for N in N_values]
    ax.plot(N_values, norms, 'o-', markersize=8)
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Substep Ratio N')
    ax.set_ylabel('Final Norm')
    ax.set_title('Norm Convergence')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)

    # Plot 2: Energy convergence
    ax = axes[0, 1]
    energies = [final_obs[N]['E_total'] for N in N_values]
    ax.plot(N_values, energies, 'o-', markersize=8, color='orange')
    ax.set_xlabel('Substep Ratio N')
    ax.set_ylabel('Final Energy')
    ax.set_title('Energy Convergence')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)

    # Plot 3: Position expectation
    ax = axes[1, 0]
    pos_x = [final_obs[N]['pos_x_re'] for N in N_values]
    pos_y = [final_obs[N]['pos_y_re'] for N in N_values]
    ax.plot(N_values, pos_x, 'o-', label='<x>', markersize=8)
    ax.plot(N_values, pos_y, 's-', label='<y>', markersize=8)
    ax.set_xlabel('Substep Ratio N')
    ax.set_ylabel('Position (grid units)')
    ax.set_title('Position Expectation Values')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Convergence error
    ax = axes[1, 1]
    if len(N_values) >= 2:
        # Compare to highest N (most accurate)
        ref_N = max(N_values)
        errors = []
        compare_N = []
        for N in N_values:
            if N != ref_N:
                norm_err = abs(final_obs[N]['norm'] - final_obs[ref_N]['norm'])
                errors.append(norm_err)
                compare_N.append(N)

        if len(errors) > 0 and any(e > 0 for e in errors):
            ax.semilogy(compare_N, errors, 'o-', markersize=8, color='red')
            ax.set_xscale('log')
        else:
            ax.plot(compare_N, errors, 'o-', markersize=8, color='red')
        ax.set_xlabel('Substep Ratio N')
        ax.set_ylabel('|Norm(N) - Norm(100)|')
        ax.set_title(f'Convergence Error (vs N={ref_N})')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('phase1_convergence.png')
    print("✓ Saved: phase1_convergence.png")
    plt.close()

def plot_wavepacket_spreading():
    """Plot wavepacket position and momentum evolution"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for N in [1, 10, 100]:
        df = load_observables(N)
        if df is None:
            continue

        # Position trajectory
        ax1.plot(df['pos_x_re'], df['pos_y_re'], label=f'N={N}', alpha=0.7)
        ax1.plot(df['pos_x_re'].iloc[0], df['pos_y_re'].iloc[0], 'o', markersize=10)

        # Momentum evolution
        p_mag = np.sqrt(df['mom_x_re']**2 + df['mom_y_re']**2)
        ax2.plot(df['time'], p_mag, label=f'N={N}', alpha=0.7)

    ax1.set_xlabel('Position <x> (grid units)')
    ax1.set_ylabel('Position <y> (grid units)')
    ax1.set_title('Wavepacket Trajectory')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')

    ax2.set_xlabel('Time (natural units)')
    ax2.set_ylabel('Momentum Magnitude |<p>|')
    ax2.set_title('Momentum Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('phase1_wavepacket.png')
    print("✓ Saved: phase1_wavepacket.png")
    plt.close()

def generate_summary():
    """Generate text summary of results"""
    print("\n" + "="*60)
    print("PHASE 1 FULL VALIDATION SUMMARY")
    print("="*60)

    for N in [1, 10, 100]:
        df = load_observables(N)
        if df is None:
            continue

        print(f"\n--- N={N} ---")
        print(f"Total steps: {len(df)}")
        print(f"Time range: {df['time'].iloc[0]:.2f} - {df['time'].iloc[-1]:.2f}")

        # Norm statistics
        norm_initial = df['norm'].iloc[0]
        norm_final = df['norm'].iloc[-1]
        norm_drift = abs(norm_final - 1.0)
        print(f"Norm: {norm_initial:.6f} → {norm_final:.6f} (drift: {norm_drift*100:.3f}%)")

        # Energy statistics
        E0 = df['E_total'].iloc[0]
        E_final = df['E_total'].iloc[-1]
        E_drift = abs(E_final - E0) / E0
        print(f"Energy: {E0:.6f} → {E_final:.6f} (drift: {E_drift*100:.3f}%)")

        # Stability check
        has_nan = df['norm'].isna().any() or df['E_total'].isna().any()
        stability = "✓ STABLE" if not has_nan else "✗ NaN DETECTED"
        print(f"Stability: {stability}")

    print("\n" + "="*60)

if __name__ == "__main__":
    print("Phase 1 Full Validation Visualization")
    print("="*60)

    # Generate all plots
    plot_long_time_stability()
    plot_convergence()
    plot_wavepacket_spreading()

    # Print summary
    generate_summary()

    print("\n✓ All visualizations complete!")
    print("  - phase1_long_time_stability.png")
    print("  - phase1_convergence.png")
    print("  - phase1_wavepacket.png")
