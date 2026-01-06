#!/usr/bin/env python3
"""
E3 Causality Analysis - Data Visualization and Verification
Analyzes dispersion relation and group velocity to verify v_g ≤ c
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Constants (natural units: c = ħ = 1)
SPEED_OF_LIGHT = 1.0
MASS_GAP = 1.0

def compute_dispersion(k):
    """Compute ω from dispersion relation ω² = k² + Δ²"""
    return np.sqrt(k**2 + MASS_GAP**2)

def compute_group_velocity(k):
    """Compute group velocity v_g = dω/dk = k/ω"""
    omega = compute_dispersion(k)
    return k / omega

def compute_phase_velocity(k):
    """Compute phase velocity v_p = ω/k"""
    omega = compute_dispersion(k)
    return omega / k

def verify_causality():
    """Verify that v_g ≤ c for all wavenumbers"""
    print("=" * 70)
    print("E3 CAUSALITY VERIFICATION - THEORETICAL ANALYSIS")
    print("=" * 70)
    print()

    # Test wavenumbers
    k_values = np.logspace(-2, 2, 100)

    v_group = compute_group_velocity(k_values)
    v_phase = compute_phase_velocity(k_values)

    # Check causality
    max_v_group = np.max(v_group)
    violations = np.sum(v_group > SPEED_OF_LIGHT)

    print(f"Wavenumber range: k ∈ [{k_values[0]:.3f}, {k_values[-1]:.3f}]")
    print(f"Mass gap: Δ = {MASS_GAP}")
    print(f"Speed of light: c = {SPEED_OF_LIGHT}")
    print()

    print("Group Velocity Analysis:")
    print(f"  Maximum v_g = {max_v_group:.6f} c")
    print(f"  Minimum v_g = {np.min(v_group):.6f} c")
    print(f"  Violations (v_g > c): {violations}")
    print()

    print("Phase Velocity Analysis:")
    print(f"  Maximum v_p = {np.max(v_phase):.3f} c")
    print(f"  Minimum v_p = {np.min(v_phase):.3f} c")
    print(f"  (Phase velocity can exceed c - no causality violation)")
    print()

    # Limiting behavior
    k_small = 0.001
    k_large = 1000.0
    v_g_small = compute_group_velocity(k_small)
    v_g_large = compute_group_velocity(k_large)

    print("Limiting Behavior:")
    print(f"  k → 0: v_g → {v_g_small:.6f} (should → 0)")
    print(f"  k → ∞: v_g → {v_g_large:.6f} (should → 1)")
    print()

    # Verdict
    print("=" * 70)
    if violations == 0 and max_v_group <= SPEED_OF_LIGHT:
        print("✅ VERDICT: GO - TRD THEORY IS CAUSAL")
        print()
        print("All group velocities satisfy v_g ≤ c")
        print("Light cone structure preserved")
        print("Theory compatible with special relativity")
        result = True
    else:
        print("❌ VERDICT: NO-GO - TRD VIOLATES CAUSALITY")
        print()
        print(f"Found {violations} superluminal modes")
        print("Theory incompatible with special relativity")
        result = False
    print("=" * 70)

    return result, k_values, v_group, v_phase

def create_plots(k_values, v_group, v_phase):
    """Create visualization plots for causality analysis"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Dispersion relation
    ax = axes[0, 0]
    omega = compute_dispersion(k_values)
    ax.plot(k_values, omega, 'b-', linewidth=2, label='ω(k)')
    ax.plot(k_values, k_values, 'r--', linewidth=1, label='ω = k (massless)')
    ax.plot(k_values, np.sqrt(k_values**2 + MASS_GAP**2), 'g:',
            linewidth=2, label=f'√(k² + Δ²), Δ={MASS_GAP}')
    ax.set_xlabel('Wavenumber k', fontsize=12)
    ax.set_ylabel('Frequency ω', fontsize=12)
    ax.set_title('Dispersion Relation: ω² = k² + Δ²', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 2: Group velocity
    ax = axes[0, 1]
    ax.plot(k_values, v_group, 'b-', linewidth=2, label='v_g(k)')
    ax.axhline(y=SPEED_OF_LIGHT, color='r', linestyle='--',
               linewidth=2, label='Speed of light (c)')
    ax.fill_between(k_values, 0, SPEED_OF_LIGHT, alpha=0.2, color='green',
                     label='Causal region (v ≤ c)')
    ax.set_xlabel('Wavenumber k', fontsize=12)
    ax.set_ylabel('Group Velocity v_g', fontsize=12)
    ax.set_title('Group Velocity: v_g = k/√(k² + Δ²)', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 3: Phase velocity
    ax = axes[1, 0]
    ax.plot(k_values, v_phase, 'g-', linewidth=2, label='v_p(k)')
    ax.axhline(y=SPEED_OF_LIGHT, color='r', linestyle='--',
               linewidth=2, label='Speed of light (c)')
    ax.set_xlabel('Wavenumber k', fontsize=12)
    ax.set_ylabel('Phase Velocity v_p', fontsize=12)
    ax.set_title('Phase Velocity: v_p = ω/k (Can Exceed c)', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 4: Comparison
    ax = axes[1, 1]
    ax.plot(k_values, v_group, 'b-', linewidth=2, label='Group velocity v_g')
    ax.plot(k_values, v_phase, 'g-', linewidth=2, label='Phase velocity v_p')
    ax.axhline(y=SPEED_OF_LIGHT, color='r', linestyle='--',
               linewidth=2, label='Speed of light (c)')
    ax.fill_between(k_values, 0, SPEED_OF_LIGHT, alpha=0.2, color='green')
    ax.set_xlabel('Wavenumber k', fontsize=12)
    ax.set_ylabel('Velocity', fontsize=12)
    ax.set_title('Causality Check: v_g ≤ c (v_p can exceed c)',
                 fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    # Save figure
    output_dir = Path('output/causality')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'causality_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nPlots saved to {output_dir / 'causality_analysis.png'}")

    return fig

def save_data(k_values, v_group, v_phase):
    """Save numerical data to CSV"""
    output_dir = Path('output/causality')
    output_dir.mkdir(parents=True, exist_ok=True)

    omega = compute_dispersion(k_values)

    data = np.column_stack([k_values, omega, v_group, v_phase])
    header = 'k,omega,v_group,v_phase'

    filepath = output_dir / 'causality_theoretical.csv'
    np.savetxt(filepath, data, delimiter=',', header=header, comments='')
    print(f"Data saved to {filepath}")

def main():
    """Main analysis function"""

    # Run causality verification
    passed, k_values, v_group, v_phase = verify_causality()

    # Create visualization plots
    print("\nGenerating plots...")
    create_plots(k_values, v_group, v_phase)

    # Save numerical data
    save_data(k_values, v_group, v_phase)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print()

    if passed:
        print("Result: TRD is CAUSAL ✅")
        print()
        print("Key findings:")
        print("  • All wave modes propagate at v_g ≤ c")
        print("  • Mass gap ensures subluminal group velocity")
        print("  • Phase velocity can exceed c (no causality issue)")
        print("  • Theory compatible with special relativity")
        return 0
    else:
        print("Result: TRD violates CAUSALITY ❌")
        print()
        print("Theory must be revised!")
        return 1

if __name__ == '__main__':
    import sys
    sys.exit(main())
