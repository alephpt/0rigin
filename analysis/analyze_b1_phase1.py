#!/usr/bin/env python3
"""
B1 Phase 1 Analysis: K-Parameter and Vortex Separation Study

Analyzes the relationship between:
1. Kuramoto coupling strength K and mass hierarchy
2. Vortex separation d and interaction energy

Goal: Identify why m₂/m₁ = 3.65 instead of target 206.768
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('b1_phase1_results.csv')

# Separate K-scan and separation-scan data
k_scan = df[df['d'] == 4.0].copy()
d_scan = df[df['K'] == 1.0].copy()

# Create figure with 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('B1 Phase 1: K-Parameter and Vortex Separation Analysis', fontsize=16, fontweight='bold')

# --- Plot 1: m₂/m₁ vs K ---
ax1 = axes[0, 0]
ax1.plot(k_scan['K'], k_scan['m2_m1'], 'o-', linewidth=2, markersize=8, label='Measured')
ax1.axhline(y=206.768, color='r', linestyle='--', linewidth=2, label='Target (μ/e)')
ax1.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax1.set_ylabel('Mass Ratio m₂/m₁', fontsize=12)
ax1.set_title('Mass Hierarchy vs Coupling Strength', fontsize=13, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# Annotate current vs target
current_ratio = k_scan.loc[k_scan['K'] == 1.0, 'm2_m1'].values[0]
ax1.text(0.05, 0.95, f'Current (K=1.0): {current_ratio:.2f}\nTarget: 206.768\nGap: {206.768/current_ratio:.1f}×',
         transform=ax1.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# --- Plot 2: Energy Scaling with K ---
ax2 = axes[0, 1]
ax2.plot(k_scan['K'], k_scan['E1'], 'o-', linewidth=2, markersize=8, label='E₁ (Q=1)')
ax2.plot(k_scan['K'], k_scan['E2'], 's-', linewidth=2, markersize=8, label='E₂ (Q=2)')
ax2.plot(k_scan['K'], k_scan['E3'], '^-', linewidth=2, markersize=8, label='E₃ (Q=3)')
ax2.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax2.set_ylabel('Vortex Energy', fontsize=12)
ax2.set_title('Energy vs Coupling Strength', fontsize=13, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

# --- Plot 3: m₂/m₁ vs Separation ---
ax3 = axes[1, 0]
ax3.plot(d_scan['d'], d_scan['m2_m1'], 'o-', linewidth=2, markersize=8, color='green', label='Measured')
ax3.axhline(y=2.0, color='b', linestyle='--', linewidth=2, label='Independent limit (E₂ = 2·E₁)')
ax3.axhline(y=206.768, color='r', linestyle='--', linewidth=2, label='Target (μ/e)')
ax3.set_xlabel('Vortex Separation d (grid units)', fontsize=12)
ax3.set_ylabel('Mass Ratio m₂/m₁', fontsize=12)
ax3.set_title('Mass Ratio vs Vortex Separation', fontsize=13, fontweight='bold')
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.set_yscale('log')

# --- Plot 4: Interaction Energy vs Separation ---
ax4 = axes[1, 1]
ax4.plot(d_scan['d'], d_scan['E_interaction'], 'o-', linewidth=2, markersize=8, color='purple')
ax4.axhline(y=0, color='k', linestyle='-', linewidth=1, alpha=0.3)
ax4.set_xlabel('Vortex Separation d (grid units)', fontsize=12)
ax4.set_ylabel('Interaction Energy E_int = E₂ - 2·E₁', fontsize=12)
ax4.set_title('Vortex Interaction Energy vs Separation', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)

# Annotate interaction type
ax4.text(0.05, 0.95, 'E_int > 0: Repulsive\nE_int decreases with d',
         transform=ax4.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

plt.tight_layout()
plt.savefig('b1_phase1_analysis.png', dpi=300, bbox_inches='tight')
print("✓ Visualization saved to: b1_phase1_analysis.png")

# Print summary statistics
print("\n" + "="*60)
print("B1 PHASE 1 ANALYSIS SUMMARY")
print("="*60)

print("\n1. K-PARAMETER DEPENDENCE (fixed d=4.0):")
print("-" * 60)
print(k_scan[['K', 'E1', 'E2', 'm2_m1', 'E_interaction']].to_string(index=False))

print("\n   Key Finding:")
print(f"   - m₂/m₁ varies from {k_scan['m2_m1'].max():.2f} to {k_scan['m2_m1'].min():.2f}")
print(f"   - Relative variation: {(k_scan['m2_m1'].max() - k_scan['m2_m1'].min()) / k_scan['m2_m1'].mean() * 100:.1f}%")
print(f"   - Conclusion: K affects overall scale, NOT hierarchy structure")

print("\n2. VORTEX SEPARATION DEPENDENCE (fixed K=1.0):")
print("-" * 60)
print(d_scan[['d', 'E2', 'm2_m1', 'E_interaction']].to_string(index=False))

print("\n   Key Finding:")
print(f"   - As d increases: m₂/m₁ decreases from {d_scan['m2_m1'].max():.2f} to {d_scan['m2_m1'].min():.2f}")
print(f"   - At d=16: m₂/m₁ = {d_scan.loc[d_scan['d'] == 16.0, 'm2_m1'].values[0]:.2f}")
print(f"   - E₂/(2·E₁) at d=16: {d_scan.loc[d_scan['d'] == 16.0, 'E2'].values[0] / (2 * d_scan['E1'].values[0]):.2f}")
print(f"   - Interaction energy decreases but remains positive (repulsive)")
print(f"   - Conclusion: Vortices repel, preventing strong binding")

print("\n3. CRITICAL PHYSICS INSIGHTS:")
print("-" * 60)
print("   ✗ Current m₂/m₁ ≈ 3.65 (target: 206.768, gap: 56.6×)")
print("   ✗ E₂/E₁ ≈ 3.65 instead of expected ≈ 2 (linear scaling)")
print("   ✗ Positive interaction energy → repulsive vortices")
print("   ✗ No binding mechanism for mass hierarchy")
print("\n   Missing Physics Identified:")
print("   1. Radial excitation modes (n,l,m quantum numbers)")
print("   2. R-field self-consistency (dynamic synchronization)")
print("   3. Topological stability constraints")
print("   4. Non-linear coupling effects")

print("\n4. NEXT STEPS (Phase 2):")
print("-" * 60)
print("   → Add radial modes to vortex configurations")
print("   → Include R-field feedback in energy functional")
print("   → Test alternative vortex geometries (knots, links)")
print("   → Investigate binding mechanisms for mass gap")

print("\n" + "="*60)
print("Analysis complete. See b1_phase1_analysis.png for plots.")
print("="*60)
