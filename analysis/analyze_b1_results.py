#!/usr/bin/env python3
"""
Analyze B1 particle mass spectrum optimization results.
Generates plots and summary statistics from parameter scans.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

# Configure plotting style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Load results
df = pd.read_csv('b1_optimization_results.csv')

# Constants
TARGET_RATIO = 206.768  # Muon/electron mass ratio
BASELINE_RATIO = 6.45   # Initial test result

# Find best configuration
best_idx = df['m2_m1'].idxmax()
best_result = df.iloc[best_idx]

print("="*60)
print("B1 PARTICLE MASS RATIO OPTIMIZATION RESULTS")
print("="*60)

print("\n--- OPTIMAL CONFIGURATION ---")
print(f"K (coupling strength): {best_result['K']:.2f}")
print(f"Δ (mass gap):          {best_result['Delta']:.2f}")
print(f"d (vortex separation): {best_result['separation']:.2f}")
print(f"\nResulting mass ratio m₂/m₁: {best_result['m2_m1']:.3f}")
print(f"Baseline ratio:              {BASELINE_RATIO:.3f}")
print(f"Target ratio (μ/e):          {TARGET_RATIO:.3f}")
print(f"\nImprovement over baseline: {best_result['m2_m1']/BASELINE_RATIO:.2f}×")
print(f"Shortfall from target:     {TARGET_RATIO/best_result['m2_m1']:.1f}×")
print(f"Error from target:         {100*(1-best_result['m2_m1']/TARGET_RATIO):.1f}%")

# Create figure with subplots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('B1 Particle Mass Ratio Optimization Analysis', fontsize=16)

# 1. K-scan (fixed Delta=1, sep=10)
ax = axes[0, 0]
k_scan = df[(df['Delta'] == 1.0) & (df['separation'] == 10.0)].sort_values('K')
ax.plot(k_scan['K'], k_scan['m2_m1'], 'o-', linewidth=2, markersize=8)
ax.axhline(BASELINE_RATIO, color='gray', linestyle='--', label='Baseline (6.45)')
ax.axhline(TARGET_RATIO/10, color='green', linestyle=':', label='Target/10 (20.7)')
ax.set_xlabel('Coupling Strength K')
ax.set_ylabel('Mass Ratio m₂/m₁')
ax.set_title('Phase 1: Coupling Strength Scan')
ax.set_xscale('log')
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Delta-scan (fixed K=10, sep=10)
ax = axes[0, 1]
delta_scan = df[(df['K'] == 10.0) & (df['separation'] == 10.0)].sort_values('Delta')
if len(delta_scan) > 0:
    ax.plot(delta_scan['Delta'], delta_scan['m2_m1'], 's-', linewidth=2, markersize=8)
    ax.axhline(BASELINE_RATIO, color='gray', linestyle='--', label='Baseline')
    ax.set_xlabel('Mass Gap Δ')
    ax.set_ylabel('Mass Ratio m₂/m₁')
    ax.set_title('Phase 2: Mass Gap Scan')
    ax.set_xscale('log')
    ax.legend()
ax.grid(True, alpha=0.3)

# 3. Separation scan (fixed K=10, Delta=5)
ax = axes[0, 2]
sep_scan = df[(df['K'] == 10.0) & (df['Delta'] == 5.0)].sort_values('separation')
if len(sep_scan) > 0:
    ax.plot(sep_scan['separation'], sep_scan['m2_m1'], '^-', linewidth=2, markersize=8)
    ax.axhline(BASELINE_RATIO, color='gray', linestyle='--', label='Baseline')
    ax.axhline(TARGET_RATIO/10, color='green', linestyle=':', label='Target/10')
    ax.set_xlabel('Vortex Separation d')
    ax.set_ylabel('Mass Ratio m₂/m₁')
    ax.set_title('Phase 3: Vortex Separation Scan')
    ax.legend()
ax.grid(True, alpha=0.3)

# 4. 2D heatmap: K vs separation
ax = axes[1, 0]
pivot_data = df[df['Delta'] == 5.0].pivot_table(values='m2_m1', index='K', columns='separation')
if not pivot_data.empty:
    im = ax.imshow(pivot_data, cmap='viridis', aspect='auto', origin='lower')
    ax.set_xlabel('Vortex Separation')
    ax.set_ylabel('Coupling Strength K')
    ax.set_title('2D Grid: m₂/m₁(K, separation)')
    plt.colorbar(im, ax=ax, label='m₂/m₁')

# 5. Mass spectrum for optimal configuration
ax = axes[1, 1]
masses = [best_result['m1'], best_result['m2'], best_result['m3']]
charges = [1, 2, 3]
ax.bar(charges, masses, color=['blue', 'green', 'red'], alpha=0.7)
ax.set_xlabel('Topological Charge Q')
ax.set_ylabel('Effective Mass')
ax.set_title(f'Mass Spectrum (K={best_result["K"]:.1f}, sep={best_result["separation"]:.1f})')
ax.set_xticks(charges)
ax.grid(True, alpha=0.3)

# 6. Scaling analysis: log-log plot
ax = axes[1, 2]
ax.loglog([1, 2, 3], masses, 'o-', linewidth=2, markersize=10)
ax.set_xlabel('Topological Charge Q')
ax.set_ylabel('Mass')
ax.set_title('Topological Scaling')
ax.grid(True, which="both", alpha=0.3)

# Add text annotations with key findings
textstr = f'Best m₂/m₁ = {best_result["m2_m1"]:.2f}\n'
textstr += f'Improvement: {best_result["m2_m1"]/BASELINE_RATIO:.1f}×\n'
textstr += f'Gap to target: {TARGET_RATIO/best_result["m2_m1"]:.0f}×'
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('b1_optimization_analysis.png', dpi=150)
plt.savefig('b1_optimization_analysis.pdf')
plt.show()

# Additional analysis: Parameter correlations
print("\n--- PARAMETER CORRELATIONS ---")
corr_matrix = df[['K', 'Delta', 'separation', 'm2_m1', 'R_std', 'grad_mag']].corr()
print("\nCorrelation with m₂/m₁:")
print(corr_matrix['m2_m1'].sort_values(ascending=False))

# Find configurations exceeding quality gates
print("\n--- QUALITY GATE ANALYSIS ---")
exceeds_baseline = df[df['m2_m1'] > BASELINE_RATIO]
within_factor_10 = df[df['m2_m1'] > TARGET_RATIO/10]
within_factor_5 = df[df['m2_m1'] > TARGET_RATIO/5]

print(f"Configurations exceeding baseline ({BASELINE_RATIO}): {len(exceeds_baseline)}/{len(df)}")
print(f"Within factor 10 of target (>{TARGET_RATIO/10:.1f}): {len(within_factor_10)}/{len(df)}")
print(f"Within factor 5 of target (>{TARGET_RATIO/5:.1f}): {len(within_factor_5)}/{len(df)}")

# Summary recommendations
print("\n--- RECOMMENDATIONS FOR PHASE 5 ---")
if best_result['m2_m1'] > TARGET_RATIO/5:
    print("✓ MAJOR SUCCESS: Achieved >5× improvement!")
    print("  → Continue with current approach")
    print("  → Fine-tune separation parameter around optimal value")
elif best_result['m2_m1'] > TARGET_RATIO/10:
    print("✓ Good progress: Within order of magnitude")
    print("  → Explore radial mode contributions")
    print("  → Test higher-order topological configurations")
    print("  → Consider beyond-mean-field corrections")
else:
    print("⚠ Further optimization needed:")
    print("  → Extend parameter ranges (K > 10, larger separations)")
    print("  → Investigate asymmetric vortex configurations")
    print("  → Consider multi-scale vortex structures")

print("\n" + "="*60)
print("Analysis complete. Plots saved to:")
print("  - b1_optimization_analysis.png")
print("  - b1_optimization_analysis.pdf")