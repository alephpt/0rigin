#!/usr/bin/env python3
"""
Visualize timesync validation results for Phase 1 sign-off
Generates required plots for N=1, 10, 100 convergence verification
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Output directory
output_dir = Path(__file__).parent

# Load data for N=1, 10, 100
data = {}
for N in [1, 10, 100]:
    csv_path = output_dir / f"N_{N}" / "observables.csv"
    data[N] = pd.read_csv(csv_path)

# Create figure with 6 subplots
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# Colors for N=1, 10, 100
colors = {1: 'blue', 10: 'green', 100: 'red'}
labels = {1: 'N=1', 10: 'N=10', 100: 'N=100'}

# Plot 1: Norm conservation
ax1 = fig.add_subplot(gs[0, 0])
for N in [1, 10, 100]:
    ax1.plot(data[N]['time'], data[N]['norm_error'],
             color=colors[N], label=labels[N], linewidth=1.5)
ax1.set_xlabel('Time')
ax1.set_ylabel('Norm Error')
ax1.set_title('Norm Conservation (Phase 1.1)')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0.01, color='k', linestyle='--', alpha=0.3, label='Threshold')

# Plot 2: Energy conservation
ax2 = fig.add_subplot(gs[0, 1])
for N in [1, 10, 100]:
    E_drift = np.abs(data[N]['E_total'] - data[N]['E_total'].iloc[0])
    ax2.plot(data[N]['time'], E_drift,
             color=colors[N], label=labels[N], linewidth=1.5)
ax2.set_xlabel('Time')
ax2.set_ylabel('|E(t) - E(0)|')
ax2.set_title('Energy Drift (Phase 1.1)')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0.01, color='k', linestyle='--', alpha=0.3, label='Threshold')

# Plot 3: Trajectory x(t)
ax3 = fig.add_subplot(gs[1, 0])
for N in [1, 10, 100]:
    ax3.plot(data[N]['time'], data[N]['pos_x_re'],
             color=colors[N], label=labels[N], linewidth=1.5)
ax3.set_xlabel('Time')
ax3.set_ylabel('x_CoM (grid points)')
ax3.set_title('Particle Trajectory X (Phase 1.4)')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Trajectory y(t)
ax4 = fig.add_subplot(gs[1, 1])
for N in [1, 10, 100]:
    ax4.plot(data[N]['time'], data[N]['pos_y_re'],
             color=colors[N], label=labels[N], linewidth=1.5)
ax4.set_xlabel('Time')
ax4.set_ylabel('y_CoM (grid points)')
ax4.set_title('Particle Trajectory Y (Phase 1.4)')
ax4.legend()
ax4.grid(True, alpha=0.3)

# Plot 5: Convergence verification (final values)
ax5 = fig.add_subplot(gs[2, 0])
final_energies = [data[N]['E_total'].iloc[-1] for N in [1, 10, 100]]
E_ref = final_energies[2]  # N=100 as reference
errors = [abs(E - E_ref) for E in final_energies]
N_values = [1, 10, 100]
ax5.semilogy(N_values, errors, 'bo-', markersize=8, linewidth=2)
ax5.set_xlabel('Substep Ratio N')
ax5.set_ylabel('|E(N) - E(100)|')
ax5.set_title('Convergence: Error vs N (Phase 1.2)')
ax5.grid(True, alpha=0.3)
ax5.set_xticks([1, 10, 100])

# Add convergence test annotation
monotonic = errors[0] > errors[1]
status = "✓ PASS" if monotonic else "✗ FAIL"
ax5.text(0.05, 0.95, f'Monotonic convergence: {status}\nError(N=1) = {errors[0]:.2e}\nError(N=10) = {errors[1]:.2e}',
         transform=ax5.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 6: Momentum N-dependence
ax6 = fig.add_subplot(gs[2, 1])
for N in [1, 10, 100]:
    mom_mag = np.sqrt(data[N]['mom_x_re']**2 + data[N]['mom_y_re']**2)
    ax6.plot(data[N]['time'], mom_mag,
             color=colors[N], label=labels[N], linewidth=1.5)
ax6.set_xlabel('Time')
ax6.set_ylabel('|p| (momentum magnitude)')
ax6.set_title('Momentum Evolution (Phase 1.3 - Coupling)')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.suptitle('Phase 1 Validation: Operator Splitting Convergence (N=1, 10, 100)',
             fontsize=14, fontweight='bold')
plt.savefig(output_dir / 'phase1_validation_complete.png', dpi=150, bbox_inches='tight')
print(f"[Phase 1] Validation plots saved to {output_dir / 'phase1_validation_complete.png'}")

# Generate trajectory comparison plot
fig2, (ax_x, ax_y) = plt.subplots(1, 2, figsize=(14, 5))

# X trajectory comparison
for N in [1, 10, 100]:
    ax_x.plot(data[N]['time'], data[N]['pos_x_re'],
              color=colors[N], label=labels[N], linewidth=2)
ax_x.set_xlabel('Time', fontsize=12)
ax_x.set_ylabel('x_CoM (grid points)', fontsize=12)
ax_x.set_title('Particle X-Position: N-Dependence', fontsize=13)
ax_x.legend(fontsize=11)
ax_x.grid(True, alpha=0.3)

# Y trajectory comparison
for N in [1, 10, 100]:
    ax_y.plot(data[N]['time'], data[N]['pos_y_re'],
              color=colors[N], label=labels[N], linewidth=2)
ax_y.set_xlabel('Time', fontsize=12)
ax_y.set_ylabel('y_CoM (grid points)', fontsize=12)
ax_y.set_title('Particle Y-Position: N-Dependence', fontsize=13)
ax_y.legend(fontsize=11)
ax_y.grid(True, alpha=0.3)

# Calculate trajectory difference
x_diff_1_100 = np.abs(data[1]['pos_x_re'] - data[100]['pos_x_re'])
y_diff_1_100 = np.abs(data[1]['pos_y_re'] - data[100]['pos_y_re'])
max_diff = np.sqrt(x_diff_1_100**2 + y_diff_1_100**2).max()

ax_x.text(0.05, 0.95, f'Max |Δr(N=1, N=100)| = {max_diff:.4f} grid points',
          transform=ax_x.transAxes, fontsize=10, verticalalignment='top',
          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

plt.suptitle('Phase 1.4: Trajectory Comparison (Operator Splitting Effect)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'phase1_trajectory_comparison.png', dpi=150, bbox_inches='tight')
print(f"[Phase 1] Trajectory comparison saved to {output_dir / 'phase1_trajectory_comparison.png'}")

# Print summary
print("\n" + "="*60)
print("PHASE 1 VALIDATION SUMMARY")
print("="*60)
print(f"Final energies:")
for N in [1, 10, 100]:
    E_final = data[N]['E_total'].iloc[-1]
    print(f"  N={N:3d}: E = {E_final:.6f}")

print(f"\nConvergence errors (vs N=100):")
for i, N in enumerate([1, 10, 100]):
    print(f"  N={N:3d}: Error = {errors[i]:.2e}")

print(f"\nMonotonic convergence: {'✓ PASS' if monotonic else '✗ FAIL'}")
print(f"Max trajectory difference: {max_diff:.4f} grid points")

# Final norm and energy status
for N in [1, 10, 100]:
    norm_err = data[N]['norm_error'].abs().max()
    E_drift = np.abs(data[N]['E_total'] - data[N]['E_total'].iloc[0]).max()
    print(f"\nN={N:3d}:")
    print(f"  Max norm error: {norm_err:.2e} ({'PASS' if norm_err < 0.01 else 'FAIL'})")
    print(f"  Max energy drift: {E_drift:.2e} ({'PASS' if E_drift < 0.01 else 'FAIL'})")

print("="*60)
print("Phase 1 validation plots generated successfully!")
print("="*60)
