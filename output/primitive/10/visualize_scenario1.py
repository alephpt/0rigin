#!/usr/bin/env python3
"""
Visualize Phase 2 Scenario 1: Static Defect Localization
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches

# Load data
trajectory = np.loadtxt('defect_localization/trajectory.dat')
R_initial = np.loadtxt('defect_localization/R_initial.dat')
density_final = np.loadtxt('defect_localization/density_final.dat')
rho_vs_R = np.loadtxt('defect_localization/rho_vs_R.dat')

# Grid parameters
NX, NY = 128, 128
DEFECT_X, DEFECT_Y = 64, 64

# Reshape R field
R_init_grid = R_initial[:, 2].reshape((NY, NX))
density_grid = density_final[:, 2].reshape((NY, NX))
R_final_grid = density_final[:, 3].reshape((NY, NX))

# ===== Figure 1: Trajectory over R field =====
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Initial state
im1 = ax1.imshow(R_init_grid, origin='lower', cmap='viridis', extent=[0, NX, 0, NY])
ax1.plot(trajectory[:, 1], trajectory[:, 2], 'r-', linewidth=2, alpha=0.8, label='Trajectory')
ax1.scatter(trajectory[0, 1], trajectory[0, 2], s=100, c='white', edgecolor='red',
           linewidth=2, marker='o', label='Start', zorder=5)
ax1.scatter(trajectory[-1, 1], trajectory[-1, 2], s=100, c='red', edgecolor='white',
           linewidth=2, marker='s', label='End', zorder=5)
ax1.scatter(DEFECT_X, DEFECT_Y, s=200, c='cyan', marker='x', linewidth=3,
           label='Defect', zorder=5)
ax1.set_xlabel('x (grid points)', fontsize=12)
ax1.set_ylabel('y (grid points)', fontsize=12)
ax1.set_title('Particle Trajectory over Synchronization Field R(x,y)', fontsize=14, fontweight='bold')
ax1.legend(loc='upper left', fontsize=10)
ax1.grid(True, alpha=0.3)
cbar1 = plt.colorbar(im1, ax=ax1)
cbar1.set_label('R (synchronization)', fontsize=11)

# Final state
im2 = ax2.imshow(density_grid, origin='lower', cmap='hot',
                extent=[0, NX, 0, NY], norm=LogNorm(vmin=1e-8, vmax=density_grid.max()))
ax2.contour(R_final_grid, levels=10, colors='cyan', alpha=0.5, linewidths=0.5, extent=[0, NX, 0, NY])
ax2.scatter(DEFECT_X, DEFECT_Y, s=200, c='cyan', marker='x', linewidth=3, label='Defect')
ax2.set_xlabel('x (grid points)', fontsize=12)
ax2.set_ylabel('y (grid points)', fontsize=12)
ax2.set_title('Final Density |Ψ|² with R Contours', fontsize=14, fontweight='bold')
ax2.legend(loc='upper left', fontsize=10)
cbar2 = plt.colorbar(im2, ax=ax2)
cbar2.set_label('|Ψ|² (log scale)', fontsize=11)

plt.tight_layout()
plt.savefig('defect_localization/trajectory_overlay.png', dpi=150, bbox_inches='tight')
print("Saved: trajectory_overlay.png")

# ===== Figure 2: Time evolution =====
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Distance to defect vs time
ax = axes[0, 0]
steps = trajectory[:, 0]
dist = np.sqrt(trajectory[:, 5]**2 + trajectory[:, 6]**2)
ax.plot(steps, dist, 'b-', linewidth=2)
ax.axhline(y=dist[0], color='r', linestyle='--', alpha=0.5, label=f'Initial: {dist[0]:.2f}')
ax.axhline(y=dist[-1], color='g', linestyle='--', alpha=0.5, label=f'Final: {dist[-1]:.2f}')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('Distance to Defect (grid points)', fontsize=12)
ax.set_title('Particle Distance to Defect vs Time', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# R at CoM vs time
ax = axes[0, 1]
ax.plot(steps, trajectory[:, 3], 'purple', linewidth=2, label='R at CoM')
ax.axhline(y=trajectory[:, 4].mean(), color='gray', linestyle='--',
          alpha=0.5, label=f'<R>_global: {trajectory[:, 4].mean():.4f}')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('R (synchronization)', fontsize=12)
ax.set_title('Synchronization at Particle Location', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# CoM trajectory (x and y separately)
ax = axes[1, 0]
ax.plot(steps, trajectory[:, 1], 'b-', linewidth=2, label='x_CoM')
ax.plot(steps, trajectory[:, 2], 'r-', linewidth=2, label='y_CoM')
ax.axhline(y=DEFECT_X, color='b', linestyle='--', alpha=0.3, label='Defect x')
ax.axhline(y=DEFECT_Y, color='r', linestyle='--', alpha=0.3, label='Defect y')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('Position (grid points)', fontsize=12)
ax.set_title('Center of Mass vs Time', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# ρ vs R correlation
ax = axes[1, 1]
ax.scatter(rho_vs_R[:, 0], rho_vs_R[:, 1], alpha=0.3, s=1, c='blue')
ax.set_xlabel('R (synchronization)', fontsize=12)
ax.set_ylabel('ρ = |Ψ|² (density)', fontsize=12)
ax.set_title('Density vs Synchronization Correlation', fontsize=13, fontweight='bold')
ax.set_yscale('log')
ax.grid(True, alpha=0.3)

# Add trend analysis
R_bins = np.linspace(rho_vs_R[:, 0].min(), rho_vs_R[:, 0].max(), 20)
rho_binned = []
R_binned = []
for i in range(len(R_bins)-1):
    mask = (rho_vs_R[:, 0] >= R_bins[i]) & (rho_vs_R[:, 0] < R_bins[i+1])
    if mask.sum() > 0:
        R_binned.append(0.5 * (R_bins[i] + R_bins[i+1]))
        rho_binned.append(rho_vs_R[mask, 1].mean())

if len(R_binned) > 0:
    ax.plot(R_binned, rho_binned, 'r-', linewidth=2, label='Binned average')
    ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig('defect_localization/evolution_analysis.png', dpi=150, bbox_inches='tight')
print("Saved: evolution_analysis.png")

# ===== Figure 3: R field comparison =====
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

# Initial R
im1 = ax1.imshow(R_init_grid, origin='lower', cmap='viridis', vmin=0.85, vmax=1.0)
ax1.scatter(DEFECT_X, DEFECT_Y, s=200, c='red', marker='x', linewidth=3)
circle = mpatches.Circle((DEFECT_X, DEFECT_Y), 10, fill=False, edgecolor='red', linewidth=2, linestyle='--')
ax1.add_patch(circle)
ax1.set_xlabel('x', fontsize=12)
ax1.set_ylabel('y', fontsize=12)
ax1.set_title('Initial R Field (t=0)', fontsize=13, fontweight='bold')
cbar1 = plt.colorbar(im1, ax=ax1)
cbar1.set_label('R', fontsize=11)

# Final R
im2 = ax2.imshow(R_final_grid, origin='lower', cmap='viridis', vmin=0.85, vmax=1.0)
ax2.scatter(DEFECT_X, DEFECT_Y, s=200, c='red', marker='x', linewidth=3)
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('y', fontsize=12)
ax2.set_title('Final R Field (t=final)', fontsize=13, fontweight='bold')
cbar2 = plt.colorbar(im2, ax=ax2)
cbar2.set_label('R', fontsize=11)

# R cross-section through defect
ax3.plot(R_init_grid[DEFECT_Y, :], 'b-', linewidth=2, label='Initial R(x, y=64)')
ax3.plot(R_final_grid[DEFECT_Y, :], 'r-', linewidth=2, label='Final R(x, y=64)')
ax3.axvline(x=DEFECT_X, color='gray', linestyle='--', alpha=0.5, label='Defect center')
ax3.set_xlabel('x (grid points)', fontsize=12)
ax3.set_ylabel('R', fontsize=12)
ax3.set_title('R Cross-Section Through Defect', fontsize=13, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('defect_localization/R_field_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: R_field_comparison.png")

print("\n=== VISUALIZATION COMPLETE ===")
print("Output:")
print("  - trajectory_overlay.png")
print("  - evolution_analysis.png")
print("  - R_field_comparison.png")
