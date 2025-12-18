#!/usr/bin/env python3
"""Visualize validation results"""

import numpy as np
import matplotlib.pyplot as plt

# Load norm vs time data
data = np.loadtxt('norm_vs_time.dat')
steps = data[:, 0]
norms = data[:, 1]
drifts = data[:, 2]

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Norm vs Time
ax1.plot(steps, norms, 'b-', linewidth=2, label='|Ψ|²')
ax1.axhline(y=1.0, color='r', linestyle='--', label='Perfect conservation')
ax1.set_xlabel('Timestep', fontsize=12)
ax1.set_ylabel('Norm |Ψ|²', fontsize=12)
ax1.set_title('Norm Conservation over 50,000 Steps', fontsize=14, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0.99, 1.001])

# Plot 2: Drift vs Time  
ax2.plot(steps, drifts * 100, 'r-', linewidth=2)
ax2.set_xlabel('Timestep', fontsize=12)
ax2.set_ylabel('Norm Drift (%)', fontsize=12)
ax2.set_title('Cumulative Drift from Initial Norm', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)

# Add linear fit to show drift rate
z = np.polyfit(steps, drifts, 1)
p = np.poly1d(z)
ax2.plot(steps, p(steps) * 100, 'k--', alpha=0.5, label=f'Linear fit: {z[0]*100:.2e}%/step')
ax2.legend()

plt.tight_layout()
plt.savefig('norm_conservation.png', dpi=150)
print("Saved: norm_conservation.png")

# Summary statistics
print(f"\n=== Norm Conservation Statistics ===")
print(f"Initial norm: {norms[0]:.6f}")
print(f"Final norm (50k): {norms[-1]:.6f}")
print(f"Total drift: {drifts[-1]:.6f} ({drifts[-1]*100:.3f}%)")
print(f"Drift rate: {z[0]:.2e} per step ({z[0]*100:.2e}% per step)")
print(f"Extrapolated to 1M steps: {z[0]*1e6*100:.2f}%")
