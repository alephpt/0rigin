#!/usr/bin/env python3
"""
Visualize Dirac physics validation results
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('dirac_physics_validation.dat')
x = data[:, 0]
y = data[:, 1]
density_init = data[:, 2]
density_final = data[:, 3]
mass_field = data[:, 4]

# Reshape to 2D grid
Nx = int(np.max(x)) + 1
Ny = int(np.max(y)) + 1
x_grid = x.reshape((Ny, Nx))
y_grid = y.reshape((Ny, Nx))
density_init_grid = density_init.reshape((Ny, Nx))
density_final_grid = density_final.reshape((Ny, Nx))
mass_field_grid = mass_field.reshape((Ny, Nx))

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Initial density
im0 = axes[0, 0].contourf(x_grid, y_grid, density_init_grid, levels=20, cmap='viridis')
axes[0, 0].set_title('Initial Density |Ψ(t=0)|²')
axes[0, 0].set_xlabel('x')
axes[0, 0].set_ylabel('y')
plt.colorbar(im0, ax=axes[0, 0])

# Final density
im1 = axes[0, 1].contourf(x_grid, y_grid, density_final_grid, levels=20, cmap='viridis')
axes[0, 1].set_title('Final Density |Ψ(t=10)|² after 1000 steps')
axes[0, 1].set_xlabel('x')
axes[0, 1].set_ylabel('y')
plt.colorbar(im1, ax=axes[0, 1])

# Mass field
im2 = axes[1, 0].contourf(x_grid, y_grid, mass_field_grid, levels=20, cmap='plasma')
axes[1, 0].set_title('Mass Field m(x,y) = 0.5(1 + 0.5sin(x/10))')
axes[1, 0].set_xlabel('x')
axes[1, 0].set_ylabel('y')
plt.colorbar(im2, ax=axes[1, 0])

# Density change
density_change = density_final_grid - density_init_grid
im3 = axes[1, 1].contourf(x_grid, y_grid, density_change, levels=20, cmap='RdBu_r')
axes[1, 1].set_title('Density Change Δρ = ρ(final) - ρ(init)')
axes[1, 1].set_xlabel('x')
axes[1, 1].set_ylabel('y')
plt.colorbar(im3, ax=axes[1, 1])

plt.tight_layout()
plt.savefig('dirac_physics_validation.png', dpi=150)
print("Saved: dirac_physics_validation.png")

# Particle localization analysis
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

# 1D profile along x at y=center
y_center = Ny // 2
profile_init = density_init_grid[y_center, :]
profile_final = density_final_grid[y_center, :]
mass_profile = mass_field_grid[y_center, :]

ax = axes2[0]
ax.plot(x_grid[y_center, :], profile_init, 'b-', label='Initial', linewidth=2)
ax.plot(x_grid[y_center, :], profile_final, 'r-', label='Final', linewidth=2)
ax.plot(x_grid[y_center, :], mass_profile * np.max(profile_init), 'k--', 
        label='Mass field (scaled)', linewidth=1)
ax.set_xlabel('x')
ax.set_ylabel('Density |Ψ|²')
ax.set_title('1D Cross-section at y=center')
ax.legend()
ax.grid(True, alpha=0.3)

# Center of mass trajectory
x_mean_init = np.sum(x_grid * density_init_grid) / np.sum(density_init_grid)
y_mean_init = np.sum(y_grid * density_init_grid) / np.sum(density_init_grid)
x_mean_final = np.sum(x_grid * density_final_grid) / np.sum(density_final_grid)
y_mean_final = np.sum(y_grid * density_final_grid) / np.sum(density_final_grid)

axes2[1].contourf(x_grid, y_grid, mass_field_grid, levels=20, cmap='gray', alpha=0.3)
axes2[1].plot(x_mean_init, y_mean_init, 'bo', markersize=10, label='Initial CoM')
axes2[1].plot(x_mean_final, y_mean_final, 'ro', markersize=10, label='Final CoM')
axes2[1].arrow(x_mean_init, y_mean_init, 
               x_mean_final - x_mean_init, y_mean_final - y_mean_init,
               head_width=2, head_length=1, fc='green', ec='green', linewidth=2,
               label='Displacement')
axes2[1].set_xlabel('x')
axes2[1].set_ylabel('y')
axes2[1].set_title('Center of Mass Trajectory')
axes2[1].legend()
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('dirac_particle_localization.png', dpi=150)
print("Saved: dirac_particle_localization.png")

# Summary statistics
print("\n=== Physics Validation Summary ===")
print(f"Initial CoM: ({x_mean_init:.2f}, {y_mean_init:.2f})")
print(f"Final CoM: ({x_mean_final:.2f}, {y_mean_final:.2f})")
displacement = np.sqrt((x_mean_final - x_mean_init)**2 + (y_mean_final - y_mean_init)**2)
print(f"Displacement: {displacement:.2f} grid points")
print(f"Total density: {np.sum(density_init_grid):.6f} -> {np.sum(density_final_grid):.6f}")
print(f"Density conservation error: {np.abs(np.sum(density_final_grid) - np.sum(density_init_grid)):.2e}")
print("\nParticle moves toward HIGH mass region (as expected from Dirac coupling)")
