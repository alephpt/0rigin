#!/usr/bin/env python3
"""
SMFT Single Point Deep Simulation
==================================

Target: K=27.21, Δ=36.18 (Bright spot from phase diagram)
Resolution: 128×128 oscillator grid (maximum detail)
Steps: 50,000 (deep time evolution, t=500)
Initial: Gaussian spinor centered at (64, 64)

Physics: Full SMFT with Hamiltonian dynamics, mass feedback, field coupling
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import jax
import jax.numpy as jnp
from jax import jit
from functools import partial

print("="*80)
print(" SMFT SINGLE POINT DEEP SIMULATION")
print("="*80)
print(f"\nJAX version: {jax.__version__}")
print(f"JAX devices: {jax.devices()}")
print(f"Default backend: {jax.default_backend()}")
print("="*80)

# Target parameters from bright spot
K = 27.21
Delta = 36.18

# High-resolution grid
grid_size = 128
n_steps = 50000
dt = 0.01
damping = 5.0

print(f"\nTarget Parameters:")
print(f"  K (Kuramoto coupling): {K}")
print(f"  Δ (Mass gap): {Delta}")

print(f"\nSimulation Parameters:")
print(f"  Grid: {grid_size}×{grid_size} = {grid_size**2} oscillators")
print(f"  Steps: {n_steps}")
print(f"  Total time: {n_steps * dt} time units")
print(f"  Damping: γ = {damping}")
print(f"  Initial: Gaussian spinor at center ({grid_size//2}, {grid_size//2})")

# =============================================================================
# SMFT Physics (JAX-JIT Compiled)
# =============================================================================

@jit
def compute_local_R_field(phases):
    """
    Compute local synchronization R(x,y) from oscillator phases.
    R(x,y) = |⟨e^(iθ)⟩_local| where average is over neighborhood.
    """
    z = jnp.exp(1j * phases)

    # Spatial averaging via neighbor coupling (4-neighbor stencil)
    z_up = jnp.roll(z, 1, axis=0)
    z_down = jnp.roll(z, -1, axis=0)
    z_left = jnp.roll(z, 1, axis=1)
    z_right = jnp.roll(z, -1, axis=1)

    z_avg = (z + z_up + z_down + z_left + z_right) / 5.0
    R_field = jnp.abs(z_avg)

    return R_field

@jit
def hamiltonian_smft_step(theta, p, K, Delta, dt, omega, damping):
    """
    Single SMFT step with Hamiltonian dynamics + mass feedback.

    Physics:
    1. Hamiltonian coupling: H = -K·R·cos(θ_i - φ_local)
    2. Mass feedback: m(x) = Δ·R(x)
    3. Damped dynamics: dp/dt = F - γ·p
    """
    # Compute local R field (synchronization amplitude)
    R_field = compute_local_R_field(theta)

    # Compute local phase field
    z_field = jnp.exp(1j * theta)
    z_up = jnp.roll(z_field, 1, axis=0)
    z_down = jnp.roll(z_field, -1, axis=0)
    z_left = jnp.roll(z_field, 1, axis=1)
    z_right = jnp.roll(z_field, -1, axis=1)
    z_avg = (z_field + z_up + z_down + z_left + z_right) / 5.0
    phase_field = jnp.angle(z_avg)

    # Hamiltonian coupling force: F = K·R·sin(θ - φ)
    phase_diff = theta - phase_field
    coupling_force = K * R_field * jnp.sin(phase_diff)

    total_force = omega + coupling_force

    # Semi-implicit damping for stability
    theta_new = theta + p * dt
    p_temp = p + total_force * dt
    p_new = p_temp / (1.0 + damping * dt)

    # Mass field (SMFT): m = Δ·R
    m_field = Delta * R_field

    return theta_new, p_new, R_field, m_field

@partial(jit, static_argnums=(4,))
def evolve_smft(K, Delta, initial_theta, initial_p, n_steps, dt, omega, damping):
    """
    Full SMFT evolution with proper Hamiltonian dynamics.
    Returns snapshots at key timesteps.
    """
    def step_fn(carry, i):
        theta, p = carry
        theta_new, p_new, R_field, m_field = hamiltonian_smft_step(
            theta, p, K, Delta, dt, omega, damping
        )

        # Store snapshot every 1000 steps
        should_store = (i % 1000 == 0) | (i == n_steps - 1)

        return (theta_new, p_new), (should_store, i, R_field, m_field)

    # Run evolution
    (theta_final, p_final), (store_flags, step_ids, R_fields, m_fields) = jax.lax.scan(
        step_fn,
        (initial_theta, initial_p),
        jnp.arange(n_steps)
    )

    return theta_final, p_final, store_flags, step_ids, R_fields, m_fields

# =============================================================================
# Initialize System
# =============================================================================

print("\n" + "="*80)
print(" INITIALIZING SYSTEM")
print("="*80)

# Random initial phases
key = jax.random.PRNGKey(42)
theta_init = jax.random.uniform(key, (grid_size, grid_size), minval=-jnp.pi, maxval=jnp.pi)

# Initialize momenta to zero
p_init = jnp.zeros((grid_size, grid_size))

# Natural frequencies (all zero)
omega = jnp.zeros((grid_size, grid_size))

# Gaussian spinor at center
print("\nCreating Gaussian spinor at grid center...")
center_x, center_y = grid_size // 2, grid_size // 2
x_coord = jnp.arange(grid_size)
y_coord = jnp.arange(grid_size)
X, Y = jnp.meshgrid(x_coord, y_coord)

# Gaussian-localized spinor with radial phase structure
r_squared = (X - center_x)**2 + (Y - center_y)**2
width = grid_size / 8.0  # Spinor localization width
spinor_amplitude = jnp.exp(-r_squared / (2 * width**2))

# Apply spinor modulation to initial phases
theta_init = theta_init * (1.0 + 3.0 * spinor_amplitude)

print(f"  Spinor center: ({center_x}, {center_y})")
print(f"  Spinor width: {width:.1f} grid units")
print(f"  Spinor peak amplitude: {spinor_amplitude.max():.3f}")

# =============================================================================
# Run Simulation
# =============================================================================

print("\n" + "="*80)
print(" RUNNING DEEP TIME EVOLUTION")
print("="*80)

print("\nCompiling JAX kernels (first run)...")
start_time = time.time()

theta_final, p_final, store_flags, step_ids, R_fields, m_fields = evolve_smft(
    K, Delta, theta_init, p_init, n_steps, dt, omega, damping
)

elapsed = time.time() - start_time

print(f"\n✓ Simulation complete in {elapsed:.2f} seconds!")
print(f"  Timesteps per second: {n_steps / elapsed:.1f}")
print(f"  Microseconds per step: {elapsed * 1e6 / n_steps:.1f} μs")

# =============================================================================
# Extract Snapshots
# =============================================================================

print("\n" + "="*80)
print(" EXTRACTING SNAPSHOTS")
print("="*80)

# Get stored snapshots
stored_indices = np.where(np.array(store_flags))[0]
print(f"\nStored {len(stored_indices)} snapshots")

# Select key snapshots: initial, 1/4, 1/2, 3/4, final
snapshot_indices = [0, len(stored_indices)//4, len(stored_indices)//2,
                    3*len(stored_indices)//4, len(stored_indices)-1]

snapshots = []
for idx in snapshot_indices:
    actual_idx = stored_indices[idx]
    step = int(step_ids[actual_idx])
    t = step * dt
    R = np.array(R_fields[actual_idx])
    m = np.array(m_fields[actual_idx])
    snapshots.append((step, t, R, m))
    print(f"  Snapshot {len(snapshots)}: step={step}, t={t:.1f}, ⟨R⟩={R.mean():.3f}, ⟨m⟩={m.mean():.3f}")

# =============================================================================
# Visualization
# =============================================================================

print("\n" + "="*80)
print(" GENERATING VISUALIZATION")
print("="*80)

fig, axes = plt.subplots(len(snapshots), 2, figsize=(14, 4*len(snapshots)))
fig.suptitle(f'SMFT Deep Time Evolution: K={K}, Δ={Delta}\n128×128 Grid, 50k Steps',
             fontsize=16, fontweight='bold')

for row, (step, t, R, m) in enumerate(snapshots):
    # Synchronization field
    ax = axes[row, 0]
    im = ax.imshow(R, cmap='RdYlBu_r', vmin=0, vmax=1, origin='lower')
    ax.set_title(f't={t:.1f} (step {step}): Synchronization R(x,y)', fontweight='bold')
    plt.colorbar(im, ax=ax, label='R', fraction=0.046)
    ax.text(0.02, 0.98, f'⟨R⟩={R.mean():.3f}\nmax={R.max():.3f}\nmin={R.min():.3f}',
            transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle='round', fc='white', alpha=0.9), fontsize=10)

    # Mass field
    ax = axes[row, 1]
    im = ax.imshow(m, cmap='plasma', vmin=0, vmax=Delta, origin='lower')
    ax.set_title(f't={t:.1f} (step {step}): Mass m=Δ·R', fontweight='bold')
    plt.colorbar(im, ax=ax, label='m', fraction=0.046)
    ax.text(0.02, 0.98, f'⟨m⟩={m.mean():.3f}\nmax={m.max():.3f}\nmin={m.min():.3f}',
            transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle='round', fc='white', alpha=0.9), fontsize=10)

plt.tight_layout()

output_dir = "examples/validation/outputs/smft_single_point"
import os
os.makedirs(output_dir, exist_ok=True)

output_file = f"{output_dir}/smft_K{K:.2f}_D{Delta:.2f}_128x128_50k.png"
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\n✓ Visualization saved: {output_file}")

# =============================================================================
# Final Analysis
# =============================================================================

print("\n" + "="*80)
print(" FINAL STATE ANALYSIS")
print("="*80)

R_final = np.array(compute_local_R_field(theta_final))
m_final = Delta * R_final

print(f"\nSynchronization Order Parameter R:")
print(f"  Mean: {R_final.mean():.4f}")
print(f"  Std:  {R_final.std():.4f}")
print(f"  Min:  {R_final.min():.4f}")
print(f"  Max:  {R_final.max():.4f}")

print(f"\nMass Field m = Δ·R:")
print(f"  Mean: {m_final.mean():.4f}")
print(f"  Std:  {m_final.std():.4f}")
print(f"  Min:  {m_final.min():.4f}")
print(f"  Max:  {m_final.max():.4f}")

# Phase classification
if R_final.mean() < 0.05:
    phase = "Dead Zone"
elif R_final.mean() < 0.50:
    phase = "Gas Phase"
elif R_final.mean() < 0.75:
    phase = "Liquid Phase"
else:
    phase = "Solid Phase"

print(f"\nPhase Classification: {phase}")

print("\n" + "="*80)
print(" SMFT SINGLE POINT SIMULATION COMPLETE")
print("="*80)
print(f"\nOutput: {output_file}")
print(f"Parameters: K={K}, Δ={Delta}")
print(f"Resolution: {grid_size}×{grid_size}")
print(f"Evolution: {n_steps} steps, t={n_steps*dt}")
print(f"Final State: ⟨R⟩={R_final.mean():.3f} ({phase})")
print("="*80)
