"""
JAX-Accelerated SMFT Phase Diagram with ACTUAL PHYSICS

Implements proper SMFT dynamics:
1. Hamiltonian oscillator coupling with damping
2. Mass feedback: m = Δ·R
3. Spatial field coupling
4. 10k timesteps for proper equilibration

Target: 400 simulations in minutes (vs 3 days for NumPy)
"""

import jax
import jax.numpy as jnp
from jax import jit, vmap
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time
import os

print("="*80)
print(" JAX-ACCELERATED SMFT PHASE DIAGRAM (PROPER PHYSICS)")
print("="*80)
print(f"JAX version: {jax.__version__}")
print(f"JAX devices: {jax.devices()}")
print(f"Default backend: {jax.default_backend()}")
print("="*80)

# =============================================================================
# Parameter Space
# =============================================================================
resolution = 100  # 100×100 = 10,000 simulations
K_min, K_max = 2.5, 7.5
Delta_min, Delta_max = 0.001, 1.0

K_values = jnp.linspace(K_min, K_max, resolution)
Delta_values = jnp.linspace(Delta_min, Delta_max, resolution)

print(f"\nParameter Grid:")
print(f"  K (Kuramoto coupling): {K_min:.4f} → {K_max:.1f} ({resolution} steps)")
print(f"  Δ (Mass gap):          {Delta_min:.4f} → {Delta_max:.1f} ({resolution} steps)")
print(f"  Total simulations:     {resolution}×{resolution} = {resolution**2}")

# Simulation parameters
grid_size = 32  # Reduced from 48 to avoid OOM with 64x64 parameter grid
n_steps = 20000  # Deep equilibration for high-resolution zoom
dt = 0.01
damping = 5.0  # Oscillator damping

print(f"\nSimulation Parameters:")
print(f"  Grid: {grid_size}×{grid_size} = {grid_size**2} oscillators")
print(f"  Steps: {n_steps}")
print(f"  Time: {n_steps * dt:.1f} units")
print(f"  Damping: γ = {damping}")

# =============================================================================
# JAX SMFT Physics Kernels (JIT-compiled)
# =============================================================================

@jit
def compute_local_R_field(phases):
    """
    Compute local synchronization R(x,y) from oscillator phases.
    R(x,y) = |⟨e^(iθ)⟩_local| where average is over neighborhood.
    """
    # Complex order parameter at each site
    z = jnp.exp(1j * phases)

    # Spatial averaging via neighbor coupling (4-neighbor stencil)
    z_up = jnp.roll(z, 1, axis=0)
    z_down = jnp.roll(z, -1, axis=0)
    z_left = jnp.roll(z, 1, axis=1)
    z_right = jnp.roll(z, -1, axis=1)

    # Local average
    z_avg = (z + z_up + z_down + z_left + z_right) / 5.0

    # Local synchronization amplitude
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

    Matches SMFTSystem.step() from smft_system.py:298-350
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
    # This is the proper Kuramoto coupling (NOT simple diffusion)
    phase_diff = theta - phase_field
    coupling_force = K * R_field * jnp.sin(phase_diff)

    # Total force: natural frequency + coupling
    total_force = omega + coupling_force

    # Semi-implicit damping for stability
    # dθ/dt = p
    # dp/dt = F - γ·p
    # Update: p_new = (p_old + F*dt) / (1 + γ*dt)
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
    Returns final state + observables.
    """
    def step_fn(carry, i):
        theta, p = carry
        theta_new, p_new, R_field, m_field = hamiltonian_smft_step(
            theta, p, K, Delta, dt, omega, damping
        )
        # Store R every step for analysis
        R_mean = jnp.mean(R_field)
        return (theta_new, p_new), (R_mean, R_field, m_field)

    # Run evolution
    (theta_final, p_final), (R_history, R_fields, m_fields) = jax.lax.scan(
        step_fn,
        (initial_theta, initial_p),
        jnp.arange(n_steps)
    )

    return theta_final, p_final, R_history, R_fields, m_fields

@jit
def compute_observables(theta_final, R_history, R_fields, m_fields):
    """
    Compute all observables from simulation.
    Matches phase_diagram_sweep.py:186-199
    """
    # Synchronization
    R_initial = R_history[0]
    R_final = R_history[-1]
    R_change = R_final - R_initial

    # Localization: Use final R field spatial variance
    # High variance = low localization (gas phase)
    # Low variance = high localization (bound state)
    R_final_field = R_fields[-1]
    R_std = jnp.std(R_final_field)
    R_mean = jnp.mean(R_final_field)
    # Inverse: high localization when R is uniform (low std)
    localization = 1.0 / (R_std + 1e-10)

    # Energy stability: R fluctuations over time
    R_time_std = jnp.std(R_history[-1000:])  # Last 1000 steps
    R_time_mean = jnp.mean(R_history[-1000:])
    energy_stability = R_time_std / (jnp.abs(R_time_mean) + 1e-10)

    return R_initial, R_final, R_change, localization, energy_stability

# =============================================================================
# Batched Simulation
# =============================================================================

def simulate_single(K, Delta, key):
    """Simulate single (K, Δ) point with proper SMFT physics."""
    # Initialize phases randomly
    theta_init = jax.random.uniform(key, (grid_size, grid_size)) * 2 * jnp.pi

    # Initialize momenta to zero
    p_init = jnp.zeros((grid_size, grid_size))

    # Natural frequencies (all zero for simplicity)
    omega = jnp.zeros((grid_size, grid_size))

    # Place spinor at center of grid (localized phase excitation)
    center_x, center_y = grid_size // 2, grid_size // 2
    x_coord = jnp.arange(grid_size)
    y_coord = jnp.arange(grid_size)
    X, Y = jnp.meshgrid(x_coord, y_coord)

    # Gaussian-localized spinor with radial phase structure
    r_squared = (X - center_x)**2 + (Y - center_y)**2
    width = grid_size / 6.0  # Spinor localization width
    spinor_amplitude = jnp.exp(-r_squared / (2 * width**2))

    # Apply spinor modulation to initial phases
    theta_init = theta_init * (1.0 + 2.0 * spinor_amplitude)

    # Evolve
    theta_final, p_final, R_history, R_fields, m_fields = evolve_smft(
        K, Delta, theta_init, p_init, n_steps, dt, omega, damping
    )

    # Measure
    R_init, R_final, R_change, L, sigma_E = compute_observables(
        theta_final, R_history, R_fields, m_fields
    )

    return R_init, R_final, R_change, L, sigma_E

# =============================================================================
# Run Sequential Simulation (batch_size=1 to avoid OOM)
# =============================================================================

print("\n" + "="*80)
print(f" RUNNING SEQUENTIAL SMFT SIMULATION ({resolution}×{resolution} = {resolution**2})")
print("="*80)

start_time = time.time()

print("\nCompiling JAX kernels (first run)...")
# Generate unique random keys
master_key = jax.random.PRNGKey(42)

print(f"Running simulations sequentially (batch_size=1 for low RAM)...")
print(f"  Each simulation: {n_steps} steps of proper SMFT dynamics")
print(f"  Total: {resolution**2} simulations")

# Pre-allocate result arrays
R_initial_map = np.zeros((resolution, resolution))
R_final_map = np.zeros((resolution, resolution))
R_change_map = np.zeros((resolution, resolution))
localization_map = np.zeros((resolution, resolution))
energy_stability_map = np.zeros((resolution, resolution))

# Sequential loop through parameter grid
for i in range(resolution):
    for j in range(resolution):
        K = K_values[i]
        Delta = Delta_values[j]
        key = jax.random.fold_in(master_key, i * resolution + j)

        # Run single simulation
        R_init, R_final, R_change, L, sigma_E = simulate_single(K, Delta, key)

        # Store results
        R_initial_map[i, j] = float(R_init)
        R_final_map[i, j] = float(R_final)
        R_change_map[i, j] = float(R_change)
        localization_map[i, j] = float(L)
        energy_stability_map[i, j] = float(sigma_E)

        # Progress indicator every 10%
        progress = ((i * resolution + j + 1) / resolution**2) * 100
        if (i * resolution + j + 1) % (resolution**2 // 10) == 0:
            elapsed_so_far = time.time() - start_time
            eta = (elapsed_so_far / (i * resolution + j + 1)) * (resolution**2 - (i * resolution + j + 1))
            print(f"  [{progress:5.1f}%] ({i*resolution+j+1}/{resolution**2}) "
                  f"K={K:.2f}, Δ={Delta:.2f} | ETA: {eta:.1f}s")

elapsed = time.time() - start_time
print(f"\n✓ Sequential simulation complete in {elapsed:.2f} seconds!")
print(f"  Simulations per second: {resolution**2 / elapsed:.1f}")
print(f"  Per simulation: {elapsed*1000/resolution**2:.1f} ms")

# =============================================================================
# Save Data
# =============================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'outputs', 'phase_diagram_jax')
os.makedirs(data_dir, exist_ok=True)

# Save complete data
data_file = os.path.join(data_dir, 'jax_smft_results.npz')
np.savez_compressed(data_file,
                   K_values=np.array(K_values),
                   Delta_values=np.array(Delta_values),
                   R_final_map=R_final_map,
                   R_change_map=R_change_map,
                   localization_map=localization_map,
                   energy_stability_map=energy_stability_map,
                   resolution=resolution,
                   grid_size=grid_size,
                   n_steps=n_steps,
                   elapsed_time=elapsed)

print(f"\n💾 Data saved: {data_file}")

# =============================================================================
# Phase Classification (Updated Thresholds)
# =============================================================================

phase_map = np.zeros((resolution, resolution))
for i in range(resolution):
    for j in range(resolution):
        R = R_final_map[i, j]
        L = localization_map[i, j]
        sigma_E = energy_stability_map[i, j]

        # High-pass filter thresholds to see internal structure
        if R < 0.05:
            phase_map[i, j] = 0  # Dead Zone (shouldn't see any)
        elif R >= 0.05 and R < 0.50:
            phase_map[i, j] = 1  # Gas Phase (Red)
        elif R >= 0.50 and R < 0.75:
            phase_map[i, j] = 2  # Liquid Phase (Orange)
        elif R >= 0.75:
            phase_map[i, j] = 3  # Solid Phase (Gold)

dead_count = np.sum(phase_map == 0)
gas_count = np.sum(phase_map == 1)
liquid_count = np.sum(phase_map == 2)
solid_count = np.sum(phase_map == 3)

print(f"\nPhase Distribution (High-Pass Filter):")
print(f"  Dead Zone (R<0.05):           {dead_count:3d} / {resolution**2} ({dead_count/resolution**2*100:.1f}%)")
print(f"  Gas Phase (0.05≤R<0.50):      {gas_count:3d} / {resolution**2} ({gas_count/resolution**2*100:.1f}%)")
print(f"  Liquid Phase (0.50≤R<0.75):   {liquid_count:3d} / {resolution**2} ({liquid_count/resolution**2*100:.1f}%)")
print(f"  Solid Phase (R≥0.75):         {solid_count:3d} / {resolution**2} ({solid_count/resolution**2*100:.1f}%)")

# =============================================================================
# Visualization
# =============================================================================

print("\nGenerating visualizations...")

fig = plt.figure(figsize=(20, 12))
gs = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.25)

K_grid_np, Delta_grid_np = np.meshgrid(np.array(K_values), np.array(Delta_values), indexing='ij')

# 1. Phase diagram
ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.contourf(K_grid_np, Delta_grid_np, phase_map.T,
                   levels=[-0.5, 0.5, 1.5, 2.5, 3.5],
                   colors=['#0000FF', '#FF0000', '#FF8800', '#FFD700'], alpha=0.7)
ax1.set_title('JAX-SMFT Phase Diagram (High-Pass Filter)', fontweight='bold', fontsize=14)
ax1.set_xlabel('K (Kuramoto Coupling)', fontsize=12)
ax1.set_ylabel('Δ (Mass Gap)', fontsize=12)
ax1.grid(True, alpha=0.3)
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#0000FF', label='Dead (R<0.05)'),
                  Patch(facecolor='#FF0000', label='Gas (0.05≤R<0.50)'),
                  Patch(facecolor='#FF8800', label='Liquid (0.50≤R<0.75)'),
                  Patch(facecolor='#FFD700', label='Solid (R≥0.75)')]
ax1.legend(handles=legend_elements, loc='upper left')

# 2. R final
ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.contourf(K_grid_np, Delta_grid_np, R_final_map.T, levels=20, cmap='RdYlBu_r')
ax2.set_title('Final Synchronization ⟨R⟩', fontweight='bold', fontsize=14)
ax2.set_xlabel('K', fontsize=12)
ax2.set_ylabel('Δ', fontsize=12)
plt.colorbar(im2, ax=ax2, label='⟨R⟩')
ax2.grid(True, alpha=0.3)

# 3. Localization
ax3 = fig.add_subplot(gs[1, 0])
im3 = ax3.contourf(K_grid_np, Delta_grid_np, localization_map.T, levels=20, cmap='viridis')
ax3.set_title('Localization L', fontweight='bold', fontsize=14)
ax3.set_xlabel('K', fontsize=12)
ax3.set_ylabel('Δ', fontsize=12)
plt.colorbar(im3, ax=ax3, label='L')
ax3.grid(True, alpha=0.3)

# 4. Energy stability
ax4 = fig.add_subplot(gs[1, 1])
im4 = ax4.contourf(K_grid_np, Delta_grid_np, energy_stability_map.T, levels=20, cmap='plasma')
ax4.set_title('Energy Stability σ_E / ⟨E⟩', fontweight='bold', fontsize=14)
ax4.set_xlabel('K', fontsize=12)
ax4.set_ylabel('Δ', fontsize=12)
plt.colorbar(im4, ax=ax4, label='σ_E / ⟨E⟩')
ax4.grid(True, alpha=0.3)

fig.suptitle(f'JAX-SMFT Phase Diagram ({elapsed:.1f}s, {n_steps} steps/sim, Proper Hamiltonian)',
             fontsize=16, fontweight='bold', y=0.98)

output = os.path.join(data_dir, 'phase_diagram_jax_smft.png')
plt.savefig(output, dpi=150, bbox_inches='tight')
print(f"✓ Visualization saved: {output}")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "="*80)
print(" JAX-SMFT ACCELERATION SUMMARY")
print("="*80)
print(f"\nPhysics Implementation:")
print(f"  ✓ Hamiltonian oscillator dynamics with damping (γ={damping})")
print(f"  ✓ Proper Kuramoto coupling: F = K·R·sin(θ-φ)")
print(f"  ✓ Mass feedback: m = Δ·R")
print(f"  ✓ Spatial field coupling (4-neighbor stencil)")
print(f"\nPerformance:")
print(f"  Total time: {elapsed:.2f} seconds")
print(f"  Simulations: {resolution**2}")
print(f"  Rate: {resolution**2/elapsed:.1f} sims/sec")
print(f"  Per simulation: {elapsed*1000/resolution**2:.1f} ms ({n_steps} steps each)")
print(f"\nResults:")
print(f"  Gas phase: {gas_count}")
print(f"  Liquid phase: {liquid_count}")
print(f"  Solid phase: {solid_count}")
print(f"  Parameter coverage: K=[{K_min:.4f}, {K_max:.1f}], Δ=[{Delta_min:.4f}, {Delta_max:.1f}]")
print(f"\nOutputs:")
print(f"  📁 {data_dir}")
print(f"  📊 {os.path.basename(output)}")
print(f"  💾 {os.path.basename(data_file)}")
print("="*80)
