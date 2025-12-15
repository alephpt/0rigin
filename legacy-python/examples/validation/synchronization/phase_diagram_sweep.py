"""
MSFT Phase Diagram: Map of Existence

Scans the (K, Δ) parameter space to find regions where stable bound states form.

The "Universe Scanner":
- K (Vacuum Coherence / Kuramoto Coupling): 0.0 → 10.0
- Δ (Mass Gap / Energy Scale): 0.0 → 5.0
- Resolution: 20×20 grid = 400 universe simulations

Expected Regions:
- Dead Zone (Blue): K too low, vacuum chaotic, no mass formation
- Gas Phase (Red): K high, Δ low, mass forms but no confinement
- Bound States (Gold): Specific K/Δ ratios where particles are trapped
  - Island 1: Low energy stable state (electron-like?)
  - Island 2: High energy stable state (muon-like?)

Output: 2D phase diagrams + 3D surface plots showing the emergence of structure
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
import time
import csv
import os

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')
sys.stdout.reconfigure(line_buffering=True)  # Force unbuffered output
from kuramoto.field_theory import MSFTSystem
from kuramoto.field_theory.dirac.gamma_matrices import get_gamma_matrices_3plus1
from kuramoto.field_theory.dirac.spinor_field import DiracSpinorField, SpatialGrid

print("="*80)
print(" MSFT PHASE DIAGRAM: MAPPING THE UNIVERSE")
print("="*80)
print("\nScanning (K, Δ) parameter space to discover regions of existence:")
print("  • Dead Zone: Vacuum fails to nucleate")
print("  • Gas Phase: Mass forms but particles scatter")
print("  • Bound States: Stable particle confinement (THE GOAL)")
print("="*80)

# =============================================================================
# Parameter Space Grid
# =============================================================================
resolution = 100
K_min, K_max = 2.5, 7.5
Delta_min, Delta_max = 0.001, 1.0

K_values = np.linspace(K_min, K_max, resolution)
Delta_values = np.linspace(Delta_min, Delta_max, resolution)

K_grid, Delta_grid = np.meshgrid(K_values, Delta_values)

print(f"\nParameter Grid:")
print(f"  K (Kuramoto coupling): {K_min:.1f} → {K_max:.1f} ({resolution} steps)")
print(f"  Δ (Mass gap):          {Delta_min:.1f} → {Delta_max:.1f} ({resolution} steps)")
print(f"  Total simulations:     {resolution}×{resolution} = {resolution**2}")

# Fixed parameters
grid_size = 32
n_steps = 10000  # Long evolution to reach equilibrium
dt_MSFT = 0.01
dt_dirac = 0.0001
chiral_angle = 0.0
momentum = (0.0, 0.0)  # Start with zero momentum for confinement test

# =============================================================================
# State Vector: What We Measure at Each (K, Δ)
# =============================================================================
# Observables that define the "state" of the universe
R_final_map = np.zeros((resolution, resolution))        # Mass formation
R_change_map = np.zeros((resolution, resolution))       # Stability
localization_map = np.zeros((resolution, resolution))   # Confinement
energy_stability_map = np.zeros((resolution, resolution))  # Energy fluctuation
norm_map = np.zeros((resolution, resolution))           # Probability conservation

# =============================================================================
# Data Logging Setup
# =============================================================================
# Get script directory
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, 'outputs', 'phase_diagram')
states_dir = os.path.join(data_dir, 'states')  # Individual state PNGs
os.makedirs(data_dir, exist_ok=True)
os.makedirs(states_dir, exist_ok=True)

csv_file = os.path.join(data_dir, 'scan_progress.csv')
npz_checkpoint = os.path.join(data_dir, 'checkpoint.npz')
log_file = os.path.join(data_dir, 'scan_log.txt')

# Setup text log file
import sys
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

log_fp = open(log_file, 'w', buffering=1)
sys.stdout = Tee(sys.stdout, log_fp)
sys.stderr = Tee(sys.stderr, log_fp)

# Create CSV header
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['index', 'K', 'Delta', 'R_initial', 'R_final', 'R_change',
                    'localization', 'energy_stability', 'norm', 'status', 'timestamp'])

print(f"\n📊 Data logging enabled:")
print(f"  Progress CSV: {csv_file}")
print(f"  Text log: {log_file}")
print(f"  State PNGs: {states_dir}/")
print(f"  Checkpoint: {npz_checkpoint}")

# =============================================================================
# The Scanner: Iterate Over All (K, Δ) Coordinates
# =============================================================================
print(f"\nScanning {resolution**2} universes...")
print(f"  Evolution time: {n_steps * dt_MSFT:.1f} time units per simulation")
print(f"  Estimated runtime: ~{resolution**2 * 0.5 / 60:.0f} minutes")
print()

start_time = time.time()
checkpoint_interval = 50  # Save checkpoint every N simulations
sim_count = 0

for i in range(resolution):
    for j in range(resolution):
        K = K_grid[i, j]
        Delta = Delta_grid[i, j]
        sim_count += 1

        progress = (i * resolution + j + 1) / resolution**2 * 100
        elapsed = time.time() - start_time
        eta = (elapsed / sim_count) * (resolution**2 - sim_count) if sim_count > 0 else 0
        print(f"[{progress:5.1f}%] ({sim_count}/{resolution**2}) K={K:5.2f}, Δ={Delta:4.2f} | ETA: {eta/60:.1f}m ... ", end='', flush=True)

        try:
            # Initialize MSFT system
            MSFT = MSFTSystem(N=grid_size**2, grid_shape=(grid_size, grid_size),
                            mass_gap=Delta)
            MSFT.oscillators.coupling_strength = K

            # Seed synchronization
            n = grid_size
            x_coord, y_coord = np.meshgrid(np.arange(n), np.arange(n))
            pattern = 0.3 + 0.4 * np.cos(2*np.pi*x_coord/n * 1.5) * np.cos(2*np.pi*y_coord/n * 1.5)
            MSFT.sync_field.values = pattern.flatten()

            # Initialize spinor
            grid = SpatialGrid(Nx=grid_size, Ny=grid_size, Lx=1.0, Ly=1.0, boundary='periodic')
            spinor = DiracSpinorField(grid, initial_state='plane_wave', momentum=momentum)

            # Record initial state
            R_initial = MSFT.sync_field.values.mean()
            energy_history = []

            # Evolve
            for step in range(n_steps):
                MSFT.step(dt=dt_MSFT)

                R_field = MSFT.sync_field.values.reshape(MSFT.grid_shape)
                m_field = MSFT.compute_effective_mass().reshape(MSFT.grid_shape)

                # Evolve spinor
                substeps = int(dt_MSFT / dt_dirac)
                for sub in range(substeps):
                    spinor.step_rk4(dt_dirac, m_field, chiral_angle=chiral_angle,
                                  Delta=Delta, normalize=True)

                # Record energy every 10 steps
                if step % 10 == 0:
                    H_psi = spinor.compute_hamiltonian_action(m_field, chiral_angle, Delta)
                    energy = np.sum(np.real(np.conj(spinor.psi) * H_psi)) * grid.dx * grid.dy
                    energy_history.append(energy)

            # Measure observables
            R_final = MSFT.sync_field.values.mean()
            R_change = R_final - R_initial

            density = spinor.compute_density()
            norm_final = np.sum(density) * grid.dx * grid.dy

            # Localization: Inverse Participation Ratio
            # L = (∫ρ²dV) / (∫ρdV)² - Higher L means more localized
            L = np.sum(density**2) * grid.dx * grid.dy / (np.sum(density) * grid.dx * grid.dy)**2
            localization = L * grid_size**2  # Normalize by grid size

            # Energy stability: σ_E / ⟨E⟩
            if len(energy_history) > 1:
                energy_stability = np.std(energy_history) / (abs(np.mean(energy_history)) + 1e-10)
            else:
                energy_stability = np.nan

            # Store results
            R_final_map[i, j] = R_final
            R_change_map[i, j] = R_change
            localization_map[i, j] = localization
            energy_stability_map[i, j] = energy_stability
            norm_map[i, j] = norm_final

            # Save individual state visualization
            fig_state, axes = plt.subplots(2, 2, figsize=(10, 10))

            # R field
            im0 = axes[0,0].imshow(R_field, cmap='RdYlBu_r', vmin=0, vmax=1, origin='lower')
            axes[0,0].set_title(f'Sync R(x,y) | ⟨R⟩={R_final:.3f}')
            plt.colorbar(im0, ax=axes[0,0])

            # Mass field
            im1 = axes[0,1].imshow(m_field, cmap='plasma', vmin=0, vmax=Delta, origin='lower')
            axes[0,1].set_title(f'Mass m=Δ·R | ⟨m⟩={Delta*R_final:.3f}')
            plt.colorbar(im1, ax=axes[0,1])

            # Spinor density
            im2 = axes[1,0].imshow(density, cmap='hot', origin='lower')
            axes[1,0].set_title(f'Spinor Density | L={localization:.1f}')
            plt.colorbar(im2, ax=axes[1,0])

            # Phase
            phase = np.angle(spinor.psi[:, :, 0])
            im3 = axes[1,1].imshow(phase, cmap='twilight', vmin=-np.pi, vmax=np.pi, origin='lower')
            axes[1,1].set_title(f'Spinor Phase arg(ψ)')
            plt.colorbar(im3, ax=axes[1,1])

            fig_state.suptitle(f'State #{sim_count:03d}: K={K:.2f}, Δ={Delta:.2f}', fontweight='bold')
            plt.tight_layout()
            state_png = os.path.join(states_dir, f'{sim_count:03d}.png')
            plt.savefig(state_png, dpi=100, bbox_inches='tight')
            plt.close(fig_state)

            # Log to CSV immediately
            with open(csv_file, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([sim_count, K, Delta, R_initial, R_final, R_change,
                               localization, energy_stability, norm_final, 'SUCCESS',
                               time.time()])

            print(f"R={R_final:.3f}, L={localization:.2f}, σ_E={energy_stability:.3f} → saved {sim_count:03d}.png")

        except Exception as e:
            print(f"FAILED: {e}")
            R_final_map[i, j] = np.nan
            R_change_map[i, j] = np.nan
            localization_map[i, j] = np.nan
            energy_stability_map[i, j] = np.nan
            norm_map[i, j] = np.nan

            # Log failure to CSV
            with open(csv_file, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([sim_count, K, Delta, np.nan, np.nan, np.nan,
                               np.nan, np.nan, np.nan, f'FAILED: {e}',
                               time.time()])

        # Save checkpoint periodically
        if sim_count % checkpoint_interval == 0:
            np.savez_compressed(npz_checkpoint,
                              K_grid=K_grid, Delta_grid=Delta_grid,
                              R_final_map=R_final_map, R_change_map=R_change_map,
                              localization_map=localization_map,
                              energy_stability_map=energy_stability_map,
                              norm_map=norm_map,
                              progress=sim_count)
            print(f"  [Checkpoint saved: {sim_count}/{resolution**2}]", flush=True)

elapsed = time.time() - start_time
print(f"\n✓ Phase diagram scan complete in {elapsed/60:.1f} minutes")

# Save final complete data
final_data_file = os.path.join(data_dir, 'complete_data.npz')
np.savez_compressed(final_data_file,
                   K_grid=K_grid, Delta_grid=Delta_grid,
                   R_final_map=R_final_map, R_change_map=R_change_map,
                   localization_map=localization_map,
                   energy_stability_map=energy_stability_map,
                   norm_map=norm_map,
                   K_values=K_values, Delta_values=Delta_values,
                   resolution=resolution)
print(f"💾 Complete data saved: {final_data_file}")

# =============================================================================
# Classification: Dead Zone, Gas Phase, Bound States
# =============================================================================
print("\n" + "="*80)
print(" PHASE CLASSIFICATION")
print("="*80)

# Define phase boundaries
phase_map = np.zeros((resolution, resolution))
# 0 = Dead Zone (blue)
# 1 = Gas Phase (red)
# 2 = Bound States (gold)

for i in range(resolution):
    for j in range(resolution):
        R = R_final_map[i, j]
        L = localization_map[i, j]
        sigma_E = energy_stability_map[i, j]

        if np.isnan(R) or np.isnan(L) or np.isnan(sigma_E):
            phase_map[i, j] = -1  # Failed
        elif R < 0.05:
            phase_map[i, j] = 0  # Dead Zone: No synchronization
        elif R >= 0.05 and R < 0.08:
            phase_map[i, j] = 1  # Gas Phase: Weak sync, no localization
        elif R >= 0.08 and L >= 5 and sigma_E < 0.1:
            phase_map[i, j] = 2  # Bound State: Strong sync + localization + stable
        else:
            phase_map[i, j] = 1  # Default to gas phase

# Count phases
dead_count = np.sum(phase_map == 0)
gas_count = np.sum(phase_map == 1)
bound_count = np.sum(phase_map == 2)
failed_count = np.sum(phase_map == -1)

print(f"\nPhase Distribution:")
print(f"  Dead Zone (Blue):   {dead_count:3d} / {resolution**2} ({dead_count/resolution**2*100:.1f}%)")
print(f"  Gas Phase (Red):    {gas_count:3d} / {resolution**2} ({gas_count/resolution**2*100:.1f}%)")
print(f"  Bound States (Gold):{bound_count:3d} / {resolution**2} ({bound_count/resolution**2*100:.1f}%)")
print(f"  Failed:             {failed_count:3d} / {resolution**2}")

# Find bound state islands
if bound_count > 0:
    print(f"\n🌟 BOUND STATE ISLANDS FOUND:")
    bound_indices = np.where(phase_map == 2)
    for idx in range(len(bound_indices[0])):
        i_idx = bound_indices[0][idx]
        j_idx = bound_indices[1][idx]
        K = K_grid[i_idx, j_idx]
        Delta = Delta_grid[i_idx, j_idx]
        R = R_final_map[i_idx, j_idx]
        L = localization_map[i_idx, j_idx]
        print(f"  K={K:5.2f}, Δ={Delta:4.2f}  →  R={R:.3f}, L={L:.1f}")

# =============================================================================
# Visualization: 2D Phase Diagrams + 3D Surface Plots
# =============================================================================
print("\nGenerating visualizations...")

fig = plt.figure(figsize=(24, 18))
gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

# Row 1: Phase diagram + R_final + Localization (2D)
ax1 = fig.add_subplot(gs[0, 0])
colors = ['blue', 'red', 'gold', 'black']
cmap_phase = plt.matplotlib.colors.ListedColormap(colors[:3])
im1 = ax1.contourf(K_grid, Delta_grid, phase_map, levels=[-0.5, 0.5, 1.5, 2.5],
                   colors=['blue', 'red', 'gold'], alpha=0.7)
ax1.set_title('Phase Diagram: Map of Existence', fontweight='bold', fontsize=14)
ax1.set_xlabel('K (Vacuum Coherence)', fontsize=12)
ax1.set_ylabel('Δ (Mass Gap)', fontsize=12)
ax1.grid(True, alpha=0.3)
# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='blue', label='Dead Zone'),
                  Patch(facecolor='red', label='Gas Phase'),
                  Patch(facecolor='gold', label='Bound States')]
ax1.legend(handles=legend_elements, loc='upper left')

ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.contourf(K_grid, Delta_grid, R_final_map, levels=20, cmap='RdYlBu_r')
ax2.set_title('Mass Formation ⟨R⟩', fontweight='bold', fontsize=14)
ax2.set_xlabel('K', fontsize=12)
ax2.set_ylabel('Δ', fontsize=12)
plt.colorbar(im2, ax=ax2, label='⟨R⟩')
ax2.grid(True, alpha=0.3)

ax3 = fig.add_subplot(gs[0, 2])
im3 = ax3.contourf(K_grid, Delta_grid, localization_map, levels=20, cmap='viridis')
ax3.set_title('Particle Localization L', fontweight='bold', fontsize=14)
ax3.set_xlabel('K', fontsize=12)
ax3.set_ylabel('Δ', fontsize=12)
plt.colorbar(im3, ax=ax3, label='L')
ax3.grid(True, alpha=0.3)

# Row 2: 3D Surface Plots
ax4 = fig.add_subplot(gs[1, 0], projection='3d')
surf1 = ax4.plot_surface(K_grid, Delta_grid, R_final_map, cmap='RdYlBu_r',
                         edgecolor='none', alpha=0.8)
ax4.set_title('3D: Mass Formation', fontweight='bold', fontsize=12)
ax4.set_xlabel('K')
ax4.set_ylabel('Δ')
ax4.set_zlabel('⟨R⟩')
ax4.view_init(elev=25, azim=45)

ax5 = fig.add_subplot(gs[1, 1], projection='3d')
surf2 = ax5.plot_surface(K_grid, Delta_grid, localization_map, cmap='viridis',
                         edgecolor='none', alpha=0.8)
ax5.set_title('3D: Localization', fontweight='bold', fontsize=12)
ax5.set_xlabel('K')
ax5.set_ylabel('Δ')
ax5.set_zlabel('L')
ax5.view_init(elev=25, azim=45)

ax6 = fig.add_subplot(gs[1, 2], projection='3d')
surf3 = ax6.plot_surface(K_grid, Delta_grid, energy_stability_map, cmap='plasma',
                         edgecolor='none', alpha=0.8)
ax6.set_title('3D: Energy Stability', fontweight='bold', fontsize=12)
ax6.set_xlabel('K')
ax6.set_ylabel('Δ')
ax6.set_zlabel('σ_E / ⟨E⟩')
ax6.view_init(elev=25, azim=45)

# Row 3: Stability metrics
ax7 = fig.add_subplot(gs[2, 0])
im7 = ax7.contourf(K_grid, Delta_grid, R_change_map, levels=20, cmap='RdYlGn',
                   vmin=-0.2, vmax=0.2)
ax7.set_title('Synchronization Change ΔR', fontweight='bold', fontsize=14)
ax7.set_xlabel('K', fontsize=12)
ax7.set_ylabel('Δ', fontsize=12)
plt.colorbar(im7, ax=ax7, label='ΔR')
ax7.grid(True, alpha=0.3)

ax8 = fig.add_subplot(gs[2, 1])
im8 = ax8.contourf(K_grid, Delta_grid, energy_stability_map, levels=20, cmap='plasma')
ax8.set_title('Energy Stability σ_E / ⟨E⟩', fontweight='bold', fontsize=14)
ax8.set_xlabel('K', fontsize=12)
ax8.set_ylabel('Δ', fontsize=12)
plt.colorbar(im8, ax=ax8, label='σ_E / ⟨E⟩')
ax8.grid(True, alpha=0.3)

ax9 = fig.add_subplot(gs[2, 2])
im9 = ax9.contourf(K_grid, Delta_grid, norm_map, levels=20, cmap='coolwarm',
                   vmin=0.95, vmax=1.05)
ax9.set_title('Norm Conservation', fontweight='bold', fontsize=14)
ax9.set_xlabel('K', fontsize=12)
ax9.set_ylabel('Δ', fontsize=12)
plt.colorbar(im9, ax=ax9, label='Norm')
ax9.grid(True, alpha=0.3)

fig.suptitle('MSFT Phase Diagram: The Map of Existence', fontsize=18, fontweight='bold', y=0.995)

output = os.path.join(data_dir, 'phase_diagram_complete.png')
plt.savefig(output, dpi=150, bbox_inches='tight')
print(f"✓ Visualization saved: {output}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "="*80)
print(" PHASE DIAGRAM COMPLETE")
print("="*80)
print(f"\nScanned {resolution**2} universes in {elapsed/60:.1f} minutes")
print(f"\nDiscovered Phases:")
print(f"  • Dead Zone:    {dead_count} regions where vacuum fails")
print(f"  • Gas Phase:    {gas_count} regions with unconfined particles")
print(f"  • Bound States: {bound_count} regions with stable confinement")
print(f"\nNext Steps:")
if bound_count > 0:
    print(f"  ✓ Found {bound_count} bound state islands!")
    print(f"  → Use these (K, Δ) values in complete_MSFT_dirac_demo.py")
    print(f"  → Investigate if multiple islands = particle spectrum")
else:
    print(f"  ✗ No bound states found in this parameter range")
    print(f"  → Try finer resolution or different ranges")
print("="*80)
print(f"\n📁 All outputs saved to: {data_dir}")
print(f"  • Visualization: phase_diagram_complete.png")
print(f"  • Full data: complete_data.npz")
print(f"  • Progress log: scan_progress.csv")
print(f"  • Checkpoint: checkpoint.npz")
print("="*80)
