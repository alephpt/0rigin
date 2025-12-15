"""
Parameter Sweep for MSFT-Dirac Coupling

Systematically explores parameter space to find stable configurations:
- K (Kuramoto coupling): [1.0, 2.0, 3.0, 5.0, 8.0, 10.0]
- Spinor momentum: [0.0, 0.2, 0.5, 1.0]

Generates phase diagrams showing which parameters lead to:
1. Stable synchronization (R > 0.7)
2. Particle confinement (density centered, not in corners)
3. Smooth phase structure (no checkerboard aliasing)

This replaces manual parameter guessing with systematic exploration.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')
from kuramoto.field_theory import MSFTSystem
from kuramoto.field_theory.dirac.gamma_matrices import get_gamma_matrices_3plus1
from kuramoto.field_theory.dirac.spinor_field import DiracSpinorField, SpatialGrid

print("="*80)
print(" PARAMETER SWEEP: Finding Stable MSFT-Dirac Configurations")
print("="*80)
print("\nGoal: Systematically explore parameter space to LEARN system behavior")
print("      and find configurations with stable synchronization + confinement")
print("="*80)

# =============================================================================
# Parameter Grid
# =============================================================================
K_values = np.array([1.0, 2.0, 3.0, 5.0, 8.0, 10.0])
momentum_values = np.array([0.0, 0.2, 0.5, 1.0])

print(f"\nParameter Grid:")
print(f"  K (Kuramoto coupling): {K_values}")
print(f"  Spinor momentum p_x:   {momentum_values}")
print(f"  Total runs: {len(K_values)} × {len(momentum_values)} = {len(K_values) * len(momentum_values)}")

# Fixed parameters
grid_size = 32
n_steps = 100  # Shorter than full demo for speed
dt_MSFT = 0.01
dt_dirac = 0.0001
chiral_angle = 0.0  # Pure scalar mass

# =============================================================================
# Storage for results
# =============================================================================
results = {
    'K': [],
    'momentum': [],
    'R_initial': [],
    'R_final': [],
    'R_change': [],
    'R_stability': [],  # |R_final - R_target|, target=0.7
    'confinement': [],  # center_density / corner_density
    'phase_smoothness': [],  # 1 / std(∇phase)
    'norm_final': [],
}

# =============================================================================
# Run parameter sweep
# =============================================================================
print(f"\nRunning {len(K_values) * len(momentum_values)} simulations...")
print(f"  Simulation length: {n_steps} steps ({n_steps * dt_MSFT:.1f} time units)")
print(f"  Estimated runtime: ~{len(K_values) * len(momentum_values) * 0.5:.0f} seconds")
print()

start_time = time.time()
run_count = 0

for K in K_values:
    for p_x in momentum_values:
        run_count += 1
        print(f"[{run_count}/{len(K_values)*len(momentum_values)}] K={K:.1f}, p_x={p_x:.1f} ... ", end='', flush=True)

        # Initialize MSFT system
        MSFT = MSFTSystem(N=grid_size**2, grid_shape=(grid_size, grid_size), mass_gap=2.5)

        # CRITICAL: Override coupling strength
        MSFT.oscillators.coupling_strength = K

        # Seed synchronization pattern
        n = grid_size
        x_coord, y_coord = np.meshgrid(np.arange(n), np.arange(n))
        pattern = 0.3 + 0.4 * np.cos(2*np.pi*x_coord/n * 1.5) * np.cos(2*np.pi*y_coord/n * 1.5)
        MSFT.sync_field.values = pattern.flatten()

        # Initialize spinor
        grid = SpatialGrid(Nx=grid_size, Ny=grid_size, Lx=1.0, Ly=1.0, boundary='periodic')
        spinor = DiracSpinorField(grid, initial_state='plane_wave', momentum=(p_x, 0.0))

        # Record initial state
        R_initial = MSFT.sync_field.values.mean()

        # Evolve coupled system
        try:
            for step in range(n_steps):
                MSFT.step(dt=dt_MSFT)

                R_field = MSFT.sync_field.values.reshape(MSFT.grid_shape)
                m_field = MSFT.compute_effective_mass().reshape(MSFT.grid_shape)

                # Evolve spinor
                substeps = int(dt_MSFT / dt_dirac)
                for sub in range(substeps):
                    spinor.step_rk4(dt_dirac, m_field, chiral_angle=chiral_angle,
                                  Delta=MSFT.Delta, normalize=True)

            # Measure observables
            R_final = MSFT.sync_field.values.mean()
            R_change = R_final - R_initial
            R_stability = abs(R_final - 0.7)  # How far from target R=0.7

            density = spinor.compute_density()
            norm_final = np.sum(density) * grid.dx * grid.dy

            # Confinement: center vs corners
            center_mask = np.zeros_like(density, dtype=bool)
            center_mask[grid_size//4:3*grid_size//4, grid_size//4:3*grid_size//4] = True
            corner_mask = np.zeros_like(density, dtype=bool)
            corner_mask[:grid_size//8, :grid_size//8] = True
            corner_mask[:grid_size//8, -grid_size//8:] = True
            corner_mask[-grid_size//8:, :grid_size//8] = True
            corner_mask[-grid_size//8:, -grid_size//8:] = True

            center_density = np.mean(density[center_mask])
            corner_density = np.mean(density[corner_mask])
            confinement = center_density / (corner_density + 1e-10)

            # Phase smoothness: inverse of gradient variance
            phase = np.angle(spinor.psi[:, :, 0])
            grad_x = np.gradient(phase, axis=0)
            grad_y = np.gradient(phase, axis=1)
            grad_mag = np.sqrt(grad_x**2 + grad_y**2)
            phase_smoothness = 1.0 / (np.std(grad_mag) + 1e-10)

            # Store results
            results['K'].append(K)
            results['momentum'].append(p_x)
            results['R_initial'].append(R_initial)
            results['R_final'].append(R_final)
            results['R_change'].append(R_change)
            results['R_stability'].append(R_stability)
            results['confinement'].append(confinement)
            results['phase_smoothness'].append(phase_smoothness)
            results['norm_final'].append(norm_final)

            print(f"R: {R_initial:.3f}→{R_final:.3f}, conf: {confinement:.1f}, norm: {norm_final:.4f}")

        except Exception as e:
            print(f"FAILED: {e}")
            # Fill with NaN
            for key in results.keys():
                if key not in ['K', 'momentum']:
                    results[key].append(np.nan)
                elif key == 'K':
                    results['K'].append(K)
                elif key == 'momentum':
                    results['momentum'].append(p_x)

elapsed = time.time() - start_time
print(f"\n✓ Sweep complete in {elapsed:.1f} seconds ({elapsed/run_count:.2f} s/run)")

# Convert to arrays
for key in results.keys():
    results[key] = np.array(results[key])

# =============================================================================
# Analysis: Find best configurations
# =============================================================================
print("\n" + "="*80)
print(" ANALYSIS: Best Configurations")
print("="*80)

# Scoring: Want R_final>0.7, confinement>5, phase_smoothness>1
score = np.zeros(len(results['K']))
for i in range(len(results['K'])):
    if not np.isnan(results['R_final'][i]):
        # R target: closer to 0.7 is better
        r_score = max(0, 1 - abs(results['R_final'][i] - 0.7)/0.7)
        # Confinement: higher is better
        c_score = min(1, results['confinement'][i] / 10.0)
        # Phase smoothness: higher is better
        p_score = min(1, results['phase_smoothness'][i] / 2.0)

        score[i] = r_score + c_score + p_score
    else:
        score[i] = -1

# Top 3 configurations
top_indices = np.argsort(-score)[:3]

print("\nTop 3 Configurations:")
for rank, idx in enumerate(top_indices, 1):
    if score[idx] > 0:
        print(f"\n{rank}. K={results['K'][idx]:.1f}, p_x={results['momentum'][idx]:.1f}")
        print(f"   Score: {score[idx]:.2f}/3.0")
        print(f"   R: {results['R_initial'][idx]:.3f} → {results['R_final'][idx]:.3f} (change: {results['R_change'][idx]:+.3f})")
        print(f"   Confinement: {results['confinement'][idx]:.2f} (want >5)")
        print(f"   Phase smoothness: {results['phase_smoothness'][idx]:.2f} (want >1)")
        print(f"   Norm: {results['norm_final'][idx]:.4f}")

# =============================================================================
# Visualization: Parameter Space Phase Diagrams
# =============================================================================
print("\nGenerating phase diagrams...")

# Reshape data for heatmaps
def reshape_for_heatmap(values):
    heatmap = np.full((len(momentum_values), len(K_values)), np.nan)
    for i in range(len(results['K'])):
        k_idx = np.where(K_values == results['K'][i])[0][0]
        p_idx = np.where(momentum_values == results['momentum'][i])[0][0]
        heatmap[p_idx, k_idx] = values[i]
    return heatmap

R_final_heatmap = reshape_for_heatmap(results['R_final'])
confinement_heatmap = reshape_for_heatmap(results['confinement'])
R_change_heatmap = reshape_for_heatmap(results['R_change'])
score_heatmap = reshape_for_heatmap(score)

fig = plt.figure(figsize=(20, 12))
gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

# 1. Final R
ax1 = fig.add_subplot(gs[0, 0])
im1 = ax1.imshow(R_final_heatmap, aspect='auto', cmap='RdYlBu_r', vmin=0, vmax=1, origin='lower')
ax1.set_title('Final Synchronization ⟨R⟩', fontweight='bold', fontsize=14)
ax1.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax1.set_ylabel('Spinor Momentum p_x', fontsize=12)
ax1.set_xticks(range(len(K_values)))
ax1.set_xticklabels([f'{k:.1f}' for k in K_values])
ax1.set_yticks(range(len(momentum_values)))
ax1.set_yticklabels([f'{p:.1f}' for p in momentum_values])
plt.colorbar(im1, ax=ax1, label='⟨R⟩')
# Mark target R=0.7
for i, k in enumerate(K_values):
    for j, p in enumerate(momentum_values):
        val = R_final_heatmap[j, i]
        if not np.isnan(val) and abs(val - 0.7) < 0.1:
            ax1.plot(i, j, 'g*', markersize=15)

# 2. Confinement
ax2 = fig.add_subplot(gs[0, 1])
im2 = ax2.imshow(confinement_heatmap, aspect='auto', cmap='viridis', vmin=0, vmax=15, origin='lower')
ax2.set_title('Particle Confinement (center/corner density)', fontweight='bold', fontsize=14)
ax2.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax2.set_ylabel('Spinor Momentum p_x', fontsize=12)
ax2.set_xticks(range(len(K_values)))
ax2.set_xticklabels([f'{k:.1f}' for k in K_values])
ax2.set_yticks(range(len(momentum_values)))
ax2.set_yticklabels([f'{p:.1f}' for p in momentum_values])
plt.colorbar(im2, ax=ax2, label='Confinement')
# Mark well-confined (>5)
for i, k in enumerate(K_values):
    for j, p in enumerate(momentum_values):
        val = confinement_heatmap[j, i]
        if not np.isnan(val) and val > 5:
            ax2.plot(i, j, 'w*', markersize=15)

# 3. R change (stability indicator)
ax3 = fig.add_subplot(gs[1, 0])
im3 = ax3.imshow(R_change_heatmap, aspect='auto', cmap='RdYlGn', vmin=-0.1, vmax=0.1, origin='lower')
ax3.set_title('Synchronization Change ΔR (positive = growing)', fontweight='bold', fontsize=14)
ax3.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax3.set_ylabel('Spinor Momentum p_x', fontsize=12)
ax3.set_xticks(range(len(K_values)))
ax3.set_xticklabels([f'{k:.1f}' for k in K_values])
ax3.set_yticks(range(len(momentum_values)))
ax3.set_yticklabels([f'{p:.1f}' for p in momentum_values])
plt.colorbar(im3, ax=ax3, label='ΔR')
# Mark stable (small |ΔR|)
for i, k in enumerate(K_values):
    for j, p in enumerate(momentum_values):
        val = R_change_heatmap[j, i]
        if not np.isnan(val) and abs(val) < 0.01 and R_final_heatmap[j, i] > 0.5:
            ax3.plot(i, j, 'k*', markersize=15)

# 4. Overall score
ax4 = fig.add_subplot(gs[1, 1])
im4 = ax4.imshow(score_heatmap, aspect='auto', cmap='hot', vmin=0, vmax=3, origin='lower')
ax4.set_title('Overall Score (R target + confinement + smoothness)', fontweight='bold', fontsize=14)
ax4.set_xlabel('Kuramoto Coupling K', fontsize=12)
ax4.set_ylabel('Spinor Momentum p_x', fontsize=12)
ax4.set_xticks(range(len(K_values)))
ax4.set_xticklabels([f'{k:.1f}' for k in K_values])
ax4.set_yticks(range(len(momentum_values)))
ax4.set_yticklabels([f'{p:.1f}' for p in momentum_values])
plt.colorbar(im4, ax=ax4, label='Score')
# Mark top 3
for rank, idx in enumerate(top_indices[:3], 1):
    if score[idx] > 0:
        k_idx = np.where(K_values == results['K'][idx])[0][0]
        p_idx = np.where(momentum_values == results['momentum'][idx])[0][0]
        ax4.text(k_idx, p_idx, str(rank), ha='center', va='center',
                color='white', fontweight='bold', fontsize=14,
                bbox=dict(boxstyle='circle', fc='blue', alpha=0.7))

fig.suptitle('MSFT-Dirac Parameter Space Exploration', fontsize=16, fontweight='bold', y=0.995)

output = '/tmp/parameter_sweep_phase_diagrams.png'
plt.savefig(output, dpi=150, bbox_inches='tight')
print(f"✓ Saved: {output}")
plt.close()

# =============================================================================
# Summary
# =============================================================================
print("\n" + "="*80)
print(" RESULTS SUMMARY")
print("="*80)
print(f"\nExplored {len(K_values)} × {len(momentum_values)} = {len(results['K'])} configurations")
print(f"Runtime: {elapsed:.1f} seconds")
print(f"\nKey Findings:")
print(f"  Best score: {score.max():.2f}/3.0")
best_idx = np.argmax(score)
print(f"  Best config: K={results['K'][best_idx]:.1f}, p_x={results['momentum'][best_idx]:.1f}")
print(f"    → R={results['R_final'][best_idx]:.3f}, confinement={results['confinement'][best_idx]:.1f}")
print(f"\nNext steps:")
print(f"  1. Use best parameters in complete_MSFT_dirac_demo.py")
print(f"  2. Run full 100-step evolution with visualization")
print(f"  3. If still unstable, refine search around best config")
print("="*80)
print(f"\nOutput: {output}")
print("="*80)
