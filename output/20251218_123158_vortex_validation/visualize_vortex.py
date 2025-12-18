#!/usr/bin/env python3
"""
Visualize vortex operator splitting results
- Spatial R(x,y) field snapshots showing vortex structure
- Comparison between N=1 and N=100
- Timeseries analysis
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

def load_timeseries(filename):
    """Load timeseries CSV data"""
    data = {
        'step': [],
        'time': [],
        'R_avg': [],
        'R_std': [],
        'norm': [],
        'energy': [],
        'max_density': []
    }

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data['step'].append(int(row['step']))
            data['time'].append(float(row['time']))
            data['R_avg'].append(float(row['R_avg']))
            data['R_std'].append(float(row['R_std']))
            data['norm'].append(float(row['norm']))
            data['energy'].append(float(row['energy']))
            data['max_density'].append(float(row['max_density']))

    return {k: np.array(v) for k, v in data.items()}

def load_R_snapshots(filename, Nx=64):
    """Load R field spatial snapshots"""
    snapshots = []
    times = []

    with open(filename, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header

        for row in reader:
            step = int(row[0])
            time = float(row[1])
            R_field = np.array([float(x) for x in row[2:]])
            R_field = R_field.reshape((Nx, Nx))

            snapshots.append(R_field)
            times.append(time)

    return times, snapshots

# Load data
print("Loading N=1 data...")
n1_ts = load_timeseries('vortex_N1/timeseries.csv')
n1_times, n1_snapshots = load_R_snapshots('vortex_N1/R_field_snapshots.csv')

print("Loading N=100 data...")
n100_ts = load_timeseries('vortex_N100/timeseries.csv')
n100_times, n100_snapshots = load_R_snapshots('vortex_N100/R_field_snapshots.csv')

print(f"Loaded {len(n1_snapshots)} snapshots for N=1")
print(f"Loaded {len(n100_snapshots)} snapshots for N=100")

# Create figure: 3x4 grid
# Row 1: R(x,y) at t=0, 50, 100 for N=1
# Row 2: R(x,y) at t=0, 50, 100 for N=100
# Row 3: Difference maps and timeseries

fig = plt.figure(figsize=(16, 12))

# Find snapshot indices closest to t=0, 50, 100
snapshot_times = [0, 50, 100]
n1_indices = []
n100_indices = []

for target_t in snapshot_times:
    n1_idx = np.argmin(np.abs(np.array(n1_times) - target_t))
    n100_idx = np.argmin(np.abs(np.array(n100_times) - target_t))
    n1_indices.append(n1_idx)
    n100_indices.append(n100_idx)
    print(f"t={target_t}: N=1 snapshot at t={n1_times[n1_idx]:.1f}, N=100 at t={n100_times[n100_idx]:.1f}")

# Row 1: N=1 snapshots
for i, (idx, t) in enumerate(zip(n1_indices, snapshot_times)):
    ax = plt.subplot(3, 4, i+1)
    im = ax.imshow(n1_snapshots[idx], cmap='viridis', origin='lower', vmin=0, vmax=1)
    ax.set_title(f'N=1: R(x,y) at t={n1_times[idx]:.1f}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# Row 2: N=100 snapshots
for i, (idx, t) in enumerate(zip(n100_indices, snapshot_times)):
    ax = plt.subplot(3, 4, i+5)
    im = ax.imshow(n100_snapshots[idx], cmap='viridis', origin='lower', vmin=0, vmax=1)
    ax.set_title(f'N=100: R(x,y) at t={n100_times[idx]:.1f}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# Row 3, Col 1: Difference map at t=100
ax = plt.subplot(3, 4, 9)
diff = n1_snapshots[n1_indices[2]] - n100_snapshots[n100_indices[2]]
im = ax.imshow(diff, cmap='RdBu_r', origin='lower', vmin=-0.05, vmax=0.05)
ax.set_title(f'Difference: N=1 - N=100 at t=100')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='ΔR')

# Row 3, Col 2: R_avg timeseries
ax = plt.subplot(3, 4, 10)
ax.plot(n1_ts['time'], n1_ts['R_avg'], 'b-', label='N=1', alpha=0.7)
ax.plot(n100_ts['time'], n100_ts['R_avg'], 'r-', label='N=100', alpha=0.7)
ax.set_xlabel('Time')
ax.set_ylabel('R_avg')
ax.set_title('Average Synchronization')
ax.legend()
ax.grid(True, alpha=0.3)

# Row 3, Col 3: R_std timeseries
ax = plt.subplot(3, 4, 11)
ax.plot(n1_ts['time'], n1_ts['R_std'], 'b-', label='N=1', alpha=0.7)
ax.plot(n100_ts['time'], n100_ts['R_std'], 'r-', label='N=100', alpha=0.7)
ax.set_xlabel('Time')
ax.set_ylabel('σ_R')
ax.set_title('R Field Fluctuations')
ax.legend()
ax.grid(True, alpha=0.3)

# Row 3, Col 4: Energy timeseries
ax = plt.subplot(3, 4, 12)
ax.plot(n1_ts['time'], n1_ts['energy'], 'b-', label='N=1', alpha=0.7)
ax.plot(n100_ts['time'], n100_ts['energy'], 'r-', label='N=100', alpha=0.7)
ax.set_xlabel('Time')
ax.set_ylabel('Energy')
ax.set_title('Dirac Energy')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('vortex_spatial_analysis.png', dpi=150, bbox_inches='tight')
print("\n✓ Saved: vortex_spatial_analysis.png")

# Create second figure: Radial profiles showing vortex structure
fig2, axes = plt.subplots(2, 2, figsize=(12, 10))

def radial_profile(data, center=(32, 32)):
    """Compute radial average of 2D field"""
    y, x = np.indices(data.shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr

    return radialprofile

# Radial profiles at different times
for i, (idx_n1, idx_n100, t) in enumerate(zip(n1_indices[:3], n100_indices[:3], snapshot_times)):
    row = i // 2
    col = i % 2
    ax = axes[row, col]

    r_n1 = radial_profile(n1_snapshots[idx_n1])
    r_n100 = radial_profile(n100_snapshots[idx_n100])

    radius = np.arange(len(r_n1))

    ax.plot(radius, r_n1, 'b-', label='N=1', linewidth=2, alpha=0.7)
    ax.plot(radius, r_n100, 'r-', label='N=100', linewidth=2, alpha=0.7)
    ax.set_xlabel('Radius from center')
    ax.set_ylabel('R(r)')
    ax.set_title(f'Radial Profile at t={t}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 32)

# Fourth plot: Time evolution of vortex core
ax = axes[1, 1]
core_n1 = [snap[32, 32] for snap in n1_snapshots]
core_n100 = [snap[32, 32] for snap in n100_snapshots]

ax.plot(n1_times, core_n1, 'b-', label='N=1 (core)', linewidth=2, alpha=0.7)
ax.plot(n100_times, core_n100, 'r-', label='N=100 (core)', linewidth=2, alpha=0.7)
ax.set_xlabel('Time')
ax.set_ylabel('R at vortex core')
ax.set_title('Vortex Core Evolution')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('vortex_radial_analysis.png', dpi=150, bbox_inches='tight')
print("✓ Saved: vortex_radial_analysis.png")

# Print statistics
print("\n" + "="*60)
print("VORTEX STRUCTURE ANALYSIS")
print("="*60)

print("\nInitial R field (t=0):")
print(f"  N=1:   R_core = {n1_snapshots[0][32,32]:.4f}, R_edge = {n1_snapshots[0][0,0]:.4f}")
print(f"  N=100: R_core = {n100_snapshots[0][32,32]:.4f}, R_edge = {n100_snapshots[0][0,0]:.4f}")

print("\nFinal R field (t=100):")
idx_final_n1 = n1_indices[2]
idx_final_n100 = n100_indices[2]
print(f"  N=1:   R_core = {n1_snapshots[idx_final_n1][32,32]:.4f}, R_edge = {n1_snapshots[idx_final_n1][0,0]:.4f}")
print(f"  N=100: R_core = {n100_snapshots[idx_final_n100][32,32]:.4f}, R_edge = {n100_snapshots[idx_final_n100][0,0]:.4f}")

print("\nVortex core difference (N=1 - N=100):")
core_diff_t0 = n1_snapshots[0][32,32] - n100_snapshots[0][32,32]
core_diff_t100 = n1_snapshots[idx_final_n1][32,32] - n100_snapshots[idx_final_n100][32,32]
print(f"  t=0:   ΔR_core = {core_diff_t0:.6f}")
print(f"  t=100: ΔR_core = {core_diff_t100:.6f}")

print("\nR field statistics (final):")
print(f"  N=1:   R_avg = {n1_ts['R_avg'][-1]:.6f} ± {n1_ts['R_std'][-1]:.6f}")
print(f"  N=100: R_avg = {n100_ts['R_avg'][-1]:.6f} ± {n100_ts['R_std'][-1]:.6f}")
print(f"  Smoothing: {100*(n1_ts['R_std'][-1] - n100_ts['R_std'][-1])/n1_ts['R_std'][-1]:.2f}% reduction")

print("\nEnergy statistics (final):")
print(f"  N=1:   E = {n1_ts['energy'][-1]:.5f}")
print(f"  N=100: E = {n100_ts['energy'][-1]:.5f}")
print(f"  ΔE = {abs(n1_ts['energy'][-1] - n100_ts['energy'][-1]):.5f} ({100*abs(n1_ts['energy'][-1] - n100_ts['energy'][-1])/n100_ts['energy'][-1]:.3f}%)")

print("\n" + "="*60)
print("✓ Visualization complete!")
print("="*60)

plt.show()
