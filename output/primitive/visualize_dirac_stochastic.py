#!/usr/bin/env python3
"""
MSFT Dirac-Kuramoto Stochastic Evolution Visualization Tool

Visualizes Dirac spinor coupling with stochastic Kuramoto dynamics:
- Timeseries: R_global, spinor_norm, particle position, drift
- Spatial evolution: density |Ψ|², mass field m(x,y), spinor components
- Particle tracking: center of mass trajectory
- Field animations: snapshots at different times

Data structure:
  dirac_evolution/
    timeseries.dat          - Global metrics vs time
    summary.dat             - Final statistics
    density_XXX.dat         - Spinor density |Ψ|² snapshots
    mass_field_XXX.dat      - Emergent mass field snapshots
    snapshot_XXX.dat        - Full field state (theta, R, mass, Ψ components)

  stochastic_evolution_*.dat - Alternative format with particle tracking
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm, TwoSlopeNorm
import argparse
import os
from pathlib import Path
import glob
import re


def load_timeseries(filepath):
    """
    Load timeseries.dat file.

    Returns:
        dict: time, R_global, spinor_norm, particle_x, particle_y, particle_drift
    """
    try:
        data = np.loadtxt(filepath, comments='#')

        if data.ndim == 1:
            data = data.reshape(1, -1)

        return {
            'time': data[:, 0],
            'R_global': data[:, 1],
            'spinor_norm': data[:, 2],
            'particle_x': data[:, 3],
            'particle_y': data[:, 4],
            'particle_drift': data[:, 5]
        }
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_stochastic_evolution(filepath):
    """
    Load stochastic_evolution_*.dat file (alternative format).

    Returns:
        dict: Time, R_global, Norm, X_pos, Y_pos
    """
    try:
        data = np.loadtxt(filepath, comments='#')

        if data.ndim == 1:
            data = data.reshape(1, -1)

        return {
            'time': data[:, 0],
            'R_global': data[:, 1],
            'norm': data[:, 2],
            'x_pos': data[:, 3],
            'y_pos': data[:, 4]
        }
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_spatial_field(filepath):
    """
    Load spatial field (density or mass_field).

    Returns:
        (x, y, field) as 2D arrays
    """
    try:
        data = np.loadtxt(filepath, comments='#')

        # Get grid size
        x = data[:, 0].astype(int)
        y = data[:, 1].astype(int)
        values = data[:, 2]

        grid_size = int(np.sqrt(len(data)))

        field = values.reshape((grid_size, grid_size))

        return field
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_snapshot(filepath):
    """
    Load full field snapshot.

    Returns:
        dict: theta, R, mass, density, psi (4 complex components)
    """
    try:
        data = np.loadtxt(filepath, comments='#')

        grid_size = int(np.sqrt(len(data)))

        result = {
            'theta': data[:, 2].reshape((grid_size, grid_size)),
            'R': data[:, 3].reshape((grid_size, grid_size)),
            'mass': data[:, 4].reshape((grid_size, grid_size)),
            'density': data[:, 5].reshape((grid_size, grid_size))
        }

        # Spinor components (4 complex values = 8 real numbers)
        if data.shape[1] >= 14:
            psi = np.zeros((4, grid_size, grid_size), dtype=complex)
            for i in range(4):
                real_col = 6 + 2*i
                imag_col = 7 + 2*i
                psi[i] = (data[:, real_col] + 1j*data[:, imag_col]).reshape((grid_size, grid_size))
            result['psi'] = psi

        return result
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def plot_timeseries_evolution(data_dir, output_file=None, dpi=150):
    """
    Plot temporal evolution from timeseries.dat.
    """
    filepath = os.path.join(data_dir, 'timeseries.dat')

    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        return

    ts = load_timeseries(filepath)

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.35)

    fig.suptitle('Dirac-Kuramoto Stochastic Evolution: Timeseries',
                 fontsize=14, fontweight='bold')

    # Plot 1: R_global vs time
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(ts['time'], ts['R_global'], 'b-', linewidth=2)
    ax1.set_xlabel('Time', fontsize=11)
    ax1.set_ylabel('R_global (Synchronization)', fontsize=11)
    ax1.set_title('Kuramoto Order Parameter', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0.999, color='red', linestyle='--', alpha=0.5, label='R = 0.999')
    ax1.legend()

    # Plot 2: Spinor norm conservation
    ax2 = fig.add_subplot(gs[0, 1])
    norm_dev = np.abs(ts['spinor_norm'] - 1.0) * 100  # Percent deviation
    ax2.plot(ts['time'], norm_dev, 'g-', linewidth=2)
    ax2.set_xlabel('Time', fontsize=11)
    ax2.set_ylabel('|Norm - 1| (%)', fontsize=11)
    ax2.set_title('Spinor Norm Conservation', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='1% tolerance')
    ax2.legend()

    # Plot 3: Particle drift
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(ts['time'], ts['particle_drift'], 'purple', linewidth=2)
    ax3.set_xlabel('Time', fontsize=11)
    ax3.set_ylabel('Particle Drift (grid units)', fontsize=11)
    ax3.set_title('Center of Mass Drift', fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(5.0, color='red', linestyle='--', alpha=0.5, label='5 unit threshold')
    ax3.legend()

    # Plot 4: Particle trajectory (2D)
    ax4 = fig.add_subplot(gs[1, 0])

    # Color by time
    scatter = ax4.scatter(ts['particle_x'], ts['particle_y'],
                         c=ts['time'], cmap='viridis', s=30, alpha=0.7)

    # Add trajectory arrows
    n_arrows = min(10, len(ts['time']) - 1)
    step = max(1, len(ts['time']) // n_arrows)
    for i in range(0, len(ts['time']) - step, step):
        ax4.annotate('', xy=(ts['particle_x'][i+step], ts['particle_y'][i+step]),
                    xytext=(ts['particle_x'][i], ts['particle_y'][i]),
                    arrowprops=dict(arrowstyle='->', color='red', alpha=0.5, lw=1.5))

    # Mark start and end
    ax4.plot(ts['particle_x'][0], ts['particle_y'][0], 'go', markersize=12,
            label=f'Start ({ts["particle_x"][0]:.1f}, {ts["particle_y"][0]:.1f})')
    ax4.plot(ts['particle_x'][-1], ts['particle_y'][-1], 'ro', markersize=12,
            label=f'End ({ts["particle_x"][-1]:.1f}, {ts["particle_y"][-1]:.1f})')

    ax4.set_xlabel('X Position (grid units)', fontsize=11)
    ax4.set_ylabel('Y Position (grid units)', fontsize=11)
    ax4.set_title('Particle Trajectory (Center of Mass)', fontweight='bold')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3)
    ax4.set_aspect('equal')

    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Time')

    # Plot 5: Velocity components
    ax5 = fig.add_subplot(gs[1, 1])

    if len(ts['time']) > 1:
        dt = np.diff(ts['time'])
        vx = np.diff(ts['particle_x']) / dt
        vy = np.diff(ts['particle_y']) / dt
        t_vel = ts['time'][:-1]

        ax5.plot(t_vel, vx, 'b-', linewidth=1.5, alpha=0.7, label='v_x')
        ax5.plot(t_vel, vy, 'r-', linewidth=1.5, alpha=0.7, label='v_y')
        ax5.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax5.set_xlabel('Time', fontsize=11)
        ax5.set_ylabel('Velocity (grid units/time)', fontsize=11)
        ax5.set_title('Particle Velocity Components', fontweight='bold')
        ax5.legend()
        ax5.grid(True, alpha=0.3)

    # Plot 6: Statistics summary
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    # Calculate statistics
    R_initial = ts['R_global'][0]
    R_final = ts['R_global'][-1]
    R_mean = ts['R_global'].mean()
    R_std = ts['R_global'].std()

    norm_initial = ts['spinor_norm'][0]
    norm_final = ts['spinor_norm'][-1]
    norm_mean = ts['spinor_norm'].mean()

    drift_final = ts['particle_drift'][-1]
    drift_max = ts['particle_drift'].max()

    stats_text = f"""
    SUMMARY STATISTICS

    Kuramoto Synchronization:
      R_initial:  {R_initial:.6f}
      R_final:    {R_final:.6f}
      R_mean:     {R_mean:.6f}
      R_std:      {R_std:.6f}

    Spinor Norm Conservation:
      Norm_initial: {norm_initial:.6f}
      Norm_final:   {norm_final:.6f}
      Norm_mean:    {norm_mean:.6f}
      Deviation:    {abs(norm_mean-1)*100:.3f}%

    Particle Drift:
      Final drift:  {drift_final:.3f} units
      Max drift:    {drift_max:.3f} units

    Time range: [{ts['time'][0]:.2f}, {ts['time'][-1]:.2f}]
    Total steps: {len(ts['time'])}
    """

    ax6.text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
            verticalalignment='center', transform=ax6.transAxes)

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def plot_spatial_evolution(data_dir, output_file=None, dpi=150):
    """
    Plot spatial field evolution (density and mass fields at different times).
    """
    # Find all density and mass field files
    density_files = sorted(glob.glob(os.path.join(data_dir, 'density_*.dat')))
    mass_files = sorted(glob.glob(os.path.join(data_dir, 'mass_field_*.dat')))

    if not density_files or not mass_files:
        print("Error: No spatial field files found")
        return

    # Select subset for visualization (first, middle, last)
    n_snapshots = min(3, len(density_files))
    indices = np.linspace(0, len(density_files)-1, n_snapshots, dtype=int)

    fig = plt.figure(figsize=(18, 6 * n_snapshots))
    gs = GridSpec(n_snapshots, 3, figure=fig, hspace=0.3, wspace=0.3)

    fig.suptitle('Dirac-Kuramoto Spatial Field Evolution',
                 fontsize=14, fontweight='bold', y=0.995)

    for row, idx in enumerate(indices):
        # Extract time index from filename
        time_idx = int(re.search(r'_(\d+)\.dat', density_files[idx]).group(1))

        # Load fields
        density = load_spatial_field(density_files[idx])
        mass = load_spatial_field(mass_files[idx])

        # Plot density
        ax1 = fig.add_subplot(gs[row, 0])

        # Use log scale if density has large dynamic range
        if density.max() / (density.min() + 1e-10) > 100:
            im1 = ax1.imshow(density, origin='lower', cmap='hot',
                           norm=LogNorm(vmin=density[density>0].min() if (density>0).any() else 1e-10,
                                       vmax=density.max()),
                           interpolation='bilinear')
        else:
            im1 = ax1.imshow(density, origin='lower', cmap='hot', interpolation='bilinear')

        ax1.set_title(f'Spinor Density |Ψ|² (t={time_idx})', fontweight='bold')
        ax1.set_xlabel('X (grid units)')
        ax1.set_ylabel('Y (grid units)')
        plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

        # Plot mass field
        ax2 = fig.add_subplot(gs[row, 1])
        im2 = ax2.imshow(mass, origin='lower', cmap='viridis', interpolation='bilinear')
        ax2.set_title(f'Mass Field m(x,y) = Δ·R(x,y) (t={time_idx})', fontweight='bold')
        ax2.set_xlabel('X (grid units)')
        ax2.set_ylabel('Y (grid units)')
        plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

        # Plot mass gradient magnitude (gravity proxy)
        ax3 = fig.add_subplot(gs[row, 2])
        grad_y, grad_x = np.gradient(mass)
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)
        im3 = ax3.imshow(grad_mag, origin='lower', cmap='plasma', interpolation='bilinear')
        ax3.set_title(f'|∇m| Gradient Magnitude (t={time_idx})', fontweight='bold')
        ax3.set_xlabel('X (grid units)')
        ax3.set_ylabel('Y (grid units)')
        plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)

    plt.tight_layout()

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def plot_spinor_components(data_dir, snapshot_idx=0, output_file=None, dpi=150):
    """
    Plot full spinor field components from snapshot.
    """
    snapshot_file = os.path.join(data_dir, f'snapshot_{snapshot_idx:03d}.dat')

    if not os.path.exists(snapshot_file):
        print(f"Error: {snapshot_file} not found")
        return

    snap = load_snapshot(snapshot_file)

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.35)

    fig.suptitle(f'Dirac Spinor Field Components (Snapshot {snapshot_idx})',
                 fontsize=14, fontweight='bold')

    # Row 1: Kuramoto fields
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(snap['theta'], origin='lower', cmap='twilight', interpolation='bilinear')
    ax1.set_title('Phase θ(x,y)', fontweight='bold')
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(snap['R'], origin='lower', cmap='viridis', interpolation='bilinear')
    ax2.set_title('Sync Order R(x,y)', fontweight='bold')
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(snap['mass'], origin='lower', cmap='plasma', interpolation='bilinear')
    ax3.set_title('Mass m = Δ·R', fontweight='bold')
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)

    ax4 = fig.add_subplot(gs[0, 3])
    im4 = ax4.imshow(snap['density'], origin='lower', cmap='hot', interpolation='bilinear')
    ax4.set_title('Density |Ψ|²', fontweight='bold')
    plt.colorbar(im4, ax=ax4, fraction=0.046, pad=0.04)

    # Rows 2-3: Spinor components (if available)
    if 'psi' in snap:
        psi = snap['psi']

        for i in range(4):
            # Real part
            ax_real = fig.add_subplot(gs[1, i])

            # Symmetric colormap around zero
            vmax_real = max(abs(psi[i].real.min()), abs(psi[i].real.max()))
            if vmax_real > 0:
                im_real = ax_real.imshow(psi[i].real, origin='lower', cmap='RdBu_r',
                                        vmin=-vmax_real, vmax=vmax_real,
                                        interpolation='bilinear')
            else:
                im_real = ax_real.imshow(psi[i].real, origin='lower', cmap='RdBu_r',
                                        interpolation='bilinear')

            ax_real.set_title(f'Re(Ψ_{i})', fontweight='bold')
            plt.colorbar(im_real, ax=ax_real, fraction=0.046, pad=0.04)

            # Imaginary part
            ax_imag = fig.add_subplot(gs[2, i])

            # Symmetric colormap around zero
            vmax_imag = max(abs(psi[i].imag.min()), abs(psi[i].imag.max()))
            if vmax_imag > 0:
                im_imag = ax_imag.imshow(psi[i].imag, origin='lower', cmap='RdBu_r',
                                        vmin=-vmax_imag, vmax=vmax_imag,
                                        interpolation='bilinear')
            else:
                im_imag = ax_imag.imshow(psi[i].imag, origin='lower', cmap='RdBu_r',
                                        interpolation='bilinear')

            ax_imag.set_title(f'Im(Ψ_{i})', fontweight='bold')
            plt.colorbar(im_imag, ax=ax_imag, fraction=0.046, pad=0.04)

    plt.tight_layout()

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def analyze_dirac_evolution(data_dir):
    """
    Comprehensive text analysis of Dirac-Kuramoto stochastic evolution.
    """
    print("\n" + "="*70)
    print("DIRAC-KURAMOTO STOCHASTIC EVOLUTION ANALYSIS")
    print("="*70)

    # Load timeseries
    ts_file = os.path.join(data_dir, 'timeseries.dat')
    if not os.path.exists(ts_file):
        print(f"Error: {ts_file} not found")
        return

    ts = load_timeseries(ts_file)

    # Load summary if available
    summary_file = os.path.join(data_dir, 'summary.dat')
    if os.path.exists(summary_file):
        with open(summary_file, 'r') as f:
            summary_text = f.read()
        print("\n--- SIMULATION PARAMETERS ---")
        for line in summary_text.split('\n'):
            if line and not line.startswith('#'):
                print(f"  {line}")

    # === KURAMOTO SYNCHRONIZATION ANALYSIS ===
    print("\n" + "-"*70)
    print("1. KURAMOTO SYNCHRONIZATION ANALYSIS")
    print("-"*70)

    R_initial = ts['R_global'][0]
    R_final = ts['R_global'][-1]
    R_mean = ts['R_global'].mean()
    R_std = ts['R_global'].std()
    R_min = ts['R_global'].min()

    print(f"\n  Order Parameter R:")
    print(f"    Initial:     {R_initial:.8f}")
    print(f"    Final:       {R_final:.8f}")
    print(f"    Mean:        {R_mean:.8f}")
    print(f"    Std Dev:     {R_std:.8f}")
    print(f"    Min:         {R_min:.8f}")
    print(f"    Change:      {R_final - R_initial:+.8f}")

    # Assess synchronization quality
    if R_min > 0.999:
        sync_quality = "EXCELLENT (R > 0.999 maintained)"
    elif R_min > 0.99:
        sync_quality = "VERY GOOD (R > 0.99 maintained)"
    elif R_min > 0.95:
        sync_quality = "GOOD (R > 0.95 maintained)"
    elif R_min > 0.8:
        sync_quality = "MODERATE (R > 0.8)"
    else:
        sync_quality = "WEAK (R < 0.8)"

    print(f"\n  Synchronization Quality: {sync_quality}")

    # Desynchronization rate
    if len(ts['time']) > 1:
        total_time = ts['time'][-1] - ts['time'][0]
        desynch_rate = (R_initial - R_final) / total_time
        print(f"  Desynchronization Rate: {desynch_rate:.3e} per time unit")

        if abs(desynch_rate) < 1e-5:
            print(f"    → Negligible desynchronization (noise well-tolerated)")
        elif desynch_rate > 0:
            print(f"    → Gradual desynchronization (noise effect)")
        else:
            print(f"    → Increasing synchronization (unexpected!)")

    # === SPINOR NORM CONSERVATION ===
    print("\n" + "-"*70)
    print("2. SPINOR NORM CONSERVATION (Dirac Equation Test)")
    print("-"*70)

    norm_initial = ts['spinor_norm'][0]
    norm_final = ts['spinor_norm'][-1]
    norm_mean = ts['spinor_norm'].mean()
    norm_std = ts['spinor_norm'].std()

    norm_dev_initial = abs(norm_initial - 1.0) * 100
    norm_dev_final = abs(norm_final - 1.0) * 100
    norm_dev_mean = abs(norm_mean - 1.0) * 100
    norm_dev_max = max(abs(ts['spinor_norm'] - 1.0)) * 100

    print(f"\n  Spinor Norm ||Ψ||:")
    print(f"    Initial:     {norm_initial:.8f}  (deviation: {norm_dev_initial:.4f}%)")
    print(f"    Final:       {norm_final:.8f}  (deviation: {norm_dev_final:.4f}%)")
    print(f"    Mean:        {norm_mean:.8f}  (deviation: {norm_dev_mean:.4f}%)")
    print(f"    Std Dev:     {norm_std:.8f}")
    print(f"    Max Dev:     {norm_dev_max:.4f}%")

    # Assess conservation quality
    if norm_dev_max < 0.01:
        conservation_quality = "EXCELLENT (<0.01% error)"
    elif norm_dev_max < 0.1:
        conservation_quality = "VERY GOOD (<0.1% error)"
    elif norm_dev_max < 1.0:
        conservation_quality = "GOOD (<1% error)"
    elif norm_dev_max < 10.0:
        conservation_quality = "ACCEPTABLE (<10% error)"
    else:
        conservation_quality = "POOR (>10% error - numerical instability!)"

    print(f"\n  Conservation Quality: {conservation_quality}")

    if norm_dev_max < 1.0:
        print(f"    ✓ Dirac equation numerics are stable")
        print(f"    ✓ Unitary evolution maintained")
    else:
        print(f"    ✗ WARNING: Numerical stability issues detected!")
        print(f"    ✗ Consider smaller timestep or better integrator")

    # === PARTICLE DYNAMICS ===
    print("\n" + "-"*70)
    print("3. PARTICLE DYNAMICS (Center of Mass)")
    print("-"*70)

    x_initial = ts['particle_x'][0]
    y_initial = ts['particle_y'][0]
    x_final = ts['particle_x'][-1]
    y_final = ts['particle_y'][-1]

    drift_final = ts['particle_drift'][-1]
    drift_max = ts['particle_drift'].max()
    drift_mean = ts['particle_drift'].mean()

    print(f"\n  Position:")
    print(f"    Initial:     ({x_initial:.3f}, {y_initial:.3f})")
    print(f"    Final:       ({x_final:.3f}, {y_final:.3f})")
    print(f"    Displacement: ({x_final - x_initial:+.3f}, {y_final - y_initial:+.3f})")

    print(f"\n  Drift Distance:")
    print(f"    Final:       {drift_final:.3f} grid units")
    print(f"    Maximum:     {drift_max:.3f} grid units")
    print(f"    Mean:        {drift_mean:.3f} grid units")

    # Assess localization
    drift_threshold = 5.0  # Common threshold from tests
    if drift_max < drift_threshold:
        localization = f"GOOD (drift < {drift_threshold} units)"
        print(f"\n  Localization Quality: {localization}")
        print(f"    ✓ Particle remains well-localized")
        print(f"    ✓ Stochastic noise does not cause runaway diffusion")
    else:
        localization = f"POOR (drift > {drift_threshold} units)"
        print(f"\n  Localization Quality: {localization}")
        print(f"    ✗ WARNING: Particle diffusing significantly")
        print(f"    ✗ May need stronger coupling or lower noise")

    # Compute velocity statistics
    if len(ts['time']) > 1:
        dt = np.diff(ts['time'])
        vx = np.diff(ts['particle_x']) / dt
        vy = np.diff(ts['particle_y']) / dt
        v_mag = np.sqrt(vx**2 + vy**2)

        print(f"\n  Velocity Statistics:")
        print(f"    Mean |v|:    {v_mag.mean():.4f} grid units/time")
        print(f"    Max |v|:     {v_mag.max():.4f} grid units/time")
        print(f"    RMS |v|:     {np.sqrt(np.mean(v_mag**2)):.4f} grid units/time")

    # === SPATIAL FIELD ANALYSIS ===
    print("\n" + "-"*70)
    print("4. SPATIAL FIELD ANALYSIS")
    print("-"*70)

    # Load initial and final snapshots
    density_files = sorted(glob.glob(os.path.join(data_dir, 'density_*.dat')))
    mass_files = sorted(glob.glob(os.path.join(data_dir, 'mass_field_*.dat')))

    if density_files:
        print(f"\n  Density Fields Found: {len(density_files)} snapshots")

        # Analyze initial density
        density_initial = load_spatial_field(density_files[0])
        density_max_initial = density_initial.max()
        density_sum_initial = density_initial.sum()

        # Analyze final density
        density_final = load_spatial_field(density_files[-1])
        density_max_final = density_final.max()
        density_sum_final = density_final.sum()

        print(f"\n  Spinor Density |Ψ|²:")
        print(f"    Initial max:   {density_max_initial:.6e}")
        print(f"    Final max:     {density_max_final:.6e}")
        print(f"    Initial sum:   {density_sum_initial:.6e}")
        print(f"    Final sum:     {density_sum_final:.6e}")
        print(f"    Sum change:    {(density_sum_final/density_sum_initial - 1)*100:+.3f}%")

        # Check normalization
        if abs(density_sum_final/density_sum_initial - 1) < 0.01:
            print(f"    ✓ Density properly normalized (conserved within 1%)")
        else:
            print(f"    ✗ WARNING: Significant density change (check normalization)")

        # Analyze spreading
        spreading_ratio = density_max_initial / (density_max_final + 1e-20)
        print(f"\n  Spreading Factor: {spreading_ratio:.3f}x")

        if spreading_ratio > 2.0:
            print(f"    → Significant spreading/diffusion observed")
        elif spreading_ratio > 1.1:
            print(f"    → Moderate spreading (expected from stochastic noise)")
        else:
            print(f"    → Particle remains highly localized")

    if mass_files:
        print(f"\n  Mass Fields Found: {len(mass_files)} snapshots")

        # Analyze initial and final mass
        mass_initial = load_spatial_field(mass_files[0])
        mass_final = load_spatial_field(mass_files[-1])

        print(f"\n  Mass Field m(x,y) = Δ·R(x,y):")
        print(f"    Initial mean:  {mass_initial.mean():.6f}")
        print(f"    Final mean:    {mass_final.mean():.6f}")
        print(f"    Initial std:   {mass_initial.std():.6f}")
        print(f"    Final std:     {mass_final.std():.6f}")

        mass_change = (mass_final.mean() / mass_initial.mean() - 1) * 100
        print(f"    Mean change:   {mass_change:+.3f}%")

        # Gradient analysis (proxy for gravitational field)
        grad_y, grad_x = np.gradient(mass_final)
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)

        print(f"\n  Mass Gradient |∇m| (gravity proxy):")
        print(f"    Mean:          {grad_mag.mean():.6e}")
        print(f"    Max:           {grad_mag.max():.6e}")
        print(f"    RMS:           {np.sqrt(np.mean(grad_mag**2)):.6e}")

    # === OVERALL ASSESSMENT ===
    print("\n" + "="*70)
    print("OVERALL ASSESSMENT")
    print("="*70)

    # Collect test results
    tests_passed = 0
    tests_total = 0

    # Test 1: R > 0.95
    tests_total += 1
    if R_min > 0.95:
        print(f"  ✓ Test 1: R_min > 0.95 (PASSED: {R_min:.6f})")
        tests_passed += 1
    else:
        print(f"  ✗ Test 1: R_min > 0.95 (FAILED: {R_min:.6f})")

    # Test 2: Norm conserved within 10%
    tests_total += 1
    if norm_dev_max < 10.0:
        print(f"  ✓ Test 2: Norm conserved <10% (PASSED: {norm_dev_max:.4f}%)")
        tests_passed += 1
    else:
        print(f"  ✗ Test 2: Norm conserved <10% (FAILED: {norm_dev_max:.4f}%)")

    # Test 3: Drift < 5 units
    tests_total += 1
    if drift_final < 5.0:
        print(f"  ✓ Test 3: Particle drift <5 units (PASSED: {drift_final:.3f})")
        tests_passed += 1
    else:
        print(f"  ✗ Test 3: Particle drift <5 units (FAILED: {drift_final:.3f})")

    print(f"\n  Tests Passed: {tests_passed}/{tests_total}")

    if tests_passed == tests_total:
        print(f"\n  🎉 ALL TESTS PASSED - Stochastic Dirac-Kuramoto coupling validated!")
        print(f"     • Synchronization maintained under noise")
        print(f"     • Spinor norm conserved (unitary evolution)")
        print(f"     • Particle localization preserved")
        print(f"     • Ready for production simulations")
    elif tests_passed >= tests_total * 0.66:
        print(f"\n  ⚠️  PARTIAL SUCCESS - Most tests passed")
        print(f"     • Core functionality working")
        print(f"     • Some issues need attention (see above)")
    else:
        print(f"\n  ❌ VALIDATION FAILED - Significant issues detected")
        print(f"     • Review failed tests above")
        print(f"     • Check parameters and numerical methods")

    print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize Dirac-Kuramoto stochastic evolution data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Timeseries evolution
  %(prog)s --timeseries

  # Spatial field evolution
  %(prog)s --spatial

  # Full spinor components
  %(prog)s --spinor --snapshot 0

  # All visualizations
  %(prog)s --all

  # Save outputs
  %(prog)s --timeseries --output dirac_timeseries.png --dpi 300
        """
    )

    parser.add_argument('-d', '--data-dir', default='dirac_evolution',
                       help='Data directory (default: dirac_evolution)')
    parser.add_argument('-t', '--timeseries', action='store_true',
                       help='Plot timeseries evolution')
    parser.add_argument('-s', '--spatial', action='store_true',
                       help='Plot spatial field evolution')
    parser.add_argument('--spinor', action='store_true',
                       help='Plot full spinor components')
    parser.add_argument('--snapshot', type=int, default=0,
                       help='Snapshot index for spinor plot (default: 0)')
    parser.add_argument('--all', action='store_true',
                       help='Generate all visualizations')
    parser.add_argument('-a', '--analyze', action='store_true',
                       help='Print comprehensive text analysis')
    parser.add_argument('-o', '--output', help='Output filename (PNG/PDF/SVG)')
    parser.add_argument('--dpi', type=int, default=150,
                       help='Output resolution (default: 150)')

    args = parser.parse_args()

    # Default to all if nothing specified
    if not (args.timeseries or args.spatial or args.spinor or args.analyze):
        args.all = True

    if args.all:
        args.timeseries = True
        args.spatial = True
        args.spinor = True

    # Check data directory
    if not os.path.isdir(args.data_dir):
        print(f"Error: Directory not found: {args.data_dir}")
        return 1

    # Text analysis
    if args.analyze or args.all:
        analyze_dirac_evolution(args.data_dir)

    # Generate visualizations
    try:
        if args.timeseries:
            output = args.output if not args.all else 'dirac_timeseries.png'
            plot_timeseries_evolution(args.data_dir, output, args.dpi)

        if args.spatial:
            output = args.output if not args.all else 'dirac_spatial.png'
            plot_spatial_evolution(args.data_dir, output, args.dpi)

        if args.spinor:
            output = args.output if not args.all else f'dirac_spinor_{args.snapshot:03d}.png'
            plot_spinor_components(args.data_dir, args.snapshot, output, args.dpi)

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
