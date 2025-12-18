#!/usr/bin/env python3
"""
MSFT Field Visualization Tool

Visualizes output from MSFT (Multi-Scale Field Theory) simulation:
- R_field: Synchronization defect field (scalar)
- theta: Phase field (scalar, radians)
- gravity_x, gravity_y: Gravitational field components (vector)
- timeseries: Temporal evolution of R statistics (min/avg/max)

Data format: 128×128 grid, space-separated floats in .dat files
Timeseries format: timestep value pairs
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import argparse
import os
from pathlib import Path


def load_field_data(filepath, grid_size=128):
    """
    Load field data from .dat file.

    Args:
        filepath: Path to .dat file
        grid_size: Grid dimension (default 128×128)

    Returns:
        2D numpy array of shape (grid_size, grid_size)
    """
    try:
        data = np.loadtxt(filepath)
        if data.size != grid_size * grid_size:
            raise ValueError(f"Expected {grid_size*grid_size} values, got {data.size}")
        return data.reshape((grid_size, grid_size))
    except Exception as e:
        raise RuntimeError(f"Failed to load {filepath}: {e}")


def load_timeseries_data(filepath):
    """
    Load timeseries data from .dat file.

    Args:
        filepath: Path to timeseries .dat file

    Returns:
        tuple: (timesteps, values) as 1D numpy arrays
    """
    try:
        data = np.loadtxt(filepath)
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError(f"Expected Nx2 array, got shape {data.shape}")
        return data[:, 0], data[:, 1]
    except Exception as e:
        raise RuntimeError(f"Failed to load timeseries {filepath}: {e}")


def plot_scalar_field(ax, field, title, cmap='viridis', show_colorbar=True):
    """Plot a scalar field as a heatmap."""
    im = ax.imshow(field, origin='lower', cmap=cmap, interpolation='bilinear')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('X (grid units)')
    ax.set_ylabel('Y (grid units)')

    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    return im


def plot_vector_field(ax, vx, vy, title, subsample=8, show_magnitude=True, cmap='plasma'):
    """
    Plot a vector field with arrows and optional magnitude background.

    Args:
        ax: Matplotlib axes
        vx, vy: Vector field components
        title: Plot title
        subsample: Arrow grid spacing (plot every Nth arrow)
        show_magnitude: Show magnitude as background color
        cmap: Colormap for magnitude
    """
    grid_size = vx.shape[0]

    # Calculate magnitude
    magnitude = np.sqrt(vx**2 + vy**2)

    # Show magnitude as background
    if show_magnitude:
        im = ax.imshow(magnitude, origin='lower', cmap=cmap, interpolation='bilinear', alpha=0.7)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Magnitude')

    # Subsample for quiver plot (too dense otherwise)
    step = subsample
    x = np.arange(0, grid_size, step)
    y = np.arange(0, grid_size, step)
    X, Y = np.meshgrid(x, y)

    # Subsample vector components
    U = vx[::step, ::step]
    V = vy[::step, ::step]

    # Normalize arrow lengths for visibility
    M = np.sqrt(U**2 + V**2)
    M[M == 0] = 1  # Avoid division by zero
    U_norm = U / M
    V_norm = V / M

    # Plot arrows
    ax.quiver(X, Y, U_norm, V_norm, M[::1, ::1],
              cmap='coolwarm', alpha=0.8, scale=25, width=0.003)

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('X (grid units)')
    ax.set_ylabel('Y (grid units)')
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)


def plot_gradient_field(ax, field, title, subsample=8):
    """
    Plot gradient of a scalar field (∇R or ∇θ).

    Args:
        ax: Matplotlib axes
        field: Scalar field
        title: Plot title
        subsample: Arrow grid spacing
    """
    # Compute gradient
    grad_y, grad_x = np.gradient(field)

    # Plot as vector field
    plot_vector_field(ax, grad_x, grad_y, title, subsample=subsample,
                     show_magnitude=True, cmap='inferno')


def plot_divergence_curl(ax1, ax2, vx, vy):
    """
    Plot divergence and curl of vector field.

    Divergence: ∇·g = ∂gx/∂x + ∂gy/∂y (sources/sinks)
    Curl: ∇×g = ∂gy/∂x - ∂gx/∂y (rotation)
    """
    # Compute derivatives
    dgx_dx, dgx_dy = np.gradient(vx)
    dgy_dx, dgy_dy = np.gradient(vy)

    # Divergence (∇·g)
    divergence = dgx_dx + dgy_dy

    # Curl (∇×g, z-component only for 2D)
    curl = dgy_dx - dgx_dy

    # Plot divergence
    plot_scalar_field(ax1, divergence, 'Divergence ∇·g\n(sources/sinks)',
                     cmap='RdBu_r', show_colorbar=True)

    # Plot curl
    plot_scalar_field(ax2, curl, 'Curl ∇×g\n(rotation)',
                     cmap='PuOr', show_colorbar=True)


def create_comprehensive_visualization(data_dir, output_file=None, dpi=150):
    """
    Create comprehensive multi-panel visualization of all MSFT fields.

    Args:
        data_dir: Directory containing .dat files
        output_file: Optional output filename (PNG/PDF)
        dpi: Resolution for saved image
    """
    # Load all fields
    print("Loading field data...")
    R_field = load_field_data(os.path.join(data_dir, 'R_field.dat'))
    theta = load_field_data(os.path.join(data_dir, 'theta.dat'))
    gx = load_field_data(os.path.join(data_dir, 'gravity_x.dat'))
    gy = load_field_data(os.path.join(data_dir, 'gravity_y.dat'))

    print(f"R_field range: [{R_field.min():.3f}, {R_field.max():.3f}]")
    print(f"theta range: [{theta.min():.3f}, {theta.max():.3f}] rad")
    print(f"gravity_x range: [{gx.min():.3f}, {gx.max():.3f}]")
    print(f"gravity_y range: [{gy.min():.3f}, {gy.max():.3f}]")

    # Create figure with subplots
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle('MSFT Field Visualization (128×128 Grid)',
                 fontsize=16, fontweight='bold', y=0.995)

    # Layout: 3 rows × 3 columns
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.35)

    # Row 1: Scalar fields
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    plot_scalar_field(ax1, R_field, 'R Field\n(Synchronization Defects)', cmap='viridis')
    plot_scalar_field(ax2, theta, 'θ Field\n(Phase, radians)', cmap='twilight')
    plot_scalar_field(ax3, np.sqrt(gx**2 + gy**2),
                     'Gravity Magnitude |g|', cmap='plasma')

    # Row 2: Vector field and gradients
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    plot_vector_field(ax4, gx, gy, 'Gravity Vector Field g = (gx, gy)',
                     subsample=6, show_magnitude=False)
    plot_gradient_field(ax5, R_field, 'Gradient ∇R', subsample=6)
    plot_gradient_field(ax6, theta, 'Gradient ∇θ', subsample=6)

    # Row 3: Divergence, curl, and relationship
    ax7 = fig.add_subplot(gs[2, 0])
    ax8 = fig.add_subplot(gs[2, 1])
    ax9 = fig.add_subplot(gs[2, 2])

    plot_divergence_curl(ax7, ax8, gx, gy)

    # Compute ∇R vs g correlation
    grad_R_y, grad_R_x = np.gradient(R_field)
    correlation = grad_R_x * gx + grad_R_y * gy
    plot_scalar_field(ax9, correlation,
                     'g · ∇R Correlation\n(Theory: g = -λ∇R)',
                     cmap='seismic', show_colorbar=True)

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def create_simple_visualization(data_dir, field_name, output_file=None, dpi=150):
    """
    Create simple single-field visualization.

    Args:
        data_dir: Directory containing .dat files
        field_name: 'R_field', 'theta', 'gravity' (for vector), or 'gravity_magnitude'
        output_file: Optional output filename
        dpi: Resolution
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    if field_name == 'R_field':
        R_field = load_field_data(os.path.join(data_dir, 'R_field.dat'))
        plot_scalar_field(ax, R_field, 'R Field (Synchronization Defects)', cmap='viridis')

    elif field_name == 'theta':
        theta = load_field_data(os.path.join(data_dir, 'theta.dat'))
        plot_scalar_field(ax, theta, 'θ Field (Phase, radians)', cmap='twilight')

    elif field_name == 'gravity':
        gx = load_field_data(os.path.join(data_dir, 'gravity_x.dat'))
        gy = load_field_data(os.path.join(data_dir, 'gravity_y.dat'))
        plot_vector_field(ax, gx, gy, 'Gravity Vector Field', subsample=6)

    elif field_name == 'gravity_magnitude':
        gx = load_field_data(os.path.join(data_dir, 'gravity_x.dat'))
        gy = load_field_data(os.path.join(data_dir, 'gravity_y.dat'))
        magnitude = np.sqrt(gx**2 + gy**2)
        plot_scalar_field(ax, magnitude, 'Gravity Magnitude |g|', cmap='plasma')

    else:
        raise ValueError(f"Unknown field: {field_name}")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def analyze_field_statistics(data_dir):
    """Print statistical analysis of all fields."""
    print("\n" + "="*60)
    print("MSFT FIELD STATISTICS")
    print("="*60)

    # Load fields
    R_field = load_field_data(os.path.join(data_dir, 'R_field.dat'))
    theta = load_field_data(os.path.join(data_dir, 'theta.dat'))
    gx = load_field_data(os.path.join(data_dir, 'gravity_x.dat'))
    gy = load_field_data(os.path.join(data_dir, 'gravity_y.dat'))

    def stats(name, field):
        print(f"\n{name}:")
        print(f"  Range: [{field.min():.6f}, {field.max():.6f}]")
        print(f"  Mean: {field.mean():.6f}")
        print(f"  Std Dev: {field.std():.6f}")
        print(f"  RMS: {np.sqrt(np.mean(field**2)):.6f}")

    stats("R_field", R_field)
    stats("theta (radians)", theta)
    stats("gravity_x", gx)
    stats("gravity_y", gy)

    # Gravity magnitude
    g_mag = np.sqrt(gx**2 + gy**2)
    stats("gravity magnitude", g_mag)

    # Divergence and curl
    dgx_dx, dgx_dy = np.gradient(gx)
    dgy_dx, dgy_dy = np.gradient(gy)
    divergence = dgx_dx + dgy_dy
    curl = dgy_dx - dgx_dy

    stats("divergence (∇·g)", divergence)
    stats("curl (∇×g)", curl)

    # Test g = -λ∇R relationship
    grad_R_y, grad_R_x = np.gradient(R_field)

    # Mask to avoid division by zero
    mask = (grad_R_x**2 + grad_R_y**2) > 1e-6

    if mask.sum() > 0:
        # Estimate λ from g = -λ∇R
        lambda_x = -gx[mask] / grad_R_x[mask]
        lambda_y = -gy[mask] / grad_R_y[mask]

        print(f"\n\nRelationship g = -λ∇R:")
        print(f"  λ (from x-component): {lambda_x.mean():.6f} ± {lambda_x.std():.6f}")
        print(f"  λ (from y-component): {lambda_y.mean():.6f} ± {lambda_y.std():.6f}")

        # Correlation coefficient
        corr_x = np.corrcoef(gx[mask], -grad_R_x[mask])[0, 1]
        corr_y = np.corrcoef(gy[mask], -grad_R_y[mask])[0, 1]
        print(f"  Correlation (x): {corr_x:.6f}")
        print(f"  Correlation (y): {corr_y:.6f}")

    print("\n" + "="*60 + "\n")


def plot_timeseries(data_dir, output_file=None, dpi=150):
    """
    Plot temporal evolution of R field statistics.

    Args:
        data_dir: Directory containing timeseries_*.dat files
        output_file: Optional output filename
        dpi: Resolution
    """
    # Load timeseries data
    print("Loading timeseries data...")

    files = {
        'max': os.path.join(data_dir, 'timeseries_R_max.dat'),
        'avg': os.path.join(data_dir, 'timeseries_R_avg.dat'),
        'min': os.path.join(data_dir, 'timeseries_R_min.dat')
    }

    # Check if files exist
    missing = [k for k, v in files.items() if not os.path.exists(v)]
    if missing:
        print(f"Warning: Missing timeseries files: {missing}")
        return

    t_max, R_max = load_timeseries_data(files['max'])
    t_avg, R_avg = load_timeseries_data(files['avg'])
    t_min, R_min = load_timeseries_data(files['min'])

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('MSFT Temporal Evolution - R Field Statistics',
                 fontsize=14, fontweight='bold')

    # Plot 1: All three timeseries
    ax1 = axes[0, 0]
    ax1.plot(t_max, R_max, 'r-', label='R_max', linewidth=1.5, alpha=0.8)
    ax1.plot(t_avg, R_avg, 'b-', label='R_avg', linewidth=1.5, alpha=0.8)
    ax1.plot(t_min, R_min, 'g-', label='R_min', linewidth=1.5, alpha=0.8)
    ax1.fill_between(t_max, R_min, R_max, alpha=0.2, color='gray', label='Range')
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('R Value')
    ax1.set_title('R Field Evolution (Min/Avg/Max)')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Range and spread
    ax2 = axes[0, 1]
    R_range = R_max - R_min
    R_spread = (R_max - R_avg) / (R_avg - R_min + 1e-10)  # Asymmetry metric

    ax2_twin = ax2.twinx()
    ln1 = ax2.plot(t_max, R_range, 'purple', label='Range (max-min)', linewidth=1.5)
    ln2 = ax2_twin.plot(t_avg, R_spread, 'orange', label='Spread asymmetry',
                        linewidth=1.5, alpha=0.7)

    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Range (max - min)', color='purple')
    ax2_twin.set_ylabel('Spread Asymmetry', color='orange')
    ax2.set_title('Field Range & Distribution')
    ax2.tick_params(axis='y', labelcolor='purple')
    ax2_twin.tick_params(axis='y', labelcolor='orange')

    # Combined legend
    lns = ln1 + ln2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='best')
    ax2.grid(True, alpha=0.3)

    # Plot 3: Rate of change (derivatives)
    ax3 = axes[1, 0]

    # Compute derivatives (forward difference)
    dt = np.diff(t_avg)
    dR_max_dt = np.diff(R_max) / dt
    dR_avg_dt = np.diff(R_avg) / dt
    dR_min_dt = np.diff(R_min) / dt

    t_deriv = t_avg[:-1]  # Derivative at midpoints

    ax3.plot(t_deriv, dR_max_dt, 'r-', label='dR_max/dt', linewidth=1.5, alpha=0.7)
    ax3.plot(t_deriv, dR_avg_dt, 'b-', label='dR_avg/dt', linewidth=1.5, alpha=0.7)
    ax3.plot(t_deriv, dR_min_dt, 'g-', label='dR_min/dt', linewidth=1.5, alpha=0.7)
    ax3.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax3.set_xlabel('Timestep')
    ax3.set_ylabel('dR/dt')
    ax3.set_title('Rate of Change (Temporal Derivative)')
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3)

    # Plot 4: Phase space (R_avg vs dR_avg/dt)
    ax4 = axes[1, 1]

    # Use gradient coloring for time
    scatter = ax4.scatter(R_avg[:-1], dR_avg_dt, c=t_deriv, cmap='viridis',
                         s=10, alpha=0.6)

    # Add trajectory arrow
    n_arrows = 20
    step = len(R_avg[:-1]) // n_arrows
    for i in range(0, len(R_avg[:-1]) - step, step):
        ax4.annotate('', xy=(R_avg[i+step], dR_avg_dt[i+step]),
                    xytext=(R_avg[i], dR_avg_dt[i]),
                    arrowprops=dict(arrowstyle='->', color='red', alpha=0.3, lw=0.5))

    ax4.set_xlabel('R_avg')
    ax4.set_ylabel('dR_avg/dt')
    ax4.set_title('Phase Space (R vs dR/dt)')
    ax4.grid(True, alpha=0.3)
    ax4.axhline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Timestep')

    plt.tight_layout()

    # Print statistics
    print(f"\nTimeseries Statistics:")
    print(f"  Total timesteps: {len(t_avg)}")
    print(f"  Final R_avg: {R_avg[-1]:.6f}")
    print(f"  Final R_max: {R_max[-1]:.6f}")
    print(f"  Final R_min: {R_min[-1]:.6f}")
    print(f"  Avg growth rate (R_avg): {dR_avg_dt.mean():.6e} per timestep")
    print(f"  Final range: {R_range[-1]:.6f}")

    # Detect equilibrium or trends
    if len(dR_avg_dt) > 100:
        recent_trend = dR_avg_dt[-100:].mean()
        if abs(recent_trend) < 1e-6:
            print(f"  Status: Near equilibrium (recent dR/dt ≈ {recent_trend:.2e})")
        elif recent_trend > 0:
            print(f"  Status: Growing (recent dR/dt = {recent_trend:.2e})")
        else:
            print(f"  Status: Decaying (recent dR/dt = {recent_trend:.2e})")

    # Save or show
    if output_file:
        print(f"Saving to {output_file}...")
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize MSFT simulation output fields',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Comprehensive visualization (all fields)
  %(prog)s --comprehensive

  # Single field
  %(prog)s --field R_field
  %(prog)s --field gravity

  # Timeseries evolution
  %(prog)s --timeseries

  # Save to file
  %(prog)s --comprehensive --output msft_fields.png --dpi 300
  %(prog)s --timeseries --output evolution.png

  # Statistics only
  %(prog)s --stats

Available fields:
  R_field           - Synchronization defect field
  theta             - Phase field
  gravity           - Gravity vector field
  gravity_magnitude - Gravity magnitude
        """
    )

    parser.add_argument('-d', '--data-dir', default='.',
                       help='Directory containing .dat files (default: current directory)')
    parser.add_argument('-f', '--field',
                       choices=['R_field', 'theta', 'gravity', 'gravity_magnitude'],
                       help='Visualize single field')
    parser.add_argument('-c', '--comprehensive', action='store_true',
                       help='Create comprehensive multi-panel visualization')
    parser.add_argument('-t', '--timeseries', action='store_true',
                       help='Plot temporal evolution of R field statistics')
    parser.add_argument('-s', '--stats', action='store_true',
                       help='Print field statistics')
    parser.add_argument('-o', '--output', help='Output filename (PNG/PDF/SVG)')
    parser.add_argument('--dpi', type=int, default=150,
                       help='Output resolution (default: 150)')

    args = parser.parse_args()

    # Default to comprehensive if nothing specified
    if not (args.field or args.comprehensive or args.stats or args.timeseries):
        args.comprehensive = True

    # Check data directory exists
    if not os.path.isdir(args.data_dir):
        print(f"Error: Directory not found: {args.data_dir}")
        return 1

    # Statistics
    if args.stats:
        analyze_field_statistics(args.data_dir)

    # Visualization
    try:
        if args.comprehensive:
            create_comprehensive_visualization(args.data_dir, args.output, args.dpi)
        elif args.field:
            create_simple_visualization(args.data_dir, args.field, args.output, args.dpi)
        elif args.timeseries:
            plot_timeseries(args.data_dir, args.output, args.dpi)

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
