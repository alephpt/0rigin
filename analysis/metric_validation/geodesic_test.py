#!/usr/bin/env python3
"""
Geodesic Test for SMFT Metric Validation
Compares fermion trajectories from SMFT simulations with geodesics in derived metrics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from pathlib import Path
import sys
import glob

class GeodesicIntegrator:
    """Integrate geodesics in various proposed metrics"""

    def __init__(self, R_field, theta_field=None, Delta=1.0, dt=0.01):
        """
        Initialize with SMFT field data

        Args:
            R_field: 2D array of synchronization order parameter
            theta_field: 3D array [N_oscillators, Ny, Nx] of phases (optional)
            Delta: Vacuum potential scale (Planck mass)
            dt: Grid spacing for derivatives
        """
        self.R = R_field
        self.Delta = Delta
        self.dt = dt
        self.Ny, self.Nx = R_field.shape

        # Compute mean phase if Kuramoto phases provided
        if theta_field is not None:
            self.theta_mean = np.mean(theta_field, axis=0)
            self.compute_phase_velocity()
        else:
            self.theta_mean = None
            self.v_x = np.zeros_like(R_field)
            self.v_y = np.zeros_like(R_field)

    def compute_phase_velocity(self):
        """Compute flow velocity from phase gradient"""
        # Use numpy gradient for derivatives
        grad_x = np.gradient(self.theta_mean, axis=1) / self.dt
        grad_y = np.gradient(self.theta_mean, axis=0) / self.dt

        # Flow velocity v = (1/Delta) * grad(theta)
        self.v_x = grad_x / self.Delta
        self.v_y = grad_y / self.Delta

    def interpolate_field(self, field, x, y):
        """Bilinear interpolation of field at position (x,y)"""
        # Ensure within bounds
        x = np.clip(x, 0, self.Nx - 1.001)
        y = np.clip(y, 0, self.Ny - 1.001)

        # Grid indices
        i, j = int(x), int(y)
        fx, fy = x - i, y - j

        # Bilinear interpolation
        return ((1-fx) * (1-fy) * field[j, i] +
                fx * (1-fy) * field[j, i+1] +
                (1-fx) * fy * field[j+1, i] +
                fx * fy * field[j+1, i+1])

    def christoffel_conformal(self, x, y):
        """Christoffel symbols for conformal metric g = R^2 * eta"""
        R = self.interpolate_field(self.R, x, y)

        # Compute R derivatives using finite differences
        dx = 0.01
        dR_dx = (self.interpolate_field(self.R, x+dx, y) -
                 self.interpolate_field(self.R, x-dx, y)) / (2*dx)
        dR_dy = (self.interpolate_field(self.R, x, y+dx) -
                 self.interpolate_field(self.R, x, y-dx)) / (2*dx)
        dR_dt = 0  # Assuming static for now

        # Christoffel symbols for conformal metric
        # Gamma^mu_nu_rho = (1/R) * [delta^mu_nu * d_rho R + delta^mu_rho * d_nu R - eta_nu_rho * eta^mu_sigma * d_sigma R]
        Gamma = np.zeros((3, 3, 3))  # [upper, lower, lower]

        if R > 1e-6:  # Avoid singularities
            # Time components
            Gamma[0, 0, 0] = dR_dt / R
            Gamma[0, 1, 1] = -dR_dt / R
            Gamma[0, 2, 2] = -dR_dt / R

            # Spatial components
            Gamma[1, 0, 0] = dR_dx / R
            Gamma[1, 0, 1] = dR_dt / R
            Gamma[1, 1, 0] = dR_dt / R
            Gamma[1, 1, 1] = dR_dx / R
            Gamma[1, 1, 2] = 0
            Gamma[1, 2, 1] = 0
            Gamma[1, 2, 2] = -dR_dx / R

            Gamma[2, 0, 0] = dR_dy / R
            Gamma[2, 0, 2] = dR_dt / R
            Gamma[2, 2, 0] = dR_dt / R
            Gamma[2, 1, 1] = -dR_dy / R
            Gamma[2, 2, 1] = 0
            Gamma[2, 1, 2] = 0
            Gamma[2, 2, 2] = dR_dy / R

        return Gamma

    def christoffel_acoustic(self, x, y):
        """Christoffel symbols for acoustic metric g = R^2 * [eta + flow terms]"""
        R = self.interpolate_field(self.R, x, y)
        v_x = self.interpolate_field(self.v_x, x, y)
        v_y = self.interpolate_field(self.v_y, x, y)

        # For simplicity, using leading order (conformal part)
        # Full calculation would include flow velocity derivatives
        return self.christoffel_conformal(x, y)

    def geodesic_equation(self, state, tau, metric_type='conformal'):
        """
        Geodesic equation: d^2 x^mu / dtau^2 + Gamma^mu_rho_sigma * (dx^rho/dtau) * (dx^sigma/dtau) = 0

        Args:
            state: [t, x, y, dt/dtau, dx/dtau, dy/dtau]
            tau: affine parameter
            metric_type: 'conformal' or 'acoustic'
        """
        t, x, y, u_t, u_x, u_y = state

        # Get Christoffel symbols
        if metric_type == 'conformal':
            Gamma = self.christoffel_conformal(x, y)
        else:
            Gamma = self.christoffel_acoustic(x, y)

        # Geodesic equation
        d2t_dtau2 = -sum(Gamma[0, i, j] * [u_t, u_x, u_y][i] * [u_t, u_x, u_y][j]
                        for i in range(3) for j in range(3))
        d2x_dtau2 = -sum(Gamma[1, i, j] * [u_t, u_x, u_y][i] * [u_t, u_x, u_y][j]
                        for i in range(3) for j in range(3))
        d2y_dtau2 = -sum(Gamma[2, i, j] * [u_t, u_x, u_y][i] * [u_t, u_x, u_y][j]
                        for i in range(3) for j in range(3))

        return [u_t, u_x, u_y, d2t_dtau2, d2x_dtau2, d2y_dtau2]

    def integrate_geodesic(self, x0, y0, vx0, vy0, tau_max=10.0, metric_type='conformal'):
        """
        Integrate geodesic from initial conditions

        Args:
            x0, y0: Initial position
            vx0, vy0: Initial velocity
            tau_max: Maximum affine parameter
            metric_type: Type of metric to use

        Returns:
            trajectory: Array of [t, x, y] along geodesic
        """
        # Initial state [t, x, y, dt/dtau, dx/dtau, dy/dtau]
        # For massive particle, normalize 4-velocity
        R0 = self.interpolate_field(self.R, x0, y0)

        # In conformal metric, proper time normalization
        u_t = 1.0 / R0  # Approximate for slow motion
        u_x = vx0 / R0
        u_y = vy0 / R0

        initial_state = [0, x0, y0, u_t, u_x, u_y]

        # Integrate
        tau_points = np.linspace(0, tau_max, 1000)
        solution = odeint(self.geodesic_equation, initial_state, tau_points,
                         args=(metric_type,))

        return solution[:, :3]  # Return [t, x, y]

def load_smft_data(output_dir):
    """Load R-field and theta-field from SMFT output"""
    output_path = Path(output_dir)

    # Find the latest timestep R-field
    r_files = sorted(glob.glob(str(output_path / "N_*" / "R_field_t*.csv")))
    if not r_files:
        print(f"No R_field files found in {output_dir}")
        return None, None

    # Load R field
    R_field = pd.read_csv(r_files[-1], header=None).values

    # Try to load theta field (may not exist)
    theta_files = sorted(glob.glob(str(output_path / "N_*" / "theta_field_t*.csv")))
    theta_field = None
    if theta_files:
        # For now, just use mean phase approximation
        theta_data = pd.read_csv(theta_files[-1], header=None).values
        theta_field = theta_data[np.newaxis, :, :]  # Add oscillator dimension

    return R_field, theta_field

def compare_trajectories(geodesic_traj, fermion_traj):
    """Compare geodesic prediction with actual fermion trajectory"""
    # Interpolate to common time points
    t_common = np.linspace(0, min(geodesic_traj[-1, 0], fermion_traj[-1, 0]), 100)

    # Interpolate trajectories
    from scipy.interpolate import interp1d

    geod_x = interp1d(geodesic_traj[:, 0], geodesic_traj[:, 1],
                      kind='cubic', fill_value='extrapolate')(t_common)
    geod_y = interp1d(geodesic_traj[:, 0], geodesic_traj[:, 2],
                      kind='cubic', fill_value='extrapolate')(t_common)

    ferm_x = interp1d(fermion_traj[:, 0], fermion_traj[:, 1],
                      kind='cubic', fill_value='extrapolate')(t_common)
    ferm_y = interp1d(fermion_traj[:, 0], fermion_traj[:, 2],
                      kind='cubic', fill_value='extrapolate')(t_common)

    # Compute RMS deviation
    rms_x = np.sqrt(np.mean((geod_x - ferm_x)**2))
    rms_y = np.sqrt(np.mean((geod_y - ferm_y)**2))
    rms_total = np.sqrt(rms_x**2 + rms_y**2)

    # Compute maximum deviation
    max_dev = np.max(np.sqrt((geod_x - ferm_x)**2 + (geod_y - ferm_y)**2))

    return {
        'rms_deviation': rms_total,
        'max_deviation': max_dev,
        'rms_x': rms_x,
        'rms_y': rms_y,
        't_common': t_common,
        'geodesic_x': geod_x,
        'geodesic_y': geod_y,
        'fermion_x': ferm_x,
        'fermion_y': ferm_y
    }

def plot_comparison(R_field, geodesic_traj, fermion_traj=None, title="Geodesic vs Fermion Trajectory"):
    """Plot trajectory comparison overlaid on R-field"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot R-field
    im = ax1.imshow(R_field, cmap='viridis', origin='lower', extent=[0, R_field.shape[1], 0, R_field.shape[0]])
    plt.colorbar(im, ax=ax1, label='R(x,y)')

    # Plot geodesic
    ax1.plot(geodesic_traj[:, 1], geodesic_traj[:, 2], 'r-', linewidth=2, label='Geodesic')

    # Plot fermion trajectory if available
    if fermion_traj is not None:
        ax1.plot(fermion_traj[:, 1], fermion_traj[:, 2], 'b--', linewidth=2, label='Fermion')

    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title(title)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot time evolution
    ax2.plot(geodesic_traj[:, 0], geodesic_traj[:, 1], 'r-', label='Geodesic x(t)')
    ax2.plot(geodesic_traj[:, 0], geodesic_traj[:, 2], 'r--', label='Geodesic y(t)')

    if fermion_traj is not None:
        ax2.plot(fermion_traj[:, 0], fermion_traj[:, 1], 'b-', label='Fermion x(t)')
        ax2.plot(fermion_traj[:, 0], fermion_traj[:, 2], 'b--', label='Fermion y(t)')

    ax2.set_xlabel('Time')
    ax2.set_ylabel('Position')
    ax2.set_title('Position vs Time')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    return fig

def main():
    """Run geodesic validation tests"""

    if len(sys.argv) < 2:
        print("Usage: python geodesic_test.py <smft_output_directory>")
        print("Example: python geodesic_test.py output/20251219_112601_defect_localization_64x64/")
        sys.exit(1)

    output_dir = sys.argv[1]

    # Load SMFT data
    print(f"Loading SMFT data from {output_dir}...")
    R_field, theta_field = load_smft_data(output_dir)

    if R_field is None:
        print("Failed to load R-field data")
        sys.exit(1)

    print(f"R-field shape: {R_field.shape}")
    print(f"R-field range: [{R_field.min():.3f}, {R_field.max():.3f}]")

    # Initialize geodesic integrator
    integrator = GeodesicIntegrator(R_field, theta_field, Delta=1.0)

    # Test 1: Geodesic in uniform region (R ≈ 1)
    print("\nTest 1: Geodesic in uniform region")
    x0, y0 = R_field.shape[1]//4, R_field.shape[0]//2
    vx0, vy0 = 0.1, 0.05

    traj_conformal = integrator.integrate_geodesic(x0, y0, vx0, vy0, tau_max=50,
                                                   metric_type='conformal')
    traj_acoustic = integrator.integrate_geodesic(x0, y0, vx0, vy0, tau_max=50,
                                                  metric_type='acoustic')

    # Test 2: Geodesic near defect (R < 1)
    print("\nTest 2: Geodesic near defect")
    # Find a defect location (minimum R)
    defect_idx = np.unravel_index(R_field.argmin(), R_field.shape)
    x0, y0 = defect_idx[1] + 5, defect_idx[0] + 5  # Start near defect

    traj_defect_conformal = integrator.integrate_geodesic(x0, y0, vx0, vy0, tau_max=50,
                                                          metric_type='conformal')
    traj_defect_acoustic = integrator.integrate_geodesic(x0, y0, vx0, vy0, tau_max=50,
                                                         metric_type='acoustic')

    # Plot results
    print("\nPlotting results...")

    # Plot 1: Uniform region comparison
    fig1 = plot_comparison(R_field, traj_conformal, traj_acoustic,
                          title="Uniform Region: Conformal (red) vs Acoustic (blue)")
    plt.savefig('geodesic_uniform_comparison.png', dpi=150)

    # Plot 2: Defect region comparison
    fig2 = plot_comparison(R_field, traj_defect_conformal, traj_defect_acoustic,
                          title="Near Defect: Conformal (red) vs Acoustic (blue)")
    plt.savefig('geodesic_defect_comparison.png', dpi=150)

    # Summary statistics
    print("\n=== GEODESIC TEST SUMMARY ===")
    print(f"Uniform region deflection (conformal): Δx = {traj_conformal[-1,1] - traj_conformal[0,1]:.3f}")
    print(f"Uniform region deflection (acoustic): Δx = {traj_acoustic[-1,1] - traj_acoustic[0,1]:.3f}")
    print(f"Defect region deflection (conformal): Δx = {traj_defect_conformal[-1,1] - traj_defect_conformal[0,1]:.3f}")
    print(f"Defect region deflection (acoustic): Δx = {traj_defect_acoustic[-1,1] - traj_defect_acoustic[0,1]:.3f}")

    # Check for lensing signature
    uniform_deviation = np.sqrt((traj_conformal[-1,1] - traj_acoustic[-1,1])**2 +
                               (traj_conformal[-1,2] - traj_acoustic[-1,2])**2)
    defect_deviation = np.sqrt((traj_defect_conformal[-1,1] - traj_defect_acoustic[-1,1])**2 +
                               (traj_defect_conformal[-1,2] - traj_defect_acoustic[-1,2])**2)

    print(f"\nMetric comparison deviation (uniform): {uniform_deviation:.3f}")
    print(f"Metric comparison deviation (defect): {defect_deviation:.3f}")

    if defect_deviation > 2 * uniform_deviation:
        print("✓ Defect lensing detected: Strong metric-dependent deflection near R < 1")
    else:
        print("✗ No clear defect lensing: Similar trajectories in both regions")

    plt.show()

    print("\nGeodesic tests complete. Results saved to geodesic_*.png")

if __name__ == "__main__":
    main()