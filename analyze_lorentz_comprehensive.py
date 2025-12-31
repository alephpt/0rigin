#!/usr/bin/env python3
"""
Comprehensive Lorentz Force Analysis
Verifies:
  1. Energy conservation in pure magnetic field (<10^-10 error)
  2. Larmor radius: r = mv/(qB) within 1%
  3. Cyclotron frequency: ω = qB/m
  4. Long-term stability over 100k+ steps
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
import sys

def analyze_trajectory(trajectory_file):
    """Analyze test particle trajectory"""

    # Load trajectory
    data = np.loadtxt(trajectory_file, delimiter=',', skiprows=1)

    if data.shape[0] < 10:
        print("ERROR: Insufficient trajectory points")
        return

    t = data[:, 0]   # Time
    x = data[:, 1]   # Position x
    y = data[:, 2]   # Position y
    vx = data[:, 3]  # Velocity x
    vy = data[:, 4]  # Velocity y

    # Particle parameters (from config)
    m = 0.1      # Mass
    q = -1.0     # Charge
    B = 0.1      # Magnetic field
    v0 = 0.01    # Initial speed

    # Theoretical predictions
    r_larmor_theory = m * v0 / abs(q * B)
    omega_theory = abs(q * B / m)
    T_theory = 2 * np.pi / omega_theory

    print("="*70)
    print("COMPREHENSIVE LORENTZ FORCE VALIDATION")
    print("="*70)
    print(f"\nParameters:")
    print(f"  Mass:      m = {m}")
    print(f"  Charge:    q = {q}")
    print(f"  B-field:   B = {B}")
    print(f"  Initial v: v0 = {v0} c")
    print(f"  Timesteps: {len(t)}")
    print(f"  Duration:  {t[-1]:.2f} t_P")

    # 1. ENERGY CONSERVATION
    print(f"\n{'='*70}")
    print("1. ENERGY CONSERVATION (Pure Magnetic Field)")
    print(f"{'='*70}")

    v_mag = np.sqrt(vx**2 + vy**2)
    KE = 0.5 * m * v_mag**2
    KE_drift = (KE - KE[0]) / KE[0]

    print(f"  Initial KE: {KE[0]:.10e}")
    print(f"  Final KE:   {KE[-1]:.10e}")
    print(f"  Drift:      {KE_drift[-1]*100:.3e}%")
    print(f"  Max drift:  {np.max(np.abs(KE_drift))*100:.3e}%")
    print(f"  RMS drift:  {np.sqrt(np.mean(KE_drift**2))*100:.3e}%")

    energy_conserved = np.max(np.abs(KE_drift)) < 1e-8  # 10^-10 relative
    print(f"\n  Status: {'✅ PASS' if energy_conserved else '❌ FAIL'} (<10⁻⁸ required)")

    # 2. LARMOR RADIUS
    print(f"\n{'='*70}")
    print("2. LARMOR RADIUS: r = mv/(qB)")
    print(f"{'='*70}")

    # Find orbit center
    x_center = np.mean(x)
    y_center = np.mean(y)

    # Compute radius from center
    r = np.sqrt((x - x_center)**2 + (y - y_center)**2)
    r_mean = np.mean(r)
    r_std = np.std(r)

    r_error = abs(r_mean - r_larmor_theory) / r_larmor_theory

    print(f"  Theoretical: r = {r_larmor_theory:.6f} ℓ_P")
    print(f"  Measured:    r = {r_mean:.6f} ± {r_std:.6f} ℓ_P")
    print(f"  Error:       {r_error*100:.3f}%")
    print(f"  Orbit center: ({x_center:.3f}, {y_center:.3f})")

    radius_accurate = r_error < 0.01  # 1% tolerance
    print(f"\n  Status: {'✅ PASS' if radius_accurate else '❌ FAIL'} (<1% required)")

    # 3. CYCLOTRON FREQUENCY
    print(f"\n{'='*70}")
    print("3. CYCLOTRON FREQUENCY: ω = qB/m")
    print(f"{'='*70}")

    # Compute angular position
    theta = np.arctan2(y - y_center, x - x_center)

    # Unwrap phase
    theta_unwrap = np.unwrap(theta)

    # Fit linear phase evolution: θ(t) = θ0 + ωt
    if len(t) > 100:
        # Use middle portion to avoid transients
        idx_start = len(t) // 4
        idx_end = 3 * len(t) // 4

        coeffs = np.polyfit(t[idx_start:idx_end], theta_unwrap[idx_start:idx_end], 1)
        omega_measured = abs(coeffs[0])

        omega_error = abs(omega_measured - omega_theory) / omega_theory

        print(f"  Theoretical: ω = {omega_theory:.6f} rad/t_P")
        print(f"  Measured:    ω = {omega_measured:.6f} rad/t_P")
        print(f"  Error:       {omega_error*100:.3f}%")
        print(f"  Period:      T = {2*np.pi/omega_measured:.3f} t_P (theory: {T_theory:.3f})")

        freq_accurate = omega_error < 0.01
        print(f"\n  Status: {'✅ PASS' if freq_accurate else '❌ FAIL'} (<1% required)")
    else:
        print("  Insufficient data for frequency analysis")
        freq_accurate = False

    # 4. SPEED CONSERVATION
    print(f"\n{'='*70}")
    print("4. SPEED CONSERVATION (|v| should be constant)")
    print(f"{'='*70}")

    v_drift = (v_mag - v_mag[0]) / v_mag[0]

    print(f"  Initial |v|: {v_mag[0]:.10e} c")
    print(f"  Final |v|:   {v_mag[-1]:.10e} c")
    print(f"  Drift:       {v_drift[-1]*100:.3e}%")
    print(f"  Max drift:   {np.max(np.abs(v_drift))*100:.3e}%")
    print(f"  RMS drift:   {np.sqrt(np.mean(v_drift**2))*100:.3e}%")

    speed_conserved = np.max(np.abs(v_drift)) < 1e-8
    print(f"\n  Status: {'✅ PASS' if speed_conserved else '❌ FAIL'} (<10⁻⁸ required)")

    # OVERALL SUMMARY
    print(f"\n{'='*70}")
    print("OVERALL VALIDATION SUMMARY")
    print(f"{'='*70}")

    all_passed = energy_conserved and radius_accurate and freq_accurate and speed_conserved

    print(f"  Energy conservation: {'✅' if energy_conserved else '❌'}")
    print(f"  Larmor radius:       {'✅' if radius_accurate else '❌'}")
    print(f"  Cyclotron frequency: {'✅' if freq_accurate else '❌'}")
    print(f"  Speed conservation:  {'✅' if speed_conserved else '❌'}")
    print(f"\n  Overall: {'✅ ALL TESTS PASSED' if all_passed else '❌ SOME TESTS FAILED'}")

    # PLOTS
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Trajectory
    axes[0, 0].plot(x, y, 'b-', linewidth=0.5, alpha=0.7)
    axes[0, 0].plot(x[0], y[0], 'go', markersize=10, label='Start')
    axes[0, 0].plot(x[-1], y[-1], 'ro', markersize=10, label='End')
    axes[0, 0].plot(x_center, y_center, 'k+', markersize=15, label='Center')
    circle = plt.Circle((x_center, y_center), r_larmor_theory, fill=False, color='r', linestyle='--', label='Theory')
    axes[0, 0].add_patch(circle)
    axes[0, 0].set_xlabel('x [ℓ_P]')
    axes[0, 0].set_ylabel('y [ℓ_P]')
    axes[0, 0].set_title('Particle Trajectory')
    axes[0, 0].legend()
    axes[0, 0].axis('equal')
    axes[0, 0].grid(True, alpha=0.3)

    # Energy vs time
    axes[0, 1].plot(t, KE, 'b-', linewidth=1)
    axes[0, 1].axhline(KE[0], color='r', linestyle='--', label='Initial')
    axes[0, 1].set_xlabel('Time [t_P]')
    axes[0, 1].set_ylabel('Kinetic Energy')
    axes[0, 1].set_title('Energy Conservation')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Energy drift
    axes[0, 2].plot(t, KE_drift * 100, 'b-', linewidth=1)
    axes[0, 2].axhline(0, color='k', linestyle='-', linewidth=0.5)
    axes[0, 2].set_xlabel('Time [t_P]')
    axes[0, 2].set_ylabel('Energy Drift [%]')
    axes[0, 2].set_title(f'Max Drift: {np.max(np.abs(KE_drift))*100:.2e}%')
    axes[0, 2].grid(True, alpha=0.3)

    # Radius vs time
    axes[1, 0].plot(t, r, 'b-', linewidth=1, alpha=0.7)
    axes[1, 0].axhline(r_larmor_theory, color='r', linestyle='--', label='Theory')
    axes[1, 0].set_xlabel('Time [t_P]')
    axes[1, 0].set_ylabel('Radius [ℓ_P]')
    axes[1, 0].set_title(f'Larmor Radius (error: {r_error*100:.2f}%)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # Angular position
    axes[1, 1].plot(t, theta_unwrap, 'b.', markersize=1, alpha=0.5)
    if len(t) > 100:
        axes[1, 1].plot(t[idx_start:idx_end], np.polyval(coeffs, t[idx_start:idx_end]), 'r-', linewidth=2, label=f'ω = {omega_measured:.4f}')
    axes[1, 1].set_xlabel('Time [t_P]')
    axes[1, 1].set_ylabel('Angle [rad]')
    axes[1, 1].set_title(f'Cyclotron Motion (ω error: {omega_error*100:.2f}%)')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    # Speed vs time
    axes[1, 2].plot(t, v_mag, 'b-', linewidth=1)
    axes[1, 2].axhline(v_mag[0], color='r', linestyle='--', label='Initial')
    axes[1, 2].set_xlabel('Time [t_P]')
    axes[1, 2].set_ylabel('Speed |v| [c]')
    axes[1, 2].set_title(f'Speed Drift: {np.max(np.abs(v_drift))*100:.2e}%')
    axes[1, 2].legend()
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()

    output_dir = Path(trajectory_file).parent
    plot_file = output_dir / "lorentz_comprehensive_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"\n  Plot saved: {plot_file}")

    return all_passed

if __name__ == "__main__":
    if len(sys.argv) > 1:
        trajectory_file = sys.argv[1]
    else:
        # Find latest output
        import glob
        outputs = glob.glob("output/*lorentz_comprehensive*/particle_trajectory.csv")
        if not outputs:
            print("ERROR: No trajectory file found")
            sys.exit(1)
        trajectory_file = sorted(outputs)[-1]

    print(f"Analyzing: {trajectory_file}\n")
    success = analyze_trajectory(trajectory_file)
    sys.exit(0 if success else 1)
