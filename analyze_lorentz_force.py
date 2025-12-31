#!/usr/bin/env python3
"""
Analyze Lorentz force test results from console output.
Extract particle positions and verify cyclotron motion.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import re

def analyze_console_log(log_file):
    """Extract particle positions from console log."""

    positions = []
    velocities = []

    with open(log_file, 'r') as f:
        for line in f:
            # Match lines like: "Particle at step 1000: pos=(54.9996,50.002), vel=(-0.00357519,0.00964126)"
            match = re.search(r'Particle at step (\d+): pos=\(([\d.-]+),([\d.-]+)\), vel=\(([\d.-]+),([\d.-]+)\)', line)
            if match:
                step = int(match.group(1))
                x = float(match.group(2))
                y = float(match.group(3))
                vx = float(match.group(4))
                vy = float(match.group(5))
                positions.append([x, y])
                velocities.append([vx, vy])

    if not positions:
        print("No particle data found in log file")
        return None

    positions = np.array(positions)
    velocities = np.array(velocities)

    # Compute center of circular motion
    cx = np.mean(positions[:, 0])
    cy = np.mean(positions[:, 1])

    # Compute radius
    radii = np.sqrt((positions[:, 0] - cx)**2 + (positions[:, 1] - cy)**2)
    r_mean = np.mean(radii)
    r_std = np.std(radii)

    # Compute velocity magnitude
    v_mag = np.sqrt(velocities[:, 0]**2 + velocities[:, 1]**2)
    v_mean = np.mean(v_mag)

    print("\n=== Lorentz Force Analysis ===")
    print(f"Number of data points: {len(positions)}")
    print(f"Center of motion: ({cx:.4f}, {cy:.4f}) ℓ_P")
    print(f"Mean radius: {r_mean:.4f} ± {r_std:.4f} ℓ_P")
    print(f"Mean velocity: {v_mean:.4f} c")

    # With B = 0.1, m = 0.1, q = -1, v = 0.01
    # Theory: r_L = mv/(qB) = 0.1 * 0.01 / (1 * 0.1) = 0.01 ℓ_P
    # Theory: ω = qB/m = 1 * 0.1 / 0.1 = 1 rad/τ_P

    B_field = 0.1  # From our test setup
    mass = 0.1
    charge = -1.0

    r_theory = mass * v_mean / (abs(charge) * B_field)
    omega_theory = abs(charge) * B_field / mass

    print(f"\nTheoretical predictions (B = {B_field}):")
    print(f"  Larmor radius: r_L = {r_theory:.4f} ℓ_P")
    print(f"  Cyclotron frequency: ω = {omega_theory:.4f} rad/τ_P")

    print(f"\nComparison:")
    r_error = abs(r_mean - r_theory) / r_theory * 100
    print(f"  Radius error: {r_error:.1f}%")

    # Check if it's a good circle
    circularity = r_std / r_mean
    print(f"  Circularity (std/mean): {circularity:.3f}")

    if r_error < 20 and circularity < 0.1:
        print("\n✓ Lorentz force validation: PASS")
        print("  Particle follows expected cyclotron motion")
    else:
        print("\n✗ Lorentz force validation: FAIL")
        if r_error >= 20:
            print(f"  Radius error too large: {r_error:.1f}% > 20%")
        if circularity >= 0.1:
            print(f"  Motion not circular enough: {circularity:.3f} > 0.1")

    # Plot trajectory
    if len(positions) > 2:
        plt.figure(figsize=(10, 5))

        plt.subplot(1, 2, 1)
        plt.plot(positions[:, 0], positions[:, 1], 'b-', alpha=0.7, label='Trajectory')
        plt.plot(positions[0, 0], positions[0, 1], 'go', markersize=8, label='Start')
        plt.plot(positions[-1, 0], positions[-1, 1], 'ro', markersize=8, label='End')
        plt.plot(cx, cy, 'kx', markersize=10, label='Center')

        # Theory circle
        theta = np.linspace(0, 2*np.pi, 100)
        x_theory = cx + r_theory * np.cos(theta)
        y_theory = cy + r_theory * np.sin(theta)
        plt.plot(x_theory, y_theory, 'r--', alpha=0.5, label=f'Theory (r={r_theory:.3f})')

        plt.xlabel('x (ℓ_P)')
        plt.ylabel('y (ℓ_P)')
        plt.title('Particle Trajectory')
        plt.legend()
        plt.axis('equal')
        plt.grid(True, alpha=0.3)

        plt.subplot(1, 2, 2)
        plt.plot(velocities[:, 0], velocities[:, 1], 'b-', alpha=0.7, label='Velocity')
        plt.plot(velocities[0, 0], velocities[0, 1], 'go', markersize=8, label='Start')
        plt.plot(velocities[-1, 0], velocities[-1, 1], 'ro', markersize=8, label='End')
        plt.xlabel('vx (c)')
        plt.ylabel('vy (c)')
        plt.title('Velocity Phase Space')
        plt.legend()
        plt.axis('equal')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('lorentz_force_trajectory.png', dpi=150)
        print("\nTrajectory plot saved to lorentz_force_trajectory.png")
        plt.show()

    return positions, velocities

if __name__ == "__main__":
    if len(sys.argv) > 1:
        log_file = sys.argv[1]
    else:
        # Find latest log file
        import glob
        import os

        logs = glob.glob("output/*/console.log")
        if logs:
            log_file = max(logs, key=os.path.getctime)
            print(f"Using latest log: {log_file}")
        else:
            print("No log files found")
            sys.exit(1)

    analyze_console_log(log_file)