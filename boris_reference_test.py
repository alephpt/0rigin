#!/usr/bin/env python3
"""
Reference Boris algorithm implementation to verify correct behavior.
The Boris algorithm should produce PERFECT circular orbits in uniform B-field.
"""

import numpy as np
import matplotlib.pyplot as plt

def boris_pusher(x, v, E, B, q, m, dt):
    """
    Standard Boris pusher for charged particle in EM field.

    Args:
        x: position [x, y, z]
        v: velocity [vx, vy, vz]
        E: electric field [Ex, Ey, Ez]
        B: magnetic field [Bx, By, Bz]
        q: charge
        m: mass
        dt: timestep

    Returns:
        x_new, v_new: updated position and velocity
    """
    # Constants
    q_over_m = q / m

    # Half acceleration from E field
    v_minus = v + 0.5 * q_over_m * E * dt

    # Rotation from B field
    t = 0.5 * q_over_m * B * dt
    s = 2.0 * t / (1.0 + np.dot(t, t))

    # Rotation step
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)

    # Half acceleration from E field
    v_new = v_plus + 0.5 * q_over_m * E * dt

    # Position update
    x_new = x + v_new * dt

    return x_new, v_new


def simulate_cyclotron():
    """Simulate charged particle in uniform B-field."""

    # Physical parameters (matching SMFT test)
    q = -1.0        # charge
    m = 0.1         # mass
    B_z = 0.1       # magnetic field
    v0 = 0.01       # initial speed
    dt = 0.0001     # timestep

    # Initial conditions
    x0 = np.array([55.0, 50.0, 0.0])  # position
    v0 = np.array([0.0, v0, 0.0])      # velocity

    # Fields
    E = np.array([0.0, 0.0, 0.0])      # no electric field
    B = np.array([0.0, 0.0, B_z])      # uniform B in z

    # Theoretical predictions
    omega_c = abs(q * B_z / m)  # cyclotron frequency
    r_L = m * np.linalg.norm(v0[:2]) / (abs(q) * B_z)  # Larmor radius
    T = 2 * np.pi / omega_c  # period

    print("Theoretical predictions:")
    print(f"  Cyclotron frequency: ω_c = {omega_c:.6f} rad/t_P")
    print(f"  Larmor radius: r_L = {r_L:.6f} ℓ_P")
    print(f"  Period: T = {T:.6f} t_P")

    # Simulate
    n_steps = 100000
    record_every = 5

    trajectory = []
    x, v = x0.copy(), v0.copy()

    for step in range(n_steps):
        if step % record_every == 0:
            trajectory.append({
                't': step * dt,
                'x': x[0], 'y': x[1],
                'vx': v[0], 'vy': v[1],
                'speed': np.linalg.norm(v[:2])
            })

        x, v = boris_pusher(x, v, E, B, q, m, dt)

    # Analyze results
    traj = trajectory
    t = np.array([p['t'] for p in traj])
    x = np.array([p['x'] for p in traj])
    y = np.array([p['y'] for p in traj])
    speed = np.array([p['speed'] for p in traj])

    # Find orbit center
    x_center = np.mean(x)
    y_center = np.mean(y)

    # Compute radius
    r = np.sqrt((x - x_center)**2 + (y - y_center)**2)
    r_mean = np.mean(r)
    r_std = np.std(r)

    print("\nNumerical results:")
    print(f"  Orbit center: ({x_center:.6f}, {y_center:.6f})")
    print(f"  Mean radius: {r_mean:.6f} ± {r_std:.6f} ℓ_P")
    print(f"  Radius variation: {r_std/r_mean*100:.6f}%")
    print(f"  Min radius: {np.min(r):.6f} ℓ_P")
    print(f"  Max radius: {np.max(r):.6f} ℓ_P")
    print(f"  Speed conservation: {(speed[-1] - speed[0])/speed[0]*100:.10f}%")

    # Count orbits
    theta = np.arctan2(y - y_center, x - x_center)
    theta_unwrap = np.unwrap(theta)
    n_orbits = (theta_unwrap[-1] - theta_unwrap[0]) / (2 * np.pi)
    print(f"  Number of orbits: {n_orbits:.3f}")

    # Plot trajectory
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Trajectory plot
    ax1.plot(x, y, 'b-', linewidth=0.5, alpha=0.7)
    ax1.plot(x[0], y[0], 'go', markersize=8, label='Start')
    ax1.plot(x[-1], y[-1], 'ro', markersize=8, label='End')
    ax1.plot(x_center, y_center, 'k+', markersize=12, label='Center')
    circle = plt.Circle((x_center, y_center), r_L, fill=False,
                        color='r', linestyle='--', label=f'Theory r={r_L:.4f}')
    ax1.add_patch(circle)
    ax1.set_xlabel('x [ℓ_P]')
    ax1.set_ylabel('y [ℓ_P]')
    ax1.set_title('Particle Trajectory')
    ax1.legend()
    ax1.axis('equal')
    ax1.grid(True, alpha=0.3)

    # Radius vs time
    ax2.plot(t, r, 'b-', linewidth=1)
    ax2.axhline(r_L, color='r', linestyle='--', label='Theory')
    ax2.set_xlabel('Time [t_P]')
    ax2.set_ylabel('Radius [ℓ_P]')
    ax2.set_title(f'Orbit Radius (variation: {r_std/r_mean*100:.2f}%)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('boris_reference_test.png', dpi=150)
    print("\nPlot saved: boris_reference_test.png")

    # Verdict
    if r_std / r_mean < 0.001:  # <0.1% variation
        print("\n✅ PASS: Boris algorithm working correctly")
    else:
        print(f"\n❌ FAIL: Boris algorithm has {r_std/r_mean*100:.2f}% radius variation")
        print("   Expected: <0.1% for proper implementation")


if __name__ == "__main__":
    simulate_cyclotron()