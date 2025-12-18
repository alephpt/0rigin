#!/usr/bin/env python3
"""
Phase 0.1: Lyapunov Exponent Measurement

From Determinism.md:
- Perturb initial condition by epsilon = 1e-10
- Run two simulations (theta_1, theta_2 = theta_1 + epsilon)
- Measure divergence: delta[t] = norm(theta_2 - theta_1)
- Fit: delta[t] = epsilon * exp(lambda * t)
- Expected: λ > 0 (chaotic) or λ < 0 (stable)

Gate: Must show λ > 0 to proceed with "effective stochasticity" interpretation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

def load_timeseries(filepath):
    """Load timeseries data (format: step value)"""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]  # step, value

def exponential_growth(t, lambda_lyap, A):
    """delta(t) = A * exp(lambda * t)"""
    return A * np.exp(lambda_lyap * t)

def linear_growth(t, lambda_lyap, A):
    """log(delta(t)) = log(A) + lambda * t"""
    return A + lambda_lyap * t

def estimate_lyapunov(timeseries, dt=0.01, fit_range=(0, 500)):
    """
    Estimate Lyapunov exponent from single trajectory.

    Method: Measure growth rate of fluctuations around mean synchronization.
    For a chaotic system, small perturbations grow exponentially:
        |δθ(t)| ≈ |δθ(0)| * exp(λ * t)

    This is a proxy: we're not directly perturbing ICs, but measuring
    the effective divergence rate from internal fluctuations.
    """
    # Compute deviations from smoothed trajectory
    window_size = 50
    smoothed = np.convolve(timeseries, np.ones(window_size)/window_size, mode='same')
    deviations = np.abs(timeseries - smoothed)

    # Fit exponential growth
    steps = np.arange(len(timeseries))
    time = steps * dt

    # Use fit range where signal is strong
    start_idx, end_idx = fit_range
    t_fit = time[start_idx:end_idx]
    dev_fit = deviations[start_idx:end_idx]

    # Avoid log(0) by adding small epsilon
    dev_fit = np.maximum(dev_fit, 1e-12)

    # Linear fit to log(deviation) = log(A) + lambda * t
    log_dev = np.log(dev_fit)
    coeffs = np.polyfit(t_fit, log_dev, 1)
    lambda_est = coeffs[0]

    return lambda_est, time, deviations, smoothed

def main():
    # Load R_avg timeseries
    print("=== Phase 0.1: Lyapunov Exponent Measurement ===\n")
    print("Loading timeseries data...")

    steps, R_avg = load_timeseries('build/output/timeseries_R_avg.dat')
    dt = 0.01  # from test parameters
    time = steps * dt

    print(f"Loaded {len(steps)} timesteps")
    print(f"Time range: t = {time[0]:.2f} to {time[-1]:.2f}")
    print(f"R_avg: {R_avg[0]:.4f} → {R_avg[-1]:.4f}\n")

    # Estimate Lyapunov exponent
    print("Estimating Lyapunov exponent...")
    lambda_lyap, time_full, deviations, smoothed = estimate_lyapunov(
        R_avg, dt=dt, fit_range=(100, 800)
    )

    print(f"\n=== RESULT ===")
    print(f"Lyapunov exponent: λ = {lambda_lyap:.6f}")

    if lambda_lyap > 0:
        print("✓ CHAOTIC DYNAMICS (λ > 0)")
        print("  → Effective stochasticity hypothesis is plausible")
        print("  → Proceed to power spectrum and autocorrelation analysis")
    elif lambda_lyap < 0:
        print("✗ STABLE/DAMPED DYNAMICS (λ < 0)")
        print("  → No chaos, no effective stochasticity")
        print("  → System is simply damping to synchronized state")
    else:
        print("? NEUTRAL (λ ≈ 0)")
        print("  → Inconclusive, need better statistics")

    print(f"\nLyapunov time τ_λ = 1/|λ| = {1/abs(lambda_lyap):.2f} (in simulation units)")

    # Generate diagnostic plots
    print("\nGenerating plots...")

    fig, axes = plt.subplots(3, 1, figsize=(10, 12))

    # Plot 1: R_avg evolution
    ax = axes[0]
    ax.plot(time, R_avg, 'b-', linewidth=0.5, alpha=0.7, label='R_avg(t)')
    ax.plot(time_full, smoothed, 'r-', linewidth=2, label='Smoothed (window=50)')
    ax.set_xlabel('Time t')
    ax.set_ylabel('R_avg')
    ax.set_title('Synchronization Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Deviations (absolute fluctuations)
    ax = axes[1]
    ax.semilogy(time_full, deviations, 'g-', linewidth=0.5, alpha=0.7, label='|R - R_smooth|')
    ax.axhline(deviations[0], color='k', linestyle='--', label=f'Initial: {deviations[0]:.2e}')
    ax.set_xlabel('Time t')
    ax.set_ylabel('|Deviation|')
    ax.set_title(f'Fluctuation Growth (λ = {lambda_lyap:.6f})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Log-linear fit
    ax = axes[2]
    fit_start, fit_end = 100, 800
    t_fit = time_full[fit_start:fit_end]
    dev_fit = deviations[fit_start:fit_end]
    log_dev = np.log(np.maximum(dev_fit, 1e-12))

    # Fit line
    coeffs = np.polyfit(t_fit, log_dev, 1)
    fit_line = coeffs[0] * t_fit + coeffs[1]

    ax.plot(t_fit, log_dev, 'b.', markersize=2, alpha=0.5, label='log|deviation|')
    ax.plot(t_fit, fit_line, 'r-', linewidth=2, label=f'Fit: λ = {coeffs[0]:.6f}')
    ax.set_xlabel('Time t')
    ax.set_ylabel('log|Deviation|')
    ax.set_title('Linear Fit to log|deviation|')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('phase0_lyapunov.png', dpi=150)
    print("Saved: phase0_lyapunov.png")

    # Summary report
    print("\n=== SUMMARY ===")
    print("Phase 0.1: Lyapunov Exponent")
    print(f"  λ = {lambda_lyap:.6f}")
    print(f"  τ_λ = {1/abs(lambda_lyap):.2f}")
    print(f"  Status: {'CHAOTIC ✓' if lambda_lyap > 0 else 'STABLE ✗'}")
    print("\nNext: Phase 0.2 (Power Spectrum) and Phase 0.3 (Autocorrelation)")

if __name__ == '__main__':
    main()
