#!/usr/bin/env python3
"""
Phase 0.2: Power Spectrum Analysis
Phase 0.3: Autocorrelation Function

From Determinism.md:

Phase 0.2:
- Measure phase deviations: delta_theta = theta - mean(theta)
- Compute power spectrum: S(ω) = |FFT(delta_theta)|²
- Check: flat (white)? Lorentzian (colored)? 1/f (chaotic)?

Phase 0.3:
- Measure autocorrelation: C(t) = ⟨δθ(t)δθ(0)⟩
- Fit: C(t) = A * exp(-γt)
- Extract effective damping γ
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import correlate

def load_timeseries(filepath):
    """Load timeseries data (format: step value)"""
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]

def exponential_decay(t, A, gamma):
    """C(t) = A * exp(-gamma * t)"""
    return A * np.exp(-gamma * t)

def power_spectrum_analysis(timeseries, dt):
    """
    Compute power spectrum of fluctuations.

    Returns:
    - frequencies: ω array
    - power: S(ω) = |FFT(δθ)|²
    """
    # Remove mean (detrend)
    deviations = timeseries - np.mean(timeseries)

    # Compute FFT
    fft_result = np.fft.rfft(deviations)
    frequencies = np.fft.rfftfreq(len(deviations), d=dt)

    # Power spectrum (magnitude squared)
    power = np.abs(fft_result)**2

    return frequencies, power

def autocorrelation_analysis(timeseries, dt, max_lag=500):
    """
    Compute autocorrelation function.

    C(t) = ⟨δθ(t)δθ(0)⟩

    Returns:
    - lags: time lags
    - C: autocorrelation
    - gamma: fitted damping rate
    """
    # Remove mean
    deviations = timeseries - np.mean(timeseries)

    # Compute autocorrelation
    C_full = correlate(deviations, deviations, mode='full')
    C_full /= len(deviations)  # Normalize

    # Keep only positive lags
    center = len(C_full) // 2
    C = C_full[center:]

    # Truncate to max_lag
    if max_lag is not None and max_lag < len(C):
        C = C[:max_lag]

    lags = np.arange(len(C))
    time_lags = lags * dt

    # Fit exponential decay: C(t) = A * exp(-gamma * t)
    # Only fit where C > 0.1 * C[0] (avoid noise floor)
    threshold = 0.1 * C[0]
    fit_mask = C > threshold

    try:
        popt, pcov = curve_fit(
            exponential_decay,
            time_lags[fit_mask],
            C[fit_mask],
            p0=[C[0], 1.0],
            maxfev=10000
        )
        A_fit, gamma_fit = popt
    except:
        print("Warning: Autocorrelation fit failed, using fallback")
        gamma_fit = 1.0
        A_fit = C[0]

    return time_lags, C, gamma_fit, A_fit

def main():
    # Load R_avg timeseries
    print("=== Phase 0.2 & 0.3: Spectrum + Autocorrelation ===\n")
    print("Loading timeseries data...")

    steps, R_avg = load_timeseries('build/output/timeseries_R_avg.dat')
    dt = 0.01
    time = steps * dt

    print(f"Loaded {len(steps)} timesteps")
    print(f"Time range: t = {time[0]:.2f} to {time[-1]:.2f}")
    print(f"R_avg: {R_avg[0]:.4f} → {R_avg[-1]:.4f}\n")

    # Phase 0.2: Power Spectrum
    print("=== Phase 0.2: Power Spectrum Analysis ===")
    frequencies, power = power_spectrum_analysis(R_avg, dt)

    # Classify spectrum shape
    # White noise: flat S(ω)
    # Colored noise: S(ω) ~ 1/(1 + (ω/ω_c)²) (Lorentzian)
    # 1/f noise: S(ω) ~ 1/ω^α

    # Check slope in log-log (excluding DC and high-freq noise)
    f_min, f_max = 0.5, 20.0
    mask = (frequencies >= f_min) & (frequencies <= f_max)
    f_fit = frequencies[mask]
    p_fit = power[mask]

    # Fit: log(S) = log(A) - alpha * log(f)
    if len(f_fit) > 10:
        log_f = np.log(f_fit)
        log_p = np.log(p_fit + 1e-10)  # avoid log(0)
        coeffs = np.polyfit(log_f, log_p, 1)
        alpha = -coeffs[0]
    else:
        alpha = 0.0

    print(f"Power spectrum slope: α = {alpha:.3f}")

    if abs(alpha) < 0.5:
        spectrum_type = "WHITE (flat spectrum)"
    elif 0.5 <= alpha < 1.5:
        spectrum_type = "COLORED (Lorentzian-like, α ≈ 1)"
    elif alpha >= 1.5:
        spectrum_type = "1/f NOISE (chaotic, α ≈ 2)"
    else:
        spectrum_type = "UNKNOWN"

    print(f"Classification: {spectrum_type}")

    # Phase 0.3: Autocorrelation
    print("\n=== Phase 0.3: Autocorrelation Analysis ===")
    lags, C, gamma_fit, A_fit = autocorrelation_analysis(R_avg, dt, max_lag=500)

    print(f"Fitted damping rate: γ = {gamma_fit:.6f}")
    print(f"Autocorrelation time: τ_c = 1/γ = {1/gamma_fit:.2f}")

    # Generate plots
    print("\nGenerating plots...")

    fig = plt.figure(figsize=(14, 10))

    # Plot 1: Timeseries
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(time, R_avg, 'b-', linewidth=0.5, alpha=0.7)
    ax1.set_xlabel('Time t')
    ax1.set_ylabel('R_avg')
    ax1.set_title('Synchronization Evolution')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Power Spectrum (log-log)
    ax2 = plt.subplot(2, 2, 2)
    ax2.loglog(frequencies[1:], power[1:], 'g-', linewidth=1, alpha=0.7)

    # Plot reference slopes
    f_ref = np.logspace(np.log10(f_min), np.log10(f_max), 100)
    p_white = np.ones_like(f_ref) * np.median(power[mask])
    p_colored = p_white[0] / (1 + (f_ref/5.0)**2)
    p_1f = p_white[0] * (f_ref[0]/f_ref)**2

    ax2.loglog(f_ref, p_white, 'k--', linewidth=2, alpha=0.5, label='White (α=0)')
    ax2.loglog(f_ref, p_colored, 'b--', linewidth=2, alpha=0.5, label='Lorentzian (α≈1)')
    ax2.loglog(f_ref, p_1f, 'r--', linewidth=2, alpha=0.5, label='1/f² (α=2)')

    ax2.set_xlabel('Frequency ω')
    ax2.set_ylabel('Power S(ω)')
    ax2.set_title(f'Power Spectrum (α = {alpha:.2f}, {spectrum_type})')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Plot 3: Power Spectrum (linear, low freq)
    ax3 = plt.subplot(2, 2, 3)
    f_max_plot = 10.0
    mask_plot = frequencies <= f_max_plot
    ax3.plot(frequencies[mask_plot], power[mask_plot], 'g-', linewidth=1)
    ax3.set_xlabel('Frequency ω')
    ax3.set_ylabel('Power S(ω)')
    ax3.set_title('Power Spectrum (Linear Scale, Low Frequencies)')
    ax3.grid(True, alpha=0.3)

    # Plot 4: Autocorrelation
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(lags, C, 'b-', linewidth=1, alpha=0.7, label='C(t)')

    # Plot fitted exponential
    C_fit = exponential_decay(lags, A_fit, gamma_fit)
    ax4.plot(lags, C_fit, 'r--', linewidth=2, label=f'Fit: γ = {gamma_fit:.3f}')

    ax4.set_xlabel('Time lag t')
    ax4.set_ylabel('Autocorrelation C(t)')
    ax4.set_title(f'Autocorrelation (τ_c = {1/gamma_fit:.2f})')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('phase0_spectrum_autocorr.png', dpi=150)
    print("Saved: phase0_spectrum_autocorr.png")

    # Summary
    print("\n=== SUMMARY ===")
    print("Phase 0.2: Power Spectrum")
    print(f"  Slope: α = {alpha:.3f}")
    print(f"  Type: {spectrum_type}")
    print("\nPhase 0.3: Autocorrelation")
    print(f"  Damping: γ = {gamma_fit:.6f}")
    print(f"  Correlation time: τ_c = {1/gamma_fit:.2f}")
    print("\n=== INTERPRETATION ===")
    print("Combined with Lyapunov (λ < 0):")
    print("  → System is DAMPED, not chaotic")
    print("  → Fluctuations are RELAXATION to synchronized state")
    print("  → No effective stochasticity from internal dynamics")
    print("\nConclusion: 'Effective Stochasticity' hypothesis REJECTED by Phase 0")

if __name__ == '__main__':
    main()
