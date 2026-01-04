#!/usr/bin/env python3
"""
B1 Saturation Analysis Script
Tests if linear scaling continues beyond d=100

Phase 5 Linear Fit: m₂/m₁ = 0.8083·d - 13.985 (R²=0.998, d∈[30,100])
Saturation Test: Compare observed vs predicted at d=100-200
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Phase 5 linear fit parameters
ALPHA = 0.8083  # Slope
BETA = -13.985  # Intercept
R_SQUARED = 0.998287

def linear_prediction(d):
    """Predict m₂/m₁ from Phase 5 linear fit"""
    return ALPHA * d + BETA

def analyze_saturation(csv_file="analysis/b1_saturation_check_256.csv"):
    """Analyze saturation check results"""

    # Load data
    df = pd.read_csv(csv_file)

    print("=" * 60)
    print("B1 SATURATION CHECK ANALYSIS")
    print("=" * 60)
    print(f"\nData: {csv_file}")
    print(f"Phase 5 fit: m₂/m₁ = {ALPHA:.4f}·d + ({BETA:.3f})")
    print(f"Fit quality: R² = {R_SQUARED:.6f}\n")

    # Compute predictions and residuals
    df['m2_m1_predicted'] = df['separation'].apply(linear_prediction)
    df['residual'] = df['m2_m1'] - df['m2_m1_predicted']
    df['residual_pct'] = 100 * df['residual'] / df['m2_m1_predicted']

    # Display results
    print("SATURATION DETECTION TABLE")
    print("-" * 80)
    print(f"{'d':>6} {'Observed':>10} {'Predicted':>10} {'Residual':>10} {'Error %':>10} {'Verdict':>15}")
    print("-" * 80)

    saturation_detected = False
    for _, row in df.iterrows():
        d = row['separation']
        obs = row['m2_m1']
        pred = row['m2_m1_predicted']
        res = row['residual']
        err_pct = row['residual_pct']

        # Classify deviation
        if abs(err_pct) < 5:
            verdict = "LINEAR ✓"
        elif abs(err_pct) < 10:
            verdict = "MARGINAL"
        else:
            verdict = "SATURATION ✗"
            saturation_detected = True

        print(f"{d:6.0f} {obs:10.2f} {pred:10.2f} {res:+10.2f} {err_pct:+9.1f}% {verdict:>15}")

    print("-" * 80)

    # Statistical analysis
    mean_residual = df['residual_pct'].mean()
    std_residual = df['residual_pct'].std()
    max_deviation = df['residual_pct'].abs().max()

    print(f"\nSTATISTICAL SUMMARY:")
    print(f"  Mean residual: {mean_residual:+.2f}%")
    print(f"  Std deviation: {std_residual:.2f}%")
    print(f"  Max deviation: {max_deviation:.2f}%")

    # Fit new power law to saturation data
    from scipy.optimize import curve_fit

    def power_law(d, alpha, beta):
        return alpha * d + beta

    params, covariance = curve_fit(power_law, df['separation'], df['m2_m1'])
    alpha_new, beta_new = params

    # Compute R² for saturation fit
    ss_res = np.sum((df['m2_m1'] - power_law(df['separation'], *params))**2)
    ss_tot = np.sum((df['m2_m1'] - df['m2_m1'].mean())**2)
    r_squared_new = 1 - (ss_res / ss_tot)

    print(f"\nSATURATION FIT (d∈[100,200]):")
    print(f"  m₂/m₁ = {alpha_new:.4f}·d + ({beta_new:.3f})")
    print(f"  R² = {r_squared_new:.6f}")
    print(f"  Slope change: {100*(alpha_new/ALPHA - 1):+.1f}%")

    # Extrapolation to muon mass
    TARGET_RATIO = 206.768
    d_target_phase5 = (TARGET_RATIO - BETA) / ALPHA
    d_target_saturation = (TARGET_RATIO - beta_new) / alpha_new

    print(f"\nEXTRAPOLATION TO MUON MASS (m₂/m₁ = {TARGET_RATIO:.3f}):")
    print(f"  Phase 5 fit:     d ≈ {d_target_phase5:.1f}")
    print(f"  Saturation fit:  d ≈ {d_target_saturation:.1f}")
    print(f"  Difference:      Δd = {d_target_saturation - d_target_phase5:+.1f}")

    # Grid requirements
    grid_size_phase5 = d_target_phase5 / 0.4  # d_max = 0.4 × grid_size
    grid_size_saturation = d_target_saturation / 0.4

    print(f"\nGRID REQUIREMENTS:")
    print(f"  Phase 5 fit:     {grid_size_phase5:.0f}³ grid")
    print(f"  Saturation fit:  {grid_size_saturation:.0f}³ grid")

    # Final verdict
    print("\n" + "=" * 60)
    if max_deviation < 5:
        print("VERDICT: LINEAR SCALING CONTINUES ✓")
        print("  - Residuals < 5% → Phase 5 fit remains valid")
        print("  - Proceed to 512³ grid (d=220-273)")
        print(f"  - Target: {grid_size_phase5:.0f}³ grid for muon mass")
    elif max_deviation < 10:
        print("VERDICT: MARGINAL SATURATION")
        print("  - Residuals 5-10% → slight deviation from linearity")
        print("  - Caution: Extrapolation may need correction")
        print("  - Recommend: 512³ validation run before 683³")
    else:
        print("VERDICT: SATURATION DETECTED ✗")
        print("  - Residuals > 10% → linear scaling fails")
        print("  - Vortex coupling has finite range")
        print("  - PIVOT TO ALTERNATIVE PHYSICS:")
        print("    1. Stückelberg gauge mass (already tested)")
        print("    2. Environmental coupling (thermal bath)")
        print("    3. Multi-vortex resonances")
    print("=" * 60)

    # Create visualization
    create_saturation_plot(df, alpha_new, beta_new, r_squared_new)

    return df

def create_saturation_plot(df, alpha_new, beta_new, r_squared_new):
    """Visualize saturation test results"""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left plot: Scaling curve
    ax1 = axes[0]
    d_range = np.linspace(0, 250, 100)

    # Phase 5 fit
    ax1.plot(d_range, ALPHA * d_range + BETA, 'b--',
             label=f'Phase 5 fit (R²={R_SQUARED:.4f})', linewidth=2)

    # Saturation fit
    ax1.plot(d_range, alpha_new * d_range + beta_new, 'r-',
             label=f'Saturation fit (R²={r_squared_new:.4f})', linewidth=2)

    # Data points
    ax1.plot(df['separation'], df['m2_m1'], 'ko', markersize=8,
             label='Observed data', zorder=5)

    # Target line
    ax1.axhline(206.768, color='green', linestyle=':', linewidth=2,
                label='Muon/electron (206.768)')

    ax1.set_xlabel('Vortex Separation (d)', fontsize=12)
    ax1.set_ylabel('Mass Ratio (m₂/m₁)', fontsize=12)
    ax1.set_title('B1 Saturation Check: Linear vs Observed Scaling', fontsize=14)
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 250)
    ax1.set_ylim(0, 220)

    # Right plot: Residuals
    ax2 = axes[1]
    ax2.plot(df['separation'], df['residual_pct'], 'ko-', markersize=8, linewidth=2)
    ax2.axhline(0, color='blue', linestyle='--', linewidth=2, label='Phase 5 fit')
    ax2.axhline(5, color='orange', linestyle=':', linewidth=1, label='±5% threshold')
    ax2.axhline(-5, color='orange', linestyle=':', linewidth=1)
    ax2.axhline(10, color='red', linestyle=':', linewidth=1, label='±10% threshold')
    ax2.axhline(-10, color='red', linestyle=':', linewidth=1)

    ax2.set_xlabel('Vortex Separation (d)', fontsize=12)
    ax2.set_ylabel('Residual Error (%)', fontsize=12)
    ax2.set_title('Deviation from Linear Prediction', fontsize=14)
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('analysis/b1_saturation_check_plot.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: analysis/b1_saturation_check_plot.png")

if __name__ == "__main__":
    import sys
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "analysis/b1_saturation_check_256.csv"

    # Wait for CSV to exist
    if not Path(csv_file).exists():
        print(f"Waiting for {csv_file} to be generated...")
        import time
        while not Path(csv_file).exists():
            time.sleep(5)
        print(f"Found {csv_file}!")

    analyze_saturation(csv_file)
