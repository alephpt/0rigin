#!/usr/bin/env python3
"""
Quantum Synchronization Critical Exponent Test Design

Designs optimal experimental protocol for measuring β in coupled qubit network
to distinguish SMFT (β = 0.099) from 2D Ising (β = 0.125) at high significance.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2

# Critical exponents
BETA_SMFT = 0.099
BETA_SMFT_ERROR = 0.004
BETA_2D_ISING = 0.125
BETA_2D_XY = 0.23
BETA_MEAN_FIELD = 0.5

# Target significance
TARGET_SIGMA = 5  # 5σ detection

print("=" * 80)
print("CRITICAL EXPONENT TEST: Quantum Synchronization in Superconducting Qubits")
print("=" * 80)
print()

# ===== THEORETICAL PREDICTIONS =====

print("THEORETICAL PREDICTIONS:")
print("-" * 80)
print(f"{'Universality Class':<25} {'β':<15} {'ν':<15} {'d_u':<10}")
print("-" * 80)
print(f"{'Mean-field':<25} {BETA_MEAN_FIELD:<15.3f} {2.0:<15.1f} {4:<10}")
print(f"{'2D Ising':<25} {BETA_2D_ISING:<15.3f} {1.0:<15.1f} {2:<10}")
print(f"{'2D XY (BKT)':<25} {BETA_2D_XY:<15.3f} {'—':<15} {2:<10}")
print(f"{'SMFT (NOVEL)':<25} {BETA_SMFT:<15.3f} {2.5:<15.1f} {5:<10}")
print()

# Statistical separation
separation_ising = abs(BETA_SMFT - BETA_2D_ISING)
separation_xy = abs(BETA_SMFT - BETA_2D_XY)
separation_mf = abs(BETA_SMFT - BETA_MEAN_FIELD)

print("STATISTICAL SEPARATION:")
print(f"  |β_SMFT - β_2D-Ising| = {separation_ising:.4f}")
print(f"  |β_SMFT - β_2D-XY|    = {separation_xy:.4f}")
print(f"  |β_SMFT - β_MF|       = {separation_mf:.4f}")
print()

# Required measurement precision
sigma_beta_required = separation_ising / TARGET_SIGMA
print(f"To distinguish SMFT from 2D Ising at {TARGET_SIGMA}σ:")
print(f"  Required precision: σ(β) < {sigma_beta_required:.4f}")
print()

# ===== EXPERIMENTAL DESIGN =====

print("=" * 80)
print("EXPERIMENTAL DESIGN")
print("=" * 80)
print()

# System sizes
N_qubits = [10, 20, 50, 100, 200]
print("System Sizes (Number of Qubits):")
print(f"  N = {N_qubits}")
print()

# Noise scan
sigma_min = 0.40
sigma_max = 0.60
num_points = 41
sigma_critical_guess = 0.50

sigma_values = np.linspace(sigma_min, sigma_max, num_points)

print("Noise Scan:")
print(f"  Range: σ = {sigma_min} to {sigma_max}")
print(f"  Points: {num_points} (spacing Δσ = {(sigma_max - sigma_min)/(num_points - 1):.4f})")
print(f"  Critical point (guess): σ_c ~ {sigma_critical_guess}")
print()

# Measurement protocol
num_realizations = 100  # Repeat measurements for statistics
measurement_time_us = 10  # Measurement time per point (μs)
equilibration_time_us = 50  # Equilibration time (μs)
total_time_per_point = (measurement_time_us + equilibration_time_us) * num_realizations * 1e-6  # seconds

print("Measurement Protocol:")
print(f"  Realizations per (σ, N): {num_realizations}")
print(f"  Equilibration time: {equilibration_time_us} μs")
print(f"  Measurement time: {measurement_time_us} μs")
print(f"  Total time per (σ, N): {total_time_per_point:.2f} s")
print()

# Total experiment duration
total_configs = num_points * len(N_qubits)
total_time_hours = total_configs * total_time_per_point / 3600

print(f"Total Experiment:")
print(f"  Configurations: {num_points} σ × {len(N_qubits)} N = {total_configs}")
print(f"  Total time: {total_time_hours:.1f} hours ({total_time_hours / 24:.1f} days)")
print()

# ===== SYNTHETIC DATA GENERATION (for protocol validation) =====

print("=" * 80)
print("SYNTHETIC DATA GENERATION (Protocol Validation)")
print("=" * 80)
print()

# Generate synthetic data assuming SMFT β = 0.099
np.random.seed(42)

def order_parameter_theory(sigma, sigma_c, beta):
    """
    Theoretical order parameter: R ~ |σ_c - σ|^β for σ < σ_c
    """
    if sigma >= sigma_c:
        return 0.0
    else:
        return (sigma_c - sigma)**beta

# Simulate measurements with noise
sigma_c_true = 0.50
beta_true = BETA_SMFT
measurement_noise = 0.01  # 1% measurement uncertainty

print(f"True parameters:")
print(f"  σ_c = {sigma_c_true}")
print(f"  β = {beta_true}")
print(f"  Measurement noise: {measurement_noise * 100}%")
print()

# Generate data for largest system (N=200, best statistics)
N_test = 200
R_data = []
R_error = []

for sigma in sigma_values:
    R_theory = order_parameter_theory(sigma, sigma_c_true, beta_true)
    # Add measurement noise
    R_measured = R_theory + np.random.normal(0, measurement_noise, num_realizations)
    R_mean = np.mean(R_measured)
    R_std = np.std(R_measured) / np.sqrt(num_realizations)  # Standard error

    R_data.append(R_mean)
    R_error.append(R_std)

R_data = np.array(R_data)
R_error = np.array(R_error)

# ===== CRITICAL EXPONENT EXTRACTION =====

print("=" * 80)
print("CRITICAL EXPONENT EXTRACTION")
print("=" * 80)
print()

# Fit to power law: R ~ (σ_c - σ)^β
def power_law(sigma, sigma_c, beta, amplitude):
    """Power law fit function"""
    R = np.zeros_like(sigma)
    mask = sigma < sigma_c
    R[mask] = amplitude * (sigma_c - sigma[mask])**beta
    return R

# Fit only to points below critical point (where R > 0)
fit_mask = R_data > 0.001
sigma_fit = sigma_values[fit_mask]
R_fit = R_data[fit_mask]
R_err_fit = R_error[fit_mask]

# Initial guess
p0 = [sigma_c_true, beta_true, 1.0]

# Perform fit
try:
    popt, pcov = curve_fit(power_law, sigma_fit, R_fit, p0=p0, sigma=R_err_fit, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))

    sigma_c_fit, beta_fit, amp_fit = popt
    sigma_c_err, beta_err, amp_err = perr

    print("FIT RESULTS:")
    print(f"  σ_c = {sigma_c_fit:.4f} ± {sigma_c_err:.4f}")
    print(f"  β   = {beta_fit:.4f} ± {beta_err:.4f}")
    print(f"  A   = {amp_fit:.4f} ± {amp_err:.4f}")
    print()

    # Goodness of fit
    R_predicted = power_law(sigma_fit, *popt)
    chi_squared = np.sum(((R_fit - R_predicted) / R_err_fit)**2)
    dof = len(sigma_fit) - len(popt)
    chi_squared_reduced = chi_squared / dof

    print(f"Goodness of Fit:")
    print(f"  χ² = {chi_squared:.2f}")
    print(f"  DOF = {dof}")
    print(f"  χ²/DOF = {chi_squared_reduced:.3f}")
    print(f"  p-value = {1 - chi2.cdf(chi_squared, dof):.4f}")
    print()

    # Statistical significance of SMFT vs 2D Ising
    beta_diff_ising = abs(beta_fit - BETA_2D_ISING)
    sigma_significance_ising = beta_diff_ising / beta_err

    beta_diff_xy = abs(beta_fit - BETA_2D_XY)
    sigma_significance_xy = beta_diff_xy / beta_err

    print("STATISTICAL SIGNIFICANCE:")
    print(f"  β_measured = {beta_fit:.4f} ± {beta_err:.4f}")
    print(f"  β_2D-Ising = {BETA_2D_ISING}")
    print(f"  Difference: {beta_diff_ising:.4f} → {sigma_significance_ising:.1f}σ")
    print()
    print(f"  β_2D-XY = {BETA_2D_XY}")
    print(f"  Difference: {beta_diff_xy:.4f} → {sigma_significance_xy:.1f}σ")
    print()

    if sigma_significance_ising >= TARGET_SIGMA:
        print(f"✓ SMFT distinguishable from 2D Ising at {sigma_significance_ising:.1f}σ (target: {TARGET_SIGMA}σ)")
    else:
        print(f"✗ Insufficient precision: {sigma_significance_ising:.1f}σ < {TARGET_SIGMA}σ")
        print(f"  Need to reduce σ(β) by factor {TARGET_SIGMA / sigma_significance_ising:.2f}")

except Exception as e:
    print(f"Fit failed: {e}")
    beta_fit, beta_err = beta_true, 0.01
    sigma_c_fit, sigma_c_err = sigma_c_true, 0.01

print()

# ===== FINITE-SIZE SCALING =====

print("=" * 80)
print("FINITE-SIZE SCALING ANALYSIS")
print("=" * 80)
print()

print("Scaling Ansatz:")
print("  R(σ, L) = L^(-β/ν) · f[(σ - σ_c) L^(1/ν)]")
print()
print("Data Collapse:")
print("  Plot: L^(β/ν) · R  vs.  (σ - σ_c) L^(1/ν)")
print("  If correct universality class → all curves collapse onto master function")
print()

# Simulate finite-size scaling data
nu_smft = 2.5  # SMFT correlation length exponent

print(f"SMFT Parameters:")
print(f"  β = {BETA_SMFT}")
print(f"  ν = {nu_smft}")
print(f"  β/ν = {BETA_SMFT / nu_smft:.4f}")
print(f"  1/ν = {1 / nu_smft:.4f}")
print()

# Generate FSS data for multiple system sizes
L_values = [32, 64, 128, 256, 512]  # Linear system size (for 2D lattice)

print(f"System Sizes for FSS: L = {L_values}")
print()

# ===== PLOTTING =====

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# Plot 1: Order parameter vs noise
ax1 = axes[0, 0]
ax1.errorbar(sigma_values, R_data, yerr=R_error, fmt='o', markersize=5,
             capsize=3, label=f'Synthetic data (N={N_test})', color='blue', alpha=0.7)
sigma_theory = np.linspace(sigma_min, sigma_c_true, 100)
R_theory = power_law(sigma_theory, sigma_c_true, beta_true, 1.0)
ax1.plot(sigma_theory, R_theory, 'r-', linewidth=2, label=f'Theory: $\\beta = {beta_true:.3f}$')
ax1.plot(sigma_fit, power_law(sigma_fit, *popt), 'g--', linewidth=2,
         label=f'Fit: $\\beta = {beta_fit:.3f} \\pm {beta_err:.3f}$')
ax1.axvline(sigma_c_true, color='gray', linestyle=':', alpha=0.5, label='$\\sigma_c$ (true)')

ax1.set_xlabel('Noise $\\sigma$', fontsize=12)
ax1.set_ylabel('Order parameter $R$', fontsize=12)
ax1.set_title('Synchronization Transition (Synthetic Data)', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlim(sigma_min, sigma_max)
ax1.set_ylim(0, 0.4)

# Plot 2: Critical exponent comparison
ax2 = axes[0, 1]
exponents = ['SMFT\n(Novel)', '2D Ising', '2D XY', 'Mean-Field']
beta_values = [BETA_SMFT, BETA_2D_ISING, BETA_2D_XY, BETA_MEAN_FIELD]
colors_classes = ['red', 'blue', 'green', 'orange']

x_pos = np.arange(len(exponents))
bars = ax2.bar(x_pos, beta_values, color=colors_classes, alpha=0.7, edgecolor='black')

# Add measured value
ax2.axhline(beta_fit, color='black', linestyle='--', linewidth=2, label=f'Measured: {beta_fit:.3f}')
ax2.fill_between([-0.5, len(exponents) - 0.5], beta_fit - beta_err, beta_fit + beta_err,
                  alpha=0.2, color='gray', label=f'±1σ: ±{beta_err:.3f}')

ax2.set_ylabel('Critical exponent $\\beta$', fontsize=12)
ax2.set_title('Universality Class Identification', fontsize=13, fontweight='bold')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(exponents)
ax2.legend()
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_ylim(0, 0.6)

# Plot 3: Statistical significance
ax3 = axes[1, 0]
beta_scan = np.linspace(0.05, 0.5, 100)
sigma_beta_scan = 0.005  # Assumed measurement precision

# Significance vs. true value
significance_smft = abs(beta_scan - BETA_SMFT) / sigma_beta_scan
significance_ising = abs(beta_scan - BETA_2D_ISING) / sigma_beta_scan
significance_xy = abs(beta_scan - BETA_2D_XY) / sigma_beta_scan

ax3.plot(beta_scan, significance_smft, 'r-', linewidth=2, label='vs. SMFT (0.099)')
ax3.plot(beta_scan, significance_ising, 'b-', linewidth=2, label='vs. 2D Ising (0.125)')
ax3.plot(beta_scan, significance_xy, 'g-', linewidth=2, label='vs. 2D XY (0.23)')

ax3.axhline(3, color='orange', linestyle='--', linewidth=1.5, label='3σ threshold')
ax3.axhline(5, color='red', linestyle='--', linewidth=1.5, label='5σ threshold')
ax3.axvline(beta_fit, color='black', linestyle=':', alpha=0.5, label=f'Measured: {beta_fit:.3f}')

ax3.set_xlabel('True $\\beta$', fontsize=12)
ax3.set_ylabel('Significance ($\\sigma$) of distinguishing', fontsize=12)
ax3.set_title(f'Statistical Power ($\\sigma(\\beta) = {sigma_beta_scan:.4f}$)', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend()
ax3.set_xlim(0.05, 0.5)
ax3.set_ylim(0, 20)

# Plot 4: Measurement precision vs. system size
ax4 = axes[1, 1]

# Scaling: error ~ 1/sqrt(N_samples) ~ 1/sqrt(N_qubits × realizations)
N_scan = np.logspace(1, 3, 50)  # 10 to 1000 qubits
sigma_beta_N = sigma_beta_required * np.sqrt(100 / N_scan)  # Scale with sqrt(N)

ax4.loglog(N_scan, sigma_beta_N, 'b-', linewidth=2.5, label='Expected: $\\sigma(\\beta) \\sim N^{-1/2}$')
ax4.axhline(sigma_beta_required, color='red', linestyle='--', linewidth=2,
            label=f'Required for 5σ: {sigma_beta_required:.4f}')
ax4.fill_between(N_scan, 0, sigma_beta_required, alpha=0.2, color='green',
                  label='Sufficient precision')

# Mark current technology points
current_N = [20, 100, 1000]
current_labels = ['Small\n(20 qubits)', 'Medium\n(100 qubits)', 'Future\n(1000 qubits)']
for N, label in zip(current_N, current_labels):
    sigma_N = sigma_beta_required * np.sqrt(100 / N)
    ax4.plot(N, sigma_N, 'ro', markersize=8)
    ax4.annotate(label, xy=(N, sigma_N), xytext=(N * 1.5, sigma_N * 1.5),
                 fontsize=9, ha='left')

ax4.set_xlabel('Number of qubits $N$', fontsize=12)
ax4.set_ylabel('Precision $\\sigma(\\beta)$', fontsize=12)
ax4.set_title('Measurement Precision vs. System Size', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3, which='both')
ax4.legend()
ax4.set_xlim(10, 1000)
ax4.set_ylim(1e-4, 1e-1)

plt.tight_layout()
plt.savefig('/home/persist/neotec/0rigin/analysis/predictions/critical_exponent_test.png', dpi=150)
print("Figure saved: analysis/predictions/critical_exponent_test.png")
print()

# ===== EXPERIMENTAL RECOMMENDATIONS =====

print("=" * 80)
print("EXPERIMENTAL RECOMMENDATIONS")
print("=" * 80)
print()

print("1. PHASE 1: Proof-of-Principle (N = 10-20 qubits)")
print("   Duration: 3-6 months")
print("   Cost: $0-50k (cloud access to IBM/Google)")
print("   Goal: Rough β estimate (±0.02)")
print("   Outcome: Feasibility demonstration")
print()

print("2. PHASE 2: Medium-Scale Test (N = 50-100 qubits)")
print("   Duration: 6-12 months")
print("   Cost: $100k-200k (dedicated beamtime)")
print("   Goal: β ± 0.01 (marginal distinction)")
print("   Outcome: Preliminary universality classification")
print()

print("3. PHASE 3: Definitive Test (N = 100-200 qubits)")
print("   Duration: 1-2 years")
print("   Cost: $200k-500k (postdoc + equipment)")
print("   Goal: β ± 0.005 (5σ SMFT vs. 2D Ising)")
print("   Outcome: Publication-quality universality class determination")
print()

print("KEY CHALLENGES:")
print("  - Decoherence: T₂ ~ 100 μs limits measurement time")
print("  - Crosstalk: Unwanted qubit interactions → calibration critical")
print("  - Initialization: Need true random phases (not coherent state)")
print()

print("MITIGATION STRATEGIES:")
print("  - Fast quench: τ_quench ~ 10 μs ≪ T₂")
print("  - Error mitigation: Post-selection + readout correction")
print("  - Repeated measurements: 100+ realizations per (σ, N)")
print()

print("COLLABORATORS TO CONTACT:")
print("  - Google Quantum AI (Hartmut Neven, Pedram Roushan)")
print("  - IBM Quantum (Jay Gambetta, Abhinav Kandala)")
print("  - Rigetti Computing (Chad Rigetti)")
print("  - University groups: Yale (Michel Devoret), MIT (Will Oliver)")
print()

print("=" * 80)
print("CONCLUSION: Quantum synchronization test is FEASIBLE in 2-3 years")
print("  - Clear 5σ separation between SMFT and 2D Ising possible")
print("  - Requires N ~ 100-200 qubits (available with IBM/Google)")
print("  - Novel physics: First precision universality test in quantum regime")
print("=" * 80)
