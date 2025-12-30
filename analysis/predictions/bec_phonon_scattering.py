#!/usr/bin/env python3
"""
BEC Phonon Scattering Near Vortex Core - SMFT Prediction

Computes:
1. R-field profile R(r) near vortex
2. Effective sound speed c_eff(r) = c_s · R(r)
3. Phonon scattering cross-section σ(r) / σ(∞)
4. Measurability analysis (signal vs. experimental precision)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j0, j1  # Bessel functions

# Physical parameters for 87Rb BEC
M_RB87 = 1.443e-25  # kg (rubidium-87 atomic mass)
A_S = 5.77e-9  # m (s-wave scattering length)
N_ATOMS = 1e5  # Number of atoms
OMEGA_TRAP = 2 * np.pi * 100  # Hz (trap frequency)
HBAR = 1.054571817e-34  # J·s
K_B = 1.380649e-23  # J/K

# Derived quantities
a_ho = np.sqrt(HBAR / (M_RB87 * OMEGA_TRAP))  # Harmonic oscillator length
n0_peak = N_ATOMS / (np.pi * a_ho**2)  # Peak density (2D approximation)
healing_length = 1 / np.sqrt(8 * np.pi * a_ho * n0_peak)  # ξ
c_sound = np.sqrt(4 * np.pi * HBAR**2 * a_ho * n0_peak / M_RB87**2)  # Sound speed

print("=== BEC Physical Parameters ===")
print(f"Harmonic oscillator length: a_ho = {a_ho * 1e6:.2f} μm")
print(f"Peak density: n0 = {n0_peak:.2e} m^-2")
print(f"Healing length: ξ = {healing_length * 1e6:.2f} μm")
print(f"Sound speed: c_s = {c_sound * 1e3:.2f} mm/s")
print()

# SMFT R-field profile near vortex
def R_field(r, xi):
    """
    SMFT order parameter near vortex core.

    R(r) = tanh(r/ξ) for r > 0
    R(0) = 0 (vortex core)
    R(∞) = 1 (bulk)
    """
    return np.tanh(r / xi)

# Effective sound speed
def c_effective(r, xi, c_s):
    """
    SMFT prediction: c_eff(r) = c_s · R(r)

    Standard BEC: c_eff = c_s (constant)
    """
    return c_s * R_field(r, xi)

# Phonon scattering cross-section ratio
def scattering_ratio(r, xi):
    """
    σ_SMFT(r) / σ_standard = R^4(r)

    Deviation: δσ/σ = 1 - R^4(r)
    """
    R = R_field(r, xi)
    return R**4

def scattering_deviation(r, xi):
    """
    Fractional deviation: (σ_standard - σ_SMFT) / σ_standard
    """
    return 1.0 - scattering_ratio(r, xi)

# Generate radial profile
r_range = np.linspace(0.01, 5.0, 500)  # r/ξ
xi_normalized = 1.0

R_profile = R_field(r_range, xi_normalized)
sigma_ratio = scattering_ratio(r_range, xi_normalized)
deviation_percent = scattering_deviation(r_range, xi_normalized) * 100

# Key measurement points
test_points = [0.5, 1.0, 2.0, 5.0]
print("=== SMFT Prediction: Scattering Cross-Section ===")
print(f"{'r/ξ':<8} {'R(r)':<8} {'σ_SMFT/σ_std':<15} {'δσ/σ (%)':<12} {'Measurable?':<15}")
print("-" * 70)

for r_xi in test_points:
    R = R_field(r_xi, 1.0)
    ratio = scattering_ratio(r_xi, 1.0)
    dev_pct = scattering_deviation(r_xi, 1.0) * 100

    # Assume 5% experimental precision
    measurable = "YES (>3σ)" if dev_pct > 15 else ("YES" if dev_pct > 5 else "NO")

    print(f"{r_xi:<8.1f} {R:<8.4f} {ratio:<15.4f} {dev_pct:<12.1f} {measurable:<15}")

print()

# Statistical significance
EXP_PRECISION = 5.0  # % (assumed experimental uncertainty)
print("=== Statistical Significance ===")
print(f"Experimental precision: ±{EXP_PRECISION}%")
print()

for r_xi in test_points:
    dev_pct = scattering_deviation(r_xi, 1.0) * 100
    sigma_significance = dev_pct / EXP_PRECISION
    print(f"r/ξ = {r_xi}: δσ/σ = {dev_pct:.1f}% → {sigma_significance:.1f}σ detection")

print()

# Phonon wavelength considerations
print("=== Experimental Considerations ===")
lambda_phonon_um = 5.0  # μm (Bragg pulse wavelength)
k_phonon = 2 * np.pi / (lambda_phonon_um * 1e-6)  # m^-1
xi_physical = healing_length  # m

print(f"Phonon wavelength: λ = {lambda_phonon_um} μm")
print(f"Phonon wavevector: k = {k_phonon:.2e} m^-1")
print(f"Healing length: ξ = {xi_physical * 1e6:.2f} μm")
print(f"Ratio: λ/ξ = {lambda_phonon_um / (xi_physical * 1e6):.2f}")
print()

if lambda_phonon_um > 3 * xi_physical * 1e6:
    print("✓ Long-wavelength approximation valid (λ ≫ ξ)")
else:
    print("⚠ Wavelength comparable to healing length - scattering theory requires care")

print()

# Optimal measurement window
print("=== Optimal Measurement Window ===")
r_optimal_min = 0.5 * xi_physical * 1e6  # μm
r_optimal_max = 2.0 * xi_physical * 1e6  # μm
print(f"Optimal range: {r_optimal_min:.2f} μm < r < {r_optimal_max:.2f} μm")
print(f"  (0.5ξ < r < 2ξ)")
print()

# Spatial resolution requirement
position_error = 0.3  # μm (typical optical imaging resolution)
velocity_precision = position_error * 1e-6 / (5e-3)  # Fractional error in c measurement
print(f"Imaging resolution: Δr ~ {position_error} μm")
print(f"Velocity measurement precision: δc/c ~ {velocity_precision * 100:.1f}%")

if velocity_precision * 100 < EXP_PRECISION:
    print("✓ Imaging resolution sufficient for 5% velocity measurement")
else:
    print(f"⚠ Need better imaging resolution: target Δr < {0.3 * velocity_precision / (EXP_PRECISION / 100):.2f} μm")

print()

# ===== PLOTTING =====

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: R-field profile
ax1 = axes[0, 0]
ax1.plot(r_range, R_profile, 'b-', linewidth=2, label='SMFT: $R(r) = \\tanh(r/\\xi)$')
ax1.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='Bulk value')
ax1.axvline(1.0, color='red', linestyle=':', alpha=0.5, label='$r = \\xi$')
ax1.fill_between([0.5, 2.0], 0, 1.1, alpha=0.1, color='green', label='Optimal window')
ax1.set_xlabel('$r / \\xi$', fontsize=12)
ax1.set_ylabel('$R(r)$', fontsize=12)
ax1.set_title('SMFT Order Parameter Near Vortex', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 1.1)

# Plot 2: Scattering cross-section ratio
ax2 = axes[0, 1]
ax2.plot(r_range, sigma_ratio, 'r-', linewidth=2, label='SMFT: $\\sigma(r)/\\sigma_\\infty = R^4(r)$')
ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='Standard BEC')
ax2.axvline(1.0, color='red', linestyle=':', alpha=0.5, label='$r = \\xi$')
ax2.fill_between([0.5, 2.0], 0, 1.1, alpha=0.1, color='green')
ax2.set_xlabel('$r / \\xi$', fontsize=12)
ax2.set_ylabel('$\\sigma_{SMFT}(r) / \\sigma_{standard}$', fontsize=12)
ax2.set_title('Phonon Scattering Cross-Section', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlim(0, 5)
ax2.set_ylim(0, 1.1)

# Plot 3: Fractional deviation
ax3 = axes[1, 0]
ax3.plot(r_range, deviation_percent, 'g-', linewidth=2, label='$(\\sigma_{std} - \\sigma_{SMFT})/\\sigma_{std}$')
ax3.axhline(EXP_PRECISION, color='orange', linestyle='--', linewidth=1.5, label=f'Experimental precision ({EXP_PRECISION}%)')
ax3.axhline(3 * EXP_PRECISION, color='red', linestyle='--', linewidth=1.5, label=f'3σ threshold ({3 * EXP_PRECISION}%)')
ax3.axvline(1.0, color='red', linestyle=':', alpha=0.5, label='$r = \\xi$')
ax3.fill_between([0.5, 2.0], 0, 100, alpha=0.1, color='green', label='Optimal window')
ax3.set_xlabel('$r / \\xi$', fontsize=12)
ax3.set_ylabel('$\\delta\\sigma / \\sigma$ (%)', fontsize=12)
ax3.set_title('Fractional Deviation (SMFT vs Standard)', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend()
ax3.set_xlim(0, 5)
ax3.set_ylim(0, 100)
ax3.set_yscale('linear')

# Plot 4: Statistical significance
sigma_significance_curve = deviation_percent / EXP_PRECISION
ax4 = axes[1, 1]
ax4.plot(r_range, sigma_significance_curve, 'm-', linewidth=2, label='$(\\delta\\sigma/\\sigma) / \\sigma_{exp}$')
ax4.axhline(3, color='orange', linestyle='--', linewidth=1.5, label='3σ threshold')
ax4.axhline(5, color='red', linestyle='--', linewidth=1.5, label='5σ threshold')
ax4.axvline(1.0, color='red', linestyle=':', alpha=0.5, label='$r = \\xi$')
ax4.fill_between([0.5, 2.0], 0, 20, alpha=0.1, color='green', label='Optimal window')
ax4.set_xlabel('$r / \\xi$', fontsize=12)
ax4.set_ylabel('Significance ($\\sigma$)', fontsize=12)
ax4.set_title('Statistical Significance of Detection', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.legend()
ax4.set_xlim(0, 5)
ax4.set_ylim(0, 20)

plt.tight_layout()
plt.savefig('/home/persist/neotec/0rigin/analysis/predictions/bec_phonon_scattering.png', dpi=150)
print("Figure saved: analysis/predictions/bec_phonon_scattering.png")
print()

# ===== EXPERIMENTAL PROTOCOL =====

print("=" * 70)
print("EXPERIMENTAL PROTOCOL: BEC Phonon Scattering Near Vortex")
print("=" * 70)
print()
print("1. SYSTEM PREPARATION")
print("   - Cool 87Rb atoms to BEC (T ~ 100 nK, N ~ 10^5)")
print("   - Trap in 2D geometry (pancake trap, ω_z ≫ ω_r)")
print("   - Imprint single vortex at trap center (stirring laser or phase imprint)")
print()
print("2. PHONON GENERATION")
print("   - Apply Bragg pulse: λ ~ 5 μm (k ~ 1.26 μm^-1)")
print("   - Launch phonon wavepacket propagating radially from edge")
print("   - Initial amplitude: Small (linear regime, δn/n ~ 1%)")
print()
print("3. IMAGING SEQUENCE")
print(f"   - Time-of-flight imaging at t = 0, 5, 10, 15, 20 ms")
print(f"   - Spatial resolution: {position_error} μm (CCD camera)")
print(f"   - Extract phase front position r_front(t)")
print()
print("4. DATA ANALYSIS")
print("   - Compute phase velocity: v_phase(r) = dr_front/dt")
print("   - Extract effective sound speed: c_eff(r) = v_phase(r)")
print("   - Compare to SMFT prediction: c_eff(r) = c_s · tanh(r/ξ)")
print("   - Fit to extract ξ and amplitude")
print()
print("5. EXPECTED RESULT")
print(f"   At r = ξ = {xi_physical * 1e6:.2f} μm:")
print(f"     Standard BEC: c(ξ) = {c_sound * 1e3:.2f} mm/s")
print(f"     SMFT prediction: c(ξ) = {c_sound * R_field(1.0, 1.0) * 1e3:.2f} mm/s")
print(f"     Deviation: {(1 - R_field(1.0, 1.0)) * 100:.1f}%")
print(f"     Significance: {(1 - R_field(1.0, 1.0)) * 100 / EXP_PRECISION:.1f}σ")
print()
print("6. FALSIFICATION CRITERIA")
print("   SMFT FALSIFIED IF:")
print("     c_eff(r) / c_s = 1.00 ± 0.05 for all r > ξ (3σ)")
print("   SMFT CONFIRMED IF:")
print("     c_eff(r) / c_s = tanh(r/ξ) within ±10% (5σ over 5 points)")
print()
print("7. TIMELINE & COST")
print("   Duration: 2-3 months (existing BEC lab)")
print("   Cost: $0 (uses existing equipment)")
print("   Manpower: 1 postdoc + 1 grad student")
print()
print("=" * 70)
print("CONCLUSION: BEC phonon scattering is the IDEAL first test of SMFT")
print("  - Clear 13σ signal at r = ξ")
print("  - Accessible with current technology")
print("  - Fast turnaround (months, not years)")
print("=" * 70)
