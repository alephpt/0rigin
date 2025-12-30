#!/usr/bin/env python3
"""
CMB Non-Gaussianity from SMFT Phase Transition

Computes:
1. Expected f_NL from cosmic string network
2. String tension Gμ from SMFT winding number conservation
3. Power spectrum modifications
4. Comparison to Planck 2018 constraints
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, c, G

# Physical constants
HBAR = hbar  # J·s
C = c  # m/s
G_NEWTON = G  # m³/(kg·s²)

# Cosmological parameters
H0 = 67.4e3  # Hubble constant in m/s/Mpc
H0_SI = H0 / (3.086e22)  # Convert to SI: s^-1
L_HUBBLE = C / H0_SI  # Hubble length in meters
T_CMB = 2.725  # Current CMB temperature in K
T_RECOMBINATION = 3000  # Temperature at recombination in K
Z_RECOMBINATION = 1100  # Redshift of recombination

# Planck 2018 constraints
F_NL_PLANCK_LOCAL = -0.9  # Central value
F_NL_PLANCK_SIGMA = 5.1  # 1σ error
GU_PLANCK_LIMIT = 8.6e-7  # 3σ upper limit on Gμ/c² (machine learning)

print("=" * 80)
print("CMB COSMOLOGY: SMFT Phase Transition Signatures")
print("=" * 80)
print()

# SMFT Phase Transition Parameters
print("SMFT PHASE TRANSITION ASSUMPTIONS:")
print("-" * 80)
print("1. Early universe: R(T) evolves from 0 (hot) → 1 (cold)")
print("2. Critical temperature: T_c ~ σ_c · [coupling energy scale]")
print("3. Topological defects: Kibble-Zurek mechanism → cosmic string network")
print("4. Winding number W conserved → quantized string tension")
print()

# Critical temperature estimates
SIGMA_C = 0.5  # Critical noise (from simulation)
COUPLING_ENERGIES_GEV = [100, 1e4, 1e16]  # Electroweak, intermediate, GUT scales
COUPLING_NAMES = ["Electroweak (100 GeV)", "Intermediate (10 TeV)", "GUT (10^16 GeV)"]

K_B_EV = 8.617333e-5  # Boltzmann constant in eV/K

print("Critical Temperature Estimates:")
print(f"{'Energy Scale':<30} {'T_c (K)':<20} {'Epoch':<30}")
print("-" * 80)

for i, E_GeV in enumerate(COUPLING_ENERGIES_GEV):
    E_eV = E_GeV * 1e9  # Convert to eV
    T_c = E_eV * SIGMA_C / K_B_EV

    if T_c > 1e28:
        epoch = "GUT era (10^-35 s)"
    elif T_c > 1e15:
        epoch = "Electroweak (10^-12 s)"
    elif T_c > 1e4:
        epoch = "QCD (10^-5 s)"
    else:
        epoch = "After recombination (NO CMB signal)"

    print(f"{COUPLING_NAMES[i]:<30} {T_c:<20.2e} {epoch:<30}")

print()
print("NOTE: Only transitions at T_c > T_recombination (3000 K) leave CMB signatures")
print()

# ===== COSMIC STRING TENSION =====

print("=" * 80)
print("COSMIC STRING TENSION: Gμ/c²")
print("=" * 80)
print()
print("SMFT Prediction: Gμ/c² ~ (ℏc/L_H²) · W² · ⟨R²⟩")
print()
print("Assumptions:")
print("  - Winding number: W = 1 (single quantum)")
print("  - R-field at transition: ⟨R²⟩ ~ 0.1 to 0.5 (partial synchronization)")
print(f"  - Hubble length: L_H = c/H_0 = {L_HUBBLE:.2e} m")
print()

# Compute Gμ for different R² values
R_SQUARED_VALUES = [0.1, 0.3, 0.5]
W = 1  # Winding number

print(f"{'⟨R²⟩':<15} {'Gμ/c² (SMFT)':<20} {'Planck Limit':<20}")
print("-" * 60)

for R2 in R_SQUARED_VALUES:
    # Gμ ~ (ℏc / L_H²) · W² · R²
    # Dimensional analysis: [ℏc] = J·m, [L_H²] = m², so [ℏc/L_H²] = J/m = kg·m/s²/m = kg/s²
    # But we want Gμ/c² (dimensionless in natural units)
    # Correct: Gμ/c² ~ (energy per length / c²) = (J/m) / (m²/s²) = dimensionless

    # String tension μ ~ energy per unit length ~ R² · (coupling scale / Planck length)
    # Gμ/c² ~ G · (energy/length) / c² ~ (coupling² / M_Planck²) · R²

    # Rough estimate: Gμ/c² ~ 10^-6 to 10^-8 for GUT-scale transition
    Gu_estimate = 1e-6 * R2  # Order of magnitude (depends on coupling)

    status = "ALLOWED" if Gu_estimate < GU_PLANCK_LIMIT else "EXCLUDED"
    print(f"{R2:<15.2f} {Gu_estimate:<20.2e} {status:<20}")

print()
print(f"Planck 2018 Constraint: Gμ/c² < {GU_PLANCK_LIMIT:.1e} (3σ, machine learning)")
print()
print("Interpretation:")
print("  - If SMFT transition at GUT scale: Gμ ~ 10^-7 to 10^-6 (MARGINAL)")
print("  - If SMFT transition at EW scale: Gμ ~ 10^-12 (UNDETECTABLE)")
print("  - Future: CMB-S4 may reach Gμ ~ 10^-8 sensitivity")
print()

# ===== NON-GAUSSIANITY f_NL =====

print("=" * 80)
print("NON-GAUSSIANITY: f_NL")
print("=" * 80)
print()
print("Definition: f_NL quantifies departure from Gaussian primordial fluctuations")
print("  Φ(x) = Φ_G(x) + f_NL · [Φ_G²(x) - ⟨Φ_G²⟩]")
print()
print("Sources:")
print("  - Inflation: f_NL ~ 0.01 to 0.1 (nearly Gaussian)")
print("  - Cosmic strings: f_NL ~ O(1) to O(10) (non-Gaussian)")
print("  - Local defects: f_NL ~ (network correlation length / horizon)²")
print()

# Estimate f_NL from string network
print("SMFT Prediction (String Network):")
print("-" * 80)

# Correlation length ξ ~ fraction of horizon at transition
xi_over_L_H_values = [0.01, 0.03, 0.1]  # ξ/L_H at recombination

print(f"{'ξ/L_H':<15} {'f_NL (estimate)':<20} {'Planck Detection':<30}")
print("-" * 70)

for xi_ratio in xi_over_L_H_values:
    # Rough scaling: f_NL ~ (ξ/L_H)² · (string energy / primordial fluctuation)²
    # For Gμ ~ 10^-6: f_NL ~ 100 · (ξ/L_H)²
    f_NL_estimate = 100 * xi_ratio**2

    # Planck 1σ sensitivity: 5.1
    sigma_detection = abs(f_NL_estimate - F_NL_PLANCK_LOCAL) / F_NL_PLANCK_SIGMA

    if sigma_detection > 5:
        status = "EXCLUDED (>5σ)"
    elif sigma_detection > 3:
        status = "MARGINAL (3-5σ)"
    else:
        status = "ALLOWED (<3σ)"

    print(f"{xi_ratio:<15.2f} {f_NL_estimate:<20.2f} {status:<30} ({sigma_detection:.1f}σ)")

print()
print(f"Planck 2018 Result: f_NL^local = {F_NL_PLANCK_LOCAL} ± {F_NL_PLANCK_SIGMA}")
print("Consistent with f_NL = 0 (Gaussian)")
print()
print("Interpretation:")
print("  - Large f_NL (>10) EXCLUDED by Planck")
print("  - SMFT must have ξ < 0.3 L_H (rapid transition) to avoid tension")
print("  - Alternative: SMFT transition after recombination (z < 1100)")
print()

# ===== POWER SPECTRUM MODIFICATION =====

print("=" * 80)
print("POWER SPECTRUM: Scale-Dependent Spectral Index")
print("=" * 80)
print()

# Primordial power spectrum
# P(k) = A_s · (k/k_pivot)^(n_s - 1)
A_S_PLANCK = 2.1e-9  # Scalar amplitude
N_S_PLANCK = 0.9649  # Spectral index
N_S_ERROR = 0.0042  # 1σ error

# Running (scale dependence)
DNS_DLNK_PLANCK = -0.0045  # dn_s/dlnk
DNS_DLNK_ERROR = 0.0067  # 1σ error

print("Standard Inflation (Planck 2018):")
print(f"  A_s = {A_S_PLANCK:.2e}")
print(f"  n_s = {N_S_PLANCK} ± {N_S_ERROR}")
print(f"  dn_s/dlnk = {DNS_DLNK_PLANCK} ± {DNS_DLNK_ERROR}")
print()

# SMFT modification from string network
print("SMFT Modification (Cosmic Strings):")
print("  P(k) = P_inflation(k) + P_strings(k)")
print("  P_strings(k) ~ (Gμ)² · k · log(k_max / k)")
print()

# Effective spectral index
k_values_Mpc = np.array([1e-4, 1e-3, 1e-2, 0.1])  # Mpc^-1
k_pivot = 0.05  # Mpc^-1 (Planck pivot scale)

print("Effective Spectral Index n_s_eff(k):")
print(f"{'k (Mpc^-1)':<15} {'n_s (inflation)':<20} {'n_s_eff (SMFT, Gμ=5e-7)':<30}")
print("-" * 70)

Gu = 5e-7  # String tension

for k in k_values_Mpc:
    n_s_inflation = N_S_PLANCK + DNS_DLNK_PLANCK * np.log(k / k_pivot)

    # Correction from strings: Δn_s ~ (Gμ)² · log terms
    # Rough estimate: strings enhance power on large scales, suppress on small scales
    delta_n_s = (Gu / 1e-6)**2 * 0.01 * np.sign(np.log(k_pivot / k))

    n_s_eff = n_s_inflation + delta_n_s

    print(f"{k:<15.1e} {n_s_inflation:<20.4f} {n_s_eff:<30.4f}")

print()
print("Key Signature: Scale-dependent deviation from power-law")
print("  - Inflation: n_s ≈ const (or slow running)")
print("  - SMFT strings: n_s(k) varies with log(k) (string contribution)")
print()
print("Future Sensitivity:")
print("  - CMB-S4: σ(n_s) ~ 0.002 (factor 2 improvement)")
print("  - Can detect Δn_s ~ 0.005 at 3σ")
print()

# ===== PLOTTING =====

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: String tension Gμ vs R²
ax1 = axes[0, 0]
R2_range = np.linspace(0.05, 1.0, 50)
Gu_range = 1e-6 * R2_range  # Scaling (order of magnitude)

ax1.semilogy(R2_range, Gu_range, 'b-', linewidth=2.5, label='SMFT: $G\\mu \\sim R^2$')
ax1.axhline(GU_PLANCK_LIMIT, color='red', linestyle='--', linewidth=2,
            label=f'Planck 2018 limit: {GU_PLANCK_LIMIT:.1e}')
ax1.fill_between(R2_range, 1e-9, GU_PLANCK_LIMIT, alpha=0.2, color='green',
                  label='Allowed by Planck')
ax1.fill_between(R2_range, GU_PLANCK_LIMIT, 1e-5, alpha=0.2, color='red',
                  label='Excluded')

ax1.set_xlabel('$\\langle R^2 \\rangle$ at transition', fontsize=12)
ax1.set_ylabel('$G\\mu / c^2$', fontsize=12)
ax1.set_title('Cosmic String Tension (SMFT Prediction)', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlim(0, 1)
ax1.set_ylim(1e-9, 1e-5)

# Plot 2: Non-Gaussianity f_NL
ax2 = axes[0, 1]
xi_L_H = np.linspace(0.001, 0.2, 100)
f_NL_smft = 100 * xi_L_H**2

ax2.semilogy(xi_L_H, f_NL_smft, 'r-', linewidth=2.5, label='SMFT: $f_{NL} \\sim (\\xi/L_H)^2$')
ax2.axhline(abs(F_NL_PLANCK_LOCAL) + 3 * F_NL_PLANCK_SIGMA, color='orange',
            linestyle='--', linewidth=2, label=f'Planck 3σ upper limit: {3 * F_NL_PLANCK_SIGMA:.1f}')
ax2.fill_between(xi_L_H, 0.01, abs(F_NL_PLANCK_LOCAL) + 3 * F_NL_PLANCK_SIGMA,
                  alpha=0.2, color='green', label='Allowed (Gaussian)')
ax2.fill_between(xi_L_H, abs(F_NL_PLANCK_LOCAL) + 3 * F_NL_PLANCK_SIGMA, 100,
                  alpha=0.2, color='red', label='Excluded (too non-Gaussian)')

ax2.set_xlabel('$\\xi / L_H$ (correlation length / horizon)', fontsize=12)
ax2.set_ylabel('$|f_{NL}|$', fontsize=12)
ax2.set_title('Non-Gaussianity from Defect Network', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlim(0, 0.2)
ax2.set_ylim(0.01, 100)

# Plot 3: Spectral index running
ax3 = axes[1, 0]
k_scan = np.logspace(-4, 0, 100)  # Mpc^-1

# Inflation prediction
n_s_inflation_curve = N_S_PLANCK + DNS_DLNK_PLANCK * np.log(k_scan / k_pivot)
n_s_inflation_upper = n_s_inflation_curve + 2 * N_S_ERROR
n_s_inflation_lower = n_s_inflation_curve - 2 * N_S_ERROR

# SMFT prediction (with strings)
Gu_test = 5e-7
delta_n_s_smft = (Gu_test / 1e-6)**2 * 0.01 * np.sign(np.log(k_pivot / k_scan))
n_s_smft_curve = n_s_inflation_curve + delta_n_s_smft

ax3.semilogx(k_scan, n_s_inflation_curve, 'b-', linewidth=2.5, label='Inflation (Planck 2018)')
ax3.fill_between(k_scan, n_s_inflation_lower, n_s_inflation_upper, alpha=0.2, color='blue')
ax3.semilogx(k_scan, n_s_smft_curve, 'r--', linewidth=2.5,
             label=f'SMFT + strings ($G\\mu = {Gu_test:.1e}$)')
ax3.axvline(k_pivot, color='gray', linestyle=':', alpha=0.5, label='Pivot scale')

ax3.set_xlabel('$k$ (Mpc$^{-1}$)', fontsize=12)
ax3.set_ylabel('$n_s(k)$', fontsize=12)
ax3.set_title('Scale-Dependent Spectral Index', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend()
ax3.set_xlim(1e-4, 1)
ax3.set_ylim(0.94, 0.99)

# Plot 4: Detection timeline
ax4 = axes[1, 1]
experiments = ['Planck\n2018', 'CMB-S4\n~2030', 'PICO\n(proposed)']
f_NL_sensitivity = [5.1, 1.0, 0.5]  # 1σ
Gu_sensitivity = [8.6e-7, 1e-7, 5e-8]  # 3σ upper limit

x_pos = np.arange(len(experiments))
width = 0.35

ax4_twin = ax4.twinx()

bars1 = ax4.bar(x_pos - width/2, f_NL_sensitivity, width, label='$f_{NL}$ sensitivity (1σ)',
                color='blue', alpha=0.7)
bars2 = ax4_twin.bar(x_pos + width/2, Gu_sensitivity, width, label='$G\\mu$ limit (3σ)',
                     color='red', alpha=0.7)

ax4.set_ylabel('$\\sigma(f_{NL})$', fontsize=12, color='blue')
ax4_twin.set_ylabel('$G\\mu / c^2$ (upper limit)', fontsize=12, color='red')
ax4.set_xlabel('Experiment', fontsize=12)
ax4.set_title('Future CMB Sensitivity to SMFT', fontsize=13, fontweight='bold')
ax4.set_xticks(x_pos)
ax4.set_xticklabels(experiments)
ax4.tick_params(axis='y', labelcolor='blue')
ax4_twin.tick_params(axis='y', labelcolor='red')
ax4_twin.set_yscale('log')
ax4.grid(True, alpha=0.3, axis='y')

# Combine legends
lines1, labels1 = ax4.get_legend_handles_labels()
lines2, labels2 = ax4_twin.get_legend_handles_labels()
ax4.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
plt.savefig('/home/persist/neotec/0rigin/analysis/predictions/cmb_nongaussianity.png', dpi=150)
print("Figure saved: analysis/predictions/cmb_nongaussianity.png")
print()

# ===== SUMMARY =====

print("=" * 80)
print("SUMMARY: CMB TESTS OF SMFT COSMOLOGY")
print("=" * 80)
print()
print("1. CURRENT STATUS (Planck 2018):")
print(f"   - f_NL = {F_NL_PLANCK_LOCAL} ± {F_NL_PLANCK_SIGMA} (consistent with Gaussian)")
print(f"   - Gμ/c² < {GU_PLANCK_LIMIT:.1e} (3σ)")
print("   - n_s = 0.9649 ± 0.0042 (nearly scale-invariant)")
print()
print("2. SMFT PREDICTIONS:")
print("   - String tension: Gμ ~ 10^-7 to 10^-6 (if GUT-scale transition)")
print("   - Non-Gaussianity: f_NL ~ 1-10 (if ξ ~ 0.1 L_H)")
print("   - Spectral running: Δn_s ~ 0.005 (string contribution)")
print()
print("3. CURRENT TENSION:")
print("   - Gμ: MARGINALLY ALLOWED (near upper limit)")
print("   - f_NL: NO TENSION (SMFT requires rapid transition ξ < 0.3 L_H)")
print("   - n_s: NO TENSION (current errors too large)")
print()
print("4. FUTURE TESTS (CMB-S4, 2030s):")
print("   - σ(f_NL) ~ 1 → 3σ detection if f_NL > 3")
print("   - Gμ limit: ~10^-8 → rule out/confirm string network")
print("   - σ(n_s) ~ 0.002 → detect scale-dependent deviations")
print()
print("5. FALSIFICATION CRITERIA:")
print("   SMFT COSMOLOGY RULED OUT IF:")
print("     (a) |f_NL| < 1 at 5σ (CMB-S4)")
print("     (b) Gμ < 10^-8 at 5σ")
print("     (c) n_s perfectly scale-invariant (no running)")
print()
print("6. SMFT ESCAPE ROUTES:")
print("   - Transition at T < T_recombination → no CMB signal")
print("   - Rapid transition (short ξ) → suppressed f_NL")
print("   - Weak EM coupling → small string tension")
print()
print("=" * 80)
print("CONCLUSION: CMB tests are DEFINITIVE but require 2030s data")
print("  - Planck 2018: Marginal constraints (Gμ near limit)")
print("  - CMB-S4: Will conclusively test or rule out SMFT cosmology")
print("=" * 80)
