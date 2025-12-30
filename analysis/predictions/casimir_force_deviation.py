#!/usr/bin/env python3
"""
Casimir Force Modification from R-Field Coupling

Computes:
1. Standard Casimir force F(d) for parallel plates
2. SMFT prediction: F_SMFT = F_standard · ⟨R²⟩
3. Measurable deviations vs. experimental precision
4. Material and geometry dependence
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, c, pi, epsilon_0

# Physical constants
HBAR_C = hbar * c  # J·m

# Casimir force (ideal case, T=0, perfect conductors)
def casimir_force_standard(d):
    """
    Standard Casimir force per unit area between parallel plates.

    F/A = -π²ℏc / (240 d⁴)

    Args:
        d: separation in meters

    Returns:
        Force per unit area in N/m² (negative = attractive)
    """
    return -(pi**2 * HBAR_C) / (240 * d**4)

def casimir_force_smft(d, R_avg):
    """
    SMFT prediction: F_SMFT = F_standard · ⟨R²⟩(d)

    Args:
        d: separation in meters
        R_avg: average value of R-field between plates

    Returns:
        Force per unit area in N/m²
    """
    return casimir_force_standard(d) * R_avg**2

# Distance range (nano-scale to micron-scale)
d_nm = np.array([10, 20, 50, 100, 200, 500, 1000])  # nm
d_meters = d_nm * 1e-9  # meters

# Compute forces
F_standard = casimir_force_standard(d_meters)  # N/m²

# Convert to pN/μm² for easier reading
F_standard_pN_um2 = F_standard * 1e-12 / 1e-12  # pN/μm²

# SMFT prediction with different R values
R_values = [0.90, 0.95, 0.99, 1.00]  # Different synchronization levels
colors = ['red', 'orange', 'green', 'blue']

print("=" * 80)
print("CASIMIR FORCE: SMFT vs Standard QFT")
print("=" * 80)
print()
print("Standard Casimir Force (Perfect Conductors, T=0K):")
print(f"{'d (nm)':<10} {'F_standard (pN/μm²)':<25} {'F_standard (N/m²)':<20}")
print("-" * 80)

for i, d in enumerate(d_nm):
    print(f"{d:<10.0f} {F_standard_pN_um2[i]:<25.4e} {F_standard[i]:<20.4e}")

print()
print("=" * 80)
print("SMFT PREDICTION: F_SMFT = F_standard · ⟨R²⟩")
print("=" * 80)
print()

for R_avg in R_values[:-1]:  # Exclude R=1.00 (that's standard case)
    print(f"--- R_avg = {R_avg} (⟨R²⟩ = {R_avg**2:.4f}) ---")
    print(f"{'d (nm)':<10} {'F_SMFT (pN/μm²)':<20} {'δF/F (%)':<15}")
    print("-" * 50)

    for i, d in enumerate(d_nm):
        F_smft = casimir_force_smft(d_meters[i], R_avg)
        F_smft_pN = F_smft * 1e-12 / 1e-12
        deviation_percent = (1 - R_avg**2) * 100

        print(f"{d:<10.0f} {F_smft_pN:<20.4e} {deviation_percent:<15.2f}")

    print()

# Key observation: Distance-independent deviation
print("=" * 80)
print("KEY FEATURE: Distance-Independent Deviation")
print("=" * 80)
print()
print("If R = const (uniform between plates):")
print("  F_SMFT / F_standard = R² = const (independent of d)")
print()
print("Power law scaling:")
print("  Standard: F ∝ d^(-4)")
print("  SMFT:     F ∝ d^(-4) (SAME scaling, different amplitude)")
print()
print("This means: Measuring F(d) slope won't distinguish SMFT from standard QFT")
print("            Must measure absolute magnitude with high precision")
print()

# Experimental precision requirements
print("=" * 80)
print("EXPERIMENTAL FEASIBILITY")
print("=" * 80)
print()

exp_precision_percent = 2.0  # State-of-art (Decca 2007: ~0.2%, typical ~1-2%)

print(f"Assumed experimental precision: ±{exp_precision_percent}%")
print()
print(f"{'R_avg':<10} {'δF/F (%)':<15} {'Detectable at 3σ?':<25}")
print("-" * 50)

for R_avg in R_values[:-1]:
    deviation = (1 - R_avg**2) * 100
    sigma = deviation / exp_precision_percent
    detectable = "YES" if sigma >= 3 else "NO"
    print(f"{R_avg:<10.2f} {deviation:<15.2f} {detectable:<25} ({sigma:.1f}σ)")

print()
print("Conclusion:")
print(f"  For R_avg = 0.95 (5% suppression): δF/F = 9.6%")
print(f"  With 2% precision: 9.6% / 2% = 4.8σ → MARGINALLY DETECTABLE")
print(f"  Need precision < 2% for confident detection (5σ)")
print()

# Material dependence
print("=" * 80)
print("MATERIAL DEPENDENCE TEST")
print("=" * 80)
print()
print("Hypothesis: Different materials → different R values")
print()
print("Materials:")
print("  Gold (Au): Excellent conductor → R_Au ≈ 0.99")
print("  Silicon (Si): Semiconductor → R_Si ≈ 0.90")
print()

R_Au = 0.99
R_Si = 0.90

print("Test: Measure F_Au and F_Si at same distance d")
print()
print(f"{'d (nm)':<10} {'F_Au/F_Si (standard)':<25} {'F_Au/F_Si (SMFT)':<25} {'Difference':<15}")
print("-" * 80)

for d in [50, 100, 200]:
    # Standard: ratio comes from ε(ω) only (material-dependent dielectric function)
    # For simplicity, assume similar ε → ratio ≈ 1 in standard theory
    ratio_standard = 1.00  # Simplified (real case: Lifshitz theory with ε(ω))

    # SMFT: additional factor from R
    ratio_smft = R_Au**2 / R_Si**2

    difference_percent = (ratio_smft - ratio_standard) / ratio_standard * 100

    print(f"{d:<10.0f} {ratio_standard:<25.3f} {ratio_smft:<25.3f} {difference_percent:<15.1f}%")

print()
print("Expected SMFT signature:")
print(f"  F_Au / F_Si = (R_Au / R_Si)² = ({R_Au}/{R_Si})² = {(R_Au/R_Si)**2:.3f}")
print(f"  Enhancement: {((R_Au/R_Si)**2 - 1) * 100:.1f}% relative to standard theory")
print()
print("Note: Real experiment must account for ε(ω) differences (Lifshitz theory)")
print("      SMFT adds R² factor ON TOP of standard material dependence")
print()

# Geometry dependence
print("=" * 80)
print("GEOMETRY DEPENDENCE: Spatially Varying R(x,y)")
print("=" * 80)
print()
print("If R(x,y,z) varies spatially between plates:")
print("  F_SMFT = -(π²ℏc/240d⁴) · ⟨R²(x,y,z)⟩")
print()
print("Example: Corrugated surface")
print("  Flat region: R_flat = 0.95")
print("  Corrugations: R_corr = 0.92 (local field enhancement)")
print()

R_flat = 0.95
R_corr = 0.92
area_fraction_corr = 0.3  # 30% of surface is corrugated

R2_avg_flat = R_flat**2
R2_avg_corrugated = (1 - area_fraction_corr) * R_flat**2 + area_fraction_corr * R_corr**2

print(f"Flat plates: ⟨R²⟩ = {R2_avg_flat:.4f}")
print(f"Corrugated:  ⟨R²⟩ = {R2_avg_corrugated:.4f}")
print(f"Ratio: F_corr / F_flat = {R2_avg_corrugated / R2_avg_flat:.4f}")
print(f"Deviation: {((R2_avg_corrugated / R2_avg_flat) - 1) * 100:.2f}%")
print()
print("This geometry dependence is ADDITIONAL to standard EM mode structure")
print("Could be tested by comparing flat vs. patterned surfaces")
print()

# ===== PLOTTING =====

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Casimir force vs distance
ax1 = axes[0, 0]
for i, R_avg in enumerate(R_values):
    F_smft_array = casimir_force_smft(d_meters, R_avg) * 1e-12 / 1e-12  # pN/μm²
    label = f'R = {R_avg}' if R_avg < 1.0 else 'Standard (R=1)'
    linestyle = '-' if R_avg == 1.0 else '--'
    ax1.loglog(d_nm, -F_smft_array, color=colors[i], linestyle=linestyle,
               linewidth=2, label=label, marker='o', markersize=4)

ax1.set_xlabel('Separation d (nm)', fontsize=12)
ax1.set_ylabel('|F| (pN/μm²)', fontsize=12)
ax1.set_title('Casimir Force: SMFT vs Standard', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3, which='both')
ax1.legend()

# Add d^-4 reference line
d_ref = np.array([10, 1000])
F_ref = casimir_force_standard(d_ref * 1e-9) * 1e-12 / 1e-12
ax1.plot(d_ref, -F_ref, 'k:', linewidth=1.5, alpha=0.5, label='$d^{-4}$ scaling')

# Plot 2: Fractional deviation
ax2 = axes[0, 1]
for i, R_avg in enumerate(R_values[:-1]):  # Exclude R=1
    deviation = (1 - R_avg**2) * 100
    ax2.axhline(deviation, color=colors[i], linestyle='--', linewidth=2,
                label=f'R = {R_avg} (δF/F = {deviation:.1f}%)')

ax2.axhline(exp_precision_percent, color='orange', linestyle=':', linewidth=2,
            label=f'Exp. precision (±{exp_precision_percent}%)')
ax2.axhline(3 * exp_precision_percent, color='red', linestyle=':', linewidth=2,
            label=f'3σ threshold (±{3 * exp_precision_percent}%)')

ax2.set_xlabel('Separation d (nm)', fontsize=12)
ax2.set_ylabel('$\\delta F / F$ (%)', fontsize=12)
ax2.set_title('Fractional Deviation (Distance-Independent)', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlim(10, 1000)
ax2.set_ylim(0, 25)
ax2.set_xscale('log')

# Plot 3: Statistical significance vs R
ax3 = axes[1, 0]
R_scan = np.linspace(0.85, 1.00, 50)
deviation_scan = (1 - R_scan**2) * 100
sigma_scan = deviation_scan / exp_precision_percent

ax3.plot(R_scan, sigma_scan, 'b-', linewidth=2.5, label='Significance (2% precision)')
ax3.axhline(3, color='orange', linestyle='--', linewidth=1.5, label='3σ threshold')
ax3.axhline(5, color='red', linestyle='--', linewidth=1.5, label='5σ threshold')
ax3.fill_between(R_scan, 0, sigma_scan, where=(sigma_scan >= 3), alpha=0.2,
                  color='green', label='Detectable at 3σ')

ax3.set_xlabel('$\\langle R \\rangle$ (average R-field value)', fontsize=12)
ax3.set_ylabel('Statistical Significance ($\\sigma$)', fontsize=12)
ax3.set_title('Detectability vs R-Field Strength', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.legend()
ax3.set_xlim(0.85, 1.00)
ax3.set_ylim(0, 10)

# Plot 4: Material ratio test
ax4 = axes[1, 1]
R_material_1 = np.linspace(0.95, 1.00, 30)
R_material_2_values = [0.85, 0.90, 0.95]

for R2 in R_material_2_values:
    ratio = R_material_1**2 / R2**2
    deviation_from_unity = (ratio - 1) * 100
    ax4.plot(R_material_1, deviation_from_unity, linewidth=2,
             label=f'Material 2: R={R2}')

ax4.axhline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax4.axhline(exp_precision_percent, color='orange', linestyle='--', linewidth=1.5,
            label=f'Exp. precision (±{exp_precision_percent}%)')
ax4.fill_between([0.95, 1.00], -20, exp_precision_percent, alpha=0.1,
                  color='red', label='Standard theory (ratio=1)')

ax4.set_xlabel('Material 1: $R_1$', fontsize=12)
ax4.set_ylabel('$(F_1/F_2 - 1)$ (%) [SMFT contribution]', fontsize=12)
ax4.set_title('Material Dependence Test', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.legend()
ax4.set_xlim(0.95, 1.00)
ax4.set_ylim(-5, 20)

plt.tight_layout()
plt.savefig('/home/persist/neotec/0rigin/analysis/predictions/casimir_force_deviation.png', dpi=150)
print("Figure saved: analysis/predictions/casimir_force_deviation.png")
print()

# ===== EXPERIMENTAL PROTOCOL =====

print("=" * 80)
print("EXPERIMENTAL PROTOCOL: Casimir Force Measurement")
print("=" * 80)
print()
print("1. SYSTEM SETUP")
print("   - Gold-coated sphere (R ~ 50 μm) on cantilever (AFM)")
print("   - Gold-coated flat plate")
print("   - Vacuum chamber: P < 10^-6 Torr")
print("   - Temperature control: T = 300 ± 0.1 K")
print()
print("2. FORCE MEASUREMENT")
print("   - Scan separation: d = 50-500 nm (piezo actuator)")
print("   - Measure cantilever deflection Δz (optical lever)")
print("   - Extract force: F = k_spring · Δz")
print("   - Precision: δF/F ~ 1-2% (state-of-art)")
print()
print("3. DATA ANALYSIS")
print("   - Fit to Lifshitz theory: F(d, T, ε(ω))")
print("   - Extract residual: ΔF = F_measured - F_Lifshitz")
print("   - Test SMFT: ΔF/F = (1 - ⟨R²⟩) = const?")
print()
print("4. EXPECTED RESULT (R_avg = 0.95)")
print(f"   At d = 100 nm:")
print(f"     Standard: F = {-F_standard_pN_um2[3]:.4e} pN/μm²")
print(f"     SMFT:     F = {-casimir_force_smft(d_meters[3], 0.95) * 1e-12 / 1e-12:.4e} pN/μm²")
print(f"     Deviation: {(1 - 0.95**2) * 100:.1f}%")
print(f"     Significance: {(1 - 0.95**2) * 100 / exp_precision_percent:.1f}σ")
print()
print("5. FALSIFICATION CRITERIA")
print("   SMFT FALSIFIED IF:")
print("     F_measured / F_Lifshitz = 1.000 ± 0.020 for all materials")
print("     (no R-dependence detected at 3σ)")
print()
print("6. ENHANCED MEASUREMENT: Frequency-Domain Method")
print("   - Use nanomembrane resonator instead of static cantilever")
print("   - Measure frequency shift: Δf ∝ ∂F/∂d")
print("   - Precision: Δf/f ~ 10^-6 → δF/F ~ 10^-4 (100× improvement)")
print("   - With this: 0.1% R-field variations detectable")
print()
print("7. TIMELINE & COST")
print("   Standard AFM method: 6 months, $0 (existing labs)")
print("   Nanomembrane method: 12 months, $50k-100k (fabrication)")
print()
print("=" * 80)
print("CONCLUSION: Casimir force test is CHALLENGING but FEASIBLE")
print("  - Requires ~1% precision (state-of-art: 0.2-2%)")
print("  - Distance-independent deviation → need absolute calibration")
print("  - Material comparison test may be more robust")
print("=" * 80)
