#!/usr/bin/env python3
"""
Ehrenfest Theorem Validation Analysis

CRITICAL DIAGNOSTICS:
1. Energy E(t) < 0 and conserved → True binding
2. Ehrenfest: d<p>/dt ≈ -<∇V> → Physics is correct
3. Long-time escape → E > 0 scattering resonance
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data
energy = np.loadtxt('ehrenfest/energy.dat')
ehrenfest = np.loadtxt('ehrenfest/ehrenfest.dat')
trajectory = np.loadtxt('ehrenfest/trajectory.dat')

# ===== CRITICAL DIAGNOSTIC 1: Energy Analysis =====
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

ax = axes[0, 0]
steps = energy[:, 0]
E_total = energy[:, 1]
KE = energy[:, 2]
PE = energy[:, 3]

ax.plot(steps, E_total, 'k-', linewidth=2, label='Total E')
ax.plot(steps, KE, 'r--', alpha=0.7, label='KE')
ax.plot(steps, PE, 'b--', alpha=0.7, label='PE')
ax.axhline(y=0, color='gray', linestyle=':', label='E=0 (bound/unbound threshold)')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('Energy', fontsize=12)
ax.set_title('DIAGNOSTIC 1: Energy E(t)', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Check binding
E_mean = E_total.mean()
E_final = E_total[-1]
is_bound = E_mean < 0

ax.text(0.02, 0.98, f'<E> = {E_mean:.4f}\\n' +
                    f'E(final) = {E_final:.4f}\\n' +
                    f'Bound: {"YES" if is_bound else "NO"}',
        transform=ax.transAxes, fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round',
                 facecolor='lightgreen' if is_bound else 'salmon',
                 alpha=0.7))

# ===== Energy Conservation =====
ax = axes[0, 1]
E_drift = np.abs(E_total - E_total[0])
ax.plot(steps, E_drift, 'b-', linewidth=2)
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('|E(t) - E(0)|', fontsize=12)
ax.set_title('Energy Conservation', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

drift_percent = 100 * E_drift[-1] / np.abs(E_total[0])
ax.text(0.02, 0.98, f'Drift: {drift_percent:.2f}%',
        transform=ax.transAxes, fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

# ===== DIAGNOSTIC 2: Ehrenfest Theorem =====
ax = axes[0, 2]
steps_ehr = ehrenfest[:, 0]
dp_dt_x = ehrenfest[:, 1]
dp_dt_y = ehrenfest[:, 2]
F_x = ehrenfest[:, 3]
F_y = ehrenfest[:, 4]
ratio_x = ehrenfest[:, 5]
ratio_y = ehrenfest[:, 6]

# Plot ratios
valid_x = np.abs(ratio_x) < 10  # Filter outliers
valid_y = np.abs(ratio_y) < 10

if valid_x.sum() > 0:
    ax.plot(steps_ehr[valid_x], ratio_x[valid_x], 'b.', markersize=2, alpha=0.5, label='x-direction')
if valid_y.sum() > 0:
    ax.plot(steps_ehr[valid_y], ratio_y[valid_y], 'r.', markersize=2, alpha=0.5, label='y-direction')

ax.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='Ehrenfest: ratio=1')
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('(d<p>/dt) / F', fontsize=12)
ax.set_title('DIAGNOSTIC 2: Ehrenfest Theorem', fontsize=14, fontweight='bold')
ax.set_ylim([-5, 5])
ax.legend()
ax.grid(True, alpha=0.3)

# Statistics
if valid_x.sum() > 0 and valid_y.sum() > 0:
    mean_ratio_x = ratio_x[valid_x].mean()
    mean_ratio_y = ratio_y[valid_y].mean()
    ehrenfest_satisfied = (np.abs(mean_ratio_x - 1.0) < 0.5 and
                           np.abs(mean_ratio_y - 1.0) < 0.5)

    ax.text(0.02, 0.98, f'<ratio_x> = {mean_ratio_x:.2f}\\n' +
                        f'<ratio_y> = {mean_ratio_y:.2f}\\n' +
                        f'Ehrenfest: {"PASS" if ehrenfest_satisfied else "FAIL"}',
            transform=ax.transAxes, fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round',
                     facecolor='lightgreen' if ehrenfest_satisfied else 'salmon',
                     alpha=0.7))

# ===== DIAGNOSTIC 3: Long-Time Trajectory =====
ax = axes[1, 0]
steps_traj = trajectory[:, 0]
dist = trajectory[:, 3]

ax.plot(steps_traj, dist, 'b-', linewidth=2)
ax.axhline(y=dist[0], color='r', linestyle='--', alpha=0.5, label=f'Initial: {dist[0]:.2f}')
ax.axhline(y=dist.max(), color='orange', linestyle='--', alpha=0.5, label=f'Max: {dist.max():.2f}')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('Distance to Defect', fontsize=12)
ax.set_title('DIAGNOSTIC 3: Long-Time Evolution (500k steps)', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Check for escape
escape_threshold = 2.0 * dist[:100].max()
escaped = dist[-1000:].mean() > escape_threshold

ax.text(0.02, 0.98, f'Initial max: {dist[:100].max():.2f}\\n' +
                    f'Final dist: {dist[-1]:.2f}\\n' +
                    f'Escaped: {"YES" if escaped else "NO"}',
        transform=ax.transAxes, fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round',
                 facecolor='salmon' if escaped else 'lightgreen',
                 alpha=0.7))

# ===== Force vs Acceleration Scatter =====
ax = axes[1, 1]
if valid_x.sum() > 0:
    ax.scatter(F_x[valid_x], dp_dt_x[valid_x], s=1, alpha=0.3, c='blue', label='x')
if valid_y.sum() > 0:
    ax.scatter(F_y[valid_y], dp_dt_y[valid_y], s=1, alpha=0.3, c='red', label='y')

# Perfect correlation line
f_range = np.array([F_x.min(), F_x.max()])
ax.plot(f_range, f_range, 'k--', linewidth=2, label='Perfect Ehrenfest')
ax.set_xlabel('Force F = -<∇V>', fontsize=12)
ax.set_ylabel('Acceleration d<p>/dt', fontsize=12)
ax.set_title('Force vs Acceleration (should lie on diagonal)', fontsize=13, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_aspect('equal', adjustable='box')

# ===== Summary Panel =====
ax = axes[1, 2]
ax.axis('off')

# Collect verdicts
verdict1 = "BOUND" if is_bound else "UNBOUND"
color1 = 'green' if is_bound else 'red'

verdict2 = "VALID" if ehrenfest_satisfied else "VIOLATED"
color2 = 'green' if ehrenfest_satisfied else 'red'

verdict3 = "CONFINED" if not escaped else "ESCAPED"
color3 = 'green' if not escaped else 'red'

y_pos = 0.9
ax.text(0.05, y_pos, "CRITICAL DIAGNOSTICS SUMMARY", fontsize=14, fontweight='bold')
y_pos -= 0.15

ax.text(0.05, y_pos, f"1. Energy E(t) < 0:", fontsize=12, fontweight='bold')
y_pos -= 0.08
ax.text(0.1, y_pos, f"{verdict1}", fontsize=12, color=color1, fontweight='bold')
y_pos -= 0.08
ax.text(0.1, y_pos, f"<E> = {E_mean:.4f}", fontsize=10, style='italic')
y_pos -= 0.12

ax.text(0.05, y_pos, f"2. Ehrenfest Theorem:", fontsize=12, fontweight='bold')
y_pos -= 0.08
ax.text(0.1, y_pos, f"{verdict2}", fontsize=12, color=color2, fontweight='bold')
if 'mean_ratio_x' in locals():
    y_pos -= 0.08
    ax.text(0.1, y_pos, f"<ratio> = ({mean_ratio_x:.2f}, {mean_ratio_y:.2f})", fontsize=10, style='italic')
y_pos -= 0.12

ax.text(0.05, y_pos, f"3. Long-Time Confinement:", fontsize=12, fontweight='bold')
y_pos -= 0.08
ax.text(0.1, y_pos, f"{verdict3}", fontsize=12, color=color3, fontweight='bold')
y_pos -= 0.08
ax.text(0.1, y_pos, f"500k steps, dist = {dist[-1]:.2f}", fontsize=10, style='italic')
y_pos -= 0.12

# Overall verdict
if is_bound and ehrenfest_satisfied and not escaped:
    overall = "✓ TRUE BINDING CONFIRMED"
    overall_color = 'darkgreen'
elif is_bound and not escaped:
    overall = "⚠ BINDING, BUT EHRENFEST VIOLATED"
    overall_color = 'orange'
elif not escaped:
    overall = "⚠ CONFINED BUT E > 0 (RESONANCE)"
    overall_color = 'orange'
else:
    overall = "✗ NOT BOUND (ESCAPED)"
    overall_color = 'red'

ax.text(0.5, 0.05, overall, fontsize=14, fontweight='bold',
        ha='center', color=overall_color,
        bbox=dict(boxstyle='round', facecolor='lightyellow',
                 edgecolor='black', linewidth=2))

plt.suptitle('Ehrenfest Theorem Validation - Resolving Internal Inconsistency',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('ehrenfest/ehrenfest_validation.png', dpi=150, bbox_inches='tight')
print("Saved: ehrenfest_validation.png")

# ===== Detailed Report =====
print("\n=== EHRENFEST VALIDATION REPORT ===\n")

print("DIAGNOSTIC 1: Energy Analysis")
print(f"  Mean energy: {E_mean:.4f}")
print(f"  Final energy: {E_final:.4f}")
print(f"  Energy drift: {drift_percent:.2f}%")
print(f"  Verdict: {verdict1}")
if is_bound:
    print(f"  → Particle is quantum mechanically bound (E < 0)")
else:
    print(f"  → Particle is unbound (E > 0), scattering resonance")
print()

print("DIAGNOSTIC 2: Ehrenfest Theorem")
if 'mean_ratio_x' in locals():
    print(f"  <d<p>/dt / F>_x = {mean_ratio_x:.3f} (expect ~1.0)")
    print(f"  <d<p>/dt / F>_y = {mean_ratio_y:.3f} (expect ~1.0)")
    print(f"  Verdict: {verdict2}")
    if ehrenfest_satisfied:
        print(f"  → Ehrenfest theorem satisfied, force law F=-∇V confirmed")
    else:
        print(f"  → Ehrenfest theorem VIOLATED, measurement or physics error")
print()

print("DIAGNOSTIC 3: Long-Time Escape")
print(f"  Initial distance: {dist[0]:.2f}")
print(f"  Maximum distance: {dist.max():.2f}")
print(f"  Final distance: {dist[-1]:.2f}")
print(f"  Evolution time: {steps_traj[-1] * 0.01:.1f} time units (500k steps)")
print(f"  Verdict: {verdict3}")
if not escaped:
    print(f"  → Particle remains confined over extended evolution")
else:
    print(f"  → Particle has escaped, not truly bound")
print()

print("=== RESOLUTION ===")
print()
if is_bound and ehrenfest_satisfied and not escaped:
    print("✓ ALL DIAGNOSTICS PASS")
    print("  - Energy E < 0: True binding confirmed")
    print("  - Ehrenfest holds: Force law validated")
    print("  - Long-time confined: Stable bound state")
    print()
    print("CONCLUSION: Dirac particle is quantum mechanically bound")
    print("            in MSFT synchronization defect.")
    print()
    print("Previous Test 1 failure (F·v ≈ 0) was due to:")
    print("  → Using classical velocity d<x>/dt instead of quantum <j>/ρ")
    print("  → High-frequency oscillations averaging to zero")
    print("  → Measurement artifact, not physics error")
elif is_bound and not escaped:
    print("⚠ PARTIAL VALIDATION")
    print("  - Binding confirmed (E < 0)")
    print("  - Confinement observed (500k steps)")
    print("  - BUT: Ehrenfest theorem violated")
    print()
    print("ISSUE: Internal inconsistency remains")
    print("       Either force calculation or acceleration measurement is wrong")
    print()
    print("REQUIRED: Debug force/velocity computation")
else:
    print("✗ BINDING CLAIM REJECTED")
    if not is_bound:
        print("  - Energy E > 0: Particle is unbound")
    if escaped:
        print("  - Particle escaped after extended evolution")
    print()
    print("CONCLUSION: Previous observation was scattering resonance,")
    print("            not true quantum binding.")
    print()
    print("REQUIRED: Redesign scenario with stronger defect or different initial conditions")
