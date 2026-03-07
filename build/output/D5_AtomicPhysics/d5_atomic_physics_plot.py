#!/usr/bin/env python3
"""D5 Atomic Physics - Publication Quality Plots"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'figure.dpi': 150,
    'font.family': 'serif',
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.color': '#cccccc',
    'lines.linewidth': 1.8,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('D5 Atomic Physics from TRD', fontsize=16, fontweight='bold', y=0.98)

# --- Panel 1: Energy Levels ---
ax = axes[0, 0]
n_vals = np.arange(1, 11)
E_trd = -13.605693 / n_vals**2  # From CSV
E_theory = -13.6057 / n_vals**2  # Bohr model
ax.plot(n_vals, E_theory, 'o--', color='#2196F3', markersize=8, label='Bohr model', zorder=3)
ax.plot(n_vals, E_trd, 's', color='#FF5722', markersize=6, label='TRD', zorder=4)
ax.set_xlabel('Principal quantum number n')
ax.set_ylabel('Energy (eV)')
ax.set_title('Hydrogen Energy Levels')
ax.legend()
ax.set_xticks(n_vals)

# --- Panel 2: Balmer Series ---
ax = axes[0, 1]
lines = ['H-α\n(3→2)', 'H-β\n(4→2)', 'H-γ\n(5→2)', 'H-δ\n(6→2)', 'H-ε\n(7→2)']
lambda_trd = [656.112, 486.009, 433.937, 410.070, 396.907]
lambda_exp = [656.281, 486.135, 434.047, 410.175, 397.008]
errors_pct = [0.0257, 0.0259, 0.0254, 0.0256, 0.0253]

x = np.arange(len(lines))
width = 0.35
bars1 = ax.bar(x - width/2, lambda_exp, width, label='Experiment', color='#2196F3', alpha=0.85)
bars2 = ax.bar(x + width/2, lambda_trd, width, label='TRD', color='#FF9800', alpha=0.85)
ax.set_ylabel('Wavelength (nm)')
ax.set_title('Balmer Series')
ax.set_xticks(x)
ax.set_xticklabels(lines, fontsize=9)
ax.legend()
# Add error annotations
for i, err in enumerate(errors_pct):
    ax.annotate(f'{err:.3f}%', xy=(x[i], min(lambda_trd[i], lambda_exp[i]) - 15),
                ha='center', fontsize=7, color='green', fontweight='bold')

# --- Panel 3: Balmer Errors ---
ax = axes[0, 2]
colors_balmer = ['#e74c3c', '#e67e22', '#f1c40f', '#2ecc71', '#3498db']
ax.bar(range(len(errors_pct)), errors_pct, color=colors_balmer, edgecolor='black', linewidth=0.5)
ax.axhline(y=0.03, color='red', linestyle='--', alpha=0.7, label='0.03% threshold')
ax.set_ylabel('Error (%)')
ax.set_title('Balmer Series Accuracy')
ax.set_xticks(range(len(lines)))
ax.set_xticklabels(['H-α', 'H-β', 'H-γ', 'H-δ', 'H-ε'], fontsize=10)
ax.legend()
ax.set_ylim(0, 0.04)
ax.text(2, 0.035, 'All lines < 0.026% error', ha='center', fontsize=10,
        fontweight='bold', color='green',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f5e9', alpha=0.8))

# --- Panel 4: Fine Structure & Lamb Shift ---
ax = axes[1, 0]
measurements = ['Fine Structure\n2P₃/₂−2P₁/₂', 'Lamb Shift\n2S₁/₂−2P₁/₂', 'Hyperfine\n21cm line']
trd_vals = [10.949, 955.321, 1421.160]
exp_vals = [10.970, 1057.800, 1420.406]
errors = [0.189, 9.688, 0.053]
units = ['GHz', 'MHz', 'MHz']

x = np.arange(len(measurements))
width = 0.35
bars1 = ax.bar(x - width/2, exp_vals, width, label='Experiment', color='#2196F3', alpha=0.85)
bars2 = ax.bar(x + width/2, trd_vals, width, label='TRD', color='#FF9800', alpha=0.85)
ax.set_ylabel('Frequency')
ax.set_title('QED Corrections')
ax.set_xticks(x)
ax.set_xticklabels(measurements, fontsize=9)
ax.legend()
ax.set_yscale('log')
for i, (err, unit) in enumerate(zip(errors, units)):
    color = 'green' if err < 1 else '#e67e22'
    ax.annotate(f'{err:.2f}%\n({unit})', xy=(x[i] + width/2 + 0.05, trd_vals[i]),
                fontsize=8, color=color, fontweight='bold')

# --- Panel 5: Rydberg Constant ---
ax = axes[1, 1]
R_trd = 1.097373e7
R_exp = 1.097373e7
error_ppb = 2.38e-8 * 1e7  # convert % to ppb: 2.38e-8 % = 0.238 ppb
ax.barh(['R∞'], [error_ppb], color='#4CAF50', height=0.4, edgecolor='black')
ax.set_xlabel('Error (parts per billion)')
ax.set_title('Rydberg Constant Precision')
ax.set_xlim(0, max(0.5, error_ppb * 1.5))
ax.text(error_ppb + 0.01, 0, f'{error_ppb:.3f} ppb\nR∞ = {R_trd:.6e} m⁻¹',
        va='center', fontsize=10, fontweight='bold', color='#2e7d32',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#e8f5e9', alpha=0.8))

# --- Panel 6: Summary scorecard ---
ax = axes[1, 2]
ax.axis('off')
summary = [
    ('Energy Levels (n=1-10)', '< 10⁻⁶ %', '#4CAF50'),
    ('Balmer Series (5 lines)', '< 0.026 %', '#4CAF50'),
    ('Fine Structure (2P split)', '0.19 %', '#4CAF50'),
    ('Hyperfine 21cm', '0.053 %', '#4CAF50'),
    ('Rydberg Constant', '0.24 ppb', '#4CAF50'),
    ('Lamb Shift', '9.7 %', '#FF9800'),
]
ax.set_title('Validation Summary', fontsize=13)
for i, (test, error, color) in enumerate(summary):
    y = 0.85 - i * 0.14
    marker = '●' if color == '#4CAF50' else '◐'
    ax.text(0.05, y, marker, fontsize=16, color=color, transform=ax.transAxes,
            va='center', fontweight='bold')
    ax.text(0.12, y, test, fontsize=11, transform=ax.transAxes, va='center')
    ax.text(0.75, y, error, fontsize=11, transform=ax.transAxes, va='center',
            fontweight='bold', color=color,
            bbox=dict(boxstyle='round,pad=0.2', facecolor=color, alpha=0.15))

ax.text(0.5, 0.02, '5/6 tests within 0.2% of experiment',
        ha='center', fontsize=12, transform=ax.transAxes,
        fontweight='bold', color='#1B5E20',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#C8E6C9', alpha=0.9))

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('d5_atomic_physics_plot.png', dpi=150)
print('Saved d5_atomic_physics_plot.png')
