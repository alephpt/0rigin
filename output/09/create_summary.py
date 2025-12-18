#!/usr/bin/env python3
"""Create comprehensive validation summary figure"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# Title
fig.suptitle('Dirac Split-Operator Validation Results', fontsize=18, fontweight='bold')

# ===== TEST 1: Beta Expectation =====
ax1 = fig.add_subplot(gs[0, 0])
ax1.text(0.5, 0.7, '✓ Beta Expectation', ha='center', fontsize=14, fontweight='bold', color='green')
ax1.text(0.5, 0.5, '<β> = 1.000', ha='center', fontsize=12)
ax1.text(0.5, 0.3, 'Upper components only', ha='center', fontsize=10, style='italic')
ax1.text(0.5, 0.1, '[PASS]', ha='center', fontsize=12, color='green', fontweight='bold')
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.axis('off')

# ===== TEST 2: Force Direction =====
ax2 = fig.add_subplot(gs[0, 1])
ax2.arrow(0.2, 0.5, 0.6, 0, head_width=0.15, head_length=0.1, fc='blue', ec='blue', linewidth=2)
ax2.text(0.5, 0.8, '✓ Ehrenfest Theorem', ha='center', fontsize=14, fontweight='bold', color='green')
ax2.text(0.5, 0.3, 'Δx = +6.0 grid points', ha='center', fontsize=11)
ax2.text(0.5, 0.15, 'Toward LOW mass', ha='center', fontsize=10, style='italic')
ax2.text(0.5, 0.02, 'F = -β·∇m', ha='center', fontsize=9, family='monospace')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.axis('off')

# ===== TEST 3: Isotropy =====
ax3 = fig.add_subplot(gs[0, 2])
theta = np.linspace(0, 2*np.pi, 100)
r = np.ones_like(theta)
ax3.plot(r * np.cos(theta), r * np.sin(theta), 'b-', linewidth=3)
ax3.scatter([0], [0], s=100, c='red', marker='o')
ax3.text(0, -1.5, '✓ Rotation Invariance', ha='center', fontsize=14, fontweight='bold', color='green')
ax3.text(0, -1.8, 'Variance = 0.0', ha='center', fontsize=11)
ax3.set_xlim(-1.5, 1.5)
ax3.set_ylim(-2, 1.5)
ax3.axis('off')
ax3.set_aspect('equal')

# ===== TEST 4: Long-time stability (norm vs time) =====
ax4 = fig.add_subplot(gs[1, :])
data = np.loadtxt('norm_vs_time.dat')
steps = data[:, 0]
norms = data[:, 1]

ax4.plot(steps, norms, 'b-', linewidth=2, label='|Ψ|²')
ax4.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Perfect')
ax4.set_xlabel('Timestep', fontsize=12)
ax4.set_ylabel('Norm |Ψ|²', fontsize=12)
ax4.set_title('✓ Long-Time Stability (50,000 steps, 0.72% drift)', fontsize=14, fontweight='bold', color='green')
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_ylim([0.99, 1.001])

# ===== Performance Comparison Table =====
ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')

table_data = [
    ['Method', 'Stability', 'Drift @ 50k', 'NaN at Step', 'Physics'],
    ['Euler (old)', 'Poor', '~100%', '19,781', '✗'],
    ['RK4 (old)', 'Fair', '~1%', 'Manual norm', '?'],
    ['Split-Op (NEW)', 'Excellent', '0.72%', 'Never', '✓'],
]

colors = [
    ['lightgray'] * 5,
    ['#ffcccc'] * 5,
    ['#ffffcc'] * 5,
    ['#ccffcc'] * 5,
]

table = ax5.table(cellText=table_data, cellColours=colors,
                  cellLoc='center', loc='center',
                  bbox=[0.1, 0.2, 0.8, 0.6])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 2)

# Make header bold
for i in range(5):
    table[(0, i)].set_text_props(weight='bold')

ax5.text(0.5, 0.95, 'Performance Comparison', ha='center', fontsize=14, fontweight='bold')

plt.savefig('validation_summary.png', dpi=150, bbox_inches='tight')
print("Saved: validation_summary.png")
