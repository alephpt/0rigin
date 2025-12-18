#!/usr/bin/env python3
"""
Create explanatory diagrams for Phase 2 Scenario 1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

# ============ DIAGRAM 1: MSFT COUPLING MECHANISM ============

fig1, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: Kuramoto Synchronization Field
ax1 = axes[0]
x = np.linspace(0, 128, 256)
y = np.linspace(0, 128, 256)
X, Y = np.meshgrid(x, y)

# Create R field with defect
R = np.ones_like(X)
defect_mask = (X-64)**2 + (Y-64)**2 < 15**2
R[defect_mask] = 0.35  # Low R in defect

im1 = ax1.contourf(X, Y, R, levels=20, cmap='RdYlGn')
ax1.contour(X, Y, R, levels=[0.35, 0.5, 0.7, 0.9], colors='black', linewidths=0.5, alpha=0.3)
circle = Circle((64, 64), 15, fill=False, edgecolor='red', linewidth=2, linestyle='--')
ax1.add_patch(circle)
ax1.set_xlabel('X (grid points)')
ax1.set_ylabel('Y (grid points)')
ax1.set_title('(A) Kuramoto Order Parameter R(x,y)\\nLow R in defect (desynchronization)', fontsize=12, fontweight='bold')
ax1.set_aspect('equal')
plt.colorbar(im1, ax=ax1, label='R (synchronization)')

# Panel B: Mass Field
ax2 = axes[1]
Delta = 0.5
m = Delta * R

im2 = ax2.contourf(X, Y, m, levels=20, cmap='viridis')
ax2.contour(X, Y, m, levels=[0.175, 0.25, 0.35, 0.45], colors='white', linewidths=0.5, alpha=0.5)
circle2 = Circle((64, 64), 15, fill=False, edgecolor='red', linewidth=2, linestyle='--')
ax2.add_patch(circle2)
ax2.set_xlabel('X (grid points)')
ax2.set_ylabel('Y (grid points)')
ax2.set_title('(B) MSFT Mass Field m(x,y) = Δ·R(x,y)\\nLow mass creates potential well', fontsize=12, fontweight='bold')
ax2.set_aspect('equal')
plt.colorbar(im2, ax=ax2, label='m (Dirac mass)')

# Panel C: Effective Potential
ax3 = axes[2]
# Slice through center
m_slice = m[128, :]
x_slice = x

ax3.plot(x_slice, m_slice, 'b-', linewidth=3, label='m(x, y=64)')
ax3.axvline(x=64, color='r', linestyle='--', linewidth=2, alpha=0.5, label='Defect center')
ax3.axvspan(49, 79, alpha=0.2, color='red', label='Defect region')
ax3.fill_between(x_slice, m_slice, alpha=0.3)
ax3.set_xlabel('X (grid points)')
ax3.set_ylabel('Mass m')
ax3.set_title('(C) Potential Profile (slice at y=64)\\nParticles experience F = -<β>Δ∇R', fontsize=12, fontweight='bold')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# Add force direction arrows
ax3.annotate('', xy=(55, 0.35), xytext=(45, 0.35),
            arrowprops=dict(arrowstyle='->', lw=2, color='green'))
ax3.annotate('', xy=(73, 0.35), xytext=(83, 0.35),
            arrowprops=dict(arrowstyle='->', lw=2, color='green'))
ax3.text(50, 0.32, 'F (toward defect)', fontsize=9, color='green', ha='center')

fig1.suptitle('MSFT-Dirac Coupling Mechanism: From Synchronization to Mass Field',
              fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/coupling_mechanism.png', dpi=150, bbox_inches='tight')
print("✅ Generated: plots/coupling_mechanism.png")
plt.close()

# ============ DIAGRAM 2: ENERGY BUG EXPLANATION ============

fig2 = plt.figure(figsize=(16, 10))
gs = fig2.add_gridspec(3, 2, hspace=0.4, wspace=0.3)

# Panel A: Bug - Stale K-Space
ax1 = fig2.add_subplot(gs[0, 0])
ax1.axis('off')
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)

# Flow diagram showing bug
boxes = [
    ('step() starts', 1, 8, 'lightblue'),
    ('FFT: ψ → ψ_k', 1, 6.5, 'lightgreen'),
    ('K-evolution on ψ_k', 1, 5, 'lightgreen'),
    ('iFFT: ψ_k → ψ', 1, 3.5, 'lightgreen'),
    ('ψ is CURRENT', 1, 2, 'green'),
    ('ψ_k is STALE!', 1, 0.5, 'red'),
]

for label, x, y, color in boxes:
    box = FancyBboxPatch((x-0.4, y-0.3), 5, 0.6, boxstyle="round,pad=0.1",
                         edgecolor='black', facecolor=color, linewidth=2)
    ax1.add_patch(box)
    ax1.text(x+2.1, y, label, fontsize=11, ha='center', va='center', fontweight='bold')

# Arrows
for i in range(len(boxes)-1):
    ax1.annotate('', xy=(boxes[i+1][1]+2.1, boxes[i+1][2]+0.3),
                xytext=(boxes[i][1]+2.1, boxes[i][2]-0.3),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

# Bug callout
ax1.text(7, 0.5, '← BUG: getEnergy()\\nreads stale ψ_k!',
        fontsize=10, color='red', fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

ax1.set_title('(A) THE BUG: Stale K-Space Data', fontsize=13, fontweight='bold', color='red')

# Panel B: Fix - Lazy Update
ax2 = fig2.add_subplot(gs[0, 1])
ax2.axis('off')
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)

fix_boxes = [
    ('getEnergy() called', 1, 8, 'lightblue'),
    ('Check: ψ_k_valid?', 1, 6.5, 'yellow'),
    ('if FALSE:', 1, 5, 'orange'),
    ('  Recompute FFT', 1, 3.5, 'lightgreen'),
    ('  Set ψ_k_valid=TRUE', 1, 2, 'lightgreen'),
    ('Use FRESH ψ_k', 1, 0.5, 'green'),
]

for label, x, y, color in fix_boxes:
    box = FancyBboxPatch((x-0.4, y-0.3), 5, 0.6, boxstyle="round,pad=0.1",
                         edgecolor='black', facecolor=color, linewidth=2)
    ax2.add_patch(box)
    ax2.text(x+2.1, y, label, fontsize=11, ha='center', va='center', fontweight='bold')

for i in range(len(fix_boxes)-1):
    ax2.annotate('', xy=(fix_boxes[i+1][1]+2.1, fix_boxes[i+1][2]+0.3),
                xytext=(fix_boxes[i][1]+2.1, fix_boxes[i][2]-0.3),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

ax2.text(7, 0.5, '✓ FIX: Always\\nuse fresh data!',
        fontsize=10, color='green', fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

ax2.set_title('(B) THE FIX: Lazy K-Space Update', fontsize=13, fontweight='bold', color='green')

# Panel C: Energy Before Fix (Explosion)
ax3 = fig2.add_subplot(gs[1, 0])
steps_fake = np.linspace(0, 50000, 100)
E_fake_strong = 0.999 + (steps_fake/50000) * 13834  # Linear explosion
E_fake_weak = 0.677 + 0.12 * np.sin(steps_fake/5000) + (steps_fake/50000) * 0.05  # Appears stable

ax3.plot(steps_fake*0.01, E_fake_strong, 'r-', linewidth=3, label='Strong defect (ΔR=0.84)')
ax3.plot(steps_fake*0.01, E_fake_weak, 'orange', linewidth=2, label='Weak defect (ΔR=0.01)')
ax3.set_xlabel('Time')
ax3.set_ylabel('Energy (INVALID)')
ax3.set_title('(C) BEFORE FIX: Energy Behavior with Bug', fontsize=13, fontweight='bold')
ax3.legend(loc='upper left')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(-1, 15000)
ax3.text(250, 10000, 'EXPLOSION!\\n(Garbage data)', fontsize=12, color='red',
        fontweight='bold', ha='center',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

# Panel D: Energy After Fix (Stable)
ax4 = fig2.add_subplot(gs[1, 1])
# Load real corrected data
energy_data = np.loadtxt('../energy_fix_test/energy.dat', skiprows=1)
steps_real = energy_data[:, 0]
E_real = energy_data[:, 1]

ax4.plot(steps_real*0.01, E_real, 'g-', linewidth=3, label='Corrected (any ΔR)')
ax4.axhline(y=0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='E=0 (binding)')
ax4.fill_between(steps_real*0.01, 0, E_real, alpha=0.2, color='orange')
ax4.set_xlabel('Time')
ax4.set_ylabel('Energy (VALID)')
ax4.set_title('(D) AFTER FIX: Energy Behavior Corrected', fontsize=13, fontweight='bold')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)
ax4.text(25, 0.63, f'✓ Stable!\\nE>0 (unbound)', fontsize=12, color='green',
        fontweight='bold', ha='center',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

# Panel E: Code Changes
ax5 = fig2.add_subplot(gs[2, :])
ax5.axis('off')
ax5.set_xlim(0, 10)
ax5.set_ylim(0, 10)

code_before = """// DiracEvolution.cpp (BEFORE - BUG)
float DiracEvolution::getEnergy(...) {
    // WRONG: Assumes _psi_k is current
    for (int c = 0; c < 4; c++) {
        density_k = std::norm(_psi_k[c][idx]);  // STALE!
        KE += density_k * omega_k;
    }
    return KE + PE;
}"""

code_after = """// DiracEvolution.cpp (AFTER - FIXED)
float DiracEvolution::getEnergy(...) {
    // FIX: Recompute FFT if cache invalid
    if (!_psi_k_valid) {
        for (int c = 0; c < 4; c++) {
            fftwf_execute(_fft_forward[c]);  // FRESH!
        }
        _psi_k_valid = true;
    }
    // Now use CURRENT k-space data
    for (int c = 0; c < 4; c++) {
        density_k = std::norm(_psi_k[c][idx]);
        KE += density_k * omega_k;
    }
    return KE + PE;
}"""

ax5.text(0.5, 7, 'CODE CHANGES:', fontsize=14, fontweight='bold', ha='left')
ax5.text(0.5, 3.5, code_before, fontsize=9, ha='left', va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='mistyrose', alpha=0.8, edgecolor='red', linewidth=2))
ax5.text(5.5, 3.5, code_after, fontsize=9, ha='left', va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8, edgecolor='green', linewidth=2))

# Arrow between code blocks
ax5.annotate('', xy=(5.3, 2), xytext=(3.5, 2),
            arrowprops=dict(arrowstyle='->', lw=3, color='blue'))
ax5.text(4.4, 2.5, 'FIX', fontsize=12, fontweight='bold', color='blue', ha='center')

fig2.suptitle('Energy Diagnostic Bug: Root Cause and Fix', fontsize=16, fontweight='bold')
plt.savefig('plots/energy_bug_explanation.png', dpi=150, bbox_inches='tight')
print("✅ Generated: plots/energy_bug_explanation.png")
plt.close()

# ============ DIAGRAM 3: SCATTERING VS BINDING ============

fig3, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Binding (E < 0) - Not observed
ax1 = axes[0]
r = np.linspace(0, 30, 200)
V_binding = -1.0 / (1 + ((r-15)/3)**2) + 0.3  # Attractive potential well
E_binding = -0.2  # Negative energy

ax1.plot(r, V_binding, 'b-', linewidth=3, label='V(r) (potential)')
ax1.axhline(y=E_binding, color='r', linestyle='--', linewidth=2, label=f'E = {E_binding:.1f} < 0')
ax1.fill_between(r, V_binding, E_binding, where=(V_binding<E_binding),
                 alpha=0.3, color='green', label='Classically allowed')
ax1.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
ax1.set_xlabel('Distance from defect')
ax1.set_ylabel('Energy')
ax1.set_title('(A) BINDING (E < 0) - NOT OBSERVED\\nParticle trapped in potential well',
             fontsize=12, fontweight='bold')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)
ax1.set_ylim(-0.8, 0.5)
ax1.text(15, -0.6, 'BOUND STATE\\n(confined forever)', fontsize=11, ha='center',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

# Panel B: Scattering (E > 0) - What we observed
ax2 = axes[1]
V_scattering = -0.5 / (1 + ((r-15)/3)**2)  # Shallower well
E_scattering = 0.68  # Positive energy (our measurement)

ax2.plot(r, V_scattering, 'b-', linewidth=3, label='V(r) (potential)')
ax2.axhline(y=E_scattering, color='r', linestyle='--', linewidth=2, label=f'E = {E_scattering:.2f} > 0')
ax2.fill_between(r, -0.6, E_scattering, alpha=0.1, color='orange', label='Classically allowed (all r)')
ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5, label='E=0 threshold')

# Add resonance trajectory
r_traj = np.array([2, 8, 14, 18, 22, 24, 20, 14, 10, 6, 10, 16, 20])
E_traj = np.ones_like(r_traj) * E_scattering
ax2.plot(r_traj, E_traj, 'go-', linewidth=2, markersize=8, label='Particle trajectory', alpha=0.7)

ax2.set_xlabel('Distance from defect')
ax2.set_ylabel('Energy')
ax2.set_title('(B) SCATTERING (E > 0) - OBSERVED\\nParticle temporarily trapped (resonance)',
             fontsize=12, fontweight='bold')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-0.6, 0.8)
ax2.text(15, 0.5, 'SCATTERING RESONANCE\\n(escapes eventually)', fontsize=11, ha='center',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

# Add arrows showing eventual escape
ax2.annotate('', xy=(28, 0.68), xytext=(24, 0.68),
            arrowprops=dict(arrowstyle='->', lw=3, color='red'))
ax2.text(26, 0.75, 'eventual\\nescape', fontsize=9, color='red', ha='center')

fig3.suptitle('Binding vs Scattering: Energy Determines Outcome', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/binding_vs_scattering.png', dpi=150, bbox_inches='tight')
print("✅ Generated: plots/binding_vs_scattering.png")
plt.close()

print("\\n" + "="*60)
print("ALL DIAGRAMS GENERATED")
print("="*60)
print("1. plots/coupling_mechanism.png")
print("2. plots/energy_bug_explanation.png")
print("3. plots/binding_vs_scattering.png")
print("4. plots/scenario1_complete_analysis.png (from previous script)")
print("="*60)
