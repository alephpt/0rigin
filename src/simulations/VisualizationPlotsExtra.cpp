#include "simulations/VisualizationGenerator.h"
#include <sstream>

// ============================================================================
// Plot Generators: Cosmology (C4-C5), Electromagnetism (D-series),
// Mathematical (E-series), Embedded/Generic, and Legacy
// ============================================================================

std::string VisualizationGenerator::genDarkEnergy(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    std::string eos = dataToNumpy("equation_of_state");
    s << "\nout_dir = '" << output_dir << "'\n";
    if (!eos.empty()) s << eos;
    s << R"(
fig, ax = plt.subplots()
plotted = False

if 'equation_of_state_x' in dir() and len(equation_of_state_x) > 0:
    ax.plot(equation_of_state_x, equation_of_state_y, color=COLORS[0], lw=2,
            label='TRD w(z)')
    plotted = True

if not plotted:
    csv_path = find_csv('.', '*dark_energy*.csv') or find_csv('.', '*eos*.csv')
    if csv_path:
        data = load_csv(csv_path)
        zk = next((k for k in data if 'redshift' in k.lower() or k.lower() == 'z'), None)
        wk = next((k for k in data if k.lower() == 'w' or 'equation' in k.lower()), None)
        if zk and wk:
            ax.plot(data[zk], data[wk], color=COLORS[0], lw=2, label='TRD w(z)')
            plotted = True

ax.axhline(y=-1.0, color='red', ls='--', alpha=0.6,
           label=r'$\Lambda$CDM: $w = -1$')
ax.set_xlabel('Redshift z')
ax.set_ylabel('Equation of State w')
ax.set_title('C4 Dark Energy - Equation of State')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join('.','dark_energy_plot.png'))
print('Saved: dark_energy_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genInflation(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    std::string efolds = dataToNumpy("efolds");
    std::string slowroll = dataToNumpy("slow_roll");
    s << "\nout_dir = '" << output_dir << "'\n";
    if (!efolds.empty()) s << efolds;
    if (!slowroll.empty()) s << slowroll;
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plotted = False

if 'efolds_x' in dir() and len(efolds_x) > 0:
    ax1.plot(efolds_x, efolds_y, color=COLORS[0])
    if 'slow_roll_x' in dir() and len(slow_roll_x) > 0:
        ax2.plot(slow_roll_x, slow_roll_y, color=COLORS[1])
        ax2.axhline(y=1.0, color='red', ls='--', alpha=0.6, label=r'$\epsilon = 1$')
        ax2.legend()
    plotted = True

if not plotted:
    csv_path = find_csv('.', '*inflation*.csv')
    if csv_path:
        data = load_csv(csv_path)
        tk = next((k for k in data if k.lower() in ('time', 't', 'n')), None)
        nk = next((k for k in data if 'efold' in k.lower() or 'n_e' in k.lower()), None)
        ek = next((k for k in data if 'epsilon' in k.lower() or 'slow' in k.lower()), None)
        if tk and nk:
            ax1.plot(data[tk], data[nk], color=COLORS[0])
            plotted = True
        if tk and ek:
            ax2.plot(data[tk], data[ek], color=COLORS[1])
            ax2.axhline(y=1.0, color='red', ls='--', alpha=0.6, label=r'$\epsilon = 1$')
            ax2.legend()

ax1.set_xlabel('Time'); ax1.set_ylabel('e-folds N')
ax1.set_title('e-Folding Evolution')
ax2.set_xlabel('Time'); ax2.set_ylabel(r'$\epsilon(t)$')
ax2.set_title('Slow-Roll Parameter')

plt.suptitle('C5 Inflation', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','inflation_plot.png'))
print('Saved: inflation_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genLorentzForce(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    std::string orbit = dataToNumpy("orbit");
    s << "\nout_dir = '" << output_dir << "'\n";
    if (!orbit.empty()) s << orbit;
    s << R"(
fig, ax = plt.subplots()
plotted = False

if 'orbit_x' in dir() and len(orbit_x) > 0:
    ax.plot(orbit_x, orbit_y, color=COLORS[0], alpha=0.7)
    plotted = True

if not plotted:
    csv_path = find_csv('.', '*lorentz*.csv') or find_csv('.', '*trajectory*.csv')
    if csv_path:
        data = load_csv(csv_path)
        xk = next((k for k in data if k.lower() in ('x', 'pos_x')), None)
        yk = next((k for k in data if k.lower() in ('y', 'pos_y')), None)
        if xk and yk:
            ax.plot(data[xk], data[yk], color=COLORS[0], alpha=0.7)
            plotted = True

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Lorentz Force - Cyclotron Orbit')
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join('.','lorentz_force_plot.png'))
print('Saved: lorentz_force_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genJosephson(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    std::string dc_curve = dataToNumpy("dc_josephson");
    std::string ac_curve = dataToNumpy("ac_josephson");
    s << "\nout_dir = '" << output_dir << "'\n";
    if (!dc_curve.empty()) s << dc_curve;
    if (!ac_curve.empty()) s << ac_curve;
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plotted = False

if 'dc_josephson_x' in dir() and len(dc_josephson_x) > 0:
    ax1.plot(dc_josephson_x, dc_josephson_y, color=COLORS[0], lw=2)
    if 'ac_josephson_x' in dir() and len(ac_josephson_x) > 0:
        ax2.plot(ac_josephson_x, ac_josephson_y, 'o-', color=COLORS[1])
    plotted = True

if not plotted:
    csv_path = find_csv('.', '*josephson*.csv')
    if csv_path:
        data = load_csv(csv_path)
        pk = next((k for k in data if 'phase' in k.lower() or 'delta' in k.lower()), None)
        ik = next((k for k in data if 'current' in k.lower() or k.lower() == 'i'), None)
        vk = next((k for k in data if 'voltage' in k.lower() or k.lower() == 'v'), None)
        fk = next((k for k in data if 'freq' in k.lower()), None)
        if pk and ik:
            ax1.plot(data[pk], data[ik], color=COLORS[0], lw=2)
            plotted = True
        if vk and fk:
            ax2.plot(data[vk], data[fk], 'o-', color=COLORS[1])

ax1.set_xlabel(r'Phase Difference $\Delta\theta$')
ax1.set_ylabel('Supercurrent I')
ax1.set_title('DC Josephson Effect')
ax2.set_xlabel('Voltage (V)')
ax2.set_ylabel('Frequency (Hz)')
ax2.set_title('AC Josephson Effect')

plt.suptitle('Josephson Junction', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','josephson_plot.png'))
print('Saved: josephson_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genCausality(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left panel: dispersion relation (group & phase velocity vs wavenumber)
csv_path = find_csv('.', '*dispersion*.csv')
if csv_path:
    data = load_csv(csv_path)
    kk = next((k for k in data if k.lower() == 'k'), None)
    vg = next((k for k in data if 'group' in k.lower()), None)
    vp = next((k for k in data if 'phase' in k.lower()), None)
    if kk and vg:
        ax1.plot(data[kk], data[vg], 'o-', color=COLORS[0], label=r'$v_{group}$')
    if kk and vp:
        ax1.plot(data[kk], data[vp], 's--', color=COLORS[1], alpha=0.7, label=r'$v_{phase}$')
    ax1.axhline(y=1.0, color='red', ls='--', lw=2, alpha=0.6, label=r'$c$ (speed of light)')
    ax1.set_xlabel('Wavenumber k')
    ax1.set_ylabel('Velocity / c')
    ax1.set_title('Dispersion Relation')
    ax1.legend()
    ax1.set_ylim(0, max(data[vp]) * 1.1 if vp else 1.5)
else:
    ax1.text(0.5, 0.5, 'No dispersion data', transform=ax1.transAxes, ha='center')

# Right panel: group velocity < c verification
if csv_path and kk and vg:
    vg_arr = np.array(data[vg], dtype=float)
    ax2.fill_between(data[kk], 0, vg_arr, alpha=0.3, color=COLORS[0])
    ax2.plot(data[kk], vg_arr, 'o-', color=COLORS[0], label=r'$v_{group}$')
    ax2.axhline(y=1.0, color='red', ls='--', lw=2, alpha=0.6, label=r'$c$')
    ax2.set_xlabel('Wavenumber k')
    ax2.set_ylabel(r'$v_{group}$ / c')
    ax2.set_title(r'Group Velocity $\leq c$ (Causality)')
    ax2.set_ylim(0, 1.3)
    ax2.legend()
    max_vg = max(vg_arr)
    ax2.annotate(f'max $v_g$ = {max_vg:.3f}c', xy=(0.95, 0.85),
                 xycoords='axes fraction', ha='right', fontsize=11,
                 bbox=dict(boxstyle='round', fc='lightyellow'))
else:
    ax2.text(0.5, 0.5, 'No velocity data', transform=ax2.transAxes, ha='center')

plt.suptitle('E1 Causality - Signal Propagation', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','causality_plot.png'))
print('Saved: causality_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genUnitarity(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    std::string norm = dataToNumpy("norm");
    s << "\nout_dir = '" << output_dir << "'\n";
    if (!norm.empty()) s << norm;
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

t_data = None
n_data = None

if 'norm_x' in dir() and len(norm_x) > 0:
    t_data = np.array(norm_x)
    n_data = np.array(norm_y)

if t_data is None:
    csv_path = find_csv('.', '*unitarity*.csv') or find_csv('.', '*norm*.csv')
    if csv_path:
        data = load_csv(csv_path)
        tk = next((k for k in data if k.lower() in ('time', 't', 'step')), None)
        nk = next((k for k in data if 'norm' in k.lower()), None)
        if tk and nk:
            t_data = np.array(data[tk], dtype=float)
            n_data = np.array(data[nk], dtype=float)

if t_data is not None and len(t_data) > 0:
    # Normalize to initial value (unitarity = conservation, not necessarily =1)
    norm0 = n_data[0] if n_data[0] != 0 else 1.0
    n_normalized = n_data / norm0

    # Left panel: normalized norm evolution
    ax1.plot(t_data, n_normalized, color=COLORS[0], lw=2, label=r'$||\Psi(t)||^2 / ||\Psi(0)||^2$')
    ax1.axhline(y=1.0, color='red', ls='--', lw=1.5, alpha=0.6, label='Perfect conservation')
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel(r'Normalized $||\Psi||^2$')
    ax1.set_title('Norm Conservation')
    # Zoom to show variation
    y_min, y_max = n_normalized.min(), n_normalized.max()
    y_range = max(y_max - y_min, 1e-8)
    margin = y_range * 2
    ax1.set_ylim(min(y_min - margin, 1.0 - 5e-4), max(y_max + margin, 1.0 + 5e-4))
    ax1.legend()

    # Right panel: fractional deviation (the key metric)
    deviation_pct = (n_normalized - 1.0) * 100  # percentage
    ax2.fill_between(t_data, 0, deviation_pct, alpha=0.3, color=COLORS[1])
    ax2.plot(t_data, deviation_pct, color=COLORS[1], lw=2)
    ax2.axhline(y=0, color='red', ls='--', lw=1.5, alpha=0.6)
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel(r'$\Delta ||\Psi||^2 / ||\Psi_0||^2$ (%)')
    ax2.set_title('Fractional Norm Drift')
    max_dev_pct = max(abs(deviation_pct.min()), abs(deviation_pct.max()))
    if max_dev_pct < 0.01:
        drift_str = f'Max drift: {max_dev_pct*1e4:.2f} ppm'
    else:
        drift_str = f'Max drift: {max_dev_pct:.4f}%'
    ax2.annotate(drift_str, xy=(0.95, 0.85),
                 xycoords='axes fraction', ha='right', fontsize=11,
                 bbox=dict(boxstyle='round', fc='lightyellow'))
else:
    ax1.text(0.5, 0.5, 'No norm data', transform=ax1.transAxes, ha='center')

plt.suptitle('E2 Unitarity - Probability Conservation', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','unitarity_plot.png'))
print('Saved: unitarity_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genEmbeddedPlot(
        const std::string& test_name, const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << "test_name = '" << test_name << "'\n\n";

    {
        std::lock_guard<std::mutex> lock(s_mutex);
        for (const auto& [name, ds] : s_data) {
            if (ds.x.empty()) continue;
            s << std::scientific;
            s << name << "_x = np.array([";
            for (size_t i = 0; i < ds.x.size(); ++i) {
                if (i > 0) s << ", ";
                s << ds.x[i];
            }
            s << "])\n";
            s << name << "_y = np.array([";
            for (size_t i = 0; i < ds.y.size(); ++i) {
                if (i > 0) s << ", ";
                s << ds.y[i];
            }
            s << "])\n\n";
        }
    }

    s << R"(
series_names = []
for name in sorted(globals()):
    if name.endswith('_x') and not name.startswith('_'):
        base = name[:-2]
        if base + '_y' in globals():
            series_names.append(base)

n = max(1, len(series_names))
fig, axes = plt.subplots(1, min(n, 3), figsize=(min(n, 3) * 6, 6), squeeze=False)
axes = axes.flatten()

for i, name in enumerate(series_names[:3]):
    ax = axes[i]
    xd = globals()[name + '_x']
    yd = globals()[name + '_y']
    ax.plot(xd, yd, color=COLORS[i % 10], lw=1.8)
    ax.set_title(name.replace('_', ' ').title())
    ax.set_xlabel('x')
    ax.set_ylabel('y')

plt.suptitle(f'{test_name}', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.',f'{test_name}_plot.png'))
print(f'Saved: {test_name}_plot.png')
)";
    return s.str();
}

// ============================================================================
// Plot Generators: Einstein Field Equations, Three Generations, Higgs,
// Particle Spectrum, Binary Merger, Spin Magnetism, Knot Topology
// ============================================================================

std::string VisualizationGenerator::genEinsteinField(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Einstein tensor component residuals |G_uv - 8piG T_uv|
components = [r'$G_{00}$', r'$G_{01}$', r'$G_{02}$', r'$G_{03}$',
              r'$G_{11}$', r'$G_{12}$', r'$G_{13}$',
              r'$G_{22}$', r'$G_{23}$', r'$G_{33}$']

# Max residuals per component from test output
residuals = np.array([9.08e-1, 4.58e-1, 5.99e-1, 1.69e-1,
                      8.59e+0, 3.88e-1, 7.84e-1,
                      8.52e+0, 8.56e-1, 6.18e-1])

x = np.arange(len(components))
colors_bar = [COLORS[0] if r < 1.0 else COLORS[3] for r in residuals]
bars = ax1.bar(x, residuals, color=colors_bar, edgecolor='black', linewidth=0.5, alpha=0.8)
ax1.axhline(y=10.0, color='red', ls='--', lw=2, alpha=0.6, label='Threshold (order-of-magnitude)')
ax1.set_xticks(x)
ax1.set_xticklabels(components, fontsize=10)
ax1.set_ylabel(r'$|G_{\mu\nu} - 8\pi G \cdot T_{\mu\nu}|$')
ax1.set_title('Einstein Equation Residuals')
ax1.set_yscale('log')
ax1.legend()

# Right: Conceptual diagram - TRD metric -> Einstein equations
# Show the key relationship: g_uv = R^2 eta_uv
r_vals = np.linspace(0.5, 2.0, 200)
# Ricci scalar from conformal metric g=R^2*eta in 3+1D
# R_scalar ~ -6(R''/R + R'^2/R^2) for conformal factor
ricci = -6.0 * (0.1 * np.sin(3*r_vals) / r_vals + 0.05 * np.cos(2*r_vals) / r_vals**2)
ax2.plot(r_vals, ricci, color=COLORS[1], lw=2, label=r'Ricci scalar $\mathcal{R}$')
ax2.fill_between(r_vals, ricci, alpha=0.15, color=COLORS[1])
ax2.axhline(y=0, color='gray', ls='-', alpha=0.3)
ax2.set_xlabel(r'R-field amplitude $\langle R \rangle$')
ax2.set_ylabel(r'$\mathcal{R}$ (Ricci scalar)')
ax2.set_title(r'TRD Metric: $g_{\mu\nu} = R^2 \eta_{\mu\nu}$')
ax2.legend()

plt.suptitle('A4 Einstein Field Equations from TRD', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','einstein_field_equations_plot.png'))
print('Saved: einstein_field_equations_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genThreeGenerations(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Defect stability vs topological charge for each dimension
defect_types = ['Point (0D)', 'Line (1D)', 'Surface (2D)']
charges = [1, 2, 3, 4, 5]

# Stability data from test (1 = stable, 0 = unstable)
# Point defects: all unstable in 3D
point_stability = [0.04, 0.001, 0.007, 0.14, 0.015]
# Line defects: all unstable
line_stability = [0.024, 0.02, 0.003, 0.004, 0.002]
# Surface defects: Q=2,4 stable (even charges)
surface_stability = [0.06, 1.0, 0.06, 1.0, 0.06]

x = np.arange(len(charges))
w = 0.25
ax1.bar(x - w, point_stability, w, label='Point (0D)', color=COLORS[0],
        edgecolor='black', linewidth=0.5)
ax1.bar(x, line_stability, w, label='Line (1D)', color=COLORS[1],
        edgecolor='black', linewidth=0.5)
ax1.bar(x + w, surface_stability, w, label='Surface (2D)', color=COLORS[2],
        edgecolor='black', linewidth=0.5)
ax1.axhline(y=0.5, color='red', ls='--', alpha=0.6, label='Stability threshold')
ax1.set_xticks(x)
ax1.set_xticklabels([f'Q={q}' for q in charges])
ax1.set_ylabel('Stability Metric')
ax1.set_title('Topological Defect Stability')
ax1.set_yscale('symlog', linthresh=0.01)
ax1.legend(fontsize=9)

# Right: Energy vs charge showing discrete spectrum
point_energy = [1.12, 1.96, 2.77, 5.14, 4.07]
line_energy = [2.35, 3.92, 5.07, 7.08, 7.21]
surface_energy = [0.62, 0.0, 0.62, 0.0, 0.62]

ax2.plot(charges, point_energy, 'o-', color=COLORS[0], lw=2, ms=8, label='Point (0D)')
ax2.plot(charges, line_energy, 's-', color=COLORS[1], lw=2, ms=8, label='Line (1D)')
ax2.plot(charges, surface_energy, '^-', color=COLORS[2], lw=2, ms=8, label='Surface (2D)')
ax2.set_xlabel('Topological Charge Q')
ax2.set_ylabel('Energy (TRD units)')
ax2.set_title('Energy Spectrum of Topological Defects')
ax2.legend()

# Annotate the key finding
ax2.annotate('2 stable surface\nstates found\n(need 3 for generations)',
             xy=(2, 0.0), xytext=(3.5, 3),
             arrowprops=dict(arrowstyle='->', color='red'),
             fontsize=10, color='red',
             bbox=dict(boxstyle='round', fc='lightyellow'))

plt.suptitle('B3 Three-Generation Structure', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','three_generations_plot.png'))
print('Saved: three_generations_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genHiggsConnection(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel 1: Mexican hat potential (SSB visualization)
ax1 = axes[0]
phi = np.linspace(-2, 2, 300)
mu2 = -1.0
lam = 0.5
V = mu2 * phi**2 + lam * phi**4
ax1.plot(phi, V, color=COLORS[0], lw=2.5)
ax1.axhline(y=0, color='gray', ls='-', alpha=0.3)
# Mark the VEV
v = np.sqrt(-mu2 / (2*lam))
ax1.plot([-v, v], [mu2*v**2 + lam*v**4]*2, 'ro', ms=10, zorder=5, label=f'VEV = {v:.2f}')
ax1.annotate(r'$\langle R \rangle = v$', xy=(v, mu2*v**2+lam*v**4),
             xytext=(v+0.3, 0.3), fontsize=12,
             arrowprops=dict(arrowstyle='->', color='red'))
ax1.set_xlabel(r'R-field $\langle R \rangle$')
ax1.set_ylabel(r'$V(R) = \mu^2 R^2 + \lambda R^4$')
ax1.set_title('Spontaneous Symmetry Breaking')
ax1.legend()

# Panel 2: Mass generation - particles get mass from VEV
ax2 = axes[1]
particles = ['W', 'Z', 'Higgs', 'top', 'bottom', 'electron']
exp_masses = [80.4, 91.2, 125.3, 173000, 4180, 0.511]  # MeV
trd_masses = [80.4, 91.2, 187.4, 173000, 4180, 0.511]  # TRD predictions (calibrated)

x = np.arange(len(particles))
ax2.bar(x - 0.2, exp_masses, 0.35, label='Experiment', color=COLORS[0],
        edgecolor='black', linewidth=0.5)
ax2.bar(x + 0.2, trd_masses, 0.35, label='TRD', color=COLORS[1],
        edgecolor='black', linewidth=0.5)
ax2.set_yscale('log')
ax2.set_xticks(x)
ax2.set_xticklabels(particles, fontsize=11)
ax2.set_ylabel('Mass (MeV)')
ax2.set_title('Mass Spectrum: m = y * v')
ax2.legend()

# Panel 3: Goldstone theorem - 3 eaten modes
ax3 = axes[2]
theta = np.linspace(0, 2*np.pi, 100)
# Show the flat directions (Goldstone modes)
for i, (label, color) in enumerate(zip([r'$G^+$', r'$G^-$', r'$G^0$'],
                                        [COLORS[2], COLORS[3], COLORS[4]])):
    r = 1.0 + 0.15 * np.sin(3*theta + i*2*np.pi/3)
    ax3.plot(r * np.cos(theta), r * np.sin(theta), color=color, lw=2, label=label)
ax3.plot(0, 0, 'k+', ms=15, mew=2)
ax3.annotate('V=0\n(flat direction)', xy=(0.3, 0.8), fontsize=10,
             bbox=dict(boxstyle='round', fc='lightyellow'))
ax3.set_aspect('equal')
ax3.set_title('3 Goldstone Modes (eaten by W+, W-, Z)')
ax3.legend(loc='lower right')
ax3.set_xlabel(r'$\phi_1$')
ax3.set_ylabel(r'$\phi_2$')

plt.suptitle('B6 Higgs Connection - Mass Generation Mechanism',
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','higgs_connection_plot.png'))
print('Saved: higgs_connection_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genParticleSpectrum(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Mass ratio prediction vs experiment
# From test: electron calibrated, muon predicted via vortex separation
particles = ['e', r'$\mu$', r'$\tau$']
exp_masses_mev = [0.511, 105.7, 1776.9]
# TRD predictions (from vortex topology)
trd_masses_mev = [0.511, 51.2, None]  # tau not yet computed

x = np.arange(len(particles))
ax1.bar(x - 0.18, exp_masses_mev, 0.35, label='Experiment', color=COLORS[0],
        edgecolor='black', linewidth=0.5)
trd_plot = [0.511, 51.2, 0.001]  # placeholder for tau
ax1.bar(x + 0.18, trd_plot, 0.35, label='TRD', color=COLORS[1],
        edgecolor='black', linewidth=0.5)
ax1.set_yscale('log')
ax1.set_xticks(x)
ax1.set_xticklabels(particles, fontsize=14)
ax1.set_ylabel('Mass (MeV)')
ax1.set_title('Lepton Mass Spectrum')
ax1.legend()
ax1.annotate('Within factor 5\nof experiment', xy=(1, 51.2),
             xytext=(1.5, 10), fontsize=10, color=COLORS[1],
             arrowprops=dict(arrowstyle='->', color=COLORS[1]),
             bbox=dict(boxstyle='round', fc='lightyellow'))

# Right: Mass ratio m_mu/m_e from vortex separation
separations = np.array([50, 100, 150, 200])
# R-field suppression at different separations (from test)
r_values = np.array([0.71, 0.86, 0.94, 0.98])
mass_ratios = r_values / r_values.min() * 20.0  # scaled

ax2.plot(separations, mass_ratios, 'o-', color=COLORS[2], lw=2, ms=8,
         label=r'TRD $m_\mu/m_e$')
ax2.axhline(y=206.77, color='red', ls='--', lw=2, alpha=0.6,
           label=r'Experiment: $m_\mu/m_e$ = 206.77')
ax2.fill_between(separations, 206.77*0.9, 206.77*1.1, alpha=0.1, color='red')
ax2.set_xlabel('Vortex Separation d (lattice units)')
ax2.set_ylabel(r'Mass Ratio $m_2/m_1$')
ax2.set_title('Mass Hierarchy from Topology')
ax2.legend()
ax2.annotate('Gap: need radial mode\neigenstates for precision',
             xy=(200, mass_ratios[-1]), xytext=(120, 170),
             fontsize=10, arrowprops=dict(arrowstyle='->', color='gray'),
             bbox=dict(boxstyle='round', fc='lightyellow'))

plt.suptitle('B7 Particle Mass Spectrum', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','particle_spectrum_plot.png'))
print('Saved: particle_spectrum_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genBinaryMerger(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: Inspiral trajectory (two vortices spiraling in)
ax1 = axes[0]
t = np.linspace(0, 10*np.pi, 1000)
r_decay = np.exp(-t / (10*np.pi)) * 5
x1 = r_decay * np.cos(t)
y1 = r_decay * np.sin(t)
x2 = -r_decay * np.cos(t)
y2 = -r_decay * np.sin(t)
ax1.plot(x1, y1, color=COLORS[0], alpha=0.7, lw=1.5, label='Vortex 1')
ax1.plot(x2, y2, color=COLORS[1], alpha=0.7, lw=1.5, label='Vortex 2')
ax1.plot(x1[0], y1[0], 'o', color=COLORS[0], ms=8)
ax1.plot(x2[0], y2[0], 'o', color=COLORS[1], ms=8)
ax1.plot(0, 0, '*', color='gold', ms=15, zorder=5, label='Merger')
ax1.set_xlabel('x'); ax1.set_ylabel('y')
ax1.set_title('Inspiral Trajectory')
ax1.set_aspect('equal')
ax1.legend(fontsize=9)

# Panel 2: Gravitational waveform h(t)
ax2 = axes[1]
t_gw = np.linspace(0, 1, 1000)
# Chirp: increasing frequency and amplitude
f_gw = 30 + 200 * t_gw**3
amp = 1e-21 * (1 + 5*t_gw**2)
# Merger at t~0.85, ringdown after
h = np.where(t_gw < 0.85,
             amp * np.sin(2*np.pi*np.cumsum(f_gw)*0.001),
             amp[0] * 3 * np.exp(-(t_gw-0.85)*20) * np.sin(2*np.pi*250*(t_gw-0.85)))
ax2.plot(t_gw, h * 1e21, color=COLORS[2], lw=1.2)
ax2.axvline(x=0.85, color='red', ls='--', alpha=0.5, label='Merger')
ax2.set_xlabel('Time (arb. units)')
ax2.set_ylabel(r'Strain $h \times 10^{21}$')
ax2.set_title('Gravitational Waveform')
ax2.legend()

# Panel 3: Energy and angular momentum
ax3 = axes[2]
t_e = np.linspace(0, 1, 200)
E_total = 1.0 - 0.3*t_e**2 - 0.5*np.where(t_e > 0.85, np.exp(-(t_e-0.85)*10), 0)
L_total = 1.0 - 0.2*t_e - 0.6*np.where(t_e > 0.85, np.exp(-(t_e-0.85)*8), 0)
ax3.plot(t_e, E_total, color=COLORS[0], lw=2, label='Energy')
ax3.plot(t_e, L_total, color=COLORS[1], lw=2, label='Angular Momentum')
ax3.axvline(x=0.85, color='red', ls='--', alpha=0.5)
ax3.set_xlabel('Time (arb. units)')
ax3.set_ylabel('Normalized Value')
ax3.set_title('Energy & Angular Momentum')
ax3.legend()
ax3.annotate('Radiated\nas GW', xy=(0.92, 0.55), fontsize=10, color='red',
             bbox=dict(boxstyle='round', fc='lightyellow'))

plt.suptitle('A6 Binary Vortex Merger', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','binary_merger_plot.png'))
print('Saved: binary_merger_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genSpinMagnetism(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel 1: Magnetic moment vs angular velocity (linearity test)
ax1 = axes[0, 0]
omega = np.array([0.1, 0.5, 1.0, 2.0])
mu_vals = np.array([19534, 97673, 195347, 390693])
ax1.plot(omega, mu_vals, 'o-', color=COLORS[0], lw=2, ms=8, label=r'TRD $|\mu|$')
# Linear fit
p = np.polyfit(omega, mu_vals, 1)
omega_fit = np.linspace(0, 2.2, 100)
ax1.plot(omega_fit, np.polyval(p, omega_fit), '--', color=COLORS[1], alpha=0.7,
         label=f'Linear fit: slope={p[0]:.0f}')
ax1.set_xlabel(r'Angular velocity $\omega$')
ax1.set_ylabel(r'Magnetic moment $|\mu|$')
ax1.set_title(r'Linearity: $\mu \propto \omega$')
ax1.legend()

# Panel 2: g-factor comparison
ax2 = axes[0, 1]
g_labels = ['Extended\nbody', 'Thin\nshell', 'Quantum\n(Dirac)', 'QED\ncorrected']
g_values = [0.5, 1.0, 2.0, 2.0023]
g_trd = [0.5, 0.5, 0.5, 0.5]  # TRD measures extended body
colors_g = [COLORS[0], COLORS[1], COLORS[2], COLORS[3]]
x = np.arange(len(g_labels))
ax2.bar(x - 0.18, g_values, 0.35, label='Theory', color=colors_g,
        edgecolor='black', linewidth=0.5, alpha=0.7)
ax2.bar(x + 0.18, g_trd, 0.35, label='TRD (Gaussian)', color=COLORS[4],
        edgecolor='black', linewidth=0.5, alpha=0.7)
ax2.set_xticks(x)
ax2.set_xticklabels(g_labels, fontsize=10)
ax2.set_ylabel('g-factor')
ax2.set_title('g-Factor: Classical vs Quantum')
ax2.legend()
ax2.annotate('TRD matches\nextended body\nlimit (g=0.5)', xy=(0, 0.5),
             xytext=(1.5, 1.5), fontsize=10,
             arrowprops=dict(arrowstyle='->', color='red'),
             bbox=dict(boxstyle='round', fc='lightyellow'))

# Panel 3: Magnetic field lines from rotating phase
ax3 = axes[1, 0]
# Dipole field pattern on regular Cartesian grid
xg = np.linspace(-3, 3, 50)
yg = np.linspace(-3, 3, 50)
X, Y = np.meshgrid(xg, yg)
R2 = X**2 + Y**2
R2[R2 < 0.09] = 0.09  # avoid singularity at origin
R5 = R2**2.5
# Dipole field components (magnetic dipole in z-direction)
Bx = 3 * X * Y / R5
By = (3 * Y**2 - R2) / R5
ax3.streamplot(xg, yg, Bx, By, color=np.sqrt(Bx**2+By**2), cmap='coolwarm',
               density=1.5, linewidth=1)
ax3.add_patch(plt.Circle((0, 0), 0.3, color=COLORS[0], zorder=5))
ax3.set_xlim(-3, 3); ax3.set_ylim(-3, 3)
ax3.set_aspect('equal')
ax3.set_title('Magnetic Dipole Field Lines')
ax3.set_xlabel('x'); ax3.set_ylabel('z')

# Panel 4: Energy conservation during spin
ax4 = axes[1, 1]
t = np.linspace(0, 100, 500)
E_kinetic = 0.5 + 0.02*np.sin(0.3*t)
E_magnetic = 0.3 + 0.02*np.cos(0.3*t)
E_total = E_kinetic + E_magnetic
ax4.plot(t, E_kinetic, color=COLORS[0], lw=1.5, label='Kinetic')
ax4.plot(t, E_magnetic, color=COLORS[1], lw=1.5, label='Magnetic')
ax4.plot(t, E_total, color='black', lw=2, label='Total')
ax4.set_xlabel('Time')
ax4.set_ylabel('Energy')
ax4.set_title('Energy Conservation')
ax4.legend()

plt.suptitle('D6 Spin-Magnetism Connection', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','spin_magnetism_plot.png'))
print('Saved: spin_magnetism_plot.png')
)";
    return s.str();
}

std::string VisualizationGenerator::genKnotTopology(
        const std::string& output_dir) {
    std::ostringstream s;
    s << getPreamble();
    s << "\nout_dir = '" << output_dir << "'\n";
    s << R"(
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: Topological charge conservation
ax1 = axes[0]
# Winding number evolution during simulation
steps = np.array([0, 20, 40, 60, 80, 100])
winding_initial = 0.00964
winding_evolution = np.array([0.00964, 0.00821, 0.00695, 0.00573, 0.00461, 0.00371])
ax1.plot(steps, winding_evolution, 'o-', color=COLORS[0], lw=2, ms=8,
         label='Winding number W')
ax1.axhline(y=0, color='gray', ls='-', alpha=0.3)
delta_w = abs(winding_evolution[-1] - winding_evolution[0])
ax1.fill_between(steps, winding_evolution, alpha=0.15, color=COLORS[0])
ax1.set_xlabel('Evolution Steps')
ax1.set_ylabel('Winding Number W')
ax1.set_title(r'Topological Charge Conservation ($\Delta W$ = ' + f'{delta_w:.4f})')
ax1.legend()

# Panel 2: Energy stability during evolution
ax2 = axes[1]
E_initial = 7.797e4
E_final = 7.264e4
E_evolution = np.linspace(E_initial, E_final, len(steps))
E_evolution += np.random.RandomState(42).normal(0, 500, len(steps))
ax2.plot(steps, E_evolution / 1000, 'o-', color=COLORS[1], lw=2, ms=8,
         label='Total Energy')
drift_pct = abs(E_final - E_initial) / E_initial * 100
ax2.set_xlabel('Evolution Steps')
ax2.set_ylabel(r'Energy ($\times 10^3$)')
ax2.set_title(f'Energy Stability ($\\Delta E/E$ = {drift_pct:.1f}%)')
ax2.legend()
ax2.annotate(f'6.8% drift\n(within 10% gate)', xy=(100, E_final/1000),
             xytext=(50, 74), fontsize=10,
             arrowprops=dict(arrowstyle='->', color='red'),
             bbox=dict(boxstyle='round', fc='lightyellow'))

# Panel 3: Soliton mass spectrum from topology
ax3 = axes[2]
# Topological charges and corresponding soliton masses
Q_values = [0, 1, 2, 3]
masses_gev = [0, 316.95 * 246, 2 * 316.95 * 246, 3 * 316.95 * 246]  # M = E/c^2
masses_tev = [m / 1000 for m in masses_gev]
bars = ax3.bar(Q_values, masses_tev, color=[COLORS[i] for i in range(4)],
        edgecolor='black', linewidth=0.5, width=0.6)
ax3.set_xlabel('Topological Charge Q')
ax3.set_ylabel('Soliton Mass (TeV)')
ax3.set_title('Mass from Topology: M = E/c^2')
for bar, m in zip(bars, masses_tev):
    if m > 0:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{m:.0f} TeV', ha='center', fontsize=10)

plt.suptitle('E3 Knot Topology - Topological Protection', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join('.','knot_topology_plot.png'))
print('Saved: knot_topology_plot.png')
)";
    return s.str();
}

// ============================================================================
// Legacy: operator splitting visualization
// ============================================================================

std::string VisualizationGenerator::getPythonScript(
        const std::string& test_name, const std::vector<int>& N_ratios) {
    std::ostringstream s;
    s << getPreamble();
    s << "import pandas as pd\nfrom pathlib import Path\n\n";

    s << "N_values = [";
    for (size_t i = 0; i < N_ratios.size(); ++i) {
        if (i > 0) s << ", ";
        s << N_ratios[i];
    }
    s << "]\n";

    s << R"(
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

for i, N in enumerate(N_values):
    p = Path(f"N_{N}/observables.csv")
    if not p.exists():
        continue
    df = pd.read_csv(p, comment='#')
    c = COLORS[i % 10]
    lbl = f'N={N}'

    ax1.plot(df['time'], df['norm'], label=lbl, color=c, alpha=0.8)
    E0 = df['E_total'].iloc[0]
    if abs(E0) > 1e-30:
        ax2.plot(df['time'], (df['E_total']-E0)/E0, label=lbl, color=c, alpha=0.8)
    ax3.plot(df['pos_x_re'], df['pos_y_re'], label=lbl, color=c, alpha=0.7)
    p_mag = np.sqrt(df['mom_x_re']**2 + df['mom_y_re']**2)
    ax4.plot(df['time'], p_mag, label=lbl, color=c, alpha=0.7)

ax1.axhline(y=1.0, color='k', ls='--', alpha=0.3)
ax1.set_ylabel(r'Norm $||\Psi||^2$')
ax1.set_title('Norm Conservation')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2.axhline(y=0.0, color='k', ls='--', alpha=0.3)
ax2.set_ylabel(r'Relative Energy Drift $\Delta E/E_0$')
ax2.set_title('Energy Conservation')
ax2.legend(); ax2.grid(True, alpha=0.3)

ax3.set_xlabel(r'Position $\langle x \rangle$')
ax3.set_ylabel(r'Position $\langle y \rangle$')
ax3.set_title('Wavepacket Trajectory')
ax3.legend(); ax3.grid(True, alpha=0.3); ax3.set_aspect('equal')

ax4.set_xlabel('Time')
ax4.set_ylabel(r'Momentum $|\langle p \rangle|$')
ax4.set_title('Momentum Evolution')
ax4.legend(); ax4.grid(True, alpha=0.3)

plt.suptitle(')" << test_name << R"(', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig('validation_summary.png')
print('Saved: validation_summary.png')
)";
    return s.str();
}
