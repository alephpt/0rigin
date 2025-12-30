#!/usr/bin/env python3
"""
Renormalization Group Flow Extraction from SMFT Simulations
Direction 2: Extract β(Δ) from multi-scale simulations and predict IR behavior

Physical Model:
- RG equation: μ dΔ/dμ = β(Δ) = -b₀Δ³ - b₁Δ⁵ - b₂Δ⁷ + ...
- Extract β coefficients from Δ_eff(grid_size) measurements
- Integrate flow to predict Δ(E → 0)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path
from scipy.integrate import odeint
from scipy.optimize import curve_fit

# Physical constants
M_PLANCK = 1.22e19  # GeV
HBAR_C = 0.1973  # GeV·fm


def extract_effective_delta(R_squared, Delta_bare):
    """
    Extract effective Δ from measured ⟨R²⟩

    Vacuum energy: ρ = (1/2) Δ_eff² ⟨R²⟩
    If we assume ρ is scale-invariant (fixed point), then:
    Δ_eff = Δ_bare / sqrt(⟨R²⟩)

    Args:
        R_squared: Measured ⟨R²⟩ from simulation
        Delta_bare: Bare mass gap (Planck mass)

    Returns:
        Delta_eff: Effective mass gap at this scale
    """
    if R_squared <= 0:
        return np.nan

    # Simple extraction (assumes vacuum energy fixed)
    Delta_eff = Delta_bare / np.sqrt(R_squared)

    return Delta_eff


def load_multiscale_data(output_base):
    """
    Load simulation data from multiple grid sizes

    Each grid size → different UV cutoff Λ = Δ·N/L
    where L = domain size, N = grid points

    Args:
        output_base: Base directory containing simulation outputs

    Returns:
        energy_scales: Array of energy scales (UV cutoffs)
        Delta_eff_values: Extracted Δ_eff at each scale
        R_squared_values: Measured ⟨R²⟩ at each scale
    """
    # Find all output directories with grid size info
    test_patterns = [
        "*relativistic_mass*",
        "*phase_2.3*",
        "*defect_localization*",
        "*scenario_2.4*",
        "*vacuum_energy*"
    ]

    all_data = []

    for pattern in test_patterns:
        dirs = glob.glob(f"{output_base}/{pattern}")
        for d in dirs:
            # Extract grid size from directory name
            basename = Path(d).name
            parts = basename.split('_')

            grid_size = None
            for part in parts:
                if 'x' in part and part.replace('x', '').isdigit():
                    grid_size = int(part.split('x')[0])
                    break

            if grid_size is None:
                continue

            # Load observables
            obs_files = glob.glob(f"{d}/**/observables.csv", recursive=True)

            for obs_file in obs_files:
                try:
                    df = pd.read_csv(obs_file)

                    # Extract N value from path
                    path_parts = obs_file.split('/')
                    N_value = None
                    for part in path_parts:
                        if part.startswith('N_'):
                            N_value = int(part.split('_')[1])
                            break

                    if N_value is None:
                        continue

                    # Late-time average of R field (approximated from norm)
                    n_points = len(df)
                    late_start = int(0.8 * n_points)
                    late_df = df.iloc[late_start:]

                    if 'norm' in df.columns:
                        R_mean = late_df['norm'].mean()
                        R_squared = (late_df['norm']**2).mean()
                    else:
                        continue

                    # Energy scale = UV cutoff
                    # Λ = π/(lattice spacing) = π·N/L
                    L_domain = 100.0  # Planck lengths
                    Lambda_UV = (np.pi * grid_size) / L_domain

                    # Effective Δ
                    Delta_eff = extract_effective_delta(R_squared, M_PLANCK)

                    all_data.append({
                        'grid_size': grid_size,
                        'N_substeps': N_value,
                        'energy_scale': Lambda_UV,
                        'R_squared': R_squared,
                        'Delta_eff': Delta_eff
                    })

                except Exception as e:
                    print(f"Could not load {obs_file}: {e}")

    if not all_data:
        return None, None, None

    # Convert to arrays, group by energy scale
    df_all = pd.DataFrame(all_data)
    grouped = df_all.groupby('energy_scale').agg({
        'R_squared': 'mean',
        'Delta_eff': 'mean'
    }).reset_index()

    energy_scales = grouped['energy_scale'].values
    R_squared_values = grouped['R_squared'].values
    Delta_eff_values = grouped['Delta_eff'].values

    return energy_scales, Delta_eff_values, R_squared_values


def fit_beta_function(energy_scales, Delta_eff_values, n_terms=3):
    """
    Fit RG beta function coefficients from data

    β(Δ) = -b₀Δ³ - b₁Δ⁵ - b₂Δ⁷ - ...

    From data:
    β(Δ) ≈ μ dΔ/dμ = d(ln Δ)/d(ln μ)

    Args:
        energy_scales: Array of UV cutoffs (proxy for μ)
        Delta_eff_values: Extracted Δ_eff at each scale
        n_terms: Number of beta function terms to fit

    Returns:
        beta_coeffs: [b₀, b₁, b₂, ...]
        fit_quality: R² goodness of fit
    """
    # Remove NaN values
    mask = ~(np.isnan(energy_scales) | np.isnan(Delta_eff_values))
    mu_clean = energy_scales[mask]
    Delta_clean = Delta_eff_values[mask]

    if len(mu_clean) < n_terms + 1:
        return None, None

    # Sort by scale
    sort_idx = np.argsort(mu_clean)
    mu_clean = mu_clean[sort_idx]
    Delta_clean = Delta_clean[sort_idx]

    # Compute β numerically: β ≈ μ dΔ/dμ
    log_mu = np.log(mu_clean)
    log_Delta = np.log(Delta_clean)

    # Numerical derivative
    beta_numerical = np.gradient(log_Delta, log_mu) * Delta_clean

    # Fit to polynomial in Δ
    # β(Δ) = -b₀Δ³ - b₁Δ⁵ - ...

    # Design matrix: [Δ³, Δ⁵, Δ⁷, ...]
    X = np.column_stack([Delta_clean**(2*i + 3) for i in range(n_terms)])
    y = -beta_numerical  # Negative sign convention

    # Least squares fit
    beta_coeffs, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)

    # R-squared
    y_pred = -X @ beta_coeffs
    ss_res = np.sum((beta_numerical - y_pred)**2)
    ss_tot = np.sum((beta_numerical - np.mean(beta_numerical))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    return beta_coeffs, r_squared


def integrate_rg_flow(Delta_UV, mu_UV, mu_IR, beta_coeffs):
    """
    Integrate RG flow equation from UV to IR

    dΔ/d(ln μ) = β(Δ)

    Args:
        Delta_UV: Initial value at UV scale
        mu_UV: UV energy scale
        mu_IR: IR energy scale (final)
        beta_coeffs: [b₀, b₁, b₂, ...]

    Returns:
        mu_values: Energy scales
        Delta_values: Δ(μ) along flow
    """
    def beta(Delta, log_mu):
        """Beta function"""
        result = 0
        for i, b in enumerate(beta_coeffs):
            result -= b * Delta**(2*i + 3)
        return result

    # Logarithmic scale grid
    log_mu_values = np.linspace(np.log(mu_UV), np.log(mu_IR), 1000)

    # Integrate ODE
    Delta_solution = odeint(beta, Delta_UV, log_mu_values)

    mu_values = np.exp(log_mu_values)
    Delta_values = Delta_solution.flatten()

    return mu_values, Delta_values


def theoretical_rg_prediction(n_oscillators=100):
    """
    Compute theoretical beta function coefficients

    From perturbative calculation:
    b₀ = N/(12π²)  where N = number of oscillators

    Args:
        n_oscillators: Number of Kuramoto oscillators per Planck volume

    Returns:
        beta_coeffs_theory: [b₀, 0, 0, ...]  (one-loop only)
    """
    b0_theory = n_oscillators / (12 * np.pi**2)
    return [b0_theory, 0, 0]


def plot_rg_analysis(energy_scales, Delta_eff, beta_coeffs, mu_flow, Delta_flow):
    """
    Create comprehensive RG flow analysis plots

    Args:
        energy_scales: Measured UV cutoffs
        Delta_eff: Extracted Δ_eff from data
        beta_coeffs: Fitted [b₀, b₁, b₂]
        mu_flow: Energy scales from RG integration
        Delta_flow: Δ(μ) from RG integration
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Plot 1: Δ_eff vs energy scale (data)
    ax = axes[0, 0]
    ax.loglog(energy_scales, Delta_eff, 'ro', markersize=8, label='Data')
    ax.axhline(y=M_PLANCK, color='k', linestyle=':', linewidth=2, label='M_Planck')

    # Theory: one-loop running
    b0_theory = theoretical_rg_prediction()[0]
    mu_theory = np.logspace(-3, np.log10(M_PLANCK), 100)
    Delta_theory = M_PLANCK / np.sqrt(1 + 2*b0_theory*M_PLANCK**2 * np.log(M_PLANCK/mu_theory))
    ax.loglog(mu_theory, Delta_theory, 'b--', linewidth=2, label=f'Theory (b₀={b0_theory:.3f})')

    ax.set_xlabel('Energy Scale μ (GeV)')
    ax.set_ylabel('Δ_eff(μ) (GeV)')
    ax.set_title('Effective Planck Mass Running (Data)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: RG flow trajectory
    ax = axes[0, 1]
    if mu_flow is not None and Delta_flow is not None:
        ax.loglog(mu_flow, Delta_flow, 'g-', linewidth=2, label='RG Flow (fitted β)')
        ax.loglog(energy_scales, Delta_eff, 'ro', markersize=6, label='Data')
        ax.axhline(y=M_PLANCK, color='k', linestyle=':', linewidth=1)

        # Mark key scales
        idx_TeV = np.argmin(np.abs(mu_flow - 1e3))
        idx_meV = np.argmin(np.abs(mu_flow - 1e-3))
        ax.plot(mu_flow[idx_TeV], Delta_flow[idx_TeV], 'bs', markersize=10, label='Δ(TeV)')
        ax.plot(mu_flow[idx_meV], Delta_flow[idx_meV], 'ms', markersize=10, label='Δ(meV)')

        ax.set_xlabel('Energy Scale μ (GeV)')
        ax.set_ylabel('Δ(μ) (GeV)')
        ax.set_title('RG Flow Integration')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 3: Beta function
    ax = axes[0, 2]
    if beta_coeffs is not None:
        Delta_range = np.logspace(-60, np.log10(M_PLANCK), 200)
        beta_values = np.zeros_like(Delta_range)

        for i, b in enumerate(beta_coeffs):
            beta_values -= b * Delta_range**(2*i + 3)

        ax.loglog(Delta_range, np.abs(beta_values), 'b-', linewidth=2, label='|β(Δ)|')
        ax.axvline(x=M_PLANCK, color='k', linestyle=':', label='M_Planck')

        # Indicate sign
        positive_mask = beta_values > 0
        negative_mask = beta_values < 0
        ax.fill_between(Delta_range[positive_mask], 1e-200, np.abs(beta_values[positive_mask]),
                        alpha=0.3, color='red', label='β > 0 (IR stable)')
        ax.fill_between(Delta_range[negative_mask], 1e-200, np.abs(beta_values[negative_mask]),
                        alpha=0.3, color='blue', label='β < 0 (UV free)')

        ax.set_xlabel('Δ (GeV)')
        ax.set_ylabel('|β(Δ)|')
        ax.set_title(f'Beta Function (b₀={beta_coeffs[0]:.3e})')
        ax.set_ylim([1e-200, 1e100])
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 4: Vacuum energy suppression
    ax = axes[1, 0]
    if mu_flow is not None and Delta_flow is not None:
        # ρ_vac(μ) ~ Δ(μ)²
        rho_vac = Delta_flow**2
        rho_observed = 1e-47  # GeV⁴

        ax.loglog(mu_flow, rho_vac, 'r-', linewidth=2, label='ρ_vac(μ) ~ Δ²(μ)')
        ax.axhline(y=M_PLANCK**2, color='k', linestyle=':', label='Planck scale')
        ax.axhline(y=rho_observed, color='g', linestyle='--', linewidth=2, label='Observed Λ')

        ax.fill_between([1e-3, 1e3], 1e-50, 1e100, alpha=0.2, color='yellow',
                       label='Low-energy regime')

        ax.set_xlabel('Energy Scale μ (GeV)')
        ax.set_ylabel('ρ_vac(μ) (GeV⁴)')
        ax.set_title('Vacuum Energy Density Evolution')
        ax.set_ylim([1e-50, 1e80])
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 5: Suppression factor
    ax = axes[1, 1]
    if mu_flow is not None and Delta_flow is not None:
        suppression = (Delta_flow / M_PLANCK)**2
        ax.loglog(mu_flow, suppression, 'b-', linewidth=2, label='(Δ/M_P)²')
        ax.axhline(y=1e-123, color='r', linestyle='--', linewidth=2, label='Required')

        # Mark achievement at key scales
        idx_meV = np.argmin(np.abs(mu_flow - 1e-3))
        supp_meV = suppression[idx_meV]
        ax.plot(mu_flow[idx_meV], supp_meV, 'ro', markersize=12,
               label=f'Achieved: {supp_meV:.2e}')

        ax.set_xlabel('Energy Scale μ (GeV)')
        ax.set_ylabel('Suppression Factor')
        ax.set_title('Vacuum Energy Suppression vs Scale')
        ax.set_ylim([1e-140, 1e10])
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Plot 6: Summary
    ax = axes[1, 2]
    ax.axis('off')

    summary_text = "RG Flow Analysis Summary\n"
    summary_text += "=" * 40 + "\n\n"

    if beta_coeffs is not None:
        summary_text += "Beta Function Coefficients:\n"
        for i, b in enumerate(beta_coeffs):
            summary_text += f"  b{i} = {b:.3e}\n"

        summary_text += "\nTheoretical One-Loop:\n"
        b0_th = theoretical_rg_prediction()[0]
        summary_text += f"  b₀^theory = {b0_th:.3e}\n"
        summary_text += f"  Ratio: {beta_coeffs[0]/b0_th:.2f}\n\n"

    if mu_flow is not None and Delta_flow is not None:
        idx_Planck = 0
        idx_meV = np.argmin(np.abs(mu_flow - 1e-3))

        summary_text += "Running Results:\n"
        summary_text += f"  Δ(M_P) = {Delta_flow[idx_Planck]:.3e} GeV\n"
        summary_text += f"  Δ(meV) = {Delta_flow[idx_meV]:.3e} GeV\n"
        summary_text += f"  Suppression = {(Delta_flow[idx_meV]/Delta_flow[idx_Planck])**2:.3e}\n\n"

        summary_text += "Required vs Achieved:\n"
        summary_text += f"  Required:  1.0e-123\n"
        summary_text += f"  Achieved:  {(Delta_flow[idx_meV]/M_PLANCK)**2:.3e}\n"
        summary_text += f"  Missing:   {(Delta_flow[idx_meV]/M_PLANCK)**2 / 1e-123:.3e}\n\n"

        # Verdict
        final_supp = (Delta_flow[idx_meV]/M_PLANCK)**2
        if final_supp < 1e-120:
            verdict = "✓ SUCCESS"
        elif final_supp < 1e-100:
            verdict = "⚠ PARTIAL"
        else:
            verdict = "✗ INSUFFICIENT"

        summary_text += f"Status: {verdict}"

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
           fontsize=9, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.tight_layout()
    return fig


def main():
    """Main RG flow extraction routine"""
    print("=" * 60)
    print("SMFT Renormalization Group Flow Extraction")
    print("Direction 2: Beta Function from Multi-Scale Simulations")
    print("=" * 60)

    # Load data from simulations
    output_base = "/home/persist/neotec/0rigin/output"

    print("\nLoading multi-scale simulation data...")
    energy_scales, Delta_eff, R_squared = load_multiscale_data(output_base)

    if energy_scales is None or len(energy_scales) < 3:
        print("\nInsufficient data for RG extraction.")
        print("Running theoretical prediction instead...\n")

        # Theoretical flow
        beta_coeffs_theory = theoretical_rg_prediction(n_oscillators=100)
        mu_UV = M_PLANCK
        mu_IR = 1e-3  # meV

        mu_flow, Delta_flow = integrate_rg_flow(M_PLANCK, mu_UV, mu_IR, beta_coeffs_theory)

        print("Theoretical Beta Function:")
        print(f"  b₀ = {beta_coeffs_theory[0]:.3e}")

        # Create plots
        fig = plot_rg_analysis(None, None, beta_coeffs_theory, mu_flow, Delta_flow)

    else:
        print(f"Found {len(energy_scales)} energy scales")
        print("\nEnergy Scale (GeV) | Δ_eff (GeV) | ⟨R²⟩")
        print("-" * 60)
        for i in range(len(energy_scales)):
            print(f"{energy_scales[i]:.3e} | {Delta_eff[i]:.3e} | {R_squared[i]:.3f}")

        # Fit beta function
        print("\nFitting beta function...")
        beta_coeffs, r_squared = fit_beta_function(energy_scales, Delta_eff, n_terms=3)

        if beta_coeffs is not None:
            print("\nBeta Function Coefficients:")
            for i, b in enumerate(beta_coeffs):
                print(f"  b{i} = {b:.3e}")
            print(f"  R² = {r_squared:.3f}")

            # Theoretical comparison
            b0_theory = theoretical_rg_prediction()[0]
            print(f"\nTheoretical b₀ = {b0_theory:.3e}")
            print(f"Fitted / Theory = {beta_coeffs[0]/b0_theory:.2f}")

            # Integrate RG flow
            print("\nIntegrating RG flow from Planck to meV...")
            mu_UV = max(energy_scales)
            Delta_UV = Delta_eff[np.argmax(energy_scales)]
            mu_IR = 1e-3  # meV

            mu_flow, Delta_flow = integrate_rg_flow(Delta_UV, mu_UV, mu_IR, beta_coeffs)

            # Results
            idx_meV = np.argmin(np.abs(mu_flow - 1e-3))
            Delta_at_meV = Delta_flow[idx_meV]
            suppression = (Delta_at_meV / M_PLANCK)**2

            print(f"\nResults:")
            print(f"  Δ(M_Planck) = {Delta_UV:.3e} GeV")
            print(f"  Δ(meV)      = {Delta_at_meV:.3e} GeV")
            print(f"  Suppression = {suppression:.3e}")
            print(f"  Required    = 1.0e-123")
            print(f"  Missing     = {suppression/1e-123:.3e}")

            # Create plots
            fig = plot_rg_analysis(energy_scales, Delta_eff, beta_coeffs, mu_flow, Delta_flow)

        else:
            print("Beta function fit failed")
            fig = None

    # Save plot
    if fig is not None:
        output_path = "/home/persist/neotec/0rigin/output/rg_flow_analysis.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {output_path}")

    # Final assessment
    print("\n" + "=" * 60)
    print("ASSESSMENT")
    print("=" * 60)

    print("\nKey Findings:")
    print("1. RG running provides logarithmic suppression")
    print("2. One-loop β insufficient for 10^-123")
    print("3. Need anomalous dimension OR asymptotic safety")
    print("4. Multi-loop corrections could enhance running")

    print("\nRecommendations:")
    print("- Extract higher-order beta coefficients from finer grid scans")
    print("- Test for fixed point (asymptotic safety)")
    print("- Combine with Direction 1 (quantum corrections)")
    print("- Consider Direction 3 (SMFT+GR coupling)")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
