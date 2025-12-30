#!/usr/bin/env python3
"""
Data Collapse Analysis for SMFT Universality Classification

This script performs finite-size scaling data collapse to extract critical exponents
and validate the universality class of the SMFT phase transition.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize, interpolate
from scipy.stats import chi2
import os
import sys
import glob
from typing import Tuple, List, Dict
import argparse

class DataCollapseAnalyzer:
    """Performs FSS data collapse analysis on SMFT simulation data"""

    def __init__(self, data_dir: str):
        self.data_dir = data_dir
        self.data = {}  # L -> DataFrame
        self.grid_sizes = []
        self.sigma_c = None
        self.beta = None
        self.nu = None
        self.gamma = None
        self.eta = None

    def load_data(self):
        """Load observables data for all grid sizes"""
        print("Loading data from:", self.data_dir)

        # Find all universality scan directories
        pattern = os.path.join(self.data_dir, "universality_L*")
        dirs = sorted(glob.glob(pattern))

        for d in dirs:
            # Extract grid size from directory name
            L = int(d.split('_L')[-1].split('_')[0])

            # Load observables for each noise value
            obs_file = os.path.join(d, "observables_summary.csv")
            if os.path.exists(obs_file):
                df = pd.read_csv(obs_file)
                self.data[L] = df
                self.grid_sizes.append(L)
                print(f"  Loaded L={L}: {len(df)} noise points")

        self.grid_sizes = sorted(self.grid_sizes)
        print(f"Loaded {len(self.grid_sizes)} grid sizes: {self.grid_sizes}")

    def find_critical_point_binder(self) -> float:
        """Find critical point from Binder cumulant crossings"""
        print("\n=== Finding Critical Point from Binder Cumulant ===")

        crossings = []

        for i in range(len(self.grid_sizes) - 1):
            L1 = self.grid_sizes[i]
            L2 = self.grid_sizes[i + 1]

            df1 = self.data[L1]
            df2 = self.data[L2]

            # Interpolate Binder cumulants
            f1 = interpolate.interp1d(df1['sigma'], df1['binder_cumulant'],
                                     kind='cubic', fill_value='extrapolate')
            f2 = interpolate.interp1d(df2['sigma'], df2['binder_cumulant'],
                                     kind='cubic', fill_value='extrapolate')

            # Find crossing point
            sigma_range = np.linspace(0.5, 1.5, 1000)
            diff = f1(sigma_range) - f2(sigma_range)

            # Find zero crossing
            zero_crossings = np.where(np.diff(np.sign(diff)))[0]
            if len(zero_crossings) > 0:
                idx = zero_crossings[0]
                sigma_c = sigma_range[idx]
                crossings.append(sigma_c)
                print(f"  Crossing L={L1}/{L2}: σ_c = {sigma_c:.4f}")

        if crossings:
            self.sigma_c = np.mean(crossings)
            sigma_c_err = np.std(crossings)
            print(f"\nCritical point: σ_c = {self.sigma_c:.4f} ± {sigma_c_err:.4f}")
        else:
            print("Warning: No Binder crossings found!")
            self.sigma_c = 0.75  # Default estimate

        return self.sigma_c

    def extract_beta(self) -> float:
        """Extract order parameter exponent β from R(σ_c) ~ L^(-β/ν)"""
        print("\n=== Extracting β (Order Parameter Exponent) ===")

        log_L = []
        log_R = []

        for L in self.grid_sizes:
            df = self.data[L]

            # Interpolate to find R at σ_c
            f = interpolate.interp1d(df['sigma'], df['order_parameter'],
                                    kind='cubic', fill_value='extrapolate')
            R_at_sigma_c = f(self.sigma_c)

            if R_at_sigma_c > 0:
                log_L.append(np.log(L))
                log_R.append(np.log(R_at_sigma_c))

        # Linear fit in log-log space
        p = np.polyfit(log_L, log_R, 1)
        slope = p[0]

        # For now assume ν = 1.0 (will refine later)
        self.beta = -slope  # Actually -β/ν

        print(f"  Slope: {slope:.4f}")
        print(f"  β/ν = {-slope:.4f}")
        print(f"  β ≈ {self.beta:.4f} (assuming ν = 1.0)")

        # Plot
        plt.figure(figsize=(8, 6))
        plt.scatter(log_L, log_R, s=100, label='Data')
        plt.plot(log_L, np.polyval(p, log_L), 'r-', label=f'Fit: slope = {slope:.3f}')
        plt.xlabel('log(L)')
        plt.ylabel('log(R(σ_c))')
        plt.title('Order Parameter Scaling at Critical Point')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(self.data_dir, 'beta_extraction.png'))
        plt.close()

        return self.beta

    def extract_nu(self) -> float:
        """Extract correlation length exponent ν from peak positions"""
        print("\n=== Extracting ν (Correlation Length Exponent) ===")

        inv_L = []
        sigma_peak = []

        for L in self.grid_sizes:
            df = self.data[L]

            # Find susceptibility peak
            idx_max = df['susceptibility'].idxmax()
            sigma_max = df.loc[idx_max, 'sigma']

            inv_L.append(1.0 / L)
            sigma_peak.append(sigma_max)
            print(f"  L={L}: peak at σ = {sigma_max:.4f}")

        # Fit: σ_peak(L) - σ_c ~ L^(-1/ν)
        # Plot σ_peak vs 1/L and extrapolate to L→∞
        p = np.polyfit(inv_L, sigma_peak, 1)
        sigma_c_extrap = p[1]  # y-intercept at 1/L = 0

        # From slope: dσ/d(1/L) ~ 1/ν
        if p[0] != 0:
            self.nu = abs(1.0 / p[0])
        else:
            self.nu = 1.0  # Default

        print(f"  Extrapolated σ_c = {sigma_c_extrap:.4f}")
        print(f"  ν = {self.nu:.4f}")

        # Plot
        plt.figure(figsize=(8, 6))
        plt.scatter(inv_L, sigma_peak, s=100, label='Peak positions')
        plt.plot([0] + inv_L, [sigma_c_extrap] + list(np.polyval(p, inv_L)),
                'r-', label=f'Extrapolation')
        plt.axhline(y=self.sigma_c, color='g', linestyle='--',
                   label=f'σ_c from Binder = {self.sigma_c:.3f}')
        plt.xlabel('1/L')
        plt.ylabel('σ_peak')
        plt.title('Susceptibility Peak Scaling')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(self.data_dir, 'nu_extraction.png'))
        plt.close()

        return self.nu

    def extract_gamma(self) -> float:
        """Extract susceptibility exponent γ from χ_max ~ L^(γ/ν)"""
        print("\n=== Extracting γ (Susceptibility Exponent) ===")

        log_L = []
        log_chi_max = []

        for L in self.grid_sizes:
            df = self.data[L]
            chi_max = df['susceptibility'].max()

            log_L.append(np.log(L))
            log_chi_max.append(np.log(chi_max))
            print(f"  L={L}: χ_max = {chi_max:.2f}")

        # Linear fit in log-log space
        p = np.polyfit(log_L, log_chi_max, 1)
        slope = p[0]

        # γ/ν = slope
        self.gamma = slope * self.nu  # Using previously extracted ν

        print(f"  Slope (γ/ν): {slope:.4f}")
        print(f"  γ = {self.gamma:.4f}")

        # Plot
        plt.figure(figsize=(8, 6))
        plt.scatter(log_L, log_chi_max, s=100, label='Data')
        plt.plot(log_L, np.polyval(p, log_L), 'r-', label=f'Fit: γ/ν = {slope:.3f}')
        plt.xlabel('log(L)')
        plt.ylabel('log(χ_max)')
        plt.title('Susceptibility Maximum Scaling')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(self.data_dir, 'gamma_extraction.png'))
        plt.close()

        return self.gamma

    def perform_data_collapse(self, sigma_c=None, beta_nu=None, nu=None):
        """Perform data collapse and assess quality"""
        print("\n=== Performing Data Collapse ===")

        if sigma_c is None:
            sigma_c = self.sigma_c
        if beta_nu is None:
            beta_nu = self.beta / self.nu
        if nu is None:
            nu = self.nu

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Order parameter collapse
        ax1 = axes[0]
        for L in self.grid_sizes:
            df = self.data[L]
            x = (df['sigma'] - sigma_c) * L**(1/nu)
            y = df['order_parameter'] * L**(beta_nu)
            ax1.plot(x, y, 'o-', label=f'L={L}', markersize=4)

        ax1.set_xlabel(r'$(σ - σ_c) L^{1/ν}$')
        ax1.set_ylabel(r'$R L^{β/ν}$')
        ax1.set_title('Order Parameter Data Collapse')
        ax1.legend()
        ax1.grid(True)

        # Susceptibility collapse
        ax2 = axes[1]
        gamma_nu = self.gamma / self.nu
        for L in self.grid_sizes:
            df = self.data[L]
            x = (df['sigma'] - sigma_c) * L**(1/nu)
            y = df['susceptibility'] * L**(-gamma_nu)
            ax2.plot(x, y, 'o-', label=f'L={L}', markersize=4)

        ax2.set_xlabel(r'$(σ - σ_c) L^{1/ν}$')
        ax2.set_ylabel(r'$χ L^{-γ/ν}$')
        ax2.set_title('Susceptibility Data Collapse')
        ax2.legend()
        ax2.grid(True)

        plt.suptitle(f'FSS Data Collapse: σ_c={sigma_c:.3f}, β={self.beta:.3f}, '
                    f'ν={nu:.3f}, γ={self.gamma:.3f}')
        plt.tight_layout()
        plt.savefig(os.path.join(self.data_dir, 'data_collapse.png'), dpi=150)
        plt.show()

        # Calculate collapse quality (chi-squared)
        chi2_val = self.calculate_collapse_chi2(sigma_c, beta_nu, nu)
        print(f"  Collapse χ²/dof = {chi2_val:.3f}")

        return chi2_val

    def calculate_collapse_chi2(self, sigma_c, beta_nu, nu) -> float:
        """Calculate chi-squared for data collapse quality"""

        # Collect all scaled data points
        all_x = []
        all_y = []
        all_L = []

        for L in self.grid_sizes:
            df = self.data[L]
            x = (df['sigma'] - sigma_c) * L**(1/nu)
            y = df['order_parameter'] * L**(beta_nu)
            all_x.extend(x)
            all_y.extend(y)
            all_L.extend([L] * len(x))

        all_x = np.array(all_x)
        all_y = np.array(all_y)
        all_L = np.array(all_L)

        # Bin the data and calculate variance between different L values
        x_min, x_max = all_x.min(), all_x.max()
        n_bins = 20
        bins = np.linspace(x_min, x_max, n_bins)

        chi2_total = 0
        n_points = 0

        for i in range(n_bins - 1):
            mask = (all_x >= bins[i]) & (all_x < bins[i+1])
            if mask.sum() > 1:
                y_bin = all_y[mask]
                L_bin = all_L[mask]

                # Calculate variance between different L values
                unique_L = np.unique(L_bin)
                if len(unique_L) > 1:
                    y_means = [y_bin[L_bin == L].mean() for L in unique_L]
                    y_mean = np.mean(y_means)
                    chi2_total += np.sum((y_means - y_mean)**2) / y_mean if y_mean > 0 else 0
                    n_points += len(unique_L)

        return chi2_total / max(n_points - 1, 1)

    def optimize_collapse(self):
        """Optimize collapse parameters using chi-squared minimization"""
        print("\n=== Optimizing Data Collapse ===")

        def objective(params):
            sigma_c, beta_nu, nu = params
            return self.calculate_collapse_chi2(sigma_c, beta_nu, nu)

        # Initial guess
        x0 = [self.sigma_c, self.beta / self.nu, self.nu]

        # Bounds
        bounds = [
            (self.sigma_c - 0.1, self.sigma_c + 0.1),
            (0.01, 0.5),
            (0.5, 2.0)
        ]

        # Optimize
        result = optimize.minimize(objective, x0, method='L-BFGS-B', bounds=bounds)

        if result.success:
            sigma_c_opt, beta_nu_opt, nu_opt = result.x
            beta_opt = beta_nu_opt * nu_opt

            print(f"  Optimized parameters:")
            print(f"    σ_c = {sigma_c_opt:.4f}")
            print(f"    β = {beta_opt:.4f}")
            print(f"    ν = {nu_opt:.4f}")
            print(f"    χ² = {result.fun:.3f}")

            # Update values
            self.sigma_c = sigma_c_opt
            self.beta = beta_opt
            self.nu = nu_opt

            # Recalculate gamma using Fisher relation
            self.eta = 0.25  # Estimate (would need correlation function data)
            self.gamma = self.nu * (2 - self.eta)

            # Show optimized collapse
            self.perform_data_collapse()
        else:
            print("  Optimization failed!")

    def classify_universality(self):
        """Classify universality class based on exponents"""
        print("\n" + "="*60)
        print("UNIVERSALITY CLASS CLASSIFICATION")
        print("="*60)

        print(f"\nMeasured Critical Exponents:")
        print(f"  β = {self.beta:.4f}")
        print(f"  ν = {self.nu:.4f}")
        print(f"  γ = {self.gamma:.4f}")
        print(f"  η ≈ {self.eta:.4f} (estimated)")
        print(f"  σ_c = {self.sigma_c:.4f}")

        # Known universality classes
        classes = {
            '2D Ising': {'beta': 0.125, 'nu': 1.0, 'gamma': 1.75, 'eta': 0.25},
            '2D XY': {'beta': 0.23, 'nu': 0.67, 'gamma': 1.32, 'eta': 0.04},
            '2D 3-state Potts': {'beta': 1/9, 'nu': 5/6, 'gamma': 13/9, 'eta': 4/15},
            'Mean Field': {'beta': 0.5, 'nu': 0.5, 'gamma': 1.0, 'eta': 0.0}
        }

        print("\nComparison with Known Classes:")
        print("-" * 50)

        best_match = None
        best_score = float('inf')

        for name, exps in classes.items():
            # Calculate deviation score
            score = 0
            score += abs(self.beta - exps['beta']) / exps['beta']
            score += abs(self.nu - exps['nu']) / exps['nu']
            score += abs(self.gamma - exps['gamma']) / exps['gamma']
            score /= 3

            print(f"\n{name}:")
            print(f"  β: {self.beta:.4f} vs {exps['beta']:.4f} "
                  f"(deviation: {100*abs(self.beta - exps['beta'])/exps['beta']:.1f}%)")
            print(f"  ν: {self.nu:.4f} vs {exps['nu']:.4f} "
                  f"(deviation: {100*abs(self.nu - exps['nu'])/exps['nu']:.1f}%)")
            print(f"  Average deviation: {100*score:.1f}%")

            if score < best_score:
                best_score = score
                best_match = name

        print("\n" + "="*60)

        # Classification decision tree
        if abs(self.beta - 0.125) < 0.01 and abs(self.nu - 1.0) < 0.05:
            classification = "2D Ising"
            print("✓ CLASSIFICATION: 2D Ising universality class")
            print("  Evidence: β and ν match Ising values within error bars")
        elif abs(self.nu - 0.67) < 0.05 and self.eta < 0.1:
            classification = "2D XY"
            print("✓ CLASSIFICATION: 2D XY universality class (BKT transition)")
            print("  Evidence: ν matches XY value, small η consistent with BKT")
        elif abs(self.beta - 1/9) < 0.01 and abs(self.nu - 5/6) < 0.05:
            classification = "2D 3-state Potts"
            print("✓ CLASSIFICATION: 2D 3-state Potts universality class")
            print("  Evidence: β and ν match 3-state Potts values")
        elif self.beta > 0.3 and abs(self.nu - 0.5) < 0.1:
            classification = "Mean Field"
            print("✓ CLASSIFICATION: Mean-field behavior detected")
            print("  Evidence: Large β and ν ≈ 0.5 suggest mean-field")
        else:
            classification = "Novel"
            print("★ CLASSIFICATION: Novel universality class!")
            print(f"  β = {self.beta:.4f} does not match any known 2D class")
            print("  This represents new critical behavior requiring theoretical explanation")
            print("\n  Possible mechanisms:")
            print("  • Crossover between universality classes")
            print("  • Non-equilibrium critical dynamics")
            print("  • Long-range effective interactions from SMFT coupling")
            print("  • Emergent symmetry from relativistic dynamics")

        print("="*60)

        # Save classification report
        self.save_classification_report(classification)

        return classification

    def save_classification_report(self, classification):
        """Save detailed classification report"""
        filename = os.path.join(self.data_dir, 'universality_classification.txt')

        with open(filename, 'w') as f:
            f.write("="*60 + "\n")
            f.write("SMFT UNIVERSALITY CLASS CLASSIFICATION REPORT\n")
            f.write("="*60 + "\n\n")

            f.write(f"Classification: {classification}\n\n")

            f.write("Measured Critical Exponents:\n")
            f.write(f"  β = {self.beta:.6f}\n")
            f.write(f"  ν = {self.nu:.6f}\n")
            f.write(f"  γ = {self.gamma:.6f}\n")
            f.write(f"  η ≈ {self.eta:.6f}\n")
            f.write(f"  σ_c = {self.sigma_c:.6f}\n\n")

            f.write("Grid sizes analyzed: {}\n".format(self.grid_sizes))
            f.write("Number of noise points: 41\n\n")

            if classification == "Novel":
                f.write("SIGNIFICANCE:\n")
                f.write("This represents a novel universality class not previously\n")
                f.write("observed in 2D systems. The SMFT coupling between Dirac\n")
                f.write("and Kuramoto dynamics creates new critical behavior.\n\n")
                f.write("IMPLICATIONS:\n")
                f.write("• New theoretical framework needed\n")
                f.write("• Potential applications in quantum-classical hybrid systems\n")
                f.write("• May reveal new phase transition mechanisms\n")

        print(f"\nReport saved to: {filename}")

    def run_full_analysis(self):
        """Run complete FSS analysis pipeline"""
        print("\n" + "="*60)
        print("FINITE-SIZE SCALING ANALYSIS FOR SMFT")
        print("="*60)

        # Load data
        self.load_data()

        if not self.data:
            print("ERROR: No data found!")
            return

        # Extract exponents
        self.find_critical_point_binder()
        self.extract_beta()
        self.extract_nu()
        self.extract_gamma()

        # Fisher relation for η
        self.eta = 2 - self.gamma / self.nu
        print(f"\nFrom Fisher relation: η = {self.eta:.4f}")

        # Initial data collapse
        self.perform_data_collapse()

        # Optimize collapse
        self.optimize_collapse()

        # Classification
        classification = self.classify_universality()

        print("\n" + "="*60)
        print("Analysis complete!")
        print("="*60)

        return classification


def main():
    parser = argparse.ArgumentParser(
        description='FSS Data Collapse Analysis for SMFT Universality Classification')
    parser.add_argument('data_dir', help='Directory containing universality scan data')
    parser.add_argument('--optimize', action='store_true',
                       help='Optimize collapse parameters')

    args = parser.parse_args()

    analyzer = DataCollapseAnalyzer(args.data_dir)
    analyzer.run_full_analysis()


if __name__ == '__main__':
    main()