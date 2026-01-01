#!/usr/bin/env python3
"""
Finite-Size Scaling Analysis Script

This script performs comprehensive FSS analysis on phase transition data
from SMFT simulations. It:
1. Loads data from multiple system sizes
2. Finds critical point from Binder cumulant crossing
3. Extracts critical exponents via data collapse
4. Identifies universality class
5. Generates publication-quality plots
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.interpolate import interp1d
import pandas as pd
from glob import glob
import os
import sys
from typing import Dict, List, Tuple, Optional

# Known 2D universality classes
UNIVERSALITY_CLASSES = {
    '2D_Ising': {'beta': 0.125, 'nu': 1.0, 'gamma': 1.75, 'eta': 0.25},
    '2D_XY': {'beta': 0.23, 'nu': 0.67, 'gamma': 1.32, 'eta': 0.04},
    '2D_3state_Potts': {'beta': 0.111, 'nu': 0.833, 'gamma': 1.44, 'eta': 0.267},
    '2D_4state_Potts': {'beta': 0.083, 'nu': 0.667, 'gamma': 1.17, 'eta': 0.25},
    'Mean_Field': {'beta': 0.5, 'nu': 0.5, 'gamma': 1.0, 'eta': 0.0}
}


class FSSAnalyzer:
    """Main FSS analysis class"""

    def __init__(self, data_dir: str):
        """
        Initialize analyzer with data directory

        Args:
            data_dir: Directory containing FSS data files
        """
        self.data_dir = data_dir
        self.data_by_L = {}
        self.sigma_c = None
        self.exponents = {}

    def load_data(self) -> None:
        """Load all FSS data files"""
        print("Loading FSS data...")

        # Find all data files
        pattern = os.path.join(self.data_dir, 'fss_data_L*.csv')
        files = glob(pattern)

        if not files:
            # Try alternative naming
            pattern = os.path.join(self.data_dir, '**/observables.csv')
            files = glob(pattern, recursive=True)

        for file in sorted(files):
            # Extract L from filename or path
            if 'L' in file:
                L = int(file.split('L')[1].split('.')[0].split('_')[0].split('/')[0])
            else:
                # Try to extract from path
                parts = file.split('/')
                for part in parts:
                    if 'x' in part and part.split('x')[0].isdigit():
                        L = int(part.split('x')[0])
                        break
                else:
                    continue

            # Load data
            df = pd.read_csv(file)

            # Store processed data
            self.data_by_L[L] = {
                'sigma': df['sigma'].values if 'sigma' in df else df.index.values,
                'R_mean': df['R_mean'].values if 'R_mean' in df else df['mean_R'].values,
                'R_variance': df['R_variance'].values if 'R_variance' in df else None,
                'binder': self._compute_binder(df),
                'susceptibility': self._compute_susceptibility(df, L)
            }

        print(f"Loaded data for L = {sorted(self.data_by_L.keys())}")

    def _compute_binder(self, df: pd.DataFrame) -> np.ndarray:
        """Compute Binder cumulant from data"""
        if 'Binder_cumulant' in df:
            return df['Binder_cumulant'].values
        elif 'R_squared' in df and 'R_fourth' in df:
            R2 = df['R_squared'].values
            R4 = df['R_fourth'].values
            return 1 - R4 / (3 * R2**2)
        else:
            return None

    def _compute_susceptibility(self, df: pd.DataFrame, L: int) -> np.ndarray:
        """Compute susceptibility from data"""
        if 'susceptibility' in df:
            return df['susceptibility'].values
        elif 'R_variance' in df:
            return L**2 * df['R_variance'].values
        else:
            return None

    def find_critical_point_binder(self) -> float:
        """Find critical point from Binder cumulant crossing"""
        print("\nFinding critical point from Binder crossing...")

        # Get two largest system sizes
        L_values = sorted(self.data_by_L.keys())
        if len(L_values) < 2:
            print("Need at least 2 system sizes for Binder crossing")
            return 0.85  # Default guess

        L1, L2 = L_values[-2], L_values[-1]
        data1 = self.data_by_L[L1]
        data2 = self.data_by_L[L2]

        if data1['binder'] is None or data2['binder'] is None:
            print("Binder data not available")
            return 0.85

        # Interpolate curves
        f1 = interp1d(data1['sigma'], data1['binder'], kind='cubic', fill_value='extrapolate')
        f2 = interp1d(data2['sigma'], data2['binder'], kind='cubic', fill_value='extrapolate')

        # Find crossing
        sigma_range = np.linspace(0.7, 1.0, 1000)
        diff = np.abs(f1(sigma_range) - f2(sigma_range))
        idx_min = np.argmin(diff)

        self.sigma_c = sigma_range[idx_min]
        print(f"Critical point: σ_c = {self.sigma_c:.4f}")

        return self.sigma_c

    def fit_critical_exponents(self) -> Dict[str, float]:
        """Extract critical exponents via data collapse optimization"""
        print("\nFitting critical exponents...")

        if self.sigma_c is None:
            self.sigma_c = self.find_critical_point_binder()

        # Define collapse quality function
        def collapse_quality(params):
            sigma_c, beta_over_nu, one_over_nu = params

            X_all = []
            Y_all = []

            for L, data in self.data_by_L.items():
                X = (data['sigma'] - sigma_c) * L**one_over_nu
                Y = data['R_mean'] * L**beta_over_nu
                X_all.extend(X)
                Y_all.extend(Y)

            # Measure spread in Y for similar X
            X_all = np.array(X_all)
            Y_all = np.array(Y_all)

            # Bin X values and compute variance in each bin
            n_bins = 20
            X_min, X_max = X_all.min(), X_all.max()
            bin_edges = np.linspace(X_min, X_max, n_bins + 1)

            total_var = 0
            for i in range(n_bins):
                mask = (X_all >= bin_edges[i]) & (X_all < bin_edges[i+1])
                if np.sum(mask) > 1:
                    Y_bin = Y_all[mask]
                    total_var += np.var(Y_bin)

            return total_var

        # Initial guess (2D Ising values)
        x0 = [self.sigma_c, 0.125, 1.0]
        bounds = [(self.sigma_c - 0.1, self.sigma_c + 0.1),
                  (0.05, 0.5),
                  (0.5, 2.0)]

        # Optimize
        result = minimize(collapse_quality, x0, bounds=bounds, method='L-BFGS-B')

        sigma_c_opt, beta_over_nu, one_over_nu = result.x

        self.exponents = {
            'sigma_c': sigma_c_opt,
            'beta': beta_over_nu / one_over_nu,
            'nu': 1.0 / one_over_nu,
            'beta_over_nu': beta_over_nu,
            'one_over_nu': one_over_nu
        }

        # Estimate gamma from susceptibility scaling
        self._fit_gamma()

        # Calculate eta from hyperscaling
        if 'gamma' in self.exponents and 'nu' in self.exponents:
            self.exponents['eta'] = 2.0 - self.exponents['gamma'] / self.exponents['nu']

        print(f"Critical exponents:")
        print(f"  σ_c = {self.exponents['sigma_c']:.4f}")
        print(f"  β = {self.exponents['beta']:.4f}")
        print(f"  ν = {self.exponents['nu']:.4f}")
        if 'gamma' in self.exponents:
            print(f"  γ = {self.exponents['gamma']:.4f}")
        if 'eta' in self.exponents:
            print(f"  η = {self.exponents['eta']:.4f}")

        return self.exponents

    def _fit_gamma(self) -> None:
        """Fit gamma exponent from susceptibility scaling"""
        log_L = []
        log_chi_max = []

        for L, data in self.data_by_L.items():
            if data['susceptibility'] is not None:
                chi_max = np.max(data['susceptibility'])
                log_L.append(np.log(L))
                log_chi_max.append(np.log(chi_max))

        if len(log_L) >= 2:
            # Linear fit in log-log space
            coeffs = np.polyfit(log_L, log_chi_max, 1)
            gamma_over_nu = coeffs[0]

            if 'nu' in self.exponents:
                self.exponents['gamma'] = gamma_over_nu * self.exponents['nu']

    def identify_universality_class(self) -> str:
        """Identify universality class from exponents"""
        print("\nIdentifying universality class...")

        if not self.exponents:
            self.fit_critical_exponents()

        best_match = None
        best_distance = float('inf')

        beta_meas = self.exponents.get('beta', 0)
        nu_meas = self.exponents.get('nu', 0)
        gamma_meas = self.exponents.get('gamma', 0)

        for name, exps in UNIVERSALITY_CLASSES.items():
            # Compute normalized distance
            d_beta = abs(beta_meas - exps['beta']) / exps['beta'] if exps['beta'] > 0 else 1
            d_nu = abs(nu_meas - exps['nu']) / exps['nu'] if exps['nu'] > 0 else 1
            d_gamma = abs(gamma_meas - exps['gamma']) / exps['gamma'] if gamma_meas > 0 and exps['gamma'] > 0 else 1

            distance = np.sqrt((d_beta**2 + d_nu**2 + d_gamma**2) / 3)

            print(f"  Distance to {name}: {distance:.3f}")

            if distance < best_distance:
                best_distance = distance
                best_match = name

        if best_distance > 0.2:
            universality = "Novel"
            print(f"\n*** NOVEL UNIVERSALITY CLASS DETECTED! ***")
            print(f"β = {beta_meas:.4f} does not match any known 2D class")
        else:
            universality = best_match
            print(f"\n✓ Best match: {universality} (distance = {best_distance:.3f})")

        return universality

    def plot_results(self, save_dir: Optional[str] = None) -> None:
        """Generate comprehensive plots"""
        print("\nGenerating plots...")

        if save_dir is None:
            save_dir = self.data_dir

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        # Plot 1: Order parameter
        ax = axes[0, 0]
        for L in sorted(self.data_by_L.keys()):
            data = self.data_by_L[L]
            ax.plot(data['sigma'], data['R_mean'], 'o-', label=f'L={L}', markersize=3)
        ax.set_xlabel(r'$\sigma$')
        ax.set_ylabel(r'$\langle R \rangle$')
        ax.set_title('Order Parameter')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 2: Binder cumulant
        ax = axes[0, 1]
        for L in sorted(self.data_by_L.keys()):
            data = self.data_by_L[L]
            if data['binder'] is not None:
                ax.plot(data['sigma'], data['binder'], 'o-', label=f'L={L}', markersize=3)
        if self.sigma_c:
            ax.axvline(self.sigma_c, color='r', linestyle='--', label=f'σ_c={self.sigma_c:.3f}')
        ax.set_xlabel(r'$\sigma$')
        ax.set_ylabel(r'$U_L$')
        ax.set_title('Binder Cumulant')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 3: Susceptibility
        ax = axes[0, 2]
        for L in sorted(self.data_by_L.keys()):
            data = self.data_by_L[L]
            if data['susceptibility'] is not None:
                ax.plot(data['sigma'], data['susceptibility']/L**2, 'o-', label=f'L={L}', markersize=3)
        ax.set_xlabel(r'$\sigma$')
        ax.set_ylabel(r'$\chi/L^2$')
        ax.set_title('Susceptibility (scaled)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 4: Data collapse
        ax = axes[1, 0]
        if self.exponents:
            colors = plt.cm.viridis(np.linspace(0, 1, len(self.data_by_L)))
            for i, (L, data) in enumerate(sorted(self.data_by_L.items())):
                X = (data['sigma'] - self.exponents['sigma_c']) * L**self.exponents['one_over_nu']
                Y = data['R_mean'] * L**self.exponents['beta_over_nu']
                ax.plot(X, Y, 'o', label=f'L={L}', markersize=3, color=colors[i])
            ax.set_xlabel(r'$(σ - σ_c)L^{1/ν}$')
            ax.set_ylabel(r'$\langle R \rangle L^{β/ν}$')
            ax.set_title('Data Collapse')
            ax.legend()
            ax.grid(True, alpha=0.3)

        # Plot 5: Critical region fit
        ax = axes[1, 1]
        if self.sigma_c:
            for L in sorted(self.data_by_L.keys()):
                data = self.data_by_L[L]
                # Select points below sigma_c
                mask = (data['sigma'] < self.sigma_c) & (data['R_mean'] > 0.01)
                if np.sum(mask) > 2:
                    x = np.log(self.sigma_c - data['sigma'][mask])
                    y = np.log(data['R_mean'][mask])
                    ax.plot(x, y, 'o', label=f'L={L}', markersize=3)

            if self.exponents and 'beta' in self.exponents:
                # Show power law fit
                x_fit = np.linspace(-4, 0, 100)
                y_fit = self.exponents['beta'] * x_fit + np.log(0.5)  # Arbitrary offset
                ax.plot(x_fit, y_fit, 'r--', label=f'β={self.exponents["beta"]:.3f}')

            ax.set_xlabel(r'$\log(σ_c - σ)$')
            ax.set_ylabel(r'$\log\langle R \rangle$')
            ax.set_title('Critical Exponent β')
            ax.legend()
            ax.grid(True, alpha=0.3)

        # Plot 6: Exponent comparison
        ax = axes[1, 2]
        if self.exponents:
            measured = [self.exponents.get('beta', 0),
                       self.exponents.get('nu', 0),
                       self.exponents.get('gamma', 0)]

            x_pos = np.arange(3)
            width = 0.15

            # Plot measured values
            ax.bar(x_pos, measured, width, label='Measured', color='blue')

            # Plot known classes
            colors = ['red', 'green', 'orange', 'purple']
            for i, (name, exps) in enumerate(list(UNIVERSALITY_CLASSES.items())[:4]):
                values = [exps['beta'], exps['nu'], exps['gamma']]
                ax.bar(x_pos + (i+1)*width, values, width, label=name, color=colors[i], alpha=0.7)

            ax.set_xlabel('Exponent')
            ax.set_ylabel('Value')
            ax.set_title('Universality Class Comparison')
            ax.set_xticks(x_pos + 2*width)
            ax.set_xticklabels(['β', 'ν', 'γ'])
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(True, alpha=0.3, axis='y')

        plt.suptitle('Finite-Size Scaling Analysis', fontsize=16)
        plt.tight_layout()

        # Save figure
        output_file = os.path.join(save_dir, 'fss_analysis.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Plots saved to {output_file}")

        plt.show()

    def generate_report(self, output_file: Optional[str] = None) -> None:
        """Generate text report"""
        if output_file is None:
            output_file = os.path.join(self.data_dir, 'fss_report.txt')

        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("                    FINITE-SIZE SCALING ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write("Data Summary:\n")
            f.write("-" * 40 + "\n")
            f.write(f"System sizes: L = {sorted(self.data_by_L.keys())}\n")
            f.write(f"Total data points: {sum(len(d['sigma']) for d in self.data_by_L.values())}\n\n")

            if self.exponents:
                f.write("Critical Exponents:\n")
                f.write("-" * 40 + "\n")
                f.write(f"σ_c = {self.exponents.get('sigma_c', 0):.4f}\n")
                f.write(f"β = {self.exponents.get('beta', 0):.4f}\n")
                f.write(f"ν = {self.exponents.get('nu', 0):.4f}\n")
                f.write(f"γ = {self.exponents.get('gamma', 0):.4f}\n")
                f.write(f"η = {self.exponents.get('eta', 0):.4f}\n\n")

            universality = self.identify_universality_class()
            f.write("Universality Class:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Identified: {universality}\n\n")

            f.write("Comparison with Known Classes:\n")
            f.write("-" * 40 + "\n")
            f.write(f"{'Class':<20} {'β':<10} {'ν':<10} {'γ':<10}\n")
            f.write("-" * 50 + "\n")

            if self.exponents:
                f.write(f"{'Measured':<20} {self.exponents.get('beta', 0):<10.4f} "
                       f"{self.exponents.get('nu', 0):<10.4f} {self.exponents.get('gamma', 0):<10.4f}\n")

            for name, exps in UNIVERSALITY_CLASSES.items():
                f.write(f"{name:<20} {exps['beta']:<10.4f} {exps['nu']:<10.4f} {exps['gamma']:<10.4f}\n")

            f.write("\n" + "=" * 80 + "\n")

        print(f"Report saved to {output_file}")


def main():
    """Main entry point"""
    if len(sys.argv) < 2:
        print("Usage: python analyze_fss.py <data_directory>")
        sys.exit(1)

    data_dir = sys.argv[1]

    # Create analyzer
    analyzer = FSSAnalyzer(data_dir)

    # Run analysis
    analyzer.load_data()
    analyzer.find_critical_point_binder()
    analyzer.fit_critical_exponents()
    universality = analyzer.identify_universality_class()

    # Generate outputs
    analyzer.plot_results()
    analyzer.generate_report()

    # Final summary
    print("\n" + "=" * 60)
    print("FINAL RESULT:")
    print(f"  Universality Class: {universality}")
    if analyzer.exponents:
        print(f"  β = {analyzer.exponents.get('beta', 0):.4f}")
        print(f"  ν = {analyzer.exponents.get('nu', 0):.4f}")
    print("=" * 60)


if __name__ == "__main__":
    main()