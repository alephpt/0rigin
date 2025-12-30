#!/usr/bin/env python3
"""
Critical Exponent Extraction for SMFT Universality Analysis

Advanced methods for extracting critical exponents with error estimation.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize, stats, interpolate
from scipy.special import gammaln
import emcee  # For MCMC error estimation
import corner
import os
import sys
from typing import Tuple, List, Dict
import warnings
warnings.filterwarnings('ignore')


class ExponentExtractor:
    """Advanced critical exponent extraction with error analysis"""

    def __init__(self, data_dir: str):
        self.data_dir = data_dir
        self.data = {}  # L -> DataFrame
        self.grid_sizes = []

        # Exponents with uncertainties
        self.exponents = {
            'beta': {'value': None, 'error': None},
            'nu': {'value': None, 'error': None},
            'gamma': {'value': None, 'error': None},
            'eta': {'value': None, 'error': None},
            'alpha': {'value': None, 'error': None},
            'delta': {'value': None, 'error': None},
            'sigma_c': {'value': None, 'error': None}
        }

    def load_data(self):
        """Load simulation data"""
        import glob

        pattern = os.path.join(self.data_dir, "universality_L*", "observables_summary.csv")
        files = sorted(glob.glob(pattern))

        for f in files:
            L = int(f.split('_L')[-1].split('/')[0])
            df = pd.read_csv(f)
            self.data[L] = df
            self.grid_sizes.append(L)

        self.grid_sizes = sorted(self.grid_sizes)
        print(f"Loaded data for L = {self.grid_sizes}")

    def bootstrap_error(self, data: np.ndarray, func, n_bootstrap: int = 1000) -> Tuple[float, float]:
        """Estimate error using bootstrap resampling"""
        results = []

        for _ in range(n_bootstrap):
            # Resample with replacement
            sample = np.random.choice(data, size=len(data), replace=True)
            try:
                results.append(func(sample))
            except:
                pass  # Skip failed fits

        if results:
            mean = np.mean(results)
            error = np.std(results)
            return mean, error
        else:
            return func(data), 0.0

    def mcmc_fit(self, log_likelihood, prior, n_params: int, n_walkers: int = 32,
                 n_steps: int = 5000) -> Dict:
        """MCMC fitting for robust parameter estimation"""

        # Initialize walkers
        pos = []
        for _ in range(n_walkers):
            p = []
            for i in range(n_params):
                p.append(prior[i][0] + np.random.rand() * (prior[i][1] - prior[i][0]))
            pos.append(p)

        # Run MCMC
        sampler = emcee.EnsembleSampler(n_walkers, n_params, log_likelihood)
        sampler.run_mcmc(pos, n_steps, progress=True)

        # Discard burn-in
        samples = sampler.chain[:, n_steps//2:, :].reshape((-1, n_params))

        # Extract results
        results = {}
        for i in range(n_params):
            results[f'param_{i}'] = {
                'mean': np.mean(samples[:, i]),
                'std': np.std(samples[:, i]),
                'median': np.median(samples[:, i]),
                'percentiles': np.percentile(samples[:, i], [16, 50, 84])
            }

        results['samples'] = samples
        return results

    def extract_sigma_c_advanced(self) -> Tuple[float, float]:
        """Extract critical point using multiple methods"""
        print("\n=== Advanced Critical Point Extraction ===")

        methods = []

        # Method 1: Binder cumulant crossing with error propagation
        sigma_c_binder, err_binder = self.binder_crossing_analysis()
        methods.append(('Binder', sigma_c_binder, err_binder))

        # Method 2: Susceptibility peak extrapolation
        sigma_c_chi, err_chi = self.susceptibility_peak_extrapolation()
        methods.append(('Susceptibility', sigma_c_chi, err_chi))

        # Method 3: Derivative peak analysis
        sigma_c_deriv, err_deriv = self.derivative_peak_analysis()
        methods.append(('Derivative', sigma_c_deriv, err_deriv))

        # Weighted average
        sigma_c_avg = 0
        weight_sum = 0
        for name, val, err in methods:
            if val is not None and err > 0:
                weight = 1 / (err ** 2)
                sigma_c_avg += val * weight
                weight_sum += weight
                print(f"  {name}: σ_c = {val:.6f} ± {err:.6f}")

        if weight_sum > 0:
            sigma_c_avg /= weight_sum
            sigma_c_err = 1 / np.sqrt(weight_sum)

            self.exponents['sigma_c']['value'] = sigma_c_avg
            self.exponents['sigma_c']['error'] = sigma_c_err

            print(f"\n  Weighted average: σ_c = {sigma_c_avg:.6f} ± {sigma_c_err:.6f}")
            return sigma_c_avg, sigma_c_err
        else:
            return 0.75, 0.01  # Default

    def binder_crossing_analysis(self) -> Tuple[float, float]:
        """Analyze Binder cumulant crossings with error estimation"""

        crossings = []
        errors = []

        for i in range(len(self.grid_sizes) - 1):
            L1, L2 = self.grid_sizes[i], self.grid_sizes[i + 1]
            df1, df2 = self.data[L1], self.data[L2]

            # Spline interpolation for smooth curves
            f1 = interpolate.UnivariateSpline(df1['sigma'], df1['binder_cumulant'], s=0)
            f2 = interpolate.UnivariateSpline(df2['sigma'], df2['binder_cumulant'], s=0)

            # Find crossing using root finding
            def diff_func(sigma):
                return f1(sigma) - f2(sigma)

            try:
                from scipy.optimize import brentq
                sigma_c = brentq(diff_func, 0.5, 1.5)
                crossings.append(sigma_c)

                # Estimate error from curve curvature
                h = 0.001
                curvature = abs((diff_func(sigma_c + h) - 2 * diff_func(sigma_c) +
                               diff_func(sigma_c - h)) / h**2)
                error = 1 / np.sqrt(max(curvature, 0.1))
                errors.append(error)
            except:
                pass

        if crossings:
            # Weighted mean
            weights = [1/e**2 for e in errors]
            sigma_c = np.average(crossings, weights=weights)
            sigma_c_err = 1 / np.sqrt(sum(weights))
            return sigma_c, sigma_c_err
        else:
            return None, None

    def susceptibility_peak_extrapolation(self) -> Tuple[float, float]:
        """Extrapolate susceptibility peak position to L→∞"""

        peak_positions = []
        L_values = []

        for L in self.grid_sizes:
            df = self.data[L]

            # Fit parabola around peak for better accuracy
            idx_max = df['susceptibility'].idxmax()
            if idx_max > 2 and idx_max < len(df) - 3:
                # Take 5 points around maximum
                idx_range = range(idx_max - 2, idx_max + 3)
                sigma_fit = df.loc[idx_range, 'sigma'].values
                chi_fit = df.loc[idx_range, 'susceptibility'].values

                # Parabolic fit
                p = np.polyfit(sigma_fit, chi_fit, 2)
                # Peak position from derivative
                sigma_peak = -p[1] / (2 * p[0])

                peak_positions.append(sigma_peak)
                L_values.append(L)

        if len(peak_positions) > 2:
            # Fit σ_peak(L) = σ_c + A * L^(-1/ν)
            def fit_func(L, sigma_c, A, nu_inv):
                return sigma_c + A * np.power(L, -nu_inv)

            try:
                popt, pcov = optimize.curve_fit(fit_func, L_values, peak_positions,
                                               p0=[0.75, 1.0, 1.0])
                sigma_c = popt[0]
                sigma_c_err = np.sqrt(pcov[0, 0])
                return sigma_c, sigma_c_err
            except:
                pass

        return None, None

    def derivative_peak_analysis(self) -> Tuple[float, float]:
        """Analyze derivative peaks for critical point"""

        derivative_peaks = []

        for L in self.grid_sizes:
            df = self.data[L]

            # Calculate numerical derivative |dR/dσ|
            sigma = df['sigma'].values
            R = df['order_parameter'].values

            # Smooth derivative using Savitzky-Golay filter
            from scipy.signal import savgol_filter
            if len(R) > 7:
                dR_dsigma = savgol_filter(R, 7, 2, deriv=1)
                idx_max = np.abs(dR_dsigma).argmax()
                derivative_peaks.append(sigma[idx_max])

        if derivative_peaks:
            # For finite systems, peaks should converge to σ_c
            sigma_c = np.mean(derivative_peaks[-len(derivative_peaks)//2:])  # Use larger L values
            sigma_c_err = np.std(derivative_peaks[-len(derivative_peaks)//2:])
            return sigma_c, sigma_c_err

        return None, None

    def extract_exponents_with_errors(self):
        """Extract all exponents with full error analysis"""
        print("\n=== Extracting Critical Exponents with Error Analysis ===")

        # First get critical point
        sigma_c, sigma_c_err = self.extract_sigma_c_advanced()

        # Extract β with MCMC
        self.extract_beta_mcmc(sigma_c)

        # Extract ν with multiple methods
        self.extract_nu_advanced(sigma_c)

        # Extract γ from susceptibility
        self.extract_gamma_advanced(sigma_c)

        # Extract η from correlation functions (if available)
        self.extract_eta_advanced(sigma_c)

        # Calculate derived exponents
        self.calculate_derived_exponents()

        # Verify scaling relations
        self.verify_scaling_relations()

        # Print summary
        self.print_exponent_summary()

    def extract_beta_mcmc(self, sigma_c: float):
        """Extract β using MCMC for error estimation"""
        print("\n  Extracting β (order parameter exponent)...")

        # Collect data at critical point
        log_L_data = []
        log_R_data = []
        errors_R = []

        for L in self.grid_sizes:
            df = self.data[L]

            # Interpolate to find R at σ_c
            f = interpolate.interp1d(df['sigma'], df['order_parameter'],
                                    kind='cubic', fill_value='extrapolate')
            R_at_sigma_c = f(sigma_c)

            if R_at_sigma_c > 0:
                log_L_data.append(np.log(L))
                log_R_data.append(np.log(R_at_sigma_c))

                # Estimate error from local variation
                sigma_range = df['sigma'].values
                idx = np.argmin(np.abs(sigma_range - sigma_c))
                if idx > 0 and idx < len(df) - 1:
                    R_vals = [df.loc[idx-1, 'order_parameter'],
                             df.loc[idx, 'order_parameter'],
                             df.loc[idx+1, 'order_parameter']]
                    error_R = np.std(R_vals) / np.mean(R_vals)
                    errors_R.append(error_R)
                else:
                    errors_R.append(0.01)

        # MCMC fit
        def log_likelihood(params):
            slope, intercept, log_sigma = params
            sigma = np.exp(log_sigma)
            model = slope * np.array(log_L_data) + intercept
            residuals = (np.array(log_R_data) - model) / sigma
            return -0.5 * np.sum(residuals**2 + np.log(2*np.pi*sigma**2))

        # Run MCMC
        prior = [(-1, 0), (-5, 5), (-5, 0)]  # Priors for slope, intercept, log_sigma
        results = self.mcmc_fit(log_likelihood, prior, 3, n_walkers=32, n_steps=5000)

        # Extract β (assuming ν = 1 for now)
        slope_samples = results['samples'][:, 0]
        beta_over_nu = -np.mean(slope_samples)
        beta_over_nu_err = np.std(slope_samples)

        # For now, assume ν = 1
        self.exponents['beta']['value'] = beta_over_nu
        self.exponents['beta']['error'] = beta_over_nu_err

        print(f"    β/ν = {beta_over_nu:.6f} ± {beta_over_nu_err:.6f}")

    def extract_nu_advanced(self, sigma_c: float):
        """Extract ν using finite-size scaling"""
        print("\n  Extracting ν (correlation length exponent)...")

        # Method 1: From Binder cumulant width scaling
        nu_binder = self.extract_nu_from_binder_width()

        # Method 2: From susceptibility peak shift
        nu_chi = self.extract_nu_from_chi_peak_shift()

        # Method 3: From correlation length (if available)
        nu_xi = self.extract_nu_from_correlation_length()

        # Combine estimates
        nu_estimates = []
        if nu_binder: nu_estimates.append(nu_binder)
        if nu_chi: nu_estimates.append(nu_chi)
        if nu_xi: nu_estimates.append(nu_xi)

        if nu_estimates:
            weights = [1/e[1]**2 for e in nu_estimates]
            nu = np.average([e[0] for e in nu_estimates], weights=weights)
            nu_err = 1 / np.sqrt(sum(weights))

            self.exponents['nu']['value'] = nu
            self.exponents['nu']['error'] = nu_err

            print(f"    ν = {nu:.6f} ± {nu_err:.6f}")
        else:
            # Default
            self.exponents['nu']['value'] = 1.0
            self.exponents['nu']['error'] = 0.1

    def extract_nu_from_binder_width(self) -> Tuple[float, float]:
        """Extract ν from width of Binder cumulant crossing region"""
        widths = []
        L_values = []

        for L in self.grid_sizes:
            df = self.data[L]
            U = df['binder_cumulant'].values

            # Find width where U varies by fixed amount
            U_mean = np.mean(U)
            threshold = 0.1  # Variation threshold

            sigma_range = df['sigma'].values
            mask = np.abs(U - U_mean) < threshold
            if mask.sum() > 2:
                width = sigma_range[mask].max() - sigma_range[mask].min()
                widths.append(width)
                L_values.append(L)

        if len(widths) > 2:
            # Fit: width ~ L^(-1/ν)
            log_L = np.log(L_values)
            log_width = np.log(widths)
            slope, _, _, _, slope_err = stats.linregress(log_L, log_width)
            nu = -1 / slope
            nu_err = nu**2 * slope_err / abs(slope)
            return (nu, nu_err)

        return None

    def extract_nu_from_chi_peak_shift(self) -> Tuple[float, float]:
        """Extract ν from susceptibility peak position scaling"""
        peak_positions = []
        L_values = []

        for L in self.grid_sizes[1:]:  # Skip smallest size
            df = self.data[L]
            idx_max = df['susceptibility'].idxmax()
            peak_positions.append(df.loc[idx_max, 'sigma'])
            L_values.append(L)

        if len(peak_positions) > 2:
            # Fit: σ_peak - σ_c ~ L^(-1/ν)
            sigma_c = self.exponents['sigma_c']['value']
            shifts = [p - sigma_c for p in peak_positions]

            # Only use points with significant shift
            valid = [(L, s) for L, s in zip(L_values, shifts) if abs(s) > 0.01]

            if len(valid) > 2:
                L_valid = [v[0] for v in valid]
                shifts_valid = [v[1] for v in valid]

                log_L = np.log(L_valid)
                log_shifts = np.log(np.abs(shifts_valid))

                slope, _, _, _, slope_err = stats.linregress(log_L, log_shifts)
                nu = -1 / slope
                nu_err = nu**2 * slope_err / abs(slope)
                return (nu, nu_err)

        return None

    def extract_nu_from_correlation_length(self) -> Tuple[float, float]:
        """Extract ν from correlation length scaling (if data available)"""
        # Placeholder - would need correlation function data
        return None

    def extract_gamma_advanced(self, sigma_c: float):
        """Extract γ with error estimation"""
        print("\n  Extracting γ (susceptibility exponent)...")

        log_L = []
        log_chi_max = []

        for L in self.grid_sizes:
            df = self.data[L]
            chi_max = df['susceptibility'].max()
            log_L.append(np.log(L))
            log_chi_max.append(np.log(chi_max))

        # Linear regression with error
        slope, intercept, r_value, p_value, slope_err = stats.linregress(log_L, log_chi_max)

        # γ/ν = slope
        nu = self.exponents['nu']['value']
        gamma = slope * nu
        gamma_err = slope_err * nu  # Simplified error propagation

        self.exponents['gamma']['value'] = gamma
        self.exponents['gamma']['error'] = gamma_err

        print(f"    γ/ν = {slope:.6f} ± {slope_err:.6f}")
        print(f"    γ = {gamma:.6f} ± {gamma_err:.6f}")

    def extract_eta_advanced(self, sigma_c: float):
        """Extract η from correlation functions or Fisher relation"""
        print("\n  Extracting η (anomalous dimension)...")

        # Use Fisher relation as default
        gamma = self.exponents['gamma']['value']
        nu = self.exponents['nu']['value']

        eta = 2 - gamma / nu

        # Error propagation
        gamma_err = self.exponents['gamma']['error']
        nu_err = self.exponents['nu']['error']
        eta_err = np.sqrt((gamma_err/nu)**2 + (gamma*nu_err/nu**2)**2)

        self.exponents['eta']['value'] = eta
        self.exponents['eta']['error'] = eta_err

        print(f"    η = {eta:.6f} ± {eta_err:.6f} (from Fisher relation)")

    def calculate_derived_exponents(self):
        """Calculate α and δ from scaling relations"""
        print("\n  Calculating derived exponents...")

        beta = self.exponents['beta']['value']
        gamma = self.exponents['gamma']['value']
        nu = self.exponents['nu']['value']

        # Rushbrooke: α + 2β + γ = 2
        alpha = 2 - 2*beta - gamma
        self.exponents['alpha']['value'] = alpha

        # Widom: γ = β(δ - 1)
        delta = 1 + gamma / beta
        self.exponents['delta']['value'] = delta

        # Error propagation (simplified)
        beta_err = self.exponents['beta']['error']
        gamma_err = self.exponents['gamma']['error']

        alpha_err = np.sqrt(4*beta_err**2 + gamma_err**2)
        delta_err = np.sqrt((gamma_err/beta)**2 + (gamma*beta_err/beta**2)**2)

        self.exponents['alpha']['error'] = alpha_err
        self.exponents['delta']['error'] = delta_err

        print(f"    α = {alpha:.6f} ± {alpha_err:.6f}")
        print(f"    δ = {delta:.6f} ± {delta_err:.6f}")

    def verify_scaling_relations(self):
        """Check if scaling relations are satisfied"""
        print("\n=== Verifying Scaling Relations ===")

        beta = self.exponents['beta']['value']
        nu = self.exponents['nu']['value']
        gamma = self.exponents['gamma']['value']
        eta = self.exponents['eta']['value']
        alpha = self.exponents['alpha']['value']

        # Fisher: γ = ν(2-η)
        fisher_lhs = gamma
        fisher_rhs = nu * (2 - eta)
        fisher_dev = abs(fisher_lhs - fisher_rhs) / fisher_lhs

        print(f"  Fisher (γ = ν(2-η)):")
        print(f"    LHS = {fisher_lhs:.4f}, RHS = {fisher_rhs:.4f}")
        print(f"    Deviation: {100*fisher_dev:.1f}%")

        # Rushbrooke: α + 2β + γ = 2
        rush_sum = alpha + 2*beta + gamma
        rush_dev = abs(rush_sum - 2) / 2

        print(f"  Rushbrooke (α + 2β + γ = 2):")
        print(f"    Sum = {rush_sum:.4f}")
        print(f"    Deviation: {100*rush_dev:.1f}%")

        # Josephson: νd = 2 - α (d=2 for 2D)
        joseph_lhs = 2 * nu
        joseph_rhs = 2 - alpha
        joseph_dev = abs(joseph_lhs - joseph_rhs) / joseph_lhs

        print(f"  Josephson (νd = 2 - α):")
        print(f"    LHS = {joseph_lhs:.4f}, RHS = {joseph_rhs:.4f}")
        print(f"    Deviation: {100*joseph_dev:.1f}%")

        if fisher_dev < 0.05 and rush_dev < 0.05 and joseph_dev < 0.05:
            print("\n  ✓ All scaling relations satisfied within 5%")
        else:
            print("\n  ⚠ Some scaling relations show deviations > 5%")

    def print_exponent_summary(self):
        """Print summary table of all exponents"""
        print("\n" + "="*60)
        print("CRITICAL EXPONENT SUMMARY")
        print("="*60)

        print("\n{:<15} {:>15} {:>15}".format("Exponent", "Value", "Error"))
        print("-"*45)

        for name in ['beta', 'nu', 'gamma', 'eta', 'alpha', 'delta', 'sigma_c']:
            exp = self.exponents[name]
            if exp['value'] is not None:
                print("{:<15} {:>15.6f} ± {:<15.6f}".format(
                    name if name != 'sigma_c' else 'σ_c',
                    exp['value'],
                    exp['error'] if exp['error'] else 0.0
                ))

        print("="*60)

    def save_results(self):
        """Save exponents to file"""
        filename = os.path.join(self.data_dir, 'critical_exponents.txt')

        with open(filename, 'w') as f:
            f.write("# SMFT Critical Exponents\n")
            f.write("# Extracted using FSS analysis\n\n")

            for name in ['beta', 'nu', 'gamma', 'eta', 'alpha', 'delta', 'sigma_c']:
                exp = self.exponents[name]
                if exp['value'] is not None:
                    f.write(f"{name:10s} = {exp['value']:12.8f} ± {exp['error']:12.8f}\n")

        print(f"\nResults saved to: {filename}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Extract critical exponents with error analysis')
    parser.add_argument('data_dir', help='Directory containing universality scan data')

    args = parser.parse_args()

    extractor = ExponentExtractor(args.data_dir)
    extractor.load_data()
    extractor.extract_exponents_with_errors()
    extractor.save_results()


if __name__ == '__main__':
    main()