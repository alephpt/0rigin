/**
 * FSSIntegration.cpp
 *
 * Implementation of FSS integration utilities
 */

#include "FSSIntegration.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>

// Extract FSS config from TestConfig
FSSIntegration::FSSTestConfig FSSIntegration::extractFSSConfig(const TestConfig& config) {
    FSSTestConfig fss_config;

    // Get system sizes
    fss_config.system_sizes = getGridSizes(config);

    // Get noise values
    fss_config.noise_values = getNoiseScanValues(config);

    // Get timing parameters
    fss_config.equilibration_steps = config.physics.total_steps * 0.7;  // Use 70% for equilibration
    fss_config.measurement_steps = config.physics.total_steps * 0.3;    // Use 30% for measurement

    // Output settings
    fss_config.save_snapshots = config.output.save_spatial_snapshots;
    fss_config.output_dir = config.output.directory;

    return fss_config;
}

// Add run result to FSS analyzer
void FSSIntegration::addRunToFSS(FiniteSizeScaling& fss, const FSSRunResult& result) {
    if (!result.R_field_snapshot.empty()) {
        // If we have the full field, use it
        fss.addDataPoint(result.L, result.sigma, result.R_field_snapshot);
    } else {
        // Otherwise use pre-computed values
        FiniteSizeScaling::FSSDataPoint point;
        point.L = result.L;
        point.sigma = result.sigma;
        point.R_mean = result.R_mean;
        point.R_variance = result.R_variance;

        // Compute moments from Binder and susceptibility if available
        if (result.susceptibility > 0) {
            point.R_squared = result.R_variance + result.R_mean * result.R_mean;
        }

        // Estimate R^4 from Binder cumulant: U_L = 1 - R^4/(3*R^2*R^2)
        if (result.binder_cumulant > -1 && point.R_squared > 0) {
            point.R_fourth = 3 * point.R_squared * point.R_squared * (1 - result.binder_cumulant);
        }

        point.num_samples = 1;
        fss.addDataPoint(point);
    }
}

// Generate FSS report
bool FSSIntegration::generateFSSReport(const FiniteSizeScaling& fss,
                                        const std::string& output_dir) {
    // Generate main report
    if (!fss.writeDataCollapseReport(output_dir)) {
        return false;
    }

    // Write data for plotting
    if (!fss.writeDataForPlotting(output_dir)) {
        std::cerr << "Warning: Failed to write plotting data" << std::endl;
    }

    // Generate visualization script
    generateVisualizationScript(output_dir);

    return true;
}

// Get grid sizes from config
std::vector<int> FSSIntegration::getGridSizes(const TestConfig& config) {
    std::vector<int> sizes;

    // Check if config has grid_sizes array (would need to be added to TestConfig)
    // For now, just use the single grid size
    sizes.push_back(config.grid.size_x);

    // If this is an FSS test, we might encode sizes in the test name
    if (config.test_name.find("L32") != std::string::npos) {
        return {32};
    } else if (config.test_name.find("L64") != std::string::npos) {
        return {64};
    } else if (config.test_name.find("L128") != std::string::npos) {
        return {128};
    } else if (config.test_name.find("L256") != std::string::npos) {
        return {256};
    } else if (config.test_name.find("L512") != std::string::npos) {
        return {512};
    }

    return sizes;
}

// Get noise scan values
std::vector<float> FSSIntegration::getNoiseScanValues(const TestConfig& config) {
    std::vector<float> noise_values;

    // Convert from double to float
    for (double sigma : config.physics.noise_scan) {
        noise_values.push_back(static_cast<float>(sigma));
    }

    return noise_values;
}

// Generate test matrix for all (L, sigma) pairs
std::vector<std::pair<TestConfig, std::pair<int, float>>>
FSSIntegration::generateFSSTestMatrix(const TestConfig& base_config) {
    std::vector<std::pair<TestConfig, std::pair<int, float>>> test_matrix;

    std::vector<int> sizes = getGridSizes(base_config);
    std::vector<float> noise_values = getNoiseScanValues(base_config);

    for (int L : sizes) {
        for (float sigma : noise_values) {
            TestConfig test = base_config;
            test.grid.size_x = L;
            test.grid.size_y = L;
            test.physics.noise_strength = sigma;

            test_matrix.push_back({test, {L, sigma}});
        }
    }

    return test_matrix;
}

// Compute FSS observables from R field
FSSIntegration::FSSRunResult FSSIntegration::computeFSSObservables(
    const std::vector<float>& R_field, int L, float sigma) {

    FSSRunResult result;
    result.L = L;
    result.sigma = sigma;
    result.R_field_snapshot = R_field;

    if (R_field.empty()) {
        result.converged = false;
        return result;
    }

    // Compute mean
    float sum = std::accumulate(R_field.begin(), R_field.end(), 0.0f);
    result.R_mean = sum / R_field.size();

    // Compute variance
    float sum_sq = 0.0f;
    for (float R : R_field) {
        sum_sq += (R - result.R_mean) * (R - result.R_mean);
    }
    result.R_variance = sum_sq / R_field.size();

    // Compute Binder cumulant
    result.binder_cumulant = FiniteSizeScaling::computeBinderCumulant(R_field);

    // Compute susceptibility
    result.susceptibility = FiniteSizeScaling::computeSusceptibility(R_field, L);

    result.converged = true;
    return result;
}

// Run complete FSS campaign (placeholder - would be integrated with SMFTTestRunner)
bool FSSIntegration::runFSSCampaign(const TestConfig& config,
                                     const std::string& output_dir) {
    std::cout << "\n===== Starting Finite-Size Scaling Campaign =====" << std::endl;

    FSSTestConfig fss_config = extractFSSConfig(config);

    std::cout << "System sizes: ";
    for (int L : fss_config.system_sizes) {
        std::cout << L << " ";
    }
    std::cout << std::endl;

    std::cout << "Noise points: " << fss_config.noise_values.size() << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;

    // This function would be called from SMFTTestRunner
    // The actual simulation runs would happen there
    // This is just the structure

    FiniteSizeScaling fss;

    // Placeholder for collecting results
    // In reality, SMFTTestRunner would call addRunToFSS after each simulation

    // Generate report
    return generateFSSReport(fss, output_dir);
}

// Validate universality class
bool FSSIntegration::validateUniversalityClass(
    const FiniteSizeScaling::FSSResult& result,
    const std::string& expected_class,
    float tolerance) {

    if (result.universality_class != expected_class) {
        std::cout << "Universality class mismatch: " << std::endl;
        std::cout << "  Expected: " << expected_class << std::endl;
        std::cout << "  Measured: " << result.universality_class << std::endl;
        return false;
    }

    // Check confidence level
    if (result.confidence_level < (1.0f - tolerance) * 100.0f) {
        std::cout << "Low confidence in universality class: "
                  << result.confidence_level << "%" << std::endl;
        return false;
    }

    std::cout << "✓ Universality class confirmed: " << expected_class
              << " (confidence: " << result.confidence_level << "%)" << std::endl;
    return true;
}

// Generate visualization script
bool FSSIntegration::generateVisualizationScript(const std::string& output_dir) {
    std::string script_path = output_dir + "/visualize_fss.py";
    std::ofstream file(script_path);

    if (!file.is_open()) {
        return false;
    }

    file << "#!/usr/bin/env python3\n";
    file << "\"\"\"\n";
    file << "Finite-Size Scaling Visualization\n";
    file << "Generates plots for FSS analysis results\n";
    file << "\"\"\"\n\n";

    file << "import numpy as np\n";
    file << "import matplotlib.pyplot as plt\n";
    file << "import pandas as pd\n";
    file << "from glob import glob\n";
    file << "import sys\n\n";

    file << "# Load data files\n";
    file << "data_files = glob('fss_data_L*.csv')\n";
    file << "collapse_file = 'fss_data_collapse.csv'\n\n";

    file << "# Create figure with subplots\n";
    file << "fig, axes = plt.subplots(2, 3, figsize=(15, 10))\n\n";

    file << "# Plot 1: Order parameter vs noise\n";
    file << "ax1 = axes[0, 0]\n";
    file << "for f in sorted(data_files):\n";
    file << "    L = int(f.split('L')[1].split('.')[0])\n";
    file << "    df = pd.read_csv(f)\n";
    file << "    ax1.plot(df['sigma'], df['R_mean'], 'o-', label=f'L={L}')\n";
    file << "ax1.set_xlabel(r'$\\sigma$')\n";
    file << "ax1.set_ylabel(r'$\\langle R \\rangle$')\n";
    file << "ax1.set_title('Order Parameter')\n";
    file << "ax1.legend()\n";
    file << "ax1.grid(True)\n\n";

    file << "# Plot 2: Binder cumulant\n";
    file << "ax2 = axes[0, 1]\n";
    file << "for f in sorted(data_files):\n";
    file << "    L = int(f.split('L')[1].split('.')[0])\n";
    file << "    df = pd.read_csv(f)\n";
    file << "    ax2.plot(df['sigma'], df['Binder_cumulant'], 'o-', label=f'L={L}')\n";
    file << "ax2.set_xlabel(r'$\\sigma$')\n";
    file << "ax2.set_ylabel(r'$U_L$')\n";
    file << "ax2.set_title('Binder Cumulant')\n";
    file << "ax2.legend()\n";
    file << "ax2.grid(True)\n\n";

    file << "# Plot 3: Susceptibility\n";
    file << "ax3 = axes[0, 2]\n";
    file << "for f in sorted(data_files):\n";
    file << "    L = int(f.split('L')[1].split('.')[0])\n";
    file << "    df = pd.read_csv(f)\n";
    file << "    ax3.plot(df['sigma'], df['susceptibility']/L**2, 'o-', label=f'L={L}')\n";
    file << "ax3.set_xlabel(r'$\\sigma$')\n";
    file << "ax3.set_ylabel(r'$\\chi/L^2$')\n";
    file << "ax3.set_title('Susceptibility (scaled)')\n";
    file << "ax3.legend()\n";
    file << "ax3.grid(True)\n\n";

    file << "# Plot 4: Data collapse\n";
    file << "ax4 = axes[1, 0]\n";
    file << "if glob(collapse_file):\n";
    file << "    df_collapse = pd.read_csv(collapse_file)\n";
    file << "    for L in df_collapse['L'].unique():\n";
    file << "        df_L = df_collapse[df_collapse['L'] == L]\n";
    file << "        ax4.plot(df_L['X_scaled'], df_L['Y_scaled'], 'o', label=f'L={L}')\n";
    file << "ax4.set_xlabel(r'$(\\sigma - \\sigma_c)L^{1/\\nu}$')\n";
    file << "ax4.set_ylabel(r'$\\langle R \\rangle L^{\\beta/\\nu}$')\n";
    file << "ax4.set_title('FSS Data Collapse')\n";
    file << "ax4.legend()\n";
    file << "ax4.grid(True)\n\n";

    file << "# Plot 5: Log-log for critical exponent\n";
    file << "ax5 = axes[1, 1]\n";
    file << "# Would plot log-log fit here\n";
    file << "ax5.set_xlabel(r'$\\log(\\sigma_c - \\sigma)$')\n";
    file << "ax5.set_ylabel(r'$\\log\\langle R \\rangle$')\n";
    file << "ax5.set_title('Critical Exponent Fit')\n";
    file << "ax5.grid(True)\n\n";

    file << "# Plot 6: Finite-size effects\n";
    file << "ax6 = axes[1, 2]\n";
    file << "# Would plot L-dependence of observables at sigma_c\n";
    file << "ax6.set_xlabel(r'$\\log L$')\n";
    file << "ax6.set_ylabel(r'$\\log$ Observable')\n";
    file << "ax6.set_title('Finite-Size Scaling')\n";
    file << "ax6.grid(True)\n\n";

    file << "plt.suptitle('Finite-Size Scaling Analysis', fontsize=16)\n";
    file << "plt.tight_layout()\n";
    file << "plt.savefig('fss_analysis.png', dpi=150)\n";
    file << "plt.show()\n";

    file.close();

    // Make script executable
    std::string chmod_cmd = "chmod +x " + script_path;
    system(chmod_cmd.c_str());

    return true;
}

// Helper: Interpolate Binder cumulant
float FSSIntegration::interpolateBinder(const std::vector<float>& sigma,
                                        const std::vector<float>& binder,
                                        float sigma_target) {
    if (sigma.size() != binder.size() || sigma.empty()) {
        return 0.0f;
    }

    // Find bracketing points
    for (size_t i = 0; i < sigma.size() - 1; ++i) {
        if (sigma_target >= sigma[i] && sigma_target <= sigma[i+1]) {
            // Linear interpolation
            float t = (sigma_target - sigma[i]) / (sigma[i+1] - sigma[i]);
            return binder[i] + t * (binder[i+1] - binder[i]);
        }
    }

    // Extrapolate if outside range
    if (sigma_target < sigma.front()) {
        return binder.front();
    } else {
        return binder.back();
    }
}

// Helper: Find crossing point
bool FSSIntegration::findCrossingPoint(const std::vector<float>& x1,
                                       const std::vector<float>& y1,
                                       const std::vector<float>& x2,
                                       const std::vector<float>& y2,
                                       float& x_cross, float& y_cross) {
    // Simple implementation - find where difference changes sign
    float min_diff = std::numeric_limits<float>::max();
    bool found = false;

    for (size_t i = 0; i < x1.size(); ++i) {
        // Find corresponding point in second curve
        float y2_interp = 0.0f;
        bool interpolated = false;

        for (size_t j = 0; j < x2.size() - 1; ++j) {
            if (x1[i] >= x2[j] && x1[i] <= x2[j+1]) {
                float t = (x1[i] - x2[j]) / (x2[j+1] - x2[j]);
                y2_interp = y2[j] + t * (y2[j+1] - y2[j]);
                interpolated = true;
                break;
            }
        }

        if (interpolated) {
            float diff = std::abs(y1[i] - y2_interp);
            if (diff < min_diff) {
                min_diff = diff;
                x_cross = x1[i];
                y_cross = (y1[i] + y2_interp) / 2.0f;
                found = true;
            }
        }
    }

    return found && (min_diff < 0.01f);  // Require curves to get close
}