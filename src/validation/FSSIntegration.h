/**
 * FSSIntegration.h
 *
 * Integration utilities for connecting FiniteSizeScaling with SMFTTestRunner
 * Provides helper functions to run FSS analysis from test configurations
 */

#pragma once

#include "FiniteSizeScaling.h"
#include "PhaseTransitionAnalyzer.h"
#include "../simulations/TestConfig.h"
#include <vector>
#include <string>
#include <map>

class FSSIntegration {
public:
    /**
     * FSS test configuration
     */
    struct FSSTestConfig {
        std::vector<int> system_sizes;        // L values to test
        std::vector<float> noise_values;      // σ values to scan
        int equilibration_steps;              // Steps to wait for equilibration
        int measurement_steps;                 // Steps to average over
        bool save_snapshots;                  // Save R field snapshots
        std::string output_dir;               // Output directory
    };

    /**
     * FSS test result for a single run
     */
    struct FSSRunResult {
        int L;                                // System size
        float sigma;                          // Noise value
        float R_mean;                         // Mean order parameter
        float R_variance;                     // Variance
        float binder_cumulant;                // U_L
        float susceptibility;                 // χ
        std::vector<float> R_field_snapshot; // Final R field (if saved)
        bool converged;                       // Equilibration achieved
    };

    /**
     * Convert TestConfig to FSS config
     */
    static FSSTestConfig extractFSSConfig(const TestConfig& config);

    /**
     * Process FSS run results and add to analyzer
     */
    static void addRunToFSS(FiniteSizeScaling& fss,
                            const FSSRunResult& result);

    /**
     * Create FSS report from collected data
     */
    static bool generateFSSReport(const FiniteSizeScaling& fss,
                                  const std::string& output_dir);

    /**
     * Check if this is an FSS scenario
     */
    static bool isFSSScenario(const std::string& scenario) {
        return scenario == "fss_analysis" ||
               scenario == "finite_size_scaling" ||
               scenario == "phase_transition_fss";
    }

    /**
     * Extract grid sizes from config
     * Returns list of L values to test (from grid_sizes array or single size)
     */
    static std::vector<int> getGridSizes(const TestConfig& config);

    /**
     * Parse noise scan values from config
     */
    static std::vector<float> getNoiseScanValues(const TestConfig& config);

    /**
     * Create individual test configs for each (L, sigma) pair
     */
    static std::vector<std::pair<TestConfig, std::pair<int, float>>>
        generateFSSTestMatrix(const TestConfig& base_config);

    /**
     * Compute FSS observables from R field
     */
    static FSSRunResult computeFSSObservables(const std::vector<float>& R_field,
                                              int L, float sigma);

    /**
     * Run complete FSS analysis campaign
     * This would be called from SMFTTestRunner when scenario is "fss_analysis"
     */
    static bool runFSSCampaign(const TestConfig& config,
                               const std::string& output_dir);

    /**
     * Validate FSS results against expected universality classes
     */
    static bool validateUniversalityClass(const FiniteSizeScaling::FSSResult& result,
                                          const std::string& expected_class,
                                          float tolerance = 0.15f);

    /**
     * Generate Python script for FSS visualization
     */
    static bool generateVisualizationScript(const std::string& output_dir);

private:
    /**
     * Helper to interpolate Binder cumulant curves
     */
    static float interpolateBinder(const std::vector<float>& sigma,
                                   const std::vector<float>& binder,
                                   float sigma_target);

    /**
     * Helper to find crossing points between curves
     */
    static bool findCrossingPoint(const std::vector<float>& x1,
                                  const std::vector<float>& y1,
                                  const std::vector<float>& x2,
                                  const std::vector<float>& y2,
                                  float& x_cross, float& y_cross);
};