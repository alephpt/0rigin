#ifndef SMFT_TEST_RUNNER_H
#define SMFT_TEST_RUNNER_H

#include "simulations/TestConfig.h"
#include "simulations/ObservableComputer.h"
#include "SMFTEngine.h"
#include "DiracEvolution.h"
#include "Nova.h"
#include <vector>
#include <string>
#include <map>
#include <fstream>

// Forward declarations
class DiracEvolution;

// Forward declare OutputManager in SMFT namespace (no longer conflicts since class renamed to SMFTCore)
namespace SMFT { class OutputManager; }

/**
 * SMFTTestRunner - Unified test execution framework for SMFT
 *
 * Features:
 * - YAML-driven configuration
 * - Multiple substep ratio testing (N=1, 10, 100, etc.)
 * - Quantitative validation (norm, energy, convergence)
 * - Automatic output generation via OutputManager
 * - Pass/fail reporting with detailed metrics
 *
 * Usage:
 *   SMFTTestRunner runner("config/timesync_validation.yaml");
 *   runner.initialize();
 *   runner.run();
 *   runner.generateReport();
 */
class SMFTTestRunner {
public:
    /**
     * Constructor - load configuration from YAML
     * @param config_path Path to YAML configuration file
     */
    SMFTTestRunner(const std::string& config_path);

    /**
     * Constructor - use programmatic configuration
     * @param config TestConfig object
     */
    SMFTTestRunner(const TestConfig& config);

    /**
     * Destructor
     */
    ~SMFTTestRunner();

    /**
     * Initialize Nova and SMFT engine
     * @return true if successful
     */
    bool initialize();

    /**
     * Run all tests specified in configuration
     * If operator splitting enabled, runs once for each substep ratio
     * @return true if all tests pass validation
     */
    bool run();

    /**
     * Generate comprehensive test report
     * Includes pass/fail status, quantitative metrics, convergence analysis
     * @param output_path Optional path for report (default: uses config output directory)
     */
    void generateReport(const std::string& output_path = "") const;

    /**
     * Get test results for specific substep ratio
     * @param N Substep ratio (1, 10, 100, etc.)
     * @return Vector of observables for each timestep
     */
    const std::vector<ObservableComputer::Observables>& getResults(int N) const;

    /**
     * Check if all tests passed validation
     * @return true if norm, energy, and convergence checks all pass
     */
    bool allTestsPassed() const;

private:
    // Configuration
    TestConfig _config;
    bool _config_loaded;

    // Nova graphics engine
    Nova* _nova;
    bool _nova_initialized;

    // SMFT engine
    SMFTEngine* _engine;
    bool _engine_initialized;

    // Test results: substep_ratio -> vector of observables
    std::map<int, std::vector<ObservableComputer::Observables>> _results;

    // Validation results
    struct ValidationResult {
        bool norm_passed = false;
        bool energy_passed = false;
        bool convergence_passed = false;

        double max_norm_error = 0.0;
        double max_energy_drift = 0.0;
        std::vector<double> convergence_errors; // Per observable

        bool passed() const {
            // For individual N tests: only check norm and energy
            // For convergence test (N=-1): check all three
            // If convergence_errors is populated, this is the convergence test
            if (!convergence_errors.empty()) {
                return convergence_passed;
            }
            // Individual N test: only check norm and energy
            return norm_passed && energy_passed;
        }
    };

    std::map<int, ValidationResult> _validation_results;

    // Output management
    SMFT::OutputManager* _output_manager;
    std::string _timestamped_output_dir;
    int _current_grid_size = 0;  // Track current grid size for output paths

    // Console logging to file
    std::ofstream _console_log;
    std::streambuf* _cout_backup;
    std::streambuf* _cerr_backup;
    class TeeStreambuf* _tee_cout;
    class TeeStreambuf* _tee_cerr;

    // Helper methods

    /**
     * Start logging console output to file
     * @param log_path Path to console log file
     */
    void startConsoleLogging(const std::string& log_path);

    /**
     * Stop logging console output and restore stdout/stderr
     */
    void stopConsoleLogging();

    /**
     * Run all tests for a specific grid size
     * @param grid_size Grid size to test
     * @return true if all tests passed
     */
    bool runForGridSize(int grid_size);

    /**
     * Run all tests for a specific grid size and velocity (Scenario 2.3)
     * @param grid_size Grid size to test
     * @param velocity Boost velocity (in units of c)
     * @return true if all tests passed
     */
    bool runForGridSizeAndVelocity(int grid_size, float velocity);

    /**
     * Run single test with given substep ratio
     * @param N Substep ratio
     * @return true if successful
     */
    bool runSingleTest(int N);

    /**
     * Initialize Kuramoto phases based on config
     * @return Vector of initial phases
     */
    std::vector<float> initializePhases() const;

    /**
     * Initialize Kuramoto natural frequencies based on config
     * @return Vector of natural frequencies
     */
    std::vector<float> initializeOmegas() const;

    /**
     * Validate norm conservation for given test results
     * @param N Substep ratio
     * @return Validation result
     */
    ValidationResult validateNorm(int N) const;

    /**
     * Validate energy conservation for given test results
     * @param N Substep ratio
     * @return Validation result
     */
    ValidationResult validateEnergy(int N) const;

    /**
     * Validate convergence between different N ratios
     * Compares observables at final timestep
     * @return Validation result
     */
    ValidationResult validateConvergence() const;

    /**
     * Save observables to CSV file
     * @param N Substep ratio
     * @param filepath Path to output file
     */
    void saveObservablesToCSV(int N, const std::string& filepath) const;

    /**
     * Create output directory structure
     */
    void createOutputDirectories();

    /**
     * Generate visualization plots automatically
     */
    void generateVisualizations() const;

    /**
     * Get output directory for specific substep ratio
     * @param N Substep ratio
     * @return Full path to output directory
     */
    std::string getOutputDirectory(int N) const;

    /**
     * Save spatial field snapshots (theta and R fields)
     * @param N Substep ratio
     * @param step Current timestep
     */
    void saveSpatialFieldSnapshot(int N, int step) const;
};

#endif // SMFT_TEST_RUNNER_H
