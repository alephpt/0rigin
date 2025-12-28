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

// Forward declare Validation namespace and types
namespace Validation {
    struct ValidationReport;
}

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

    // EM field observables (Phase 5): substep_ratio -> vector of EM observables
    std::map<int, std::vector<ObservableComputer::EMObservables>> _em_results;

    // Previous theta field for temporal derivatives (needed for EM field computation)
    std::vector<float> _theta_previous;

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

    // Validation framework reports (stored for report generation)
    std::map<int, Validation::ValidationReport> _global_validation_reports;
    std::map<int, Validation::ValidationReport> _scenario_validation_reports;
    std::string _detected_scenario_type;

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
     * Save EM field observables to CSV file
     * Format: timestep,field_energy,max_E,max_B,total_flux,avg_lorentz_force,charge_rms,current_rms,maxwell_violation
     */
    void saveEMObservablesToCSV(int N, const std::string& filepath) const;

    /**
     * Generate plots from observables using PlottingManager
     * @param N Substep ratio for this test run
     */
    void generatePlots(int N) const;

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
     * Get the base output directory for the current test run.
     * @return Full path to the timestamped output directory.
     */
    std::string getOutputDirectory() const;

    /**
     * Run the dispersion relation analysis.
     * @return true if successful
     */
    bool runDispersionAnalysis();

    /**
     * Save spatial field snapshots (theta and R fields)
     * @param N Substep ratio
     * @param step Current timestep
     */
    void saveSpatialFieldSnapshot(int N, int step) const;

    /**
     * Detect scenario type from test configuration
     * @return Scenario identifier ("2.1", "2.2", "2.3", etc.) or empty string
     */
    std::string detectScenarioType() const;

    /**
     * Analyze velocity sweep results and find critical velocity
     * @param velocities Vector of tested velocities
     * @param grid_size Grid size used for the tests
     */
    void analyzeVelocitySweep(const std::vector<float>& velocities, int grid_size) const;

    /**
     * Find critical velocity where error exceeds threshold
     * @param velocities Vector of velocities tested
     * @param errors Vector of corresponding errors
     * @param threshold Error threshold (default 5%)
     * @return Interpolated critical velocity
     */
    float findCriticalVelocity(const std::vector<float>& velocities,
                              const std::vector<double>& errors,
                              double threshold = 5.0) const;

    /**
     * Save velocity sweep report with aggregate analysis
     * @param velocities Vector of velocities tested
     * @param errors Vector of corresponding errors
     * @param v_critical Critical velocity found
     * @param grid_size Grid size used for tests
     */
    void saveVelocitySweepReport(const std::vector<float>& velocities,
                                const std::vector<double>& errors,
                                float v_critical,
                                int grid_size) const;

    /**
     * Read observables from CSV file for analysis
     * @param filepath Path to observables.csv file
     * @return Vector of observables (empty if file not found)
     */
    std::vector<ObservableComputer::Observables> readObservablesFromFile(
        const std::string& filepath) const;

    /**
     * Phase 3 Analysis: Casimir Force (Test 3.1)
     * Analyze force vs separation data and fit power law
     */
    void analyzeCasimirForce() const;

    /**
     * Phase 3 Analysis: Vacuum Energy (Test 3.2)
     * Analyze energy density vs order parameter and fit power law
     */
    void analyzeVacuumEnergy() const;

    /**
     * Phase 3 Analysis: Phase Transition (Test 3.3)
     * Find critical noise and fit critical exponent
     */
    void analyzePhaseTransition() const;

    // Store velocity sweep results for aggregate analysis
    mutable std::map<float, ValidationResult> _velocity_results;
    mutable std::map<float, std::string> _velocity_output_dirs;

    // Store Phase 3 analysis data
    struct CasimirRun {
        double defect_separation;
        double average_force;
    };
    struct VacuumEnergyRun {
        double R_average;
        double energy_density;
    };
    struct PhaseTransitionRun {
        double noise_sigma;
        double R_equilibrium;
    };
    mutable std::vector<CasimirRun> _casimir_runs;
    mutable std::vector<VacuumEnergyRun> _vacuum_energy_runs;
    mutable std::vector<PhaseTransitionRun> _phase_transition_runs;
};

#endif // SMFT_TEST_RUNNER_H
