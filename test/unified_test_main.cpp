/**
 * unified_test_main.cpp
 *
 * Unified test framework entry point for SMFT
 *
 * Usage:
 *   ./smft_test <config_file.yaml>
 *   ./smft_test --help
 *
 * Example:
 *   ./smft_test ../config/timesync_validation.yaml
 *   ./smft_test ../config/quick_validation.yaml
 */

#include "SMFTTestRunner.h"
#include <iostream>
#include <string>

void printUsage(const char* program_name) {
    std::cout << "\n===== SMFT Unified Test Framework =====" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << " <config_file.yaml>" << std::endl;
    std::cout << "  " << program_name << " --help" << std::endl;
    std::cout << "\nDescription:" << std::endl;
    std::cout << "  Runs SMFT tests based on YAML configuration" << std::endl;
    std::cout << "  Supports operator splitting with multiple N ratios" << std::endl;
    std::cout << "  Performs quantitative validation (norm, energy, convergence)" << std::endl;
    std::cout << "  Generates comprehensive test reports" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << program_name << " ../config/timesync_validation.yaml" << std::endl;
    std::cout << "  " << program_name << " ../config/quick_validation.yaml" << std::endl;
    std::cout << "\nConfiguration Files:" << std::endl;
    std::cout << "  timesync_validation.yaml - Full validation with N=1,10,100" << std::endl;
    std::cout << "  quick_validation.yaml    - Quick test with N=1,10" << std::endl;
    std::cout << "\nOutput:" << std::endl;
    std::cout << "  Test results saved to output directory specified in config" << std::endl;
    std::cout << "  Observables saved as CSV files" << std::endl;
    std::cout << "  Test report with pass/fail status generated" << std::endl;
    std::cout << "\n======================================\n" << std::endl;
}

int main(int argc, char** argv) {
    // Check arguments
    if (argc < 2) {
        std::cerr << "Error: No configuration file specified" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    std::string arg = argv[1];

    // Handle help
    if (arg == "--help" || arg == "-h") {
        printUsage(argv[0]);
        return 0;
    }

    // Load configuration
    std::string config_path = arg;

    std::cout << "\n===== SMFT Unified Test Framework =====" << std::endl;
    std::cout << "Configuration: " << config_path << std::endl;

    // Create test runner
    SMFTTestRunner runner(config_path);

    // Initialize
    std::cout << "\n[1/3] Initializing..." << std::endl;
    if (!runner.initialize()) {
        std::cerr << "\n✗ Initialization failed" << std::endl;
        return 1;
    }

    // Run tests
    std::cout << "\n[2/3] Running tests..." << std::endl;
    bool success = runner.run();

    // Generate report
    std::cout << "\n[3/3] Generating report..." << std::endl;
    runner.generateReport();

    // Summary
    std::cout << "\n===== Test Summary =====" << std::endl;
    if (success && runner.allTestsPassed()) {
        std::cout << "✓ ALL TESTS PASSED" << std::endl;
        std::cout << "\nValidation criteria met:" << std::endl;
        std::cout << "  ✓ Norm conservation ||Ψ||² - 1 < 10⁻⁴" << std::endl;
        std::cout << "  ✓ Energy conservation |ΔE/E₀| < 1%" << std::endl;
        std::cout << "  ✓ Convergence between N ratios < 5%" << std::endl;
        return 0;
    } else {
        std::cout << "✗ TESTS FAILED" << std::endl;
        std::cout << "\nSome validation criteria not met." << std::endl;
        std::cout << "Check test report for details." << std::endl;
        return 1;
    }
}
