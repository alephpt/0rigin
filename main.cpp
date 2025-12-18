/**
 * main.cpp - Unified SMFT Application Entry Point
 *
 * Usage:
 *   ./smft                              - Run interactive visualization
 *   ./smft --test <config.yaml>         - Run test simulation
 *   ./smft --help                       - Show help
 */

#include "SMFTCore.h"
#include "simulations/SMFTTestRunner.h"
#include <iostream>
#include <string>
#include <cstring>

void printUsage(const char* program_name) {
    std::cout << "\n===== SMFT - Synchronization Mass Field Theory =====" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << "                         - Run interactive visualization" << std::endl;
    std::cout << "  " << program_name << " --test <config.yaml>    - Run test simulation" << std::endl;
    std::cout << "  " << program_name << " --help                  - Show this help message" << std::endl;

    std::cout << "\nModes:" << std::endl;
    std::cout << "  Interactive Mode:" << std::endl;
    std::cout << "    Launches the full SMFT visualization with GPU rendering" << std::endl;
    std::cout << "    Supports real-time parameter tuning and visual analysis" << std::endl;

    std::cout << "\n  Test Mode:" << std::endl;
    std::cout << "    Runs automated tests with quantitative validation" << std::endl;
    std::cout << "    Configuration-driven via YAML files" << std::endl;
    std::cout << "    Validates norm conservation, energy conservation, convergence" << std::endl;
    std::cout << "    Generates CSV output and test reports" << std::endl;

    std::cout << "\nTest Configuration Files:" << std::endl;
    std::cout << "  config/timesync_validation.yaml - Full validation (N=1,10,100, 64x64 grid)" << std::endl;
    std::cout << "  config/quick_validation.yaml    - Quick test (N=1,10, 32x32 grid)" << std::endl;

    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << program_name << std::endl;
    std::cout << "      → Launch interactive visualization" << std::endl;
    std::cout << "\n  " << program_name << " --test config/timesync_validation.yaml" << std::endl;
    std::cout << "      → Run timesync validation test" << std::endl;
    std::cout << "\n  " << program_name << " --test config/quick_validation.yaml" << std::endl;
    std::cout << "      → Run quick validation test" << std::endl;

    std::cout << "\n==================================================\n" << std::endl;
}

int runTestMode(const std::string& config_path) {
    std::cout << "\n===== SMFT Test Mode =====" << std::endl;
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

int runInteractiveMode() {
    std::cout << "\n===== SMFT Interactive Mode =====" << std::endl;
    std::cout << "Launching visualization...\n" << std::endl;

    SMFTCore* creation = SMFTCore::manifest();
    creation->actualize();
    delete creation;

    return 0;
}

int main(int argc, char* argv[]) {
    // No arguments - interactive mode
    if (argc == 1) {
        return runInteractiveMode();
    }

    // Parse arguments
    std::string arg1 = argv[1];

    // Help mode
    if (arg1 == "--help" || arg1 == "-h") {
        printUsage(argv[0]);
        return 0;
    }

    // Test mode
    if (arg1 == "--test" || arg1 == "-t") {
        if (argc < 3) {
            std::cerr << "Error: --test requires a configuration file" << std::endl;
            std::cerr << "Usage: " << argv[0] << " --test <config.yaml>" << std::endl;
            return 1;
        }
        std::string config_path = argv[2];
        return runTestMode(config_path);
    }

    // Unknown argument
    std::cerr << "Error: Unknown argument '" << arg1 << "'" << std::endl;
    printUsage(argv[0]);
    return 1;
}
