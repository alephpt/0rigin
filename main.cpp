/**
 * main.cpp - Unified TRD Application Entry Point
 *
 * Usage:
 *   ./smft                              - Run interactive visualization
 *   ./smft --test <config.yaml>         - Run test simulation
 *   ./smft --help                       - Show help
 */

#include "TRDCore.h"
#include "simulations/TRDTestRunner.h"
#include <iostream>
#include <string>
#include <cstring>

void printUsage(const char* program_name) {
    std::cout << "\n===== TRD - Topological Resonance Dynamics =====" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << "                         - Run interactive visualization" << std::endl;
    std::cout << "  " << program_name << " --test <config.yaml>    - Run test simulation" << std::endl;
    std::cout << "  " << program_name << " --help                  - Show this help message" << std::endl;

    std::cout << "\nModes:" << std::endl;
    std::cout << "  Interactive Mode:" << std::endl;
    std::cout << "    Launches the full TRD visualization with GPU rendering" << std::endl;
    std::cout << "    Supports real-time parameter tuning and visual analysis" << std::endl;

    std::cout << "\n  Test Mode:" << std::endl;
    std::cout << "    Runs automated tests with quantitative validation" << std::endl;
    std::cout << "    Configuration-driven via YAML files" << std::endl;
    std::cout << "    Validates norm conservation, energy conservation, convergence" << std::endl;
    std::cout << "    Generates CSV output and test reports" << std::endl;

    std::cout << "\nTest Configuration Files:" << std::endl;
    std::cout << "  config/timesync_validation.yaml - Full validation (N=1,10,100, 64x64 grid)" << std::endl;
    std::cout << "  config/quick_validation.yaml    - Quick test (N=1,10, 32x32 grid)" << std::endl;
    std::cout << "  config/dark_matter.yaml         - C3 dark matter prediction test" << std::endl;
    std::cout << "  config/dark_energy.yaml         - C4 dark energy mechanism test" << std::endl;
    std::cout << "  config/inflation.yaml           - C5 primordial inflation test" << std::endl;
    std::cout << "  config/weak_field_3d.yaml       - Weak field gravity validation" << std::endl;

    std::cout << "\nExamples:" << std::endl;
    std::cout << "  " << program_name << std::endl;
    std::cout << "      → Launch interactive visualization" << std::endl;
    std::cout << "\n  " << program_name << " --test config/dark_matter.yaml" << std::endl;
    std::cout << "      → Test if TRD explains flat galaxy rotation curves" << std::endl;
    std::cout << "\n  " << program_name << " --test config/dark_energy.yaml" << std::endl;
    std::cout << "      → Test dark energy mechanism (accelerating expansion)" << std::endl;
    std::cout << "\n  " << program_name << " --test config/inflation.yaml" << std::endl;
    std::cout << "      → Validate primordial inflation (e-foldings, spectral index)" << std::endl;
    std::cout << "\n  " << program_name << " --test config/weak_field_3d.yaml" << std::endl;
    std::cout << "      → Validate weak field gravity (Newton's law)" << std::endl;

    std::cout << "\n==================================================\n" << std::endl;
}

// Forward declarations for 3D particle physics tests
int runLorentzForce3DTest();
int runStuckelbergVortex3DTest();
int runGeodesic3DTest();
int runWeakField3DTest();
int runThreeBodyEM3DTest();
int runEMGravityCoupling3DTest();
int runEinsteinFieldEquationsTest();
int runLightDeflection3DTest();
int runParticleSpectrumUnifiedTest();
int runTimeDilation3DTest();
int runCosmologicalConstantTest();
int runFriedmannEquationsTest();
int runDarkMatterTest();
int runUnitarityTest();
int runScaleInvarianceTest();
int runSymmetryAnalysisTest();
int runDarkEnergyTest();
int runInflationTest();
// int runExperimentalPredictionsTest();  // Compiled as separate executable

int runTestMode(const std::string& config_path) {
    std::cout << "\n===== TRD Test Mode =====" << std::endl;
    std::cout << "Configuration: " << config_path << std::endl;

    // Detect test type from config path
    std::string test_type;
    if (config_path.find("lorentz_force_3d") != std::string::npos) {
        return runLorentzForce3DTest();
    } else if (config_path.find("stuckelberg_vortex_3d") != std::string::npos) {
        return runStuckelbergVortex3DTest();
    } else if (config_path.find("geodesic_3d") != std::string::npos) {
        return runGeodesic3DTest();
    } else if (config_path.find("weak_field_3d") != std::string::npos) {
        return runWeakField3DTest();
    } else if (config_path.find("time_dilation_3d") != std::string::npos) {
        return runTimeDilation3DTest();
    } else if (config_path.find("three_body_em_3d") != std::string::npos) {
        return runThreeBodyEM3DTest();
    } else if (config_path.find("em_gravity_coupling_3d") != std::string::npos) {
        return runEMGravityCoupling3DTest();
    } else if (config_path.find("einstein_field_equations") != std::string::npos) {
        return runEinsteinFieldEquationsTest();
    } else if (config_path.find("light_deflection_3d") != std::string::npos) {
        return runLightDeflection3DTest();
    } else if (config_path.find("particle_spectrum") != std::string::npos) {
        return runParticleSpectrumUnifiedTest();
    } else if (config_path.find("cosmological_constant") != std::string::npos) {
        return runCosmologicalConstantTest();
    } else if (config_path.find("friedmann_equations") != std::string::npos) {
        return runFriedmannEquationsTest();
    } else if (config_path.find("dark_matter") != std::string::npos) {
        return runDarkMatterTest();
    } else if (config_path.find("dark_energy") != std::string::npos) {
        return runDarkEnergyTest();
    } else if (config_path.find("inflation") != std::string::npos) {
        return runInflationTest();
    } else if (config_path.find("unitarity") != std::string::npos) {
        return runUnitarityTest();
    } else if (config_path.find("scale_invariance") != std::string::npos) {
        return runScaleInvarianceTest();
    } else if (config_path.find("symmetry_analysis") != std::string::npos) {
        return runSymmetryAnalysisTest();
    // } else if (config_path.find("experimental_predictions") != std::string::npos) {
    //     return runExperimentalPredictionsTest();  // Compiled as separate executable
    }

    // Default: TRD field theory test (timesync, etc.)
    // Create test runner
    TRDTestRunner runner(config_path);

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
    std::cout << "\n===== TRD Interactive Mode =====" << std::endl;
    std::cout << "Interactive visualization mode is currently under development.\n" << std::endl;
    std::cout << "GPU-accelerated Kuramoto dynamics are fully restored!" << std::endl;
    std::cout << "Use --test mode to run physics simulations with test framework.\n" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  ./smft --test config/timesync_validation.yaml" << std::endl;
    std::cout << "  ./smft --test config/quick_validation.yaml" << std::endl;

    // TODO: Integrate Nova visualization once TRDEngine is connected
    // This will show real-time Kuramoto phase evolution + Dirac wavepacket motion

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
