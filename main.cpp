/**
 * main.cpp - Unified TRD Application Entry Point
 *
 * Usage:
 *   ./trd                              - Run interactive visualization
 *   ./trd --test <config.yaml>         - Run test simulation
 *   ./trd --help                       - Show help
 */

#include "simulations/TRDTestRunner.h"
#include <iostream>
#include <string>

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

    std::cout << "\nTest Configuration Files (config/*.yaml):" << std::endl;
    std::cout << "  Particle physics:   particle_spectrum_unified, particle_scattering, three_generations" << std::endl;
    std::cout << "  Standard Model:     electroweak, strong_force, higgs_connection, fine_structure_constant" << std::endl;
    std::cout << "  Cosmology:          dark_matter, dark_energy, inflation, friedmann_equations" << std::endl;
    std::cout << "  General relativity: einstein_field_equations, gravitational_waves, geodesic_3d" << std::endl;
    std::cout << "  Electromagnetism:   lorentz_force_3d, em_gravity_coupling_3d, stuckelberg_vortex_3d" << std::endl;
    std::cout << "  Condensed matter:   josephson_junction, spin_magnetism, knot_topology" << std::endl;
    std::cout << "  Mathematical rigor: unitarity, renormalizability, causality, symmetry_analysis" << std::endl;
    std::cout << "  Astrophysics:       binary_merger, solar_system, astrophysical_observations" << std::endl;

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
int runParticleSpectrumTest(int argc, char* argv[]);
int runTimeDilation3DTest();
int runCosmologicalConstantTest();
int runFriedmannEquationsTest();
int runDarkMatterTest();
int runUnitarityTest();
int runScaleInvarianceTest();
int runSymmetryAnalysisTest();
int runDarkEnergyTest();
int runInflationTest();
int runRenormalizabilityTest();
int runCausalityTest();
int runThreeGenerationsTest();
int runElectroweakTest();
int runStrongForceTest();
int runHiggsConnectionTest();
int runSpinMagnetismTest();
int runKnotTopologyTest();
int runParticleScatteringTest();
int runJosephsonJunctionTest();
int runBinaryMergerTest();
int runMultiScaleTest();
int runQuantumFluctuationsTest();
int runFiniteTemperatureTest();
int runHPCScalingTest();
int runFineStructureConstantTest();
int runLHCPredictionsTest();
int runAstrophysicalObservationsTest();
int runAtomicPhysicsTest();
int runLaboratoryScaleTest();

int runTestMode(const std::string& config_path) {
    std::cout << "\n===== TRD Test Mode =====" << std::endl;
    std::cout << "Configuration: " << config_path << std::endl;

    // Set global config path for test functions
    extern std::string g_test_config_path;
    g_test_config_path = config_path;

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
        // Create argv array with config path for particle spectrum test
        char* ps_argv[2];
        ps_argv[0] = const_cast<char*>("trd");
        ps_argv[1] = const_cast<char*>(config_path.c_str());
        return runParticleSpectrumTest(2, ps_argv);
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
    } else if (config_path.find("renormalizability") != std::string::npos) {
        return runRenormalizabilityTest();
    } else if (config_path.find("causality") != std::string::npos) {
        return runCausalityTest();
    } else if (config_path.find("three_generations") != std::string::npos) {
        return runThreeGenerationsTest();
    } else if (config_path.find("electroweak") != std::string::npos) {
        return runElectroweakTest();
    } else if (config_path.find("strong_force") != std::string::npos) {
        return runStrongForceTest();
    } else if (config_path.find("higgs_connection") != std::string::npos) {
        return runHiggsConnectionTest();
    } else if (config_path.find("particle_scattering") != std::string::npos) {
        return runParticleScatteringTest();
    } else if (config_path.find("josephson_junction") != std::string::npos) {
        return runJosephsonJunctionTest();
    } else if (config_path.find("binary_merger") != std::string::npos) {
        return runBinaryMergerTest();
    } else if (config_path.find("spin_magnetism") != std::string::npos) {
        return runSpinMagnetismTest();
    } else if (config_path.find("knot_topology") != std::string::npos) {
        return runKnotTopologyTest();
    } else if (config_path.find("multiscale") != std::string::npos) {
        return runMultiScaleTest();
    } else if (config_path.find("quantum_fluctuations") != std::string::npos) {
        return runQuantumFluctuationsTest();
    } else if (config_path.find("finite_temperature") != std::string::npos) {
        return runFiniteTemperatureTest();
    } else if (config_path.find("hpc_scaling") != std::string::npos) {
        return runHPCScalingTest();
    } else if (config_path.find("fine_structure_constant") != std::string::npos) {
        return runFineStructureConstantTest();
    } else if (config_path.find("atomic_physics") != std::string::npos) {
        return runAtomicPhysicsTest();
    } else if (config_path.find("lhc_predictions") != std::string::npos) {
        return runLHCPredictionsTest();
    } else if (config_path.find("astrophysical_observations") != std::string::npos) {
        return runAstrophysicalObservationsTest();
    } else if (config_path.find("laboratory_scale") != std::string::npos) {
        return runLaboratoryScaleTest();
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

// Global configuration path for test functions
std::string g_test_config_path;

int runInteractiveMode() {
    std::cout << "\n===== TRD Interactive Mode =====" << std::endl;
    std::cout << "GPU-accelerated visualization is available when Nova is connected.\n" << std::endl;
    std::cout << "Use --test mode to run physics simulations with test framework.\n" << std::endl;
    std::cout << "\nExamples:" << std::endl;
    std::cout << "  ./trd --test config/josephson_junction.yaml" << std::endl;
    std::cout << "  ./trd --test config/unitarity.yaml" << std::endl;

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
