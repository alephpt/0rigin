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
    std::cout << "  config/three_generations.yaml   - B3 fermion generation structure" << std::endl;
    std::cout << "  config/electroweak.yaml         - B4 W/Z boson masses" << std::endl;
    std::cout << "  config/strong_force.yaml        - B5 QCD emergence and confinement" << std::endl;
    std::cout << "  config/higgs_connection.yaml    - B6 Higgs mechanism via R-field" << std::endl;
    std::cout << "  config/laboratory_scale.yaml    - D2 Laboratory-scale tests (BEC, atomic clocks, superfluid, decoherence)" << std::endl;
    std::cout << "  config/josephson_junction.yaml  - D2 AC/DC Josephson effects" << std::endl;
    std::cout << "  config/binary_merger.yaml       - D3 Gravitational wave emission" << std::endl;
    std::cout << "  config/atomic_physics.yaml      - D5 Atomic physics and spectroscopy" << std::endl;
    std::cout << "  config/spin_magnetism.yaml      - H3 Spin-magnetism connection" << std::endl;
    std::cout << "  config/knot_topology.yaml       - H1 Topological excitations (knots)" << std::endl;
    std::cout << "  config/multiscale.yaml          - F2 Multi-scale RG flow validation" << std::endl;
    std::cout << "  config/finite_temperature.yaml  - F3 Finite temperature effects" << std::endl;
    std::cout << "  config/hpc_scaling.yaml         - F5 HPC scaling (OpenMP parallelization)" << std::endl;
    std::cout << "  config/renormalizability.yaml   - E1 Renormalizability analysis" << std::endl;
    std::cout << "  config/causality.yaml           - E3 Causality validation" << std::endl;
    std::cout << "  config/symmetry_analysis.yaml   - E5 Symmetry analysis (Noether)" << std::endl;
    std::cout << "  config/fine_structure_constant.yaml - B2 Fine structure constant (α ≈ 1/137)" << std::endl;

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

// E1/E3/E5 Wave 1 Mathematical Rigor tests
int runRenormalizabilityTest();
int runCausalityTest();
int runInflationTest();
// int runExperimentalPredictionsTest();  // Compiled as separate executable

// B3-B6 Standard Model tests
int runThreeGenerationsTest();
int runElectroweakTest();
int runStrongForceTest();
int runHiggsConnectionTest();

// H3 Spin-Magnetism test
int runSpinMagnetismTest();

// H1 Knot Topology test
int runKnotTopologyTest();

// D4: Particle scattering test
int runParticleScatteringTest();

// D2 Hardware experimental tests
int runJosephsonJunctionTest();

// D3 Gravitational wave tests
int runBinaryMergerTest();

// F2: Multi-Scale Validation test
int runMultiScaleTest();

// F4: Quantum Fluctuation Incorporation test
int runQuantumFluctuationsTest();

// F3: Finite Temperature Effects test
int runFiniteTemperatureTest();

// F5: High-Performance Computing Scaling test
int runHPCScalingTest();

// B2: Fine Structure Constant test
int runFineStructureConstantTest();

// D4: LHC Predictions test
int runLHCPredictionsTest();

// D3: Astrophysical Observations test
int runAstrophysicalObservationsTest();

// D5: Atomic Physics and Precision Spectroscopy test
int runAtomicPhysicsTest();

// D2: Laboratory-Scale Tests
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
    // } else if (config_path.find("experimental_predictions") != std::string::npos) {
    //     return runExperimentalPredictionsTest();  // Compiled as separate executable
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
