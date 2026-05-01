/**
 * main.cpp - Unified TRD Application Entry Point
 *
 * Usage:
 *   ./trd                              - Run interactive visualization
 *   ./trd --test <config.yaml>         - Run test simulation
 *   ./trd --help                       - Show help
 */

#include "simulations/TRDTestRunner.h"
#include "simulations/VisualizationGenerator.h"
#include <iostream>
#include <string>

void printUsage(const char* program_name) {
    std::cout << "\n===== TRD - Topological Resonance Dynamics =====" << std::endl;
    std::cout << "\nUsage:" << std::endl;
    std::cout << "  " << program_name << "                         - Run interactive visualization" << std::endl;
    std::cout << "  " << program_name << " --test <config.yaml>    - Run test simulation" << std::endl;
    std::cout << "  " << program_name << " --test <config.yaml> --plot  - Run test + generate plots" << std::endl;
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
int runChiralOperatorComparisonTest();
int runKuramotoPhaseDiagramTest();
int runChiralChannelSelectorTest();
int runNJLGapCurveTest();

int runTestMode(const std::string& config_path, bool generate_plots) {
    std::cout << "\n===== TRD Test Mode =====" << std::endl;
    std::cout << "Configuration: " << config_path << std::endl;
    if (generate_plots) {
        std::cout << "Plot generation: ENABLED" << std::endl;
    }

    // Set global config path for test functions
    extern std::string g_test_config_path;
    g_test_config_path = config_path;

    // Extract test name from config path for visualization
    std::string test_name = config_path;
    auto slash = test_name.rfind('/');
    if (slash != std::string::npos) test_name = test_name.substr(slash + 1);
    auto dot = test_name.rfind('.');
    if (dot != std::string::npos) test_name = test_name.substr(0, dot);

    // Clear any previous visualization data
    VisualizationGenerator::clearData();

    int result = 0;

    // Detect test type from config path
    std::string test_type;
    if (config_path.find("lorentz_force_3d") != std::string::npos) {
        result = runLorentzForce3DTest();
    } else if (config_path.find("stuckelberg_vortex_3d") != std::string::npos) {
        result = runStuckelbergVortex3DTest();
    } else if (config_path.find("geodesic_3d") != std::string::npos) {
        result = runGeodesic3DTest();
    } else if (config_path.find("weak_field_3d") != std::string::npos) {
        result = runWeakField3DTest();
    } else if (config_path.find("time_dilation_3d") != std::string::npos) {
        result = runTimeDilation3DTest();
    } else if (config_path.find("three_body_em_3d") != std::string::npos) {
        result = runThreeBodyEM3DTest();
    } else if (config_path.find("em_gravity_coupling_3d") != std::string::npos) {
        result = runEMGravityCoupling3DTest();
    } else if (config_path.find("einstein_field_equations") != std::string::npos) {
        result = runEinsteinFieldEquationsTest();
    } else if (config_path.find("light_deflection_3d") != std::string::npos) {
        result = runLightDeflection3DTest();
    } else if (config_path.find("particle_spectrum") != std::string::npos) {
        char* ps_argv[2];
        ps_argv[0] = const_cast<char*>("trd");
        ps_argv[1] = const_cast<char*>(config_path.c_str());
        result = runParticleSpectrumTest(2, ps_argv);
    } else if (config_path.find("cosmological_constant") != std::string::npos) {
        result = runCosmologicalConstantTest();
    } else if (config_path.find("friedmann_equations") != std::string::npos) {
        result = runFriedmannEquationsTest();
    } else if (config_path.find("dark_matter") != std::string::npos) {
        result = runDarkMatterTest();
    } else if (config_path.find("dark_energy") != std::string::npos) {
        result = runDarkEnergyTest();
    } else if (config_path.find("inflation") != std::string::npos) {
        result = runInflationTest();
    } else if (config_path.find("unitarity") != std::string::npos) {
        result = runUnitarityTest();
    } else if (config_path.find("scale_invariance") != std::string::npos) {
        result = runScaleInvarianceTest();
    } else if (config_path.find("symmetry_analysis") != std::string::npos) {
        result = runSymmetryAnalysisTest();
    } else if (config_path.find("renormalizability") != std::string::npos) {
        result = runRenormalizabilityTest();
    } else if (config_path.find("causality") != std::string::npos) {
        result = runCausalityTest();
    } else if (config_path.find("three_generations") != std::string::npos) {
        result = runThreeGenerationsTest();
    } else if (config_path.find("electroweak") != std::string::npos) {
        result = runElectroweakTest();
    } else if (config_path.find("strong_force") != std::string::npos) {
        result = runStrongForceTest();
    } else if (config_path.find("higgs_connection") != std::string::npos) {
        result = runHiggsConnectionTest();
    } else if (config_path.find("particle_scattering") != std::string::npos) {
        result = runParticleScatteringTest();
    } else if (config_path.find("josephson_junction") != std::string::npos) {
        result = runJosephsonJunctionTest();
    } else if (config_path.find("binary_merger") != std::string::npos) {
        result = runBinaryMergerTest();
    } else if (config_path.find("spin_magnetism") != std::string::npos) {
        result = runSpinMagnetismTest();
    } else if (config_path.find("knot_topology") != std::string::npos) {
        result = runKnotTopologyTest();
    } else if (config_path.find("multiscale") != std::string::npos) {
        result = runMultiScaleTest();
    } else if (config_path.find("quantum_fluctuations") != std::string::npos) {
        result = runQuantumFluctuationsTest();
    } else if (config_path.find("finite_temperature") != std::string::npos) {
        result = runFiniteTemperatureTest();
    } else if (config_path.find("hpc_scaling") != std::string::npos) {
        result = runHPCScalingTest();
    } else if (config_path.find("fine_structure_constant") != std::string::npos) {
        result = runFineStructureConstantTest();
    } else if (config_path.find("atomic_physics") != std::string::npos) {
        result = runAtomicPhysicsTest();
    } else if (config_path.find("lhc_predictions") != std::string::npos) {
        result = runLHCPredictionsTest();
    } else if (config_path.find("astrophysical_observations") != std::string::npos) {
        result = runAstrophysicalObservationsTest();
    } else if (config_path.find("laboratory_scale") != std::string::npos) {
        result = runLaboratoryScaleTest();
    } else if (config_path.find("chiral_operator_comparison") != std::string::npos) {
        result = runChiralOperatorComparisonTest();
    } else if (config_path.find("kuramoto_phase_diagram") != std::string::npos) {
        result = runKuramotoPhaseDiagramTest();
    } else if (config_path.find("chiral_channel_selector") != std::string::npos) {
        result = runChiralChannelSelectorTest();
    } else if (config_path.find("njl_gap_curve") != std::string::npos) {
        result = runNJLGapCurveTest();
    } else {
        // Default: TRD field theory test (timesync, etc.)
        TRDTestRunner runner(config_path);

        std::cout << "\n[1/3] Initializing..." << std::endl;
        if (!runner.initialize()) {
            std::cerr << "\n✗ Initialization failed" << std::endl;
            return 1;
        }

        std::cout << "\n[2/3] Running tests..." << std::endl;
        bool success = runner.run();

        std::cout << "\n[3/3] Generating report..." << std::endl;
        runner.generateReport();

        std::cout << "\n===== Test Summary =====" << std::endl;
        if (success && runner.allTestsPassed()) {
            std::cout << "✓ ALL TESTS PASSED" << std::endl;
            result = 0;
        } else {
            std::cout << "✗ TESTS FAILED" << std::endl;
            result = 1;
        }
    }

    // Generate visualization if --plot flag was passed
    if (generate_plots) {
        std::cout << "\n[Visualization] Generating plots for " << test_name << "..." << std::endl;
        std::string output_dir = "output/" + test_name;
        VisualizationGenerator::generateTestPlot(test_name, output_dir);
    }

    return result;
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
        bool generate_plots = false;
        for (int i = 3; i < argc; ++i) {
            if (std::string(argv[i]) == "--plot" || std::string(argv[i]) == "-p") {
                generate_plots = true;
            }
        }
        return runTestMode(config_path, generate_plots);
    }

    // Unknown argument
    std::cerr << "Error: Unknown argument '" << arg1 << "'" << std::endl;
    printUsage(argv[0]);
    return 1;
}
