#ifndef VISUALIZATION_GENERATOR_H
#define VISUALIZATION_GENERATOR_H

#include <string>
#include <vector>
#include <map>
#include <mutex>

/**
 * VisualizationGenerator - Generate Python visualization scripts for TRD tests
 *
 * Creates publication-quality matplotlib plots from test results:
 * - Reads existing CSV data from output/<test>/ directories
 * - Generates self-contained Python scripts with embedded data for stdout-only tests
 * - Supports all test categories (GR, particle physics, cosmology, EM, mathematical)
 */
class VisualizationGenerator {
public:
    /**
     * Generate visualization script for operator splitting test results
     * @param output_dir Directory containing test results
     * @param test_name Name of the test
     * @param N_ratios Vector of N ratios tested
     * @return true if successful
     */
    static bool generateScript(const std::string& output_dir,
                               const std::string& test_name,
                               const std::vector<int>& N_ratios);

    /**
     * Execute a Python visualization script
     * @param script_path Path to Python script
     * @return true if successful
     */
    static bool executeScript(const std::string& script_path);

    /**
     * Generate a test-specific plot based on test_name
     * Main entry point for visualization generation.
     * @param test_name Name of the test (matches config filename stem)
     * @param output_dir Output directory for script and PNG
     * @return true if successful
     */
    static bool generateTestPlot(const std::string& test_name,
                                 const std::string& output_dir);

    // --- Static data collection API ---
    // Tests call these during execution; generateTestPlot() uses collected data.

    /**
     * Add a single data point to a named series
     * @param series_name Identifier for the data series
     * @param x X-axis value
     * @param y Y-axis value
     */
    static void addDataPoint(const std::string& series_name, float x, float y);

    /**
     * Add a complete data series (bulk)
     * @param series_name Identifier for the data series
     * @param x X-axis values
     * @param y Y-axis values
     */
    static void addDataSeries(const std::string& series_name,
                              const std::vector<float>& x,
                              const std::vector<float>& y);

    /**
     * Clear all collected data series
     */
    static void clearData();

private:
    // Collected data: series_name -> {x_values, y_values}
    struct DataSeries {
        std::vector<float> x;
        std::vector<float> y;
    };
    static std::map<std::string, DataSeries> s_data;
    static std::mutex s_mutex;

    // Python script preamble with publication-quality style
    static std::string getPreamble();

    // Category-specific script generators (return full Python script)
    static std::string genEnergyConservation(const std::string& test_name,
                                             const std::string& output_dir);
    static std::string genWeakField(const std::string& output_dir);
    static std::string genGravitationalWaves(const std::string& output_dir);
    static std::string genFineStructure(const std::string& output_dir);
    static std::string genElectroweak(const std::string& output_dir);
    static std::string genStrongForce(const std::string& output_dir);
    static std::string genFriedmann(const std::string& output_dir);
    static std::string genDarkMatter(const std::string& output_dir);
    static std::string genDarkEnergy(const std::string& output_dir);
    static std::string genInflation(const std::string& output_dir);
    static std::string genLorentzForce(const std::string& output_dir);
    static std::string genJosephson(const std::string& output_dir);
    static std::string genCausality(const std::string& output_dir);
    static std::string genUnitarity(const std::string& output_dir);
    static std::string genEinsteinField(const std::string& output_dir);
    static std::string genThreeGenerations(const std::string& output_dir);
    static std::string genHiggsConnection(const std::string& output_dir);
    static std::string genParticleSpectrum(const std::string& output_dir);
    static std::string genBinaryMerger(const std::string& output_dir);
    static std::string genSpinMagnetism(const std::string& output_dir);
    static std::string genKnotTopology(const std::string& output_dir);
    static std::string genEmbeddedPlot(const std::string& test_name,
                                       const std::string& output_dir);

    // Helper to write script and optionally execute
    static bool writeAndRun(const std::string& script,
                            const std::string& output_dir,
                            const std::string& filename);

    // Format collected data as Python arrays
    static std::string dataToNumpy(const std::string& series_name);

    // Legacy helper
    static std::string getPythonScript(const std::string& test_name,
                                       const std::vector<int>& N_ratios);
};

#endif // VISUALIZATION_GENERATOR_H
