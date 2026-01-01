#ifndef VISUALIZATION_GENERATOR_H
#define VISUALIZATION_GENERATOR_H

#include <string>
#include <vector>
#include "simulations/ObservableComputer.h"

/**
 * VisualizationGenerator - Generate Python visualization scripts on-the-fly
 * 
 * Creates publication-quality plots from test results:
 * - Long-time stability (norm & energy evolution)
 * - Operator splitting convergence
 * - Wavepacket dynamics
 */
class VisualizationGenerator {
public:
    /**
     * Generate visualization script for test results
     * @param output_dir Directory containing test results
     * @param test_name Name of the test
     * @param N_ratios Vector of N ratios tested
     * @return true if successful
     */
    static bool generateScript(const std::string& output_dir,
                               const std::string& test_name,
                               const std::vector<int>& N_ratios);

    /**
     * Execute visualization script to generate plots
     * @param script_path Path to Python script
     * @return true if successful
     */
    static bool executeScript(const std::string& script_path);

private:
    static std::string getPythonScript(const std::string& test_name,
                                       const std::vector<int>& N_ratios);
};

#endif // VISUALIZATION_GENERATOR_H
