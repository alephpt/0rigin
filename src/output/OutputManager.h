/**
 * @file OutputManager.h
 * @brief Systematic output generation for MSFT experiments
 *
 * Provides centralized management for experiment output:
 * - Timestamped directory creation
 * - Consistent file formats (.dat with headers)
 * - Metadata tracking (simulation parameters)
 * - Compatible with existing visualization scripts
 *
 * USAGE EXAMPLE:
 * @code
 * #include "output/OutputManager.h"
 *
 * // Initialize manager
 * MSFT::OutputManager output_manager;
 *
 * // Create experiment directory with timestamp
 * std::string exp_dir = output_manager.createExperimentDirectory("my_experiment");
 * // Creates: /home/persist/neotec/0rigin/output/YYYYMMDD_HHMMSS_my_experiment/
 *
 * // Prepare metadata
 * std::map<std::string, std::string> metadata;
 * metadata["grid_size"] = "256x256";
 * metadata["dt"] = "0.01";
 * metadata["steps"] = "1000";
 *
 * // Write timeseries data
 * std::vector<float> time_points = {0.0, 0.01, 0.02, ...};
 * std::vector<float> R_values = {0.99, 0.98, 0.97, ...};
 * output_manager.writeTimeseries(exp_dir + "/R_evolution.dat",
 *                                time_points, R_values, metadata, "R_global");
 *
 * // Write 2D field snapshot
 * std::vector<float> field_data(256*256);  // Your field data
 * output_manager.writeSnapshot(exp_dir + "/field_snapshot.dat",
 *                              field_data, 256, 256, metadata);
 *
 * // Create subdirectories as needed
 * std::string snapshots_dir = output_manager.createSubdirectory(exp_dir, "snapshots");
 * @endcode
 *
 * OUTPUT FORMAT:
 * - All .dat files include metadata headers (# prefixed lines)
 * - Compatible with existing Python visualization scripts
 * - CSV files for structured data
 * - Metadata files in key=value format
 */

#ifndef MSFT_OUTPUT_MANAGER_H
#define MSFT_OUTPUT_MANAGER_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <filesystem>

namespace MSFT {

/**
 * @brief Manages systematic output generation for experiments
 */
class OutputManager {
public:
    /**
     * @brief Constructor
     * @param base_path Base directory for all outputs (default: /home/persist/neotec/0rigin/output)
     */
    explicit OutputManager(const std::string& base_path = "/home/persist/neotec/0rigin/output");

    /**
     * @brief Create timestamped experiment directory
     * @param experiment_name Name of the experiment (e.g., "noise_sweep", "dirac_evolution")
     * @return Full path to the created directory
     */
    std::string createExperimentDirectory(const std::string& experiment_name);

    /**
     * @brief Write timeseries data with metadata header
     * @param filepath Full path to output file
     * @param time_data Time points
     * @param observable_data Observable values at each time
     * @param metadata Key-value pairs for simulation parameters
     * @param observable_name Name of the observable (for header)
     */
    void writeTimeseries(const std::string& filepath,
                        const std::vector<float>& time_data,
                        const std::vector<float>& observable_data,
                        const std::map<std::string, std::string>& metadata,
                        const std::string& observable_name = "observable");

    /**
     * @brief Write 2D field snapshot with metadata
     * @param filepath Full path to output file
     * @param field_data Flattened 2D field data (row-major)
     * @param nx Grid width
     * @param ny Grid height
     * @param metadata Key-value pairs for simulation parameters
     */
    void writeSnapshot(const std::string& filepath,
                      const std::vector<float>& field_data,
                      int nx, int ny,
                      const std::map<std::string, std::string>& metadata);

    /**
     * @brief Write multiple observables timeseries (for convenience)
     * @param filepath Full path to output file
     * @param time_data Time points
     * @param observables Map of observable names to data vectors
     * @param metadata Key-value pairs for simulation parameters
     */
    void writeMultipleTimeseries(const std::string& filepath,
                                 const std::vector<float>& time_data,
                                 const std::map<std::string, std::vector<float>>& observables,
                                 const std::map<std::string, std::string>& metadata);

    /**
     * @brief Write simulation metadata to JSON-like format
     * @param filepath Full path to output file
     * @param parameters Map of parameter names to values
     */
    void writeMetadata(const std::string& filepath,
                      const std::map<std::string, std::string>& parameters);

    /**
     * @brief Create subdirectory within experiment directory
     * @param parent_dir Parent directory path
     * @param subdir_name Name of subdirectory
     * @return Full path to created subdirectory
     */
    std::string createSubdirectory(const std::string& parent_dir,
                                   const std::string& subdir_name);

    /**
     * @brief Get timestamp string (for filenames)
     * @return Timestamp in format YYYYMMDD_HHMMSS
     */
    static std::string getTimestamp();

    /**
     * @brief Write CSV file with headers
     * @param filepath Full path to output file
     * @param headers Column headers
     * @param data 2D data (rows x columns)
     */
    void writeCSV(const std::string& filepath,
                  const std::vector<std::string>& headers,
                  const std::vector<std::vector<float>>& data);

    /**
     * @brief Append a single row to existing CSV
     * @param filepath Full path to output file
     * @param row Data row to append
     * @param create_if_missing If true, create file with headers if it doesn't exist
     * @param headers Headers to use if creating new file
     */
    void appendCSV(const std::string& filepath,
                   const std::vector<float>& row,
                   bool create_if_missing = true,
                   const std::vector<std::string>& headers = {});

private:
    std::string base_path_;
    std::string current_experiment_dir_;

    /**
     * @brief Write metadata header to stream
     */
    void writeMetadataHeader(std::ofstream& stream,
                            const std::map<std::string, std::string>& metadata);

    /**
     * @brief Ensure directory exists, create if necessary
     */
    void ensureDirectoryExists(const std::string& path);

    /**
     * @brief Validate file can be written
     */
    bool validateFileWrite(const std::string& filepath);
};

} // namespace MSFT

#endif // MSFT_OUTPUT_MANAGER_H