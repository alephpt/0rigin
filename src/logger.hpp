#pragma once

#include <string>
#include <fstream>
#include <memory>
#include <chrono>
#include <map>
#include "observables.hpp"

namespace MSFT {

/**
 * @file logger.hpp
 * @brief Production logging system for MSFT simulations
 *
 * FEATURES:
 * - Multiple output formats: CSV, JSON, Binary
 * - Configurable logging frequency
 * - Real-time progress updates
 * - Automatic file management
 * - Performance metrics
 */

/**
 * @brief Logging configuration
 */
struct LogConfig {
    std::string output_dir = "./output";      // Output directory
    std::string run_name = "MSFT_run";        // Run identifier
    bool enable_csv = true;                   // Enable CSV output
    bool enable_json = true;                  // Enable JSON output
    bool enable_binary = false;               // Enable binary output
    bool enable_progress = true;              // Enable progress updates
    uint32_t log_frequency = 1;               // Log every N steps
    uint32_t progress_frequency = 100;        // Progress update every N steps
    bool save_fields = false;                 // Save full field data
    uint32_t field_save_frequency = 1000;     // Save fields every N steps
};

/**
 * @brief Performance metrics
 */
struct PerformanceMetrics {
    double total_time_ms = 0.0;              // Total simulation time
    double avg_step_time_ms = 0.0;           // Average time per step
    double steps_per_second = 0.0;           // Throughput
    uint64_t total_steps = 0;                // Total steps executed
    double gpu_compute_time_ms = 0.0;        // GPU compute time
    double io_time_ms = 0.0;                 // I/O time
};

/**
 * @brief MSFT simulation logger
 */
class MSFTLogger {
public:
    /**
     * @brief Constructor
     */
    MSFTLogger(const LogConfig& config = LogConfig());

    /**
     * @brief Destructor - flushes all buffers
     */
    ~MSFTLogger();

    /**
     * @brief Initialize logger (create directories, open files)
     */
    void initialize();

    /**
     * @brief Log observables at current timestep
     *
     * @param obs Observables to log
     * @param force Force logging even if not at log frequency
     */
    void log(const AllObservables& obs, bool force = false);

    /**
     * @brief Log full field data
     *
     * @param obs Observables with field data
     * @param step Current timestep
     */
    void logFields(const AllObservables& obs, uint32_t step);

    /**
     * @brief Update progress display
     *
     * @param current_step Current timestep
     * @param total_steps Total steps to execute
     * @param time Current simulation time
     */
    void updateProgress(uint32_t current_step, uint32_t total_steps, float time);

    /**
     * @brief Log performance metrics
     *
     * @param metrics Performance data
     */
    void logPerformance(const PerformanceMetrics& metrics);

    /**
     * @brief Log simulation parameters
     *
     * @param params Parameter map
     */
    void logParameters(const std::map<std::string, float>& params);

    /**
     * @brief Finalize logging (write summaries, close files)
     */
    void finalize();

    /**
     * @brief Get output directory
     */
    const std::string& getOutputDir() const { return config_.output_dir; }

    /**
     * @brief Get run name
     */
    const std::string& getRunName() const { return config_.run_name; }

    /**
     * @brief Flush all output streams
     */
    void flush();

private:
    LogConfig config_;
    ObservableTimeSeries time_series_;

    // File streams
    std::ofstream csv_file_;
    std::ofstream json_file_;
    std::ofstream performance_file_;
    std::ofstream params_file_;

    // State tracking
    uint32_t log_counter_ = 0;
    uint32_t progress_counter_ = 0;
    bool initialized_ = false;

    // Timing
    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point last_progress_time_;

    // Helper methods
    void createOutputDirectory();
    void openFiles();
    void closeFiles();
    void writeCSVHeader();
    void writeJSONHeader();
    std::string getTimestamp() const;
    std::string formatDuration(double seconds) const;
};

/**
 * @brief Simple progress bar for console output
 */
class ProgressBar {
public:
    ProgressBar(uint32_t total, uint32_t width = 50);

    void update(uint32_t current, const std::string& message = "");
    void finish();

private:
    uint32_t total_;
    uint32_t width_;
    uint32_t last_printed_;
};

/**
 * @brief Summary statistics generator
 */
class SummaryStats {
public:
    /**
     * @brief Generate summary from time series
     */
    static std::string generate(const ObservableTimeSeries& series);

    /**
     * @brief Generate markdown report
     */
    static std::string generateMarkdown(
        const ObservableTimeSeries& series,
        const PerformanceMetrics& metrics,
        const std::map<std::string, float>& params
    );
};

} // namespace MSFT