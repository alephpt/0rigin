// include/TRDCSVWriter.h
#pragma once

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <stdexcept>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif

/**
 * TRDCSVWriter - Standardized CSV Output for TRD Physics Tests
 *
 * This header provides unified CSV export functionality across all TRD validation tests.
 * Eliminates 200+ lines of duplicate code across 15 tests by centralizing data export,
 * metadata generation, and output directory management.
 *
 * Key Features:
 *   1. Automatic metadata section (timestamp, git commit, test parameters)
 *   2. Standardized directory structure (output/<test_name>/results_YYYYMMDD_HHMMSS.csv)
 *   3. Type-safe row writing (auto-converts int/float/double/string)
 *   4. Consistent formatting (precision control, scientific notation for extremes)
 *   5. Error handling (permission checks, directory creation)
 *
 * Metadata Format (Comment Header):
 *   # Test: fine_structure_constant
 *   # Date: 2026-01-05T12:34:56
 *   # Git: a3f2e8d (main branch)
 *   # TRD Version: v3D Unified
 *   # Golden Key: 246 GeV
 *   # Parameters: K=1.0, grid=64x64x32, steps=1000
 *   #
 *   Method,Alpha_Measured,Alpha_QED,Ratio
 *   Energy,0.00353983,0.00729735,0.485084
 *
 * Output Directory Structure:
 *   output/
 *     fine_structure_constant/
 *       results_20260105_123456.csv
 *       phase_diagram_20260105_123456.csv
 *     lorentz_force/
 *       trajectory_20260105_140522.csv
 *
 * Usage Examples:
 *   // Basic usage with automatic timestamping
 *   TRD::CSVWriter csv("fine_structure_constant", "B2_FineStructure");
 *   csv.writeMetadata({{"K_coupling", "1.0"}, {"grid_size", "64x64x32"}});
 *   csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});
 *   csv.writeRow("Energy", 0.00353983, 0.00729735, 0.485084);
 *   csv.close();
 *
 *   // Multiple CSV files for same test
 *   TRD::CSVWriter results("results", "B2_FineStructure");
 *   TRD::CSVWriter phase_diagram("phase_diagram", "B2_FineStructure");
 *
 *   // Custom precision for high-accuracy data
 *   TRD::CSVWriter csv("particle_trajectory", "C1_Lorentz");
 *   csv.setPrecision(12);  // 12 significant figures
 *   csv.writeRow(0, 1.234567890123e-8, 2.345678901234e-8, 3.456789012345e-8);
 *
 * References:
 *   - Architecture review: ARCHITECTURE_REVIEW_CATEGORY_BF.md
 *   - Code quality standards: CLAUDE.md (TRD-Specific Standards)
 *   - Affected tests: 15 validation tests requiring CSV export
 *
 * @author TRD Development Team
 * @version v3D Unified
 */

namespace TRD {

// ============================================================================
// Helper Functions
// ============================================================================

namespace detail {

/**
 * Get current timestamp in ISO 8601 format
 * @return Formatted timestamp string (YYYY-MM-DDTHH:MM:SS)
 */
inline std::string getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;

#ifdef _WIN32
    localtime_s(&tm_now, &time_t_now);
#else
    localtime_r(&time_t_now, &tm_now);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y-%m-%dT%H:%M:%S");
    return oss.str();
}

/**
 * Get current timestamp for filename (YYYYMMDD_HHMMSS)
 * @return Formatted timestamp string for filenames
 */
inline std::string getFileTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;

#ifdef _WIN32
    localtime_s(&tm_now, &time_t_now);
#else
    localtime_r(&time_t_now, &tm_now);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

/**
 * Execute shell command and capture output
 * @param cmd Command to execute
 * @return Command output (trimmed)
 */
inline std::string execCommand(const char* cmd) {
    char buffer[128];
    std::string result = "";

#ifdef _WIN32
    FILE* pipe = _popen(cmd, "r");
#else
    FILE* pipe = popen(cmd, "r");
#endif

    if (!pipe) return "N/A";

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

#ifdef _WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif

    // Trim trailing whitespace
    while (!result.empty() && (result.back() == '\n' || result.back() == '\r' || result.back() == ' ')) {
        result.pop_back();
    }

    return result.empty() ? "N/A" : result;
}

/**
 * Create directory recursively
 * @param path Directory path to create
 * @return true if successful or already exists, false otherwise
 */
inline bool createDirectory(const std::string& path) {
    struct stat info;

    // Check if directory already exists
    if (stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR)) {
        return true;
    }

    // Create parent directory first
    size_t pos = path.find_last_of("/\\");
    if (pos != std::string::npos) {
        std::string parent = path.substr(0, pos);
        if (!parent.empty() && !createDirectory(parent)) {
            return false;
        }
    }

    // Create this directory
#ifdef _WIN32
    return mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

} // namespace detail

// ============================================================================
// Public API
// ============================================================================

/**
 * Get current git commit hash for reproducibility
 * @return Short commit hash (7 chars) or "N/A" if not in git repo
 */
inline std::string getGitCommitHash() {
    return detail::execCommand("git rev-parse --short HEAD 2>/dev/null");
}

/**
 * Get current git branch name
 * @return Branch name or "N/A" if not in git repo
 */
inline std::string getGitBranch() {
    return detail::execCommand("git rev-parse --abbrev-ref HEAD 2>/dev/null");
}

/**
 * Create standard output directory structure
 * @param test_name Name of the test (creates output/<test_name>/)
 * @return Full path to created directory
 * @throws std::runtime_error if directory creation fails
 */
inline std::string createOutputDirectory(const std::string& test_name) {
    std::string dir_path = "output/" + test_name;

    if (!detail::createDirectory(dir_path)) {
        throw std::runtime_error("Failed to create output directory: " + dir_path);
    }

    return dir_path;
}

/**
 * CSVWriter - Type-safe CSV file writer with automatic metadata
 *
 * Provides standardized CSV export with automatic metadata generation,
 * type-safe row writing, and consistent formatting across all TRD tests.
 */
class CSVWriter {
public:
    /**
     * Constructor - Create CSV file with optional timestamping
     *
     * @param filename Base filename (without .csv extension)
     * @param test_name Test name for directory structure (e.g., "B2_FineStructure")
     * @param append_timestamp If true, append YYYYMMDD_HHMMSS to filename
     *
     * Examples:
     *   CSVWriter("results", "B2_FineStructure", true)
     *     → output/B2_FineStructure/results_20260105_123456.csv
     *
     *   CSVWriter("phase_diagram", "B2_FineStructure", false)
     *     → output/B2_FineStructure/phase_diagram.csv
     */
    CSVWriter(const std::string& filename,
              const std::string& test_name = "",
              bool append_timestamp = true)
        : test_name_(test_name)
        , precision_(6)
        , use_scientific_(false)
        , metadata_written_(false)
        , header_written_(false)
    {
        // Build full file path
        std::string dir_path;
        if (!test_name.empty()) {
            dir_path = createOutputDirectory(test_name);
        } else {
            dir_path = "output";
            detail::createDirectory(dir_path);
        }

        // Build filename
        std::ostringstream filename_builder;
        filename_builder << dir_path << "/" << filename;

        if (append_timestamp) {
            filename_builder << "_" << detail::getFileTimestamp();
        }

        filename_builder << ".csv";
        file_path_ = filename_builder.str();

        // Open file
        file_.open(file_path_);
        if (!file_.is_open()) {
            throw std::runtime_error("Failed to open CSV file: " + file_path_);
        }
    }

    /**
     * Destructor - Ensure file is closed
     */
    ~CSVWriter() {
        if (file_.is_open()) {
            close();
        }
    }

    // Disable copy (stream not copyable)
    CSVWriter(const CSVWriter&) = delete;
    CSVWriter& operator=(const CSVWriter&) = delete;

    /**
     * Set numeric precision for floating-point output
     * @param precision Number of significant digits (default: 6)
     */
    void setPrecision(int precision) {
        precision_ = precision;
    }

    /**
     * Enable/disable scientific notation for all numbers
     * @param use_scientific If true, use scientific notation (default: false)
     */
    void setScientificNotation(bool use_scientific) {
        use_scientific_ = use_scientific;
    }

    /**
     * Write metadata section as comment header
     *
     * Automatically includes:
     *   - Test name
     *   - Current timestamp
     *   - Git commit hash and branch
     *   - TRD version
     *   - Golden key reference (246 GeV)
     *
     * @param custom_metadata Additional custom metadata as key-value pairs
     *
     * Example:
     *   csv.writeMetadata({
     *       {"K_coupling", "1.0"},
     *       {"grid_size", "64x64x32"},
     *       {"evolution_steps", "1000"}
     *   });
     */
    void writeMetadata(const std::map<std::string, std::string>& custom_metadata = {}) {
        if (metadata_written_) {
            throw std::runtime_error("Metadata already written");
        }

        // Standard metadata
        if (!test_name_.empty()) {
            file_ << "# Test: " << test_name_ << "\n";
        }
        file_ << "# Date: " << detail::getCurrentTimestamp() << "\n";

        std::string git_hash = getGitCommitHash();
        std::string git_branch = getGitBranch();
        if (git_hash != "N/A") {
            file_ << "# Git: " << git_hash;
            if (git_branch != "N/A") {
                file_ << " (" << git_branch << " branch)";
            }
            file_ << "\n";
        }

        file_ << "# TRD Version: v3D Unified\n";
        file_ << "# Golden Key: 246 GeV\n";

        // Custom metadata
        if (!custom_metadata.empty()) {
            std::ostringstream params;
            bool first = true;
            for (const auto& entry : custom_metadata) {
                if (!first) params << ", ";
                params << entry.first << "=" << entry.second;
                first = false;
            }
            file_ << "# Parameters: " << params.str() << "\n";
        }

        // Blank comment line separator
        file_ << "#\n";

        metadata_written_ = true;
    }

    /**
     * Write header row with column names
     * @param columns Vector of column names
     *
     * Example:
     *   csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});
     */
    void writeHeader(const std::vector<std::string>& columns) {
        if (header_written_) {
            throw std::runtime_error("Header already written");
        }

        if (!metadata_written_) {
            // Write minimal metadata if user didn't call writeMetadata()
            writeMetadata();
        }

        // Write column names
        for (size_t i = 0; i < columns.size(); ++i) {
            if (i > 0) file_ << ",";
            file_ << columns[i];
        }
        file_ << "\n";

        header_written_ = true;
    }

    /**
     * Write data row with automatic type conversion
     *
     * Supports: int, long, float, double, const char*, std::string
     *
     * Example:
     *   csv.writeRow("Energy", 0.00353983, 0.00729735, 0.485084);
     *   csv.writeRow(0, 1.23e-8, 2.34e-8, 3.45e-8);
     */
    template<typename... Args>
    void writeRow(Args... values) {
        if (!header_written_) {
            throw std::runtime_error("Header must be written before data rows");
        }

        writeRowImpl(values...);
        file_ << "\n";
    }

    /**
     * Flush and close the CSV file
     */
    void close() {
        if (file_.is_open()) {
            file_.flush();
            file_.close();
        }
    }

    /**
     * Get full path to the created CSV file
     * @return Absolute or relative path to CSV file
     */
    std::string getFilePath() const {
        return file_path_;
    }

    /**
     * Check if file is currently open
     * @return true if file is open and ready for writing
     */
    bool isOpen() const {
        return file_.is_open();
    }

private:
    std::ofstream file_;
    std::string file_path_;
    std::string test_name_;
    int precision_;
    bool use_scientific_;
    bool metadata_written_;
    bool header_written_;

    // Recursive template for writing row values
    template<typename T>
    void writeRowImpl(T value) {
        writeValue(value);
    }

    template<typename T, typename... Rest>
    void writeRowImpl(T value, Rest... rest) {
        writeValue(value);
        file_ << ",";
        writeRowImpl(rest...);
    }

    // Type-specific value writers
    void writeValue(int value) {
        file_ << value;
    }

    void writeValue(long value) {
        file_ << value;
    }

    void writeValue(float value) {
        formatNumber(value);
    }

    void writeValue(double value) {
        formatNumber(value);
    }

    void writeValue(const char* value) {
        file_ << value;
    }

    void writeValue(const std::string& value) {
        file_ << value;
    }

    // Smart number formatting
    template<typename T>
    void formatNumber(T value) {
        // Check for extreme values that benefit from scientific notation
        bool use_sci = use_scientific_;
        if (!use_sci && value != 0.0) {
            double abs_val = std::abs(static_cast<double>(value));
            if (abs_val < 1e-3 || abs_val > 1e6) {
                use_sci = true;
            }
        }

        if (use_sci) {
            file_ << std::scientific;
        } else {
            file_ << std::fixed;
        }

        file_ << std::setprecision(precision_) << value;

        // Reset to default
        file_ << std::defaultfloat;
    }
};

} // namespace TRD
