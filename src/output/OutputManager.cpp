/**
 * @file OutputManager.cpp
 * @brief Implementation of systematic output generation for SMFT experiments
 */

#include "OutputManager.h"
#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>

namespace SMFT {

OutputManager::OutputManager(const std::string& base_path)
    : base_path_(base_path) {
    ensureDirectoryExists(base_path_);
}

std::string OutputManager::createExperimentDirectory(const std::string& experiment_name) {
    // Create directory with format: YYYYMMDD_HHMMSS_experiment_name
    std::string timestamp = getTimestamp();
    std::string dir_name = timestamp + "_" + experiment_name;
    std::string full_path = base_path_ + "/" + dir_name;

    ensureDirectoryExists(full_path);
    current_experiment_dir_ = full_path;

    std::cout << "Created experiment directory: " << full_path << std::endl;
    return full_path;
}

void OutputManager::writeTimeseries(const std::string& filepath,
                                   const std::vector<float>& time_data,
                                   const std::vector<float>& observable_data,
                                   const std::map<std::string, std::string>& metadata,
                                   const std::string& observable_name) {
    if (time_data.size() != observable_data.size()) {
        std::cerr << "Error: time_data and observable_data sizes mismatch" << std::endl;
        return;
    }

    std::ofstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write metadata header
    writeMetadataHeader(file, metadata);

    // Write column headers
    file << "# time " << observable_name << "\n";

    // Write data
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < time_data.size(); ++i) {
        file << time_data[i] << " " << observable_data[i] << "\n";
    }

    file.close();
    std::cout << "Wrote timeseries to: " << filepath << std::endl;
}

void OutputManager::writeSnapshot(const std::string& filepath,
                                 const std::vector<float>& field_data,
                                 int nx, int ny,
                                 const std::map<std::string, std::string>& metadata) {
    if (field_data.size() != static_cast<size_t>(nx * ny)) {
        std::cerr << "Error: field_data size doesn't match nx*ny" << std::endl;
        return;
    }

    std::ofstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write metadata header
    writeMetadataHeader(file, metadata);

    // Add grid dimensions to header
    file << "# Grid: " << nx << " x " << ny << "\n";
    file << "# Format: Row-major order, space-separated values\n";

    // Write field data
    file << std::fixed << std::setprecision(6);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            file << field_data[y * nx + x];
            if (x < nx - 1) file << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Wrote snapshot to: " << filepath << std::endl;
}

void OutputManager::writeMultipleTimeseries(const std::string& filepath,
                                           const std::vector<float>& time_data,
                                           const std::map<std::string, std::vector<float>>& observables,
                                           const std::map<std::string, std::string>& metadata) {
    // Validate all observables have same length as time_data
    for (const auto& [name, data] : observables) {
        if (data.size() != time_data.size()) {
            std::cerr << "Error: Observable '" << name << "' size mismatch with time_data" << std::endl;
            return;
        }
    }

    std::ofstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write metadata header
    writeMetadataHeader(file, metadata);

    // Write column headers
    file << "# time";
    for (const auto& [name, data] : observables) {
        file << " " << name;
    }
    file << "\n";

    // Write data
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < time_data.size(); ++i) {
        file << time_data[i];
        for (const auto& [name, data] : observables) {
            file << " " << data[i];
        }
        file << "\n";
    }

    file.close();
    std::cout << "Wrote multiple timeseries to: " << filepath << std::endl;
}

void OutputManager::writeMetadata(const std::string& filepath,
                                 const std::map<std::string, std::string>& parameters) {
    std::ofstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write in simple key=value format (easy to parse)
    file << "# SMFT Experiment Metadata\n";
    file << "# Generated: " << getTimestamp() << "\n";
    file << "#\n";

    for (const auto& [key, value] : parameters) {
        file << key << "=" << value << "\n";
    }

    file.close();
    std::cout << "Wrote metadata to: " << filepath << std::endl;
}

std::string OutputManager::createSubdirectory(const std::string& parent_dir,
                                             const std::string& subdir_name) {
    std::string full_path = parent_dir + "/" + subdir_name;
    ensureDirectoryExists(full_path);
    return full_path;
}

std::string OutputManager::getTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S");
    return ss.str();
}

void OutputManager::writeCSV(const std::string& filepath,
                            const std::vector<std::string>& headers,
                            const std::vector<std::vector<float>>& data) {
    std::ofstream file(filepath);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write headers
    for (size_t i = 0; i < headers.size(); ++i) {
        file << headers[i];
        if (i < headers.size() - 1) file << ",";
    }
    file << "\n";

    // Write data
    file << std::fixed << std::setprecision(6);
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Wrote CSV to: " << filepath << std::endl;
}

void OutputManager::appendCSV(const std::string& filepath,
                             const std::vector<float>& row,
                             bool create_if_missing,
                             const std::vector<std::string>& headers) {
    bool file_exists = std::filesystem::exists(filepath);

    if (!file_exists && !create_if_missing) {
        std::cerr << "File doesn't exist and create_if_missing=false: " << filepath << std::endl;
        return;
    }

    // Open in append mode
    std::ios_base::openmode mode = std::ios::app;
    if (!file_exists && create_if_missing) {
        mode = std::ios::out;  // Create new file
    }

    std::ofstream file(filepath, mode);
    if (!file) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return;
    }

    // Write headers if creating new file
    if (!file_exists && create_if_missing && !headers.empty()) {
        for (size_t i = 0; i < headers.size(); ++i) {
            file << headers[i];
            if (i < headers.size() - 1) file << ",";
        }
        file << "\n";
    }

    // Write data row
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < row.size(); ++i) {
        file << row[i];
        if (i < row.size() - 1) file << ",";
    }
    file << "\n";

    file.close();
}

// Private methods

void OutputManager::writeMetadataHeader(std::ofstream& stream,
                                       const std::map<std::string, std::string>& metadata) {
    stream << "# SMFT Experiment Output\n";
    stream << "# Generated: " << getTimestamp() << "\n";

    if (!metadata.empty()) {
        stream << "# Parameters:\n";
        for (const auto& [key, value] : metadata) {
            stream << "#   " << key << ": " << value << "\n";
        }
    }
    stream << "#\n";
}

void OutputManager::ensureDirectoryExists(const std::string& path) {
    try {
        std::filesystem::create_directories(path);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create directory " << path << ": " << e.what() << std::endl;
    }
}

bool OutputManager::validateFileWrite(const std::string& filepath) {
    std::ofstream test(filepath);
    bool can_write = test.good();
    test.close();
    if (can_write) {
        std::filesystem::remove(filepath);  // Remove test file
    }
    return can_write;
}

} // namespace SMFT