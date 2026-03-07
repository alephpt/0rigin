#include "simulations/VisualizationGenerator.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

// Static member definitions
std::map<std::string, VisualizationGenerator::DataSeries> VisualizationGenerator::s_data;
std::mutex VisualizationGenerator::s_mutex;

// ============================================================================
// Public API
// ============================================================================

void VisualizationGenerator::addDataPoint(const std::string& series_name,
                                          float x, float y) {
    std::lock_guard<std::mutex> lock(s_mutex);
    s_data[series_name].x.push_back(x);
    s_data[series_name].y.push_back(y);
}

void VisualizationGenerator::addDataSeries(const std::string& series_name,
                                           const std::vector<float>& x,
                                           const std::vector<float>& y) {
    std::lock_guard<std::mutex> lock(s_mutex);
    auto& ds = s_data[series_name];
    ds.x.insert(ds.x.end(), x.begin(), x.end());
    ds.y.insert(ds.y.end(), y.begin(), y.end());
}

void VisualizationGenerator::clearData() {
    std::lock_guard<std::mutex> lock(s_mutex);
    s_data.clear();
}

bool VisualizationGenerator::generateTestPlot(const std::string& test_name,
                                              const std::string& output_dir) {
    std::string script;

    if (test_name.find("weak_field") != std::string::npos) {
        script = genWeakField(output_dir);
    } else if (test_name.find("gravitational_waves") != std::string::npos) {
        script = genGravitationalWaves(output_dir);
    } else if (test_name.find("fine_structure") != std::string::npos) {
        script = genFineStructure(output_dir);
    } else if (test_name.find("electroweak") != std::string::npos) {
        script = genElectroweak(output_dir);
    } else if (test_name.find("strong_force") != std::string::npos) {
        script = genStrongForce(output_dir);
    } else if (test_name.find("friedmann") != std::string::npos) {
        script = genFriedmann(output_dir);
    } else if (test_name.find("dark_matter") != std::string::npos) {
        script = genDarkMatter(output_dir);
    } else if (test_name.find("dark_energy") != std::string::npos) {
        script = genDarkEnergy(output_dir);
    } else if (test_name.find("inflation") != std::string::npos) {
        script = genInflation(output_dir);
    } else if (test_name.find("lorentz_force") != std::string::npos) {
        script = genLorentzForce(output_dir);
    } else if (test_name.find("josephson") != std::string::npos) {
        script = genJosephson(output_dir);
    } else if (test_name.find("causality") != std::string::npos) {
        script = genCausality(output_dir);
    } else if (test_name.find("unitarity") != std::string::npos) {
        script = genUnitarity(output_dir);
    } else if (test_name.find("einstein_field") != std::string::npos) {
        script = genEinsteinField(output_dir);
    } else if (test_name.find("three_generation") != std::string::npos) {
        script = genThreeGenerations(output_dir);
    } else if (test_name.find("higgs_connection") != std::string::npos) {
        script = genHiggsConnection(output_dir);
    } else if (test_name.find("particle_spectrum") != std::string::npos) {
        script = genParticleSpectrum(output_dir);
    } else if (test_name.find("binary_merger") != std::string::npos) {
        script = genBinaryMerger(output_dir);
    } else if (test_name.find("spin_magnetism") != std::string::npos) {
        script = genSpinMagnetism(output_dir);
    } else if (test_name.find("knot_topology") != std::string::npos ||
               test_name.find("knot_stability") != std::string::npos) {
        script = genKnotTopology(output_dir);
    } else {
        // Generic: use embedded data if available, else energy conservation
        if (!s_data.empty()) {
            script = genEmbeddedPlot(test_name, output_dir);
        } else {
            script = genEnergyConservation(test_name, output_dir);
        }
    }

    if (script.empty()) {
        std::cerr << "Warning: No visualization generated for " << test_name
                  << std::endl;
        return false;
    }

    return writeAndRun(script, output_dir, test_name + "_plot.py");
}

bool VisualizationGenerator::generateScript(const std::string& output_dir,
                                            const std::string& test_name,
                                            const std::vector<int>& N_ratios) {
    std::string script = getPythonScript(test_name, N_ratios);
    return writeAndRun(script, output_dir, "visualize.py");
}

bool VisualizationGenerator::executeScript(const std::string& script_path) {
    std::string cmd = "python3 " + script_path + " 2>&1";
    int result = system(cmd.c_str());
    if (result != 0) {
        std::cerr << "Warning: Script execution failed: " << script_path
                  << std::endl;
        return false;
    }
    return true;
}

// ============================================================================
// Helpers
// ============================================================================

std::string VisualizationGenerator::getPreamble() {
    return R"(#!/usr/bin/env python3
"""TRD Visualization - Publication Quality Plots"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, glob

plt.rcParams.update({
    'figure.figsize': (10, 7),
    'figure.dpi': 150,
    'font.family': 'serif',
    'font.size': 12,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.color': '#cccccc',
    'lines.linewidth': 1.8,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

COLORS = plt.cm.tab10.colors

def find_csv(directory, pattern):
    """Find CSV matching pattern. Search given dir, then current dir, then sibling output dirs."""
    for d in [directory, '.', '..']:
        matches = sorted(glob.glob(os.path.join(d, pattern)))
        if matches:
            return matches[-1]
    # Search all sibling directories under parent
    parent = os.path.dirname(directory) if directory != '.' else '..'
    if os.path.isdir(parent):
        for sub in sorted(os.listdir(parent)):
            sub_path = os.path.join(parent, sub)
            if os.path.isdir(sub_path):
                matches = sorted(glob.glob(os.path.join(sub_path, pattern)))
                if matches:
                    return matches[-1]
    return None

def load_csv(path):
    """Load CSV skipping comment lines (# prefix)."""
    rows = []
    headers = None
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if headers is None:
                headers = line.split(',')
            else:
                rows.append(line.split(','))
    data = {}
    for i, h in enumerate(headers):
        vals = []
        for row in rows:
            try:
                vals.append(float(row[i]))
            except (ValueError, IndexError):
                vals.append(row[i] if i < len(row) else '')
        data[h.strip()] = vals
    return data

)";
}

bool VisualizationGenerator::writeAndRun(const std::string& script,
                                         const std::string& output_dir,
                                         const std::string& filename) {
    // Ensure directory exists
    struct stat info;
    if (stat(output_dir.c_str(), &info) != 0) {
        std::string cmd = "mkdir -p " + output_dir;
        system(cmd.c_str());
    }

    std::string path = output_dir + "/" + filename;
    std::ofstream f(path);
    if (!f.is_open()) {
        std::cerr << "Failed to write script: " << path << std::endl;
        return false;
    }
    f << script;
    f.close();

    // Execute from output_dir — scripts use '.' as working directory
    std::string cmd = "cd " + output_dir + " && python3 " + filename + " 2>&1";
    int result = system(cmd.c_str());
    if (result == 0) {
        std::cout << "  Plots saved in " << output_dir << "/" << std::endl;
    } else {
        std::cerr << "  Warning: Plot generation failed for " << filename
                  << std::endl;
    }
    return result == 0;
}

std::string VisualizationGenerator::dataToNumpy(const std::string& series_name) {
    std::lock_guard<std::mutex> lock(s_mutex);
    auto it = s_data.find(series_name);
    if (it == s_data.end() || it->second.x.empty()) {
        return "";
    }

    std::ostringstream oss;
    oss << std::scientific;

    oss << series_name << "_x = np.array([";
    for (size_t i = 0; i < it->second.x.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << it->second.x[i];
    }
    oss << "])\n";

    oss << series_name << "_y = np.array([";
    for (size_t i = 0; i < it->second.y.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << it->second.y[i];
    }
    oss << "])\n";

    return oss.str();
}
