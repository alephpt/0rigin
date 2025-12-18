/**
 * Mass Quantization Analysis
 *
 * Purpose: Check if defect masses are quantized
 *
 * Theory (from topological field theory):
 * - If defects are topological solitons, mass may be quantized
 * - m_defect = n · m_quantum (n = winding number)
 * - Histogram of m should show peaks at discrete values
 *
 * Test:
 * 1. Collect large sample of defects
 * 2. Measure m = Δ·R at each defect core
 * 3. Create histogram
 * 4. Look for peaks (quantization)
 * 5. Fit peaks to check if spacing is uniform
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <random>
#include <map>

struct Defect {
    int x, y;
    float winding;
    float R_local;
    float m_local;
};

float phase_diff(float theta1, float theta2) {
    float diff = theta1 - theta2;
    while (diff > M_PI) diff -= 2.0f * M_PI;
    while (diff < -M_PI) diff += 2.0f * M_PI;
    return diff;
}

std::vector<Defect> detect_defects(const std::vector<float>& theta,
                                     const std::vector<float>& R_field,
                                     uint32_t Nx, uint32_t Ny,
                                     float Delta) {
    std::vector<Defect> defects;
    
    for (uint32_t y = 0; y < Ny - 1; y++) {
        for (uint32_t x = 0; x < Nx - 1; x++) {
            uint32_t idx00 = y * Nx + x;
            uint32_t idx10 = y * Nx + (x + 1);
            uint32_t idx01 = (y + 1) * Nx + x;
            uint32_t idx11 = (y + 1) * Nx + (x + 1);
            
            float circ = 0.0f;
            circ += phase_diff(theta[idx10], theta[idx00]);
            circ += phase_diff(theta[idx11], theta[idx10]);
            circ += phase_diff(theta[idx01], theta[idx11]);
            circ += phase_diff(theta[idx00], theta[idx01]);
            
            float winding = circ / (2.0f * M_PI);
            
            if (std::abs(winding) > 0.5f) {
                Defect d;
                d.x = x;
                d.y = y;
                d.winding = winding;
                d.R_local = (R_field[idx00] + R_field[idx10] + 
                             R_field[idx01] + R_field[idx11]) / 4.0f;
                d.m_local = Delta * d.R_local;
                defects.push_back(d);
            }
        }
    }
    
    return defects;
}

std::vector<float> compute_R_field(const std::vector<float>& theta, 
                                     uint32_t Nx, uint32_t Ny) {
    std::vector<float> R_field(Nx * Ny);
    
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            
            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;
            
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int nx = static_cast<int>(x) + dx;
                    int ny = static_cast<int>(y) + dy;
                    
                    if (nx < 0) nx += Nx;
                    if (nx >= static_cast<int>(Nx)) nx -= Nx;
                    if (ny < 0) ny += Ny;
                    if (ny >= static_cast<int>(Ny)) ny -= Ny;
                    
                    uint32_t neighbor_idx = ny * Nx + nx;
                    sum_cos += std::cos(theta[neighbor_idx]);
                    sum_sin += std::sin(theta[neighbor_idx]);
                    count++;
                }
            }
            
            float avg_cos = sum_cos / count;
            float avg_sin = sum_sin / count;
            R_field[idx] = std::sqrt(avg_cos * avg_cos + avg_sin * avg_sin);
        }
    }
    
    return R_field;
}

void cpu_kuramoto_step(std::vector<float>& theta,
                       const std::vector<float>& omega,
                       uint32_t Nx, uint32_t Ny,
                       float K, float dt, float damping) {
    std::vector<float> theta_new(theta.size());
    
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            float theta_i = theta[idx];
            float omega_i = omega[idx];
            
            float coupling = 0.0f;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    
                    int nx = static_cast<int>(x) + dx;
                    int ny = static_cast<int>(y) + dy;
                    
                    if (nx < 0) nx += Nx;
                    if (nx >= static_cast<int>(Nx)) nx -= Nx;
                    if (ny < 0) ny += Ny;
                    if (ny >= static_cast<int>(Ny)) ny -= Ny;
                    
                    uint32_t neighbor_idx = ny * Nx + nx;
                    coupling += std::sin(theta[neighbor_idx] - theta_i);
                }
            }
            coupling *= (K / 8.0f);
            
            float damping_force = -damping * std::sin(theta_i);
            float total_force = omega_i + coupling + damping_force;
            
            theta_new[idx] = theta_i + dt * total_force;
            
            float temp = theta_new[idx] + M_PI;
            float wrapped = temp - 2.0f * M_PI * std::floor(temp / (2.0f * M_PI));
            theta_new[idx] = wrapped - M_PI;
        }
    }
    
    theta = theta_new;
}

int main() {
    std::cout << "=== Mass Quantization Analysis ===" << std::endl;
    std::cout << "Measuring defect mass distribution\n" << std::endl;
    
    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const int N_samples = 20;  // Collect defects from 20 timesteps
    const float dt = 0.01f;
    const float K = 1.0f;
    const float damping = 0.1f;
    const float Delta = 2.5f;
    
    std::cout << "Parameters:" << std::endl;
    std::cout << "  Grid: " << Nx << " × " << Ny << std::endl;
    std::cout << "  Sample timesteps: " << N_samples << std::endl;
    std::cout << "  Delta: " << Delta << std::endl;
    std::cout << std::endl;
    
    // Initialize random state
    std::vector<float> theta(Nx * Ny);
    std::vector<float> omega(Nx * Ny);
    
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> phase_dist(-M_PI, M_PI);
    std::uniform_real_distribution<float> omega_dist(-0.25f, 0.25f);
    
    for (size_t i = 0; i < theta.size(); i++) {
        theta[i] = phase_dist(rng);
        omega[i] = omega_dist(rng);
    }
    
    // Collect defects from multiple timesteps
    std::vector<float> all_masses;
    std::vector<float> all_R_values;
    std::vector<float> all_windings;
    
    std::cout << "Collecting defect samples..." << std::endl;
    
    for (int sample = 0; sample < N_samples; sample++) {
        // Evolve a bit between samples
        for (int step = 0; step < 5; step++) {
            cpu_kuramoto_step(theta, omega, Nx, Ny, K, dt, damping);
        }
        
        // Detect defects
        std::vector<float> R_field = compute_R_field(theta, Nx, Ny);
        std::vector<Defect> defects = detect_defects(theta, R_field, Nx, Ny, Delta);
        
        // Collect properties
        for (const auto& d : defects) {
            all_masses.push_back(d.m_local);
            all_R_values.push_back(d.R_local);
            all_windings.push_back(d.winding);
        }
        
        if (sample % 5 == 0) {
            std::cout << "  Sample " << sample << ": " << defects.size() << " defects" << std::endl;
        }
    }
    
    std::cout << "\nTotal defects collected: " << all_masses.size() << std::endl;
    
    if (all_masses.empty()) {
        std::cout << "⚠ No defects found" << std::endl;
        return 1;
    }
    
    // Statistical analysis
    float m_min = *std::min_element(all_masses.begin(), all_masses.end());
    float m_max = *std::max_element(all_masses.begin(), all_masses.end());
    
    float m_sum = 0.0f;
    for (float m : all_masses) m_sum += m;
    float m_mean = m_sum / all_masses.size();
    
    float m_var = 0.0f;
    for (float m : all_masses) {
        float diff = m - m_mean;
        m_var += diff * diff;
    }
    m_var /= all_masses.size();
    float m_std = std::sqrt(m_var);
    
    std::cout << "\n=== MASS DISTRIBUTION STATISTICS ===" << std::endl;
    std::cout << "Mass m = Δ·R at defect cores:" << std::endl;
    std::cout << "  Min: " << m_min << std::endl;
    std::cout << "  Max: " << m_max << std::endl;
    std::cout << "  Mean: " << m_mean << " ± " << m_std << std::endl;
    std::cout << "  Range: " << (m_max - m_min) << std::endl;
    std::cout << "  Coefficient of variation: " << (m_std / m_mean * 100) << "%" << std::endl;
    
    // Create histogram
    const int n_bins = 50;
    std::vector<int> histogram(n_bins, 0);
    float bin_width = (m_max - m_min) / n_bins;
    
    for (float m : all_masses) {
        int bin = static_cast<int>((m - m_min) / bin_width);
        if (bin >= n_bins) bin = n_bins - 1;
        if (bin < 0) bin = 0;
        histogram[bin]++;
    }
    
    // Save histogram
    std::ofstream f_hist("output/mass_histogram.dat");
    f_hist << "# bin_center count\n";
    for (int i = 0; i < n_bins; i++) {
        float bin_center = m_min + (i + 0.5f) * bin_width;
        f_hist << bin_center << " " << histogram[i] << "\n";
    }
    f_hist.close();
    std::cout << "\nSaved: output/mass_histogram.dat" << std::endl;
    
    // Save all masses
    std::ofstream f_masses("output/defect_masses.dat");
    f_masses << "# m R winding\n";
    for (size_t i = 0; i < all_masses.size(); i++) {
        f_masses << all_masses[i] << " " << all_R_values[i] << " " << all_windings[i] << "\n";
    }
    f_masses.close();
    std::cout << "Saved: output/defect_masses.dat" << std::endl;
    
    // Look for peaks (simple peak detection)
    std::cout << "\n=== PEAK DETECTION ===" << std::endl;
    std::cout << "Searching for quantization peaks..." << std::endl;
    
    std::vector<int> peaks;
    for (int i = 1; i < n_bins - 1; i++) {
        if (histogram[i] > histogram[i-1] && histogram[i] > histogram[i+1]) {
            if (histogram[i] > 20) {  // Significant peak threshold
                peaks.push_back(i);
            }
        }
    }
    
    if (peaks.size() >= 2) {
        std::cout << "Found " << peaks.size() << " peaks:" << std::endl;
        for (size_t i = 0; i < peaks.size(); i++) {
            int bin = peaks[i];
            float mass_peak = m_min + (bin + 0.5f) * bin_width;
            std::cout << "  Peak " << (i+1) << ": m = " << mass_peak 
                      << " (count: " << histogram[bin] << ")" << std::endl;
        }
        
        // Check if peaks are evenly spaced
        if (peaks.size() >= 3) {
            std::vector<float> spacings;
            for (size_t i = 1; i < peaks.size(); i++) {
                float spacing = (peaks[i] - peaks[i-1]) * bin_width;
                spacings.push_back(spacing);
            }
            
            float spacing_mean = 0.0f;
            for (float s : spacings) spacing_mean += s;
            spacing_mean /= spacings.size();
            
            float spacing_var = 0.0f;
            for (float s : spacings) {
                float diff = s - spacing_mean;
                spacing_var += diff * diff;
            }
            spacing_var /= spacings.size();
            float spacing_std = std::sqrt(spacing_var);
            
            std::cout << "\nPeak spacing: " << spacing_mean << " ± " << spacing_std << std::endl;
            
            if (spacing_std / spacing_mean < 0.2f) {
                std::cout << "✓ QUANTIZATION DETECTED" << std::endl;
                std::cout << "  → Masses are quantized with Δm ≈ " << spacing_mean << std::endl;
            } else {
                std::cout << "✗ No clear quantization (spacing is irregular)" << std::endl;
            }
        }
    } else {
        std::cout << "Insufficient peaks for quantization analysis" << std::endl;
        std::cout << "  → Distribution appears continuous (not quantized)" << std::endl;
    }
    
    std::cout << "\n=== CONCLUSION ===" << std::endl;
    std::cout << "Defect mass distribution analyzed" << std::endl;
    std::cout << "See output/mass_histogram.dat for full distribution" << std::endl;
    
    return 0;
}
