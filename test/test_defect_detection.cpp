/**
 * Defect Detection and Analysis
 *
 * Purpose: Detect topological defects in phase field θ(x,y)
 *
 * Theory:
 * - Defects are vortices where phase winds by ±2π around a point
 * - Winding number: w = (1/2π) ∮ ∇θ·dl
 * - w = +1: counterclockwise vortex
 * - w = -1: clockwise vortex
 * - w = 0: no defect
 *
 * Detection method:
 * - For each plaquette (2×2 cell), compute phase circulation
 * - Circulation = Σ Δθ around edges (unwrapped)
 * - If |circulation| > π: defect detected
 *
 * Analysis:
 * 1. Count defects vs time
 * 2. Measure defect spacing
 * 3. Compute local mass m(defect) = Δ·R(defect)
 * 4. Check for mass quantization
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdint>

struct Defect {
    int x, y;           // Grid position
    float winding;      // Winding number (should be ±1)
    float R_local;      // Local synchronization at defect
    float m_local;      // Local mass m = Δ·R
};

// Unwrap phase difference to [-π, π]
float phase_diff(float theta1, float theta2) {
    float diff = theta1 - theta2;
    // Wrap to [-π, π]
    while (diff > M_PI) diff -= 2.0f * M_PI;
    while (diff < -M_PI) diff += 2.0f * M_PI;
    return diff;
}

// Detect defects using plaquette circulation method
std::vector<Defect> detect_defects(const std::vector<float>& theta,
                                     const std::vector<float>& R_field,
                                     uint32_t Nx, uint32_t Ny,
                                     float Delta) {
    std::vector<Defect> defects;
    
    // Scan all plaquettes (2×2 cells)
    for (uint32_t y = 0; y < Ny - 1; y++) {
        for (uint32_t x = 0; x < Nx - 1; x++) {
            // Get corners of plaquette
            //  (x,y) ---> (x+1,y)
            //    |           |
            //    v           v
            // (x,y+1) -> (x+1,y+1)
            
            uint32_t idx00 = y * Nx + x;
            uint32_t idx10 = y * Nx + (x + 1);
            uint32_t idx01 = (y + 1) * Nx + x;
            uint32_t idx11 = (y + 1) * Nx + (x + 1);
            
            float theta00 = theta[idx00];
            float theta10 = theta[idx10];
            float theta01 = theta[idx01];
            float theta11 = theta[idx11];
            
            // Compute circulation around plaquette
            // Path: (0,0) → (1,0) → (1,1) → (0,1) → (0,0)
            float circ = 0.0f;
            circ += phase_diff(theta10, theta00);  // Right edge
            circ += phase_diff(theta11, theta10);  // Top edge
            circ += phase_diff(theta01, theta11);  // Left edge
            circ += phase_diff(theta00, theta01);  // Bottom edge
            
            // Winding number = circulation / 2π
            float winding = circ / (2.0f * M_PI);
            
            // Detect defect if |winding| > 0.5 (should be ±1)
            if (std::abs(winding) > 0.5f) {
                Defect d;
                d.x = x;
                d.y = y;
                d.winding = winding;
                
                // Compute local R at defect center
                d.R_local = (R_field[idx00] + R_field[idx10] + 
                             R_field[idx01] + R_field[idx11]) / 4.0f;
                
                d.m_local = Delta * d.R_local;
                
                defects.push_back(d);
            }
        }
    }
    
    return defects;
}

// Compute R field (sync field)
std::vector<float> compute_R_field(const std::vector<float>& theta, 
                                     uint32_t Nx, uint32_t Ny) {
    std::vector<float> R_field(Nx * Ny);
    
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            
            float sum_cos = 0.0f;
            float sum_sin = 0.0f;
            int count = 0;
            
            // 3×3 neighborhood
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int nx = static_cast<int>(x) + dx;
                    int ny = static_cast<int>(y) + dy;
                    
                    // Periodic boundaries
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

int main() {
    std::cout << "=== Defect Detection and Analysis ===" << std::endl;
    std::cout << "Detecting topological defects in phase field\n" << std::endl;
    
    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const float Delta = 2.5f;
    
    std::cout << "Loading final phase field from 10k-step simulation..." << std::endl;
    
    // Load theta field from previous simulation
    std::vector<float> theta(Nx * Ny);
    std::ifstream f_theta("output/theta.dat", std::ios::binary);
    if (!f_theta.is_open()) {
        std::cerr << "Error: Could not open output/theta.dat" << std::endl;
        std::cerr << "Run test_smft_compute_only first to generate data." << std::endl;
        return 1;
    }
    
    f_theta.read(reinterpret_cast<char*>(theta.data()), theta.size() * sizeof(float));
    f_theta.close();
    
    std::cout << "Loaded " << theta.size() << " phase values" << std::endl;
    
    // Compute R field
    std::cout << "Computing synchronization field R(x,y)..." << std::endl;
    std::vector<float> R_field = compute_R_field(theta, Nx, Ny);
    
    // Detect defects
    std::cout << "Scanning for topological defects..." << std::endl;
    std::vector<Defect> defects = detect_defects(theta, R_field, Nx, Ny, Delta);
    
    std::cout << "\n=== DETECTION RESULTS ===" << std::endl;
    std::cout << "Total defects found: " << defects.size() << std::endl;
    
    if (defects.empty()) {
        std::cout << "\n✓ No defects detected" << std::endl;
        std::cout << "  → System is highly synchronized (R → 1)" << std::endl;
        std::cout << "  → Phase field is smooth" << std::endl;
        return 0;
    }
    
    // Classify defects by winding number
    int n_plus = 0;   // w = +1
    int n_minus = 0;  // w = -1
    
    for (const auto& d : defects) {
        if (d.winding > 0) n_plus++;
        else n_minus++;
    }
    
    std::cout << "  Positive vortices (w=+1): " << n_plus << std::endl;
    std::cout << "  Negative vortices (w=-1): " << n_minus << std::endl;
    std::cout << "  Net topological charge: " << (n_plus - n_minus) << std::endl;
    
    // Analyze defect properties
    std::vector<float> R_at_defects;
    std::vector<float> m_at_defects;
    
    for (const auto& d : defects) {
        R_at_defects.push_back(d.R_local);
        m_at_defects.push_back(d.m_local);
    }
    
    float R_min = *std::min_element(R_at_defects.begin(), R_at_defects.end());
    float R_max = *std::max_element(R_at_defects.begin(), R_at_defects.end());
    float R_mean = 0.0f;
    for (float r : R_at_defects) R_mean += r;
    R_mean /= R_at_defects.size();
    
    std::cout << "\n=== DEFECT PROPERTIES ===" << std::endl;
    std::cout << "R at defect cores:" << std::endl;
    std::cout << "  Min: " << R_min << std::endl;
    std::cout << "  Max: " << R_max << std::endl;
    std::cout << "  Mean: " << R_mean << std::endl;
    
    // Check if defects are mass holes (R → 0)
    if (R_mean < 0.3f) {
        std::cout << "  → Defects are MASS HOLES (low R)" << std::endl;
    } else if (R_mean < 0.7f) {
        std::cout << "  → Defects have REDUCED MASS" << std::endl;
    } else {
        std::cout << "  → Defects are NOT mass singularities (high R)" << std::endl;
    }
    
    // Mass distribution at defects
    float m_min = *std::min_element(m_at_defects.begin(), m_at_defects.end());
    float m_max = *std::max_element(m_at_defects.begin(), m_at_defects.end());
    float m_mean = 0.0f;
    for (float m : m_at_defects) m_mean += m;
    m_mean /= m_at_defects.size();
    
    std::cout << "\nMass m = Δ·R at defects:" << std::endl;
    std::cout << "  Min: " << m_min << std::endl;
    std::cout << "  Max: " << m_max << std::endl;
    std::cout << "  Mean: " << m_mean << std::endl;
    
    // Save defect positions and properties
    std::ofstream f_defects("output/defects.dat");
    f_defects << "# x y winding R_local m_local\n";
    for (const auto& d : defects) {
        f_defects << d.x << " " << d.y << " " 
                  << d.winding << " " << d.R_local << " " << d.m_local << "\n";
    }
    f_defects.close();
    std::cout << "\nSaved: output/defects.dat" << std::endl;
    
    // Defect density
    float defect_density = float(defects.size()) / (Nx * Ny);
    float avg_spacing = std::sqrt(1.0f / defect_density);
    
    std::cout << "\n=== SPATIAL STATISTICS ===" << std::endl;
    std::cout << "Defect density: " << defect_density << " per site" << std::endl;
    std::cout << "Average spacing: " << avg_spacing << " lattice units" << std::endl;
    
    // Print first 10 defects
    std::cout << "\nFirst 10 defects:" << std::endl;
    std::cout << std::setw(4) << "x" << std::setw(4) << "y" 
              << std::setw(10) << "winding" << std::setw(10) << "R" 
              << std::setw(10) << "m" << std::endl;
    for (size_t i = 0; i < std::min(size_t(10), defects.size()); i++) {
        const auto& d = defects[i];
        std::cout << std::setw(4) << d.x << std::setw(4) << d.y
                  << std::fixed << std::setprecision(3)
                  << std::setw(10) << d.winding
                  << std::setw(10) << d.R_local
                  << std::setw(10) << d.m_local << std::endl;
    }
    
    std::cout << "\n=== CONCLUSION ===" << std::endl;
    std::cout << "✓ Defect detection complete" << std::endl;
    std::cout << "  Found " << defects.size() << " topological defects" << std::endl;
    std::cout << "  Next: Analyze mass quantization" << std::endl;
    
    return 0;
}
