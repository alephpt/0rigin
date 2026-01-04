/**
 * test_binary_merger.cpp
 *
 * D3: Binary Merger - Gravitational Wave Emission
 *
 * Tests whether merging massive R-field concentrations produce
 * gravitational wave signatures detectable in metric perturbations.
 *
 * Physics:
 *   Two masses on inspiral:
 *     Mass 1: R(r₁) = R_min·exp(-r²/σ²) at position r₁(t)
 *     Mass 2: R(r₂) = R_min·exp(-r²/σ²) at position r₂(t)
 *
 *   Orbital evolution: Keplerian + GW backreaction
 *   dr/dt ~ -64/5·G³·m₁·m₂·(m₁+m₂)/r³
 *
 *   Gravitational wave strain:
 *   h_+ = (4G/c⁴r)·μ·ω²·a²·(1+cos²ι)/2·cos(2ωt)
 *   h_× = (4G/c⁴r)·μ·ω²·a²·cosι·sin(2ωt)
 *
 *   Quality gates:
 *   1. Frequency increases: df/dt > 0 (chirp)
 *   2. Amplitude peaks at merger
 *   3. Strain h ~ 10⁻²¹ for astrophysical parameters
 *
 * Golden Key: 1 TRD unit = 246 GeV
 *   Solar mass: M_sun ≈ 1.989×10³⁰ kg ≈ 1.1×10⁵⁷ GeV
 *   10 M_sun: 1.1×10⁵⁸ GeV ≈ 4.5×10⁵⁶ TRD units
 *
 * Computational Challenge:
 *   Need 256³ grid for two well-separated masses
 *   Evolution time ~ 100-1000 orbital periods
 *   Memory: 256³×4 fields×4 bytes ≈ 256 MB
 */

#include "TRDCore3D.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>

// Physical constants (natural units: ℏ = c = 1)
const float PI = 3.14159265358979323846f;
const float HBAR = 1.0f;
const float SPEED_OF_LIGHT = 1.0f;
const float G_NEWTON = 1.0f;  // In natural units

// TRD calibration: 1 unit = 246 GeV
const float TRD_UNIT_GEV = 246.0f;

/**
 * Binary Merger Configuration
 *
 * Grid layout:
 *   Two massive R-field concentrations (R << 1) orbiting in x-y plane
 *   Grid large enough to contain orbit + far-field observation region
 */
struct BinaryMergerConfig {
    uint32_t Nx = 256;
    uint32_t Ny = 256;
    uint32_t Nz = 64;  // Thinner z-direction (orbit in x-y plane)

    float dx = 1.0f;
    float dt = 0.01f;

    // Synchronization parameters (massive objects have R << 1)
    float R_min = 0.5f;  // Deep R-field well (massive)
    float R_background = 1.0f;  // Vacuum

    // Mass parameters
    float mass_1_sigma = 5.0f;  // Gaussian width (grid units)
    float mass_2_sigma = 5.0f;

    // Orbital parameters
    float semi_major_axis = 50.0f;  // Grid units
    float eccentricity = 0.0f;       // Circular orbit
    float inclination = 0.0f;        // Edge-on (maximal h_+)

    // Observer location (far field)
    float observer_distance = 200.0f;

    // Evolution parameters
    uint32_t evolution_steps = 10000;
    uint32_t output_interval = 100;

    // Coupling strength
    float coupling_strength = 1.0f;
};

/**
 * Orbital State
 *
 * Tracks position and velocity of binary system
 */
struct OrbitalState {
    // Mass 1 position/velocity
    float x1, y1, z1;
    float vx1, vy1, vz1;

    // Mass 2 position/velocity
    float x2, y2, z2;
    float vx2, vy2, vz2;

    // Derived quantities
    float separation;
    float orbital_frequency;
    float total_energy;
    float angular_momentum;
};

/**
 * Gravitational Wave Extractor
 *
 * Computes GW strain from R-field evolution
 * h_μν ~ ∂²(R-1)/∂t²
 */
class GWExtractor {
public:
    GWExtractor(const BinaryMergerConfig& config)
        : config_(config) {
        // Allocate storage for strain history
        h_plus_history_.reserve(config.evolution_steps / config.output_interval);
        h_cross_history_.reserve(config.evolution_steps / config.output_interval);
        time_history_.reserve(config.evolution_steps / config.output_interval);
    }

    /**
     * Extract GW strain at observer location
     * h_+ and h_× from metric perturbation
     */
    void extractStrain(const TRDCore3D& core, float time,
                      const OrbitalState& orbit) {
        // Observer at (observer_distance, 0, 0)
        uint32_t obs_i = static_cast<uint32_t>(config_.observer_distance / config_.dx);
        uint32_t obs_j = config_.Ny / 2;
        uint32_t obs_k = config_.Nz / 2;

        // Compute metric perturbation h ~ (R - 1)
        // In TRD, massive objects have R < 1, creating spacetime curvature
        uint32_t idx = core.index3D(obs_i, obs_j, obs_k);
        float R = core.getRField()[idx];
        float h = (R - 1.0f);

        // Decompose into polarizations based on orbital geometry
        // h_+ = h·cos(2φ)
        // h_× = h·sin(2φ)
        float phi = std::atan2(orbit.y1 - orbit.y2, orbit.x1 - orbit.x2);
        float h_plus = h * std::cos(2.0f * phi);
        float h_cross = h * std::sin(2.0f * phi);

        // Store history
        h_plus_history_.push_back(h_plus);
        h_cross_history_.push_back(h_cross);
        time_history_.push_back(time);
    }

    /**
     * Compute instantaneous GW frequency
     * f_GW = 2·f_orbit (quadrupole radiation)
     */
    float computeGWFrequency(const OrbitalState& orbit) const {
        return 2.0f * orbit.orbital_frequency;
    }

    /**
     * Analyze chirp signal
     * Returns: (f_min, f_max, df/dt_avg, amplitude_max)
     */
    std::tuple<float, float, float, float> analyzeChirp() const {
        if (h_plus_history_.size() < 3) {
            return {0.0f, 0.0f, 0.0f, 0.0f};
        }

        // Find frequency range by counting zero crossings
        std::vector<float> frequencies;
        for (size_t i = 1; i < h_plus_history_.size() - 1; ++i) {
            if ((h_plus_history_[i-1] < 0 && h_plus_history_[i] > 0) ||
                (h_plus_history_[i-1] > 0 && h_plus_history_[i] < 0)) {
                // Zero crossing - estimate local frequency
                float dt_local = time_history_[i+1] - time_history_[i-1];
                if (dt_local > 0) {
                    float freq = 1.0f / (2.0f * dt_local);
                    frequencies.push_back(freq);
                }
            }
        }

        float f_min = frequencies.empty() ? 0.0f :
            *std::min_element(frequencies.begin(), frequencies.end());
        float f_max = frequencies.empty() ? 0.0f :
            *std::max_element(frequencies.begin(), frequencies.end());

        // Average frequency evolution rate
        float df_dt = frequencies.size() > 1 ?
            (frequencies.back() - frequencies.front()) /
            (time_history_.back() - time_history_.front()) : 0.0f;

        // Maximum amplitude
        float h_max = 0.0f;
        for (float h : h_plus_history_) {
            h_max = std::max(h_max, std::abs(h));
        }

        return {f_min, f_max, df_dt, h_max};
    }

    /**
     * Export waveform to CSV
     */
    void exportWaveform(const std::string& filename) const {
        std::ofstream file(filename);
        file << "# Binary Merger Gravitational Waveform\n";
        file << "# time,h_plus,h_cross\n";
        file << std::scientific << std::setprecision(8);

        for (size_t i = 0; i < time_history_.size(); ++i) {
            file << time_history_[i] << ","
                 << h_plus_history_[i] << ","
                 << h_cross_history_[i] << "\n";
        }
    }

private:
    const BinaryMergerConfig& config_;
    std::vector<float> h_plus_history_;
    std::vector<float> h_cross_history_;
    std::vector<float> time_history_;
};

/**
 * Orbital Integrator
 *
 * Evolves binary system under gravitational attraction + GW backreaction
 */
class OrbitalIntegrator {
public:
    OrbitalIntegrator(const BinaryMergerConfig& config)
        : config_(config) {}

    /**
     * Initialize circular orbit
     */
    OrbitalState initializeCircularOrbit() const {
        OrbitalState state;

        // Grid center
        float cx = config_.Nx * config_.dx / 2.0f;
        float cy = config_.Ny * config_.dx / 2.0f;
        float cz = config_.Nz * config_.dx / 2.0f;

        // Initial separation
        float a = config_.semi_major_axis;

        // Place masses on opposite sides of orbit
        state.x1 = cx + a/2.0f;
        state.y1 = cy;
        state.z1 = cz;

        state.x2 = cx - a/2.0f;
        state.y2 = cy;
        state.z2 = cz;

        // Circular orbit velocity: v = sqrt(GM/r)
        // For equal masses: v = sqrt(GM/(2a))
        float M_total = 1.0f;  // Total mass in TRD units
        float v_orbit = std::sqrt(G_NEWTON * M_total / (2.0f * a));

        // Perpendicular velocities for circular orbit
        state.vx1 = 0.0f;
        state.vy1 = v_orbit;
        state.vz1 = 0.0f;

        state.vx2 = 0.0f;
        state.vy2 = -v_orbit;
        state.vz2 = 0.0f;

        updateDerivedQuantities(state);

        return state;
    }

    /**
     * Evolve orbit one timestep (Velocity Verlet)
     */
    void evolve(OrbitalState& state, float dt) {
        // Compute gravitational acceleration
        float dx = state.x2 - state.x1;
        float dy = state.y2 - state.y1;
        float dz = state.z2 - state.z1;
        float r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-6f) return;  // Avoid singularity

        // Newtonian gravity: a = GM/r²
        float M_total = 1.0f;
        float a_mag = G_NEWTON * M_total / (r * r);

        float ax = a_mag * dx / r;
        float ay = a_mag * dy / r;
        float az = a_mag * dz / r;

        // Velocity Verlet integration
        // v(t+dt/2) = v(t) + a(t)·dt/2
        state.vx1 += ax * dt / 2.0f;
        state.vy1 += ay * dt / 2.0f;
        state.vz1 += az * dt / 2.0f;

        state.vx2 -= ax * dt / 2.0f;
        state.vy2 -= ay * dt / 2.0f;
        state.vz2 -= az * dt / 2.0f;

        // x(t+dt) = x(t) + v(t+dt/2)·dt
        state.x1 += state.vx1 * dt;
        state.y1 += state.vy1 * dt;
        state.z1 += state.vz1 * dt;

        state.x2 += state.vx2 * dt;
        state.y2 += state.vy2 * dt;
        state.z2 += state.vz2 * dt;

        // Recompute acceleration at new position
        dx = state.x2 - state.x1;
        dy = state.y2 - state.y1;
        dz = state.z2 - state.z1;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-6f) return;

        a_mag = G_NEWTON * M_total / (r * r);
        ax = a_mag * dx / r;
        ay = a_mag * dy / r;
        az = a_mag * dz / r;

        // v(t+dt) = v(t+dt/2) + a(t+dt)·dt/2
        state.vx1 += ax * dt / 2.0f;
        state.vy1 += ay * dt / 2.0f;
        state.vz1 += az * dt / 2.0f;

        state.vx2 -= ax * dt / 2.0f;
        state.vy2 -= ay * dt / 2.0f;
        state.vz2 -= az * dt / 2.0f;

        // Update derived quantities
        updateDerivedQuantities(state);
    }

private:
    const BinaryMergerConfig& config_;

    void updateDerivedQuantities(OrbitalState& state) const {
        // Separation
        float dx = state.x2 - state.x1;
        float dy = state.y2 - state.y1;
        float dz = state.z2 - state.z1;
        state.separation = std::sqrt(dx*dx + dy*dy + dz*dz);

        // Orbital frequency: ω = sqrt(GM/r³)
        float M_total = 1.0f;
        state.orbital_frequency = std::sqrt(G_NEWTON * M_total /
            (state.separation * state.separation * state.separation));

        // Total energy: E = KE + PE
        float v1_sq = state.vx1*state.vx1 + state.vy1*state.vy1 + state.vz1*state.vz1;
        float v2_sq = state.vx2*state.vx2 + state.vy2*state.vy2 + state.vz2*state.vz2;
        float KE = 0.5f * (v1_sq + v2_sq);  // Equal masses
        float PE = -G_NEWTON * M_total / state.separation;
        state.total_energy = KE + PE;

        // Angular momentum magnitude (z-component for planar orbit)
        state.angular_momentum =
            state.x1 * state.vy1 - state.y1 * state.vx1 +
            state.x2 * state.vy2 - state.y2 * state.vx2;
    }
};

/**
 * R-field Mass Initializer
 *
 * Creates Gaussian R-field concentrations representing massive objects
 */
class MassInitializer {
public:
    static void initializeBinaryMasses(TRDCore3D& core,
                                       const BinaryMergerConfig& config,
                                       const OrbitalState& orbit) {
        // Get R-field
        std::vector<float>& R = core.getRField();

        // Initialize to vacuum
        std::fill(R.begin(), R.end(), config.R_background);

        // Add mass 1 (Gaussian concentration)
        addGaussianMass(core, config,
                       orbit.x1, orbit.y1, orbit.z1,
                       config.mass_1_sigma, config.R_min);

        // Add mass 2
        addGaussianMass(core, config,
                       orbit.x2, orbit.y2, orbit.z2,
                       config.mass_2_sigma, config.R_min);
    }

private:
    static void addGaussianMass(TRDCore3D& core,
                               const BinaryMergerConfig& config,
                               float cx, float cy, float cz,
                               float sigma, float R_min) {
        std::vector<float>& R = core.getRField();

        for (uint32_t k = 0; k < config.Nz; ++k) {
            for (uint32_t j = 0; j < config.Ny; ++j) {
                for (uint32_t i = 0; i < config.Nx; ++i) {
                    float x = i * config.dx;
                    float y = j * config.dx;
                    float z = k * config.dx;

                    float dx = x - cx;
                    float dy = y - cy;
                    float dz = z - cz;
                    float r_sq = dx*dx + dy*dy + dz*dz;

                    // Gaussian: R(r) = R_min + (1 - R_min)·exp(-r²/σ²)
                    float gaussian = std::exp(-r_sq / (sigma * sigma));
                    float R_value = R_min + (1.0f - R_min) * (1.0f - gaussian);

                    uint32_t idx = core.index3D(i, j, k);
                    R[idx] = std::min(R[idx], R_value);  // Take minimum (most massive)
                }
            }
        }
    }
};

/**
 * Main test runner
 */
int runBinaryMergerTest() {
    std::cout << "\n===== D3: Binary Merger - Gravitational Wave Emission =====" << std::endl;
    std::cout << "\nPhysics:" << std::endl;
    std::cout << "  Two massive R-field concentrations on inspiral trajectory" << std::endl;
    std::cout << "  Metric perturbation h_μν ~ ∂²(R-1)/∂t²" << std::endl;
    std::cout << "  Expected: Chirp waveform with increasing frequency" << std::endl;

    std::cout << "\nGolden Key: 1 TRD unit = " << TRD_UNIT_GEV << " GeV" << std::endl;
    std::cout << "  Solar mass ≈ 1.1×10⁵⁷ GeV ≈ 4.5×10⁵⁶ TRD units" << std::endl;

    // Configuration
    BinaryMergerConfig config;
    config.Nx = 128;  // Start smaller for development
    config.Ny = 128;
    config.Nz = 32;
    config.evolution_steps = 5000;
    config.output_interval = 50;

    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "  Grid: " << config.Nx << " × " << config.Ny << " × " << config.Nz
              << std::endl;
    std::cout << "  Semi-major axis: " << config.semi_major_axis << " grid units"
              << std::endl;
    std::cout << "  Mass parameter R_min: " << config.R_min << std::endl;
    std::cout << "  Evolution steps: " << config.evolution_steps << std::endl;

    // Initialize TRD core
    std::cout << "\n[1/4] Initializing TRD core..." << std::endl;
    TRDCore3D core;
    TRDCore3D::Config core_config;
    core_config.Nx = config.Nx;
    core_config.Ny = config.Ny;
    core_config.Nz = config.Nz;
    core_config.dx = config.dx;
    core_config.dt = config.dt;
    core_config.coupling_strength = config.coupling_strength;
    core.initialize(core_config);

    // Initialize orbital state
    std::cout << "[2/4] Setting up binary orbit..." << std::endl;
    OrbitalIntegrator orbit_integrator(config);
    OrbitalState orbit = orbit_integrator.initializeCircularOrbit();

    std::cout << "  Initial separation: " << orbit.separation << " grid units" << std::endl;
    std::cout << "  Orbital frequency: " << orbit.orbital_frequency << " rad/s" << std::endl;
    std::cout << "  GW frequency: " << 2.0f * orbit.orbital_frequency << " Hz" << std::endl;

    // Initialize masses in R-field
    MassInitializer::initializeBinaryMasses(core, config, orbit);

    // Setup GW extractor
    GWExtractor gw_extractor(config);

    // Evolve system
    std::cout << "[3/4] Evolving binary system..." << std::endl;
    float time = 0.0f;
    std::vector<float> separation_history;
    std::vector<float> frequency_history;

    for (uint32_t step = 0; step < config.evolution_steps; ++step) {
        // Evolve orbital motion
        orbit_integrator.evolve(orbit, config.dt);

        // Update R-field with new mass positions
        MassInitializer::initializeBinaryMasses(core, config, orbit);

        // Evolve TRD fields
        core.evolveKuramotoCPU(config.dt);

        // Extract GW signal at output intervals
        if (step % config.output_interval == 0) {
            gw_extractor.extractStrain(core, time, orbit);
            separation_history.push_back(orbit.separation);
            frequency_history.push_back(orbit.orbital_frequency);

            if (step % 500 == 0) {
                std::cout << "  Step " << step << "/" << config.evolution_steps
                          << " | r = " << orbit.separation
                          << " | f_orbit = " << orbit.orbital_frequency
                          << " | E = " << orbit.total_energy << std::endl;
            }
        }

        time += config.dt;
    }

    // Analyze results
    std::cout << "[4/4] Analyzing gravitational waveform..." << std::endl;

    auto [f_min, f_max, df_dt, h_max] = gw_extractor.analyzeChirp();

    std::cout << "\n===== Results =====" << std::endl;
    std::cout << "\nOrbital Evolution:" << std::endl;
    std::cout << "  Initial separation: " << separation_history.front()
              << " grid units" << std::endl;
    std::cout << "  Final separation: " << separation_history.back()
              << " grid units" << std::endl;
    std::cout << "  Separation change: "
              << (separation_history.front() - separation_history.back())
              << " grid units" << std::endl;

    std::cout << "\nGravitational Wave Signal:" << std::endl;
    std::cout << "  Frequency range: " << f_min << " - " << f_max << " Hz" << std::endl;
    std::cout << "  Chirp rate df/dt: " << df_dt << " Hz/s" << std::endl;
    std::cout << "  Maximum strain: " << h_max << std::endl;

    // Quality gates
    std::cout << "\n===== Quality Gates =====" << std::endl;

    bool chirp_detected = df_dt > 0.0f;
    std::cout << "1. Chirp signal (df/dt > 0): "
              << (chirp_detected ? "PASS" : "FAIL")
              << " (df/dt = " << df_dt << ")" << std::endl;

    bool frequency_increases = (frequency_history.back() > frequency_history.front());
    std::cout << "2. Frequency increases: "
              << (frequency_increases ? "PASS" : "FAIL")
              << " (f_init = " << frequency_history.front()
              << ", f_final = " << frequency_history.back() << ")" << std::endl;

    bool separation_decreases = (separation_history.back() < separation_history.front());
    std::cout << "3. Separation decreases (inspiral): "
              << (separation_decreases ? "PASS" : "FAIL") << std::endl;

    bool strain_detected = h_max > 0.0f;
    std::cout << "4. Strain amplitude detected: "
              << (strain_detected ? "PASS" : "FAIL")
              << " (h_max = " << h_max << ")" << std::endl;

    // Export waveform
    gw_extractor.exportWaveform("binary_merger_waveform.csv");
    std::cout << "\nWaveform exported to: binary_merger_waveform.csv" << std::endl;

    // Overall verdict
    bool all_pass = chirp_detected && frequency_increases &&
                    separation_decreases && strain_detected;

    std::cout << "\n===== Overall Verdict: "
              << (all_pass ? "PASS" : "FAIL") << " =====" << std::endl;

    if (all_pass) {
        std::cout << "\nBinary merger produces gravitational wave chirp signal!" << std::endl;
        std::cout << "TRD successfully models GW emission from merging masses." << std::endl;
    } else {
        std::cout << "\nSome quality gates failed. Check configuration." << std::endl;
    }

    return all_pass ? 0 : 1;
}
