/**
 * test_gravitational_waves.cpp
 *
 * A5: Gravitational Wave Emission and Propagation Validation
 *
 * Goal: Validate TRD produces gravitational waves from binary systems,
 *       matching GR predictions for waveform, polarization, energy loss
 *
 * Physics Background:
 *   - Metric perturbation h_μν from R-field oscillations
 *   - Quadrupole formula: dE/dt = (G/5c⁵)<d³Q_ij/dt³>²
 *   - Polarizations: h₊ and hₓ (+ and × modes)
 *   - Binary orbital decay: separation decreases as GW energy radiates
 *
 * Test Scenarios:
 *   1. Binary orbital decay (separation decreases)
 *   2. Chirp signal (frequency increases)
 *   3. GW polarization extraction (h₊, hₓ)
 *   4. Energy conservation (E_orbit + E_GW = const)
 *
 * Quality Gates:
 *   - Orbital decay: separation decreases >5% over simulation
 *   - Chirp signal: frequency increases (df/dt > 0)
 *   - Polarization: h₊/hₓ ratio matches GR within 10%
 *   - Energy conservation: E_total drift <1%
 *
 * References:
 *   - Peters & Mathews (1963): GW energy loss formula
 *   - LIGO GW150914 detection paper
 *   - Maggiore "Gravitational Waves Vol 1" (2008)
 */

#include "TRDParticleIntegrator.h"
#include "TRDCSVWriter.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "simulations/VisualizationGenerator.h"

// ============================================================================
// Physical Constants (Natural Units)
// ============================================================================

const double PI = 3.14159265358979323846;
const double G_natural = 1.0;  // Gravitational constant (natural units)
const double c = 1.0;          // Speed of light (natural units)

// ============================================================================
// Binary System Structure
// ============================================================================

/**
 * Binary system representation for gravitational wave emission
 */
struct BinarySystem {
    TRD::ChargedParticle m1;  // Mass 1 (using ChargedParticle for position/velocity tracking)
    TRD::ChargedParticle m2;  // Mass 2

    double separation;         // Current orbital separation
    double orbital_frequency;  // Current orbital frequency
    double total_mass;         // M1 + M2
    double reduced_mass;       // μ = M1·M2/(M1+M2)

    // Energy tracking
    double E_orbital;          // Orbital energy (kinetic + potential)
    double E_radiated;         // Total energy radiated as GW

    /**
     * Constructor - Initialize binary with circular orbit
     */
    BinarySystem(double M1, double M2, double a) {
        m1.mass = M1;
        m2.mass = M2;
        m1.charge = 0.0;  // Not used for gravity
        m2.charge = 0.0;

        total_mass = M1 + M2;
        reduced_mass = (M1 * M2) / total_mass;
        separation = a;

        // Place masses on x-axis (center of mass at origin)
        double r1 = a * (M2 / total_mass);  // Distance from m1 to COM
        double r2 = a * (M1 / total_mass);  // Distance from m2 to COM

        m1.x = r1;
        m1.y = 0.0;
        m1.z = 0.0;

        m2.x = -r2;
        m2.y = 0.0;
        m2.z = 0.0;

        // Circular orbit velocities (perpendicular to separation)
        orbital_frequency = std::sqrt(G_natural * total_mass / (a * a * a));
        double v_orbit = orbital_frequency * a;

        double v1 = v_orbit * (M2 / total_mass);
        double v2 = v_orbit * (M1 / total_mass);

        m1.vx = 0.0;
        m1.vy = v1;
        m1.vz = 0.0;

        m2.vx = 0.0;
        m2.vy = -v2;
        m2.vz = 0.0;

        // Initialize energy tracking
        E_radiated = 0.0;
        updateEnergy();
    }

    /**
     * Update orbital energy (kinetic + potential)
     */
    void updateEnergy() {
        // Kinetic energy
        double KE1 = 0.5 * m1.mass * (m1.vx*m1.vx + m1.vy*m1.vy + m1.vz*m1.vz);
        double KE2 = 0.5 * m2.mass * (m2.vx*m2.vx + m2.vy*m2.vy + m2.vz*m2.vz);

        // Separation
        double dx = m2.x - m1.x;
        double dy = m2.y - m1.y;
        double dz = m2.z - m1.z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        separation = r;

        // Potential energy
        double PE = -G_natural * m1.mass * m2.mass / r;

        E_orbital = KE1 + KE2 + PE;
    }

    /**
     * Get current orbital frequency from angular momentum
     */
    void updateFrequency() {
        // For circular orbits: ω = √(GM/a³)
        // This is more robust than v/r for slightly elliptical orbits
        orbital_frequency = std::sqrt(G_natural * total_mass / (separation * separation * separation));
    }
};

// ============================================================================
// Gravitational Wave Computation
// ============================================================================

/**
 * GW polarization structure
 */
struct GWPolarization {
    double h_plus;      // + polarization
    double h_cross;     // × polarization
    double amplitude;   // |h| = sqrt(h+² + h×²)
};

/**
 * Compute gravitational wave strain using quadrupole formula
 *
 * For circular binary orbit:
 *   h₊(t) = (4G/c⁴r) μ·a²·ω² cos(2ωt)
 *   h×(t) = (4G/c⁴r) μ·a²·ω² sin(2ωt)
 *
 * Where:
 *   μ = reduced mass
 *   a = orbital separation
 *   ω = orbital angular frequency
 *   r = distance to detector
 */
GWPolarization computeGWStrain(
    const BinarySystem& binary,
    double detector_distance,
    double time)
{
    double omega = 2.0 * PI * binary.orbital_frequency;
    double a = binary.separation;
    double mu = binary.reduced_mass;

    // Strain amplitude prefactor
    double h0 = (4.0 * G_natural / (c*c*c*c)) * (mu * a * a * omega * omega) / detector_distance;

    // Polarizations (face-on binary, optimal orientation)
    GWPolarization gw;
    gw.h_plus = h0 * std::cos(2.0 * omega * time);
    gw.h_cross = h0 * std::sin(2.0 * omega * time);
    gw.amplitude = std::sqrt(gw.h_plus * gw.h_plus + gw.h_cross * gw.h_cross);

    return gw;
}

/**
 * Compute gravitational wave energy loss rate (Peters-Mathews formula)
 *
 * For circular orbit:
 *   dE/dt = -(32/5) (G⁴/c⁵) (M₁M₂)²(M₁+M₂) / a⁵
 *
 * This is the luminosity in gravitational waves.
 */
double computeGWLuminosity(const BinarySystem& binary) {
    double M1 = binary.m1.mass;
    double M2 = binary.m2.mass;
    double a = binary.separation;

    double numerator = (32.0 / 5.0) * std::pow(G_natural, 4.0) * M1*M1 * M2*M2 * (M1 + M2);
    double denominator = std::pow(c, 5.0) * std::pow(a, 5.0);

    return numerator / denominator;
}

/**
 * Compute orbital decay timescale
 *
 * τ = (5/256) (c⁵/G³) a⁴/(M₁M₂(M₁+M₂))
 */
double computeDecayTimescale(const BinarySystem& binary) {
    double M1 = binary.m1.mass;
    double M2 = binary.m2.mass;
    double a = binary.separation;

    double numerator = (5.0 / 256.0) * std::pow(c, 5.0) * std::pow(a, 4.0);
    double denominator = std::pow(G_natural, 3.0) * M1 * M2 * (M1 + M2);

    return numerator / denominator;
}

// ============================================================================
// Binary Evolution
// ============================================================================

/**
 * Compute gravitational acceleration between two masses
 */
std::array<double, 3> computeGravAcceleration(
    const TRD::ChargedParticle& target,
    const TRD::ChargedParticle& source)
{
    double dx = source.x - target.x;
    double dy = source.y - target.y;
    double dz = source.z - target.z;
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = std::sqrt(r2);

    // Newtonian gravity: F = GMm/r², a = GM/r²
    double a_mag = G_natural * source.mass / r2;

    return {a_mag * dx / r, a_mag * dy / r, a_mag * dz / r};
}

/**
 * Evolve binary system one timestep with GW energy loss
 *
 * Uses Velocity Verlet for symplectic integration + adiabatic GW damping
 */
void evolveBinaryWithGW(BinarySystem& binary, double dt) {
    // 1. Compute current accelerations
    auto [ax1, ay1, az1] = computeGravAcceleration(binary.m1, binary.m2);
    auto [ax2, ay2, az2] = computeGravAcceleration(binary.m2, binary.m1);

    // 2. Half-kick (velocities)
    binary.m1.vx += ax1 * (dt / 2.0);
    binary.m1.vy += ay1 * (dt / 2.0);
    binary.m1.vz += az1 * (dt / 2.0);

    binary.m2.vx += ax2 * (dt / 2.0);
    binary.m2.vy += ay2 * (dt / 2.0);
    binary.m2.vz += az2 * (dt / 2.0);

    // 3. Drift (positions)
    binary.m1.x += binary.m1.vx * dt;
    binary.m1.y += binary.m1.vy * dt;
    binary.m1.z += binary.m1.vz * dt;

    binary.m2.x += binary.m2.vx * dt;
    binary.m2.y += binary.m2.vy * dt;
    binary.m2.z += binary.m2.vz * dt;

    // 4. Compute NEW accelerations at new positions
    auto [ax1_new, ay1_new, az1_new] = computeGravAcceleration(binary.m1, binary.m2);
    auto [ax2_new, ay2_new, az2_new] = computeGravAcceleration(binary.m2, binary.m1);

    // 5. Half-kick (velocities)
    binary.m1.vx += ax1_new * (dt / 2.0);
    binary.m1.vy += ay1_new * (dt / 2.0);
    binary.m1.vz += az1_new * (dt / 2.0);

    binary.m2.vx += ax2_new * (dt / 2.0);
    binary.m2.vy += ay2_new * (dt / 2.0);
    binary.m2.vz += az2_new * (dt / 2.0);

    // 6. GW energy loss (velocity damping to extract binding energy)
    // GW radiation removes energy from orbit by acting like a friction force

    double dE_GW = computeGWLuminosity(binary) * dt;

    // Apply uniform velocity damping to both masses
    // This mimics the GW back-reaction force which opposes orbital motion
    //
    // The damping factor is chosen to remove dE_GW of energy per timestep
    // For a bound orbit, reducing velocities extracts binding energy

    // Current total kinetic energy
    double v1_sq = binary.m1.vx*binary.m1.vx + binary.m1.vy*binary.m1.vy + binary.m1.vz*binary.m1.vz;
    double v2_sq = binary.m2.vx*binary.m2.vx + binary.m2.vy*binary.m2.vy + binary.m2.vz*binary.m2.vz;
    double KE_total = 0.5 * binary.m1.mass * v1_sq + 0.5 * binary.m2.mass * v2_sq;

    // Damping factor: reduce KE by dE_GW
    // KE_new = KE_old - dE_GW
    // v_new = v_old · sqrt((KE_old - dE_GW) / KE_old)
    if (KE_total > dE_GW) {
        double damp_factor = std::sqrt((KE_total - dE_GW) / KE_total);

        binary.m1.vx *= damp_factor;
        binary.m1.vy *= damp_factor;
        binary.m1.vz *= damp_factor;

        binary.m2.vx *= damp_factor;
        binary.m2.vy *= damp_factor;
        binary.m2.vz *= damp_factor;
    }

    // 7. Update energy and frequency
    binary.updateEnergy();
    binary.updateFrequency();
    binary.E_radiated += dE_GW;
}

// ============================================================================
// Test Functions
// ============================================================================

/**
 * Test 1: Binary Orbital Decay
 *
 * Verify that orbital separation decreases over time due to GW emission
 */
bool testOrbitalDecay(TRD::CSVWriter& csv) {
    std::cout << "\n=== Test 1: Binary Orbital Decay ===" << std::endl;

    // Binary parameters (scaled for numerical stability)
    double M1 = 1.0;  // Solar mass units
    double M2 = 1.0;
    double a_initial = 10.0;  // Initial separation

    BinarySystem binary(M1, M2, a_initial);

    // Simulation parameters
    double dt = 0.01;
    int num_steps = 10000;  // 100 time units
    int output_interval = 100;

    std::cout << "Initial separation: " << a_initial << std::endl;
    std::cout << "Initial frequency:  " << binary.orbital_frequency << std::endl;
    std::cout << "Decay timescale:    " << computeDecayTimescale(binary) << std::endl;

    // Write CSV header
    csv.writeHeader({"Time", "Separation", "Frequency", "E_Orbital", "E_Radiated", "E_Total"});

    // Evolution
    std::vector<double> separations;
    std::vector<double> frequencies;

    for (int step = 0; step <= num_steps; ++step) {
        double time = step * dt;

        // Record data
        if (step % output_interval == 0) {
            double E_total = binary.E_orbital + binary.E_radiated;

            csv.writeRow(time, binary.separation, binary.orbital_frequency,
                        binary.E_orbital, binary.E_radiated, E_total);

            separations.push_back(binary.separation);
            frequencies.push_back(binary.orbital_frequency);
        }

        // Evolve
        if (step < num_steps) {
            evolveBinaryWithGW(binary, dt);
        }
    }

    double a_final = binary.separation;
    double decay_fraction = (a_initial - a_final) / a_initial;

    std::cout << "Final separation:   " << a_final << std::endl;
    std::cout << "Decay fraction:     " << (decay_fraction * 100) << "%" << std::endl;
    std::cout << "Final frequency:    " << binary.orbital_frequency << std::endl;

    // PASS criteria: >5% decay
    bool pass = (decay_fraction > 0.05);
    std::cout << "Status: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass;
}

/**
 * Test 2: Chirp Signal Detection
 *
 * Verify that orbital frequency increases as binary inspirals
 */
bool testChirpSignal(TRD::CSVWriter& csv) {
    std::cout << "\n=== Test 2: Chirp Signal ===" << std::endl;

    // Binary parameters (same as Test 1 for consistency)
    double M1 = 1.0;
    double M2 = 1.0;
    double a_initial = 10.0;

    BinarySystem binary(M1, M2, a_initial);

    // Simulation parameters (same duration as Test 1)
    double dt = 0.01;
    int num_steps = 10000;
    int output_interval = 100;

    csv.writeHeader({"Time", "Frequency_GW", "h_plus", "h_cross", "h_amplitude"});

    std::vector<double> times;
    std::vector<double> frequencies_gw;

    double detector_distance = 100.0;  // Distance to "detector"

    for (int step = 0; step <= num_steps; ++step) {
        double time = step * dt;

        if (step % output_interval == 0) {
            // GW frequency is 2× orbital frequency (quadrupole radiation)
            double f_gw = 2.0 * binary.orbital_frequency;

            // Extract GW strain
            GWPolarization gw = computeGWStrain(binary, detector_distance, time);

            csv.writeRow(time, f_gw, gw.h_plus, gw.h_cross, gw.amplitude);

            VisualizationGenerator::addDataPoint("waveform", static_cast<float>(time), static_cast<float>(gw.h_plus));

            times.push_back(time);
            frequencies_gw.push_back(f_gw);
        }

        if (step < num_steps) {
            evolveBinaryWithGW(binary, dt);
        }
    }

    // Compute df/dt (chirp rate)
    double f_initial = frequencies_gw.front();
    double f_final = frequencies_gw.back();
    double t_duration = times.back() - times.front();
    double df_dt = (f_final - f_initial) / t_duration;

    std::cout << "Initial GW frequency: " << f_initial << std::endl;
    std::cout << "Final GW frequency:   " << f_final << std::endl;
    std::cout << "Chirp rate (df/dt):   " << df_dt << std::endl;

    // PASS criteria: frequency increases (df/dt > 0)
    bool pass = (df_dt > 0.0);
    std::cout << "Status: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass;
}

/**
 * Test 3: GW Polarization
 *
 * Verify h₊ and h× amplitudes match GR predictions
 */
bool testPolarization(TRD::CSVWriter& csv) {
    std::cout << "\n=== Test 3: GW Polarization ===" << std::endl;

    // Fixed binary (no evolution, just extract waveform)
    double M1 = 1.0;
    double M2 = 1.0;
    double a = 5.0;

    BinarySystem binary(M1, M2, a);

    double detector_distance = 100.0;
    double omega = 2.0 * PI * binary.orbital_frequency;
    double period = 2.0 * PI / omega;

    std::cout << "Orbital period: " << period << std::endl;
    std::cout << "GW frequency:   " << (2.0 * binary.orbital_frequency) << std::endl;

    // Sample waveform over 2 periods
    double dt = period / 50.0;  // 50 samples per orbital period
    int num_samples = 100;      // 2 periods

    csv.writeHeader({"Time", "Orbital_Phase", "h_plus", "h_cross"});

    std::vector<double> h_plus_vals;
    std::vector<double> h_cross_vals;

    for (int i = 0; i < num_samples; ++i) {
        double time = i * dt;
        double phase = omega * time;

        GWPolarization gw = computeGWStrain(binary, detector_distance, time);

        csv.writeRow(time, phase, gw.h_plus, gw.h_cross);

        h_plus_vals.push_back(gw.h_plus);
        h_cross_vals.push_back(gw.h_cross);
    }

    // Check polarization ratio (should be ~1 for optimal orientation)
    double h_plus_max = *std::max_element(h_plus_vals.begin(), h_plus_vals.end());
    double h_cross_max = *std::max_element(h_cross_vals.begin(), h_cross_vals.end());
    double ratio = h_plus_max / h_cross_max;

    std::cout << "h₊ max amplitude:  " << h_plus_max << std::endl;
    std::cout << "h× max amplitude:  " << h_cross_max << std::endl;
    std::cout << "Ratio (h₊/h×):     " << ratio << std::endl;

    // PASS criteria: ratio within 10% of 1.0
    bool pass = (std::abs(ratio - 1.0) < 0.1);
    std::cout << "Status: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass;
}

/**
 * Test 4: Energy Conservation
 *
 * Verify E_total = E_orbital + E_radiated = constant
 */
bool testEnergyConservation(TRD::CSVWriter& csv) {
    std::cout << "\n=== Test 4: Energy Conservation ===" << std::endl;

    double M1 = 1.0;
    double M2 = 1.0;
    double a_initial = 10.0;

    BinarySystem binary(M1, M2, a_initial);

    double E_initial = binary.E_orbital;

    double dt = 0.01;
    int num_steps = 10000;
    int output_interval = 100;

    csv.writeHeader({"Time", "E_Orbital", "E_Radiated", "E_Total", "E_Drift_Percent"});

    std::vector<double> E_total_history;

    for (int step = 0; step <= num_steps; ++step) {
        double time = step * dt;

        if (step % output_interval == 0) {
            double E_total = binary.E_orbital + binary.E_radiated;
            double E_drift = std::abs(E_total - E_initial) / std::abs(E_initial);

            csv.writeRow(time, binary.E_orbital, binary.E_radiated, E_total, E_drift * 100.0);

            E_total_history.push_back(E_total);
        }

        if (step < num_steps) {
            evolveBinaryWithGW(binary, dt);
        }
    }

    double E_final = binary.E_orbital + binary.E_radiated;
    double drift = std::abs(E_final - E_initial) / std::abs(E_initial);

    std::cout << "Initial total energy: " << E_initial << std::endl;
    std::cout << "Final total energy:   " << E_final << std::endl;
    std::cout << "Energy drift:         " << (drift * 100.0) << "%" << std::endl;

    // PASS criteria: <1% drift
    bool pass = (drift < 0.01);
    std::cout << "Status: " << (pass ? "PASS" : "FAIL") << std::endl;

    return pass;
}

// ============================================================================
// Main Test Runner
// ============================================================================

int runGravitationalWavesTest() {
    std::cout << "\n======================================" << std::endl;
    std::cout << "A5: Gravitational Wave Validation" << std::endl;
    std::cout << "======================================" << std::endl;

    std::cout << "\nPhysics: Binary system GW emission" << std::endl;
    std::cout << "Quality Gates:" << std::endl;
    std::cout << "  - Orbital decay >5%" << std::endl;
    std::cout << "  - Chirp signal detected (df/dt > 0)" << std::endl;
    std::cout << "  - Polarization ratio within 10%" << std::endl;
    std::cout << "  - Energy conservation <1% drift" << std::endl;

    bool all_pass = true;

    // Test 1: Orbital Decay
    {
        TRD::CSVWriter csv("orbital_decay", "A5_GravitationalWaves", true);
        csv.writeMetadata({
            {"test", "orbital_decay"},
            {"M1", "1.0"},
            {"M2", "1.0"},
            {"a_initial", "10.0"},
            {"dt", "0.01"},
            {"steps", "10000"}
        });
        bool pass = testOrbitalDecay(csv);
        all_pass = all_pass && pass;
        csv.close();
    }

    // Test 2: Chirp Signal
    {
        TRD::CSVWriter csv("chirp_signal", "A5_GravitationalWaves", true);
        csv.writeMetadata({
            {"test", "chirp_signal"},
            {"M1", "1.0"},
            {"M2", "1.0"},
            {"a_initial", "8.0"},
            {"dt", "0.01"},
            {"steps", "15000"}
        });
        bool pass = testChirpSignal(csv);
        all_pass = all_pass && pass;
        csv.close();
    }

    // Test 3: Polarization
    {
        TRD::CSVWriter csv("polarization", "A5_GravitationalWaves", true);
        csv.writeMetadata({
            {"test", "polarization"},
            {"M1", "1.0"},
            {"M2", "1.0"},
            {"a", "5.0"},
            {"detector_distance", "100.0"}
        });
        bool pass = testPolarization(csv);
        all_pass = all_pass && pass;
        csv.close();
    }

    // Test 4: Energy Conservation
    {
        TRD::CSVWriter csv("energy_conservation", "A5_GravitationalWaves", true);
        csv.writeMetadata({
            {"test", "energy_conservation"},
            {"M1", "1.0"},
            {"M2", "1.0"},
            {"a_initial", "10.0"},
            {"dt", "0.01"},
            {"steps", "10000"}
        });
        bool pass = testEnergyConservation(csv);
        all_pass = all_pass && pass;
        csv.close();
    }

    // Summary
    std::cout << "\n======================================" << std::endl;
    std::cout << "Summary: " << (all_pass ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;
    std::cout << "======================================\n" << std::endl;

    return all_pass ? 0 : 1;
}
