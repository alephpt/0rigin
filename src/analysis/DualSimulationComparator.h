/**
 * DualSimulationComparator.h
 *
 * Phase 4 Test 4.1: Time Dilation Measurement
 *
 * Runs two parallel SMFT simulations to compare:
 * - Standard evolution: dτ = dt (normal time flow)
 * - Time-dilation evolution: dτ = R(x)·dt (R-modulated time)
 *
 * Physics Hypothesis:
 * R-field modulates local time flow rate. Particles in desynchronized
 * regions (R < 1) experience slower proper time, leading to phase
 * accumulation differences.
 *
 * Prediction: Δφ(t) = φ_standard(t) - φ_timedilation(t) grows linearly
 * in regions where R < 1.
 */

#pragma once

#include "SMFTEngine.h"
#include "DiracEvolution.h"
#include <memory>
#include <vector>
#include <string>

class DualSimulationComparator {
public:
    struct DualObservables {
        double time;

        // Standard evolution observables
        double norm_standard;
        double energy_standard;
        double phase_standard;
        double pos_x_standard;
        double pos_y_standard;

        // Time-dilation evolution observables
        double norm_timedilation;
        double energy_timedilation;
        double phase_timedilation;
        double pos_x_timedilation;
        double pos_y_timedilation;

        // Phase accumulation difference
        double delta_phase;        // φ_standard - φ_timedilation

        // R-field at particle center
        double R_center_standard;
        double R_center_timedilation;
    };

    /**
     * Constructor - initializes dual SMFT engines
     * @param Nx Grid width
     * @param Ny Grid height
     * @param Delta Vacuum potential parameter
     */
    DualSimulationComparator(uint32_t Nx, uint32_t Ny, float Delta);

    /**
     * Initialize both simulations with identical Kuramoto and Dirac states
     * @param theta_field Initial phase field
     * @param omega_field Natural frequency field
     * @param x0 Dirac wavepacket center x
     * @param y0 Dirac wavepacket center y
     * @param sigma Gaussian width
     */
    void initialize(const std::vector<float>& theta_field,
                   const std::vector<float>& omega_field,
                   float x0, float y0, float sigma);

    /**
     * Execute one dual timestep
     * - Engine A: Standard evolution (R-independent time)
     * - Engine B: Time-dilation evolution (dτ = R·dt)
     *
     * @param dt Time step size
     * @param K Kuramoto coupling
     * @param damping Phase damping
     * @param lambda_coupling Dirac-Kuramoto feedback
     */
    void step(float dt, float K, float damping, float lambda_coupling);

    /**
     * Compute all observables for both simulations
     * @param time Current simulation time
     * @return DualObservables struct with comparison data
     */
    DualObservables computeObservables(double time) const;

    /**
     * Write CSV header for dual observables
     */
    static std::string getCSVHeader();

    /**
     * Convert DualObservables to CSV line
     */
    static std::string toCSVLine(const DualObservables& obs);

private:
    uint32_t _Nx, _Ny;
    float _Delta;

    // Dual SMFT engines
    std::unique_ptr<SMFTEngine> _engine_standard;
    std::unique_ptr<SMFTEngine> _engine_timedilation;

    // Nova instance (required by SMFTEngine - CPU mode only)
    std::unique_ptr<class Nova> _nova;

    // Compute phase accumulation from wavepacket
    // φ = arg(⟨Ψ|Ψ⟩) = arg(∫Ψ*Ψ dx)
    double computePhaseAccumulation(const DiracEvolution* dirac) const;
};
