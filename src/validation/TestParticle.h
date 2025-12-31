/**
 * TestParticle.h
 *
 * Test particle dynamics for electromagnetic force validation.
 *
 * Purpose: Validate that EM fields extracted from Kuramoto phase
 * exert physically correct Lorentz forces on charged particles.
 *
 * Physics:
 *   - Lorentz force: F = q(E + v×B)
 *   - In 2D with B = B_z ẑ:
 *     F_x = q(E_x + v_y B_z)
 *     F_y = q(E_y - v_x B_z)
 *
 * Expected behaviors:
 *   - Cyclotron motion in B field with ω_c = qB/m
 *   - Larmor radius: r_L = mv⊥/(qB)
 *   - Electric field acceleration along field lines
 */

#pragma once

#include "../physics/EMFieldComputer.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

class TestParticle {
public:
    /**
     * Particle state in 2D
     */
    struct State {
        double x, y;      // Position [Planck lengths]
        double vx, vy;    // Velocity [c = 1]
        double t;         // Time [Planck time]

        // Constructor
        State(double x0 = 0, double y0 = 0,
              double vx0 = 0, double vy0 = 0,
              double t0 = 0)
            : x(x0), y(y0), vx(vx0), vy(vy0), t(t0) {}
    };

    /**
     * Trajectory point for history recording
     */
    struct TrajectoryPoint {
        double t;           // Time
        double x, y;        // Position
        double vx, vy;      // Velocity
        double E_x, E_y;    // Electric field at position
        double B_z;         // Magnetic field at position
        double F_x, F_y;    // Lorentz force at position
        double KE;          // Kinetic energy

        TrajectoryPoint() = default;
    };

private:
    // Physical parameters
    double charge_;     // Charge q [natural units]
    double mass_;       // Mass m [natural units]

    // Current state
    State state_;

    // Grid parameters for field interpolation
    int Nx_, Ny_;
    double dx_, dy_;
    double L_x_, L_y_;  // Domain size

    // Trajectory history
    std::vector<TrajectoryPoint> trajectory_;

    // Energy conservation tracking
    double initial_speed_;  // Speed at initialization (for energy conservation check)

    // Helper: Bilinear interpolation of field at continuous position
    double interpolateField(const Eigen::MatrixXd& field, double x, double y) const;

    // Helper: Apply periodic boundary conditions
    void applyPeriodicBC();

public:
    /**
     * Constructor
     * @param charge: Particle charge (default = 1.0)
     * @param mass: Particle mass (default = 1.0)
     * @param Nx, Ny: Grid dimensions
     * @param dx, dy: Grid spacing
     */
    TestParticle(double charge = 1.0, double mass = 1.0,
                 int Nx = 128, int Ny = 128,
                 double dx = 0.1, double dy = 0.1);

    /**
     * Initialize particle state
     * @param x0, y0: Initial position
     * @param vx0, vy0: Initial velocity
     * @param t0: Initial time
     */
    void initialize(double x0, double y0, double vx0, double vy0, double t0 = 0);

    /**
     * Evolve particle under Lorentz force using symplectic Velocity Verlet integration
     * Preserves phase space volume (Liouville's theorem) → energy conserved
     * @param fields: Electromagnetic fields
     * @param dt: Timestep
     * @param record: Whether to record trajectory point
     */
    void evolveLorentzForce(const EMFieldComputer::EMFields& fields,
                           double dt, bool record = true);

    /**
     * Compute Lorentz force at current position
     * @param fields: EM fields
     * @return Force vector (F_x, F_y)
     */
    Eigen::Vector2d computeLorentzForce(const EMFieldComputer::EMFields& fields) const;

    /**
     * Compute cyclotron frequency from current motion
     * Estimates ω_c from circular motion parameters
     * @return Cyclotron frequency [rad/Planck time]
     */
    double computeCyclotronFrequency() const;

    /**
     * Compute Larmor radius from trajectory
     * @return Radius of circular motion [Planck lengths]
     */
    double computeLarmorRadius() const;

    /**
     * Compute kinetic energy
     * @return KE = (1/2)m(v_x^2 + v_y^2)
     */
    double computeKineticEnergy() const;

    /**
     * Check if particle motion is stable (bounded orbit)
     * @param max_radius: Maximum allowed radius
     * @return true if orbit is bounded
     */
    bool isOrbitStable(double max_radius) const;

    /**
     * Record current state to trajectory history
     * @param fields: EM fields for field values at position
     */
    void recordTrajectory(const EMFieldComputer::EMFields& fields);

    /**
     * Write trajectory history to file
     * @param output_file: Output filename
     */
    void writeTrajectory(const std::string& output_file) const;

    /**
     * Clear trajectory history
     */
    void clearTrajectory() { trajectory_.clear(); }

    /**
     * Get current state
     */
    const State& getState() const { return state_; }

    /**
     * Get trajectory history
     */
    const std::vector<TrajectoryPoint>& getTrajectory() const { return trajectory_; }

    /**
     * Get physical parameters
     */
    double getCharge() const { return charge_; }
    double getMass() const { return mass_; }
    double getSpeed() const { return std::sqrt(state_.vx * state_.vx + state_.vy * state_.vy); }

    /**
     * Analyze orbital motion parameters
     * Computes period, radius, frequency from trajectory
     * @param n_periods: Number of periods to analyze
     * @return (period, radius, frequency)
     */
    std::tuple<double, double, double> analyzeOrbit(int n_periods = 1) const;

    /**
     * Validate cyclotron motion against theory
     * @param B_field: Magnetic field strength
     * @param tolerance: Relative error tolerance
     * @return true if motion matches theoretical predictions
     */
    bool validateCyclotronMotion(double B_field, double tolerance = 0.05) const;
};