/**
 * TestParticle.cpp
 *
 * Implementation of test particle dynamics for EM force validation.
 */

#include "TestParticle.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <filesystem>

TestParticle::TestParticle(double charge, double mass,
                           int Nx, int Ny,
                           double dx, double dy)
    : charge_(charge), mass_(mass),
      Nx_(Nx), Ny_(Ny),
      dx_(dx), dy_(dy),
      L_x_(Nx * dx), L_y_(Ny * dy),
      state_(0, 0, 0, 0, 0) {
    trajectory_.reserve(10000);  // Pre-allocate for efficiency
}

void TestParticle::initialize(double x0, double y0,
                              double vx0, double vy0,
                              double t0) {
    state_.x = x0;
    state_.y = y0;
    state_.vx = vx0;
    state_.vy = vy0;
    state_.t = t0;

    // Store initial speed for energy conservation check
    initial_speed_ = std::sqrt(vx0 * vx0 + vy0 * vy0);

    // Apply periodic boundary conditions
    applyPeriodicBC();

    // Clear previous trajectory
    trajectory_.clear();
}

void TestParticle::applyPeriodicBC() {
    // Wrap position to [0, L_x) and [0, L_y)
    while (state_.x < 0) state_.x += L_x_;
    while (state_.x >= L_x_) state_.x -= L_x_;
    while (state_.y < 0) state_.y += L_y_;
    while (state_.y >= L_y_) state_.y -= L_y_;
}

double TestParticle::interpolateField(const Eigen::MatrixXd& field,
                                      double x, double y) const {
    // Convert continuous position to grid indices
    double ix_cont = x / dx_;
    double iy_cont = y / dy_;

    // Get integer indices (with periodic BC)
    int ix0 = static_cast<int>(std::floor(ix_cont)) % Nx_;
    int iy0 = static_cast<int>(std::floor(iy_cont)) % Ny_;
    int ix1 = (ix0 + 1) % Nx_;
    int iy1 = (iy0 + 1) % Ny_;

    // Handle negative indices
    if (ix0 < 0) ix0 += Nx_;
    if (iy0 < 0) iy0 += Ny_;

    // Compute interpolation weights
    double wx = ix_cont - std::floor(ix_cont);
    double wy = iy_cont - std::floor(iy_cont);

    // Bilinear interpolation
    double f00 = field(ix0, iy0);
    double f10 = field(ix1, iy0);
    double f01 = field(ix0, iy1);
    double f11 = field(ix1, iy1);

    return (1 - wx) * (1 - wy) * f00 +
           wx * (1 - wy) * f10 +
           (1 - wx) * wy * f01 +
           wx * wy * f11;
}

Eigen::Vector2d TestParticle::computeLorentzForce(
    const EMFieldComputer::EMFields& fields) const {

    // Interpolate fields at particle position
    double E_x = interpolateField(fields.E_x, state_.x, state_.y);
    double E_y = interpolateField(fields.E_y, state_.x, state_.y);
    double B_z = interpolateField(fields.B_z, state_.x, state_.y);

    // Lorentz force: F = q(E + v×B)
    // In 2D: v×B = (v_y B_z, -v_x B_z, 0)
    double F_x = charge_ * (E_x + state_.vy * B_z);
    double F_y = charge_ * (E_y - state_.vx * B_z);

    return Eigen::Vector2d(F_x, F_y);
}

void TestParticle::evolveLorentzForce(const EMFieldComputer::EMFields& fields,
                                      double dt, bool record) {
    // Record initial state if requested
    if (record && trajectory_.empty()) {
        recordTrajectory(fields);
    }

    // BORIS ALGORITHM - Standard for charged particle dynamics
    // Reference: Boris, J.P. (1970) "Relativistic plasma simulation-optimization of a hybrid code"
    // Properties:
    //   - Exact energy conservation for pure magnetic field (E=0)
    //   - Exact circular orbits (no breathing modes)
    //   - Symplectic structure preserved
    //   - Second-order accurate in dt

    // Interpolate fields at current position
    double E_x = interpolateField(fields.E_x, state_.x, state_.y);
    double E_y = interpolateField(fields.E_y, state_.x, state_.y);
    double B_z = interpolateField(fields.B_z, state_.x, state_.y);

    // Compute q/m and q*B/(2m) for Boris rotation
    double q_over_m = charge_ / mass_;
    double half_qB_over_m = 0.5 * q_over_m * B_z;

    // Step 1: Half electric impulse (v^- = v^n + (qE/m) * dt/2)
    double v_minus_x = state_.vx + 0.5 * q_over_m * E_x * dt;
    double v_minus_y = state_.vy + 0.5 * q_over_m * E_y * dt;

    // Step 2: Magnetic rotation
    // t = (qB/2m) * dt (rotation vector)
    double t_x = 0.0;
    double t_y = 0.0;
    double t_z = half_qB_over_m * dt;
    double t_mag_sq = t_z * t_z;

    // s = 2t / (1 + |t|²) (scaled rotation vector)
    double s_factor = 2.0 / (1.0 + t_mag_sq);
    double s_x = 0.0;
    double s_y = 0.0;
    double s_z = s_factor * t_z;

    // v' = v^- + v^- × t (first rotation)
    // In 2D: v×(0,0,t_z) = (v_y*t_z, -v_x*t_z, 0)
    double v_prime_x = v_minus_x + v_minus_y * t_z;
    double v_prime_y = v_minus_y - v_minus_x * t_z;

    // v^+ = v^- + v' × s (second rotation)
    double v_plus_x = v_minus_x + v_prime_y * s_z;
    double v_plus_y = v_minus_y - v_prime_x * s_z;

    // Step 3: Half electric impulse (v^{n+1} = v^+ + (qE/m) * dt/2)
    state_.vx = v_plus_x + 0.5 * q_over_m * E_x * dt;
    state_.vy = v_plus_y + 0.5 * q_over_m * E_y * dt;

    // Step 4: Position update (x^{n+1} = x^n + v^{n+1} * dt)
    state_.x += state_.vx * dt;
    state_.y += state_.vy * dt;
    state_.t += dt;

    // Step 5: Apply periodic boundary conditions
    applyPeriodicBC();

    // Step 6: Energy conservation check (CRITICAL for validation)
    // For pure magnetic field, speed should be conserved exactly
    double v_current = std::sqrt(state_.vx * state_.vx + state_.vy * state_.vy);
    if (initial_speed_ > 1e-10) {
        double speed_error = std::abs(v_current - initial_speed_) / initial_speed_;
        if (speed_error > 1e-6) {  // Tightened threshold for Boris algorithm
            std::cerr << "WARNING: Particle speed changed by "
                      << speed_error * 100.0 << "% "
                      << "(Boris algorithm should conserve speed exactly)\n"
                      << "  Initial speed: " << initial_speed_ << "\n"
                      << "  Current speed: " << v_current << "\n";
        }
    }

    // Step 7: Record new state if requested
    if (record) {
        recordTrajectory(fields);
    }
}

void TestParticle::recordTrajectory(const EMFieldComputer::EMFields& fields) {
    TrajectoryPoint point;

    // State variables
    point.t = state_.t;
    point.x = state_.x;
    point.y = state_.y;
    point.vx = state_.vx;
    point.vy = state_.vy;

    // Interpolate fields at position
    point.E_x = interpolateField(fields.E_x, state_.x, state_.y);
    point.E_y = interpolateField(fields.E_y, state_.x, state_.y);
    point.B_z = interpolateField(fields.B_z, state_.x, state_.y);

    // Compute force
    Eigen::Vector2d F = computeLorentzForce(fields);
    point.F_x = F(0);
    point.F_y = F(1);

    // Kinetic energy
    point.KE = computeKineticEnergy();

    trajectory_.push_back(point);
}

void TestParticle::writeTrajectory(const std::string& output_file) const {
    // Ensure parent directory exists
    std::filesystem::path file_path(output_file);
    std::filesystem::path dir_path = file_path.parent_path();

    if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
        std::filesystem::create_directories(dir_path);
    }

    // Check if trajectory has data
    if (trajectory_.empty()) {
        std::cerr << "Warning: Trajectory is empty (no points recorded)\n";
        return;
    }

    std::ofstream file(output_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << output_file << " for writing\n";
        std::cerr << "  Directory exists: " << std::filesystem::exists(dir_path) << "\n";
        std::cerr << "  Trajectory size: " << trajectory_.size() << "\n";
        return;
    }

    // Header
    file << "# Test Particle Trajectory\n";
    file << "# charge = " << charge_ << ", mass = " << mass_ << "\n";
    file << "# t,x,y,vx,vy,E_x,E_y,B_z,F_x,F_y,KE\n";

    // Data
    for (const auto& point : trajectory_) {
        file << std::fixed << std::setprecision(6)
             << point.t << ","
             << point.x << ","
             << point.y << ","
             << point.vx << ","
             << point.vy << ","
             << point.E_x << ","
             << point.E_y << ","
             << point.B_z << ","
             << point.F_x << ","
             << point.F_y << ","
             << point.KE << "\n";
    }

    file.close();
    std::cout << "  Trajectory written to " << output_file
              << " (" << trajectory_.size() << " points)\n";
}

double TestParticle::computeKineticEnergy() const {
    return 0.5 * mass_ * (state_.vx * state_.vx + state_.vy * state_.vy);
}

double TestParticle::computeLarmorRadius() const {
    if (trajectory_.size() < 10) {
        return 0.0;  // Not enough data
    }

    // Compute average radius from center of motion
    // First find approximate center by averaging positions
    double cx = 0, cy = 0;
    for (const auto& point : trajectory_) {
        cx += point.x;
        cy += point.y;
    }
    cx /= trajectory_.size();
    cy /= trajectory_.size();

    // Compute average radius
    double r_sum = 0;
    for (const auto& point : trajectory_) {
        double dx = point.x - cx;
        double dy = point.y - cy;
        r_sum += std::sqrt(dx*dx + dy*dy);
    }

    return r_sum / trajectory_.size();
}

double TestParticle::computeCyclotronFrequency() const {
    if (trajectory_.size() < 20) {
        return 0.0;  // Not enough data
    }

    // Estimate frequency from velocity direction changes
    // Count sign changes in vx to detect periods
    int crossings = 0;
    for (size_t i = 1; i < trajectory_.size(); ++i) {
        if (trajectory_[i-1].vx * trajectory_[i].vx < 0) {
            crossings++;
        }
    }

    if (crossings < 2) {
        return 0.0;  // Not enough crossings
    }

    // Period = 2 * time / crossings
    double total_time = trajectory_.back().t - trajectory_.front().t;
    double period = 2.0 * total_time / crossings;

    return 2.0 * M_PI / period;  // ω = 2π/T
}

bool TestParticle::isOrbitStable(double max_radius) const {
    if (trajectory_.empty()) {
        return true;  // No data yet
    }

    // Check if particle stays within max_radius
    for (const auto& point : trajectory_) {
        double r = std::sqrt(point.x * point.x + point.y * point.y);
        if (r > max_radius) {
            return false;
        }
    }

    return true;
}

std::tuple<double, double, double> TestParticle::analyzeOrbit(int n_periods) const {
    if (trajectory_.size() < 20) {
        return std::make_tuple(0.0, 0.0, 0.0);
    }

    // Find orbit center
    double cx = 0, cy = 0;
    for (const auto& point : trajectory_) {
        cx += point.x;
        cy += point.y;
    }
    cx /= trajectory_.size();
    cy /= trajectory_.size();

    // Compute average radius
    double r_avg = 0;
    for (const auto& point : trajectory_) {
        double dx = point.x - cx;
        double dy = point.y - cy;
        r_avg += std::sqrt(dx*dx + dy*dy);
    }
    r_avg /= trajectory_.size();

    // Compute orbital frequency from angular velocity
    // Method: Track angle θ = atan2(y-cy, x-cx) and measure dθ/dt
    std::vector<double> angles;
    for (const auto& point : trajectory_) {
        double dx = point.x - cx;
        double dy = point.y - cy;
        angles.push_back(std::atan2(dy, dx));
    }

    // Unwrap angles to handle 2π discontinuity
    std::vector<double> unwrapped_angles;
    unwrapped_angles.push_back(angles[0]);
    for (size_t i = 1; i < angles.size(); ++i) {
        double delta = angles[i] - angles[i-1];
        // Wrap delta to [-π, π]
        while (delta > M_PI) delta -= 2.0 * M_PI;
        while (delta < -M_PI) delta += 2.0 * M_PI;
        unwrapped_angles.push_back(unwrapped_angles[i-1] + delta);
    }

    // Total angular displacement
    double total_angle = unwrapped_angles.back() - unwrapped_angles.front();
    double total_time = trajectory_.back().t - trajectory_.front().t;

    // Average angular velocity (frequency in rad/s)
    double frequency = std::abs(total_angle / total_time);

    // Period
    double period = (frequency > 1e-10) ? 2.0 * M_PI / frequency : 0.0;

    return std::make_tuple(period, r_avg, frequency);
}

bool TestParticle::validateCyclotronMotion(double B_field, double tolerance) const {
    // Theoretical predictions
    double omega_c_theory = std::abs(charge_ * B_field / mass_);  // Cyclotron frequency
    double v_perp = std::sqrt(state_.vx * state_.vx + state_.vy * state_.vy);
    double r_L_theory = mass_ * v_perp / std::abs(charge_ * B_field);  // Larmor radius

    // Measured values
    auto [period, radius, frequency] = analyzeOrbit();

    if (frequency == 0 || radius == 0) {
        return false;  // Not enough data
    }

    // Compare with theory
    double freq_error = std::abs(frequency - omega_c_theory) / omega_c_theory;
    double radius_error = std::abs(radius - r_L_theory) / r_L_theory;

    bool freq_valid = freq_error < tolerance;
    bool radius_valid = radius_error < tolerance;

    // Energy conservation check
    if (trajectory_.size() > 1) {
        double KE_initial = trajectory_.front().KE;
        double KE_final = trajectory_.back().KE;
        double energy_error = std::abs(KE_final - KE_initial) / KE_initial;
        bool energy_valid = energy_error < tolerance;

        return freq_valid && radius_valid && energy_valid;
    }

    return freq_valid && radius_valid;
}