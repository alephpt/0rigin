#include "TrajectoryComparator.h"
#include "GeometryAnalyzer.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iomanip>

// ============================================================================
// Trajectory Building from Observables
// ============================================================================

TrajectoryComparator::TrajectoryData TrajectoryComparator::buildFromObservables(
    const std::string& csv_path)
{
    TrajectoryData traj;

    std::ifstream file(csv_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open observables CSV: " + csv_path);
    }

    std::string line;
    std::getline(file, line);  // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        // Parse CSV line
        // Format: time,norm,norm_error,E_total,E_kin,E_pot,<x_re>,<x_im>,<y_re>,<y_im>,...
        std::vector<double> values;
        while (std::getline(ss, token, ',')) {
            values.push_back(std::stod(token));
        }

        if (values.size() < 10) continue;  // Need at least time + position data

        double t = values[0];
        double x_re = values[6];  // <x> real part
        double y_re = values[8];  // <y> real part

        traj.time.push_back(t);
        traj.position.push_back(Eigen::Vector2d(x_re, y_re));
    }

    file.close();

    // Compute velocity and observed acceleration from position
    if (traj.time.size() > 2) {
        traj.velocity = computeVelocity(traj.position, traj.time);
        traj.acceleration_observed = computeAcceleration(traj.velocity, traj.time);
    }

    return traj;
}

// ============================================================================
// Geodesic Acceleration Computation
// ============================================================================

void TrajectoryComparator::computeGeodesicAccelerations(
    TrajectoryData& traj,
    const std::vector<std::vector<double>>& R_field_timeseries,
    int Nx, int Ny,
    double L_domain)
{
    size_t n_steps = traj.position.size();
    traj.acceleration_geodesic.resize(n_steps);

    for (size_t i = 0; i < n_steps; ++i) {
        // Get R-field at this timestep
        if (i >= R_field_timeseries.size()) {
            // No R-field data available, use zero acceleration
            traj.acceleration_geodesic[i] = Eigen::Vector2d::Zero();
            continue;
        }

        const auto& R_field = R_field_timeseries[i];
        const auto& pos = traj.position[i];
        Eigen::Vector2d vel = (i < traj.velocity.size()) ?
            traj.velocity[i] : Eigen::Vector2d::Zero();

        // Compute geodesic acceleration at this position and velocity
        traj.acceleration_geodesic[i] = GeometryAnalyzer::computeGeodesicAcceleration(
            R_field, Nx, Ny, pos(0), pos(1), vel(0), vel(1), L_domain);
    }
}

// ============================================================================
// Error Metrics Computation
// ============================================================================

void TrajectoryComparator::computeErrorMetrics(TrajectoryData& traj) {
    size_t n_steps = traj.position.size();
    traj.acceleration_error.resize(n_steps);
    traj.acceleration_magnitude_obs.resize(n_steps);
    traj.acceleration_magnitude_geo.resize(n_steps);

    for (size_t i = 0; i < n_steps; ++i) {
        if (i >= traj.acceleration_observed.size() ||
            i >= traj.acceleration_geodesic.size()) {
            traj.acceleration_error[i] = 0.0;
            traj.acceleration_magnitude_obs[i] = 0.0;
            traj.acceleration_magnitude_geo[i] = 0.0;
            continue;
        }

        const auto& a_obs = traj.acceleration_observed[i];
        const auto& a_geo = traj.acceleration_geodesic[i];

        double mag_obs = a_obs.norm();
        double mag_geo = a_geo.norm();

        traj.acceleration_magnitude_obs[i] = mag_obs;
        traj.acceleration_magnitude_geo[i] = mag_geo;

        // Relative error: ε = |a_obs - a_geo| / |a_geo|
        if (mag_geo > 1e-12) {
            double error = (a_obs - a_geo).norm() / mag_geo;
            traj.acceleration_error[i] = error;
        } else {
            // Geodesic acceleration near zero: use absolute error
            traj.acceleration_error[i] = mag_obs;
        }
    }

    // Statistical summary
    if (!traj.acceleration_error.empty()) {
        traj.mean_error = computeMean(traj.acceleration_error);
        traj.std_error = computeStdDev(traj.acceleration_error);
        traj.max_error = computeMax(traj.acceleration_error);

        // Correlation coefficients (component-wise)
        std::vector<double> a_obs_x = extractComponent(traj.acceleration_observed, 0);
        std::vector<double> a_obs_y = extractComponent(traj.acceleration_observed, 1);
        std::vector<double> a_geo_x = extractComponent(traj.acceleration_geodesic, 0);
        std::vector<double> a_geo_y = extractComponent(traj.acceleration_geodesic, 1);

        if (!a_obs_x.empty() && !a_geo_x.empty()) {
            traj.correlation_x = computeCorrelation(a_obs_x, a_geo_x);
            traj.correlation_y = computeCorrelation(a_obs_y, a_geo_y);
        } else {
            traj.correlation_x = 0.0;
            traj.correlation_y = 0.0;
        }

        // Validation: mean error < 10%
        traj.passes_validation = (traj.mean_error < ERROR_THRESHOLD);
    } else {
        traj.mean_error = 0.0;
        traj.std_error = 0.0;
        traj.max_error = 0.0;
        traj.correlation_x = 0.0;
        traj.correlation_y = 0.0;
        traj.passes_validation = false;
    }
}

// ============================================================================
// Error Analysis Methods
// ============================================================================

std::vector<double> TrajectoryComparator::computePercentDeviation(
    const TrajectoryData& traj)
{
    std::vector<double> percent_deviation;
    percent_deviation.reserve(traj.acceleration_error.size());

    for (double error : traj.acceleration_error) {
        percent_deviation.push_back(error * 100.0);  // Convert to percentage
    }

    return percent_deviation;
}

double TrajectoryComparator::computeAverageError(const TrajectoryData& traj) {
    if (traj.acceleration_error.empty()) return 0.0;
    return traj.mean_error;
}

// ============================================================================
// Derivative Computation (Finite Differences)
// ============================================================================

std::vector<Eigen::Vector2d> TrajectoryComparator::computeVelocity(
    const std::vector<Eigen::Vector2d>& position,
    const std::vector<double>& time)
{
    size_t n = position.size();
    std::vector<Eigen::Vector2d> velocity(n);

    if (n < 2) return velocity;

    // Forward difference for first point
    double dt = time[1] - time[0];
    velocity[0] = (position[1] - position[0]) / dt;

    // Centered differences for interior points
    for (size_t i = 1; i < n - 1; ++i) {
        dt = time[i + 1] - time[i - 1];
        velocity[i] = (position[i + 1] - position[i - 1]) / dt;
    }

    // Backward difference for last point
    dt = time[n - 1] - time[n - 2];
    velocity[n - 1] = (position[n - 1] - position[n - 2]) / dt;

    return velocity;
}

std::vector<Eigen::Vector2d> TrajectoryComparator::computeAcceleration(
    const std::vector<Eigen::Vector2d>& velocity,
    const std::vector<double>& time)
{
    size_t n = velocity.size();
    std::vector<Eigen::Vector2d> acceleration(n);

    if (n < 2) return acceleration;

    // Forward difference for first point
    double dt = time[1] - time[0];
    acceleration[0] = (velocity[1] - velocity[0]) / dt;

    // Centered differences for interior points
    for (size_t i = 1; i < n - 1; ++i) {
        dt = time[i + 1] - time[i - 1];
        acceleration[i] = (velocity[i + 1] - velocity[i - 1]) / dt;
    }

    // Backward difference for last point
    dt = time[n - 1] - time[n - 2];
    acceleration[n - 1] = (velocity[n - 1] - velocity[n - 2]) / dt;

    return acceleration;
}

// ============================================================================
// Statistical Utilities
// ============================================================================

double TrajectoryComparator::computeMean(const std::vector<double>& series) {
    if (series.empty()) return 0.0;
    double sum = std::accumulate(series.begin(), series.end(), 0.0);
    return sum / static_cast<double>(series.size());
}

double TrajectoryComparator::computeStdDev(const std::vector<double>& series) {
    if (series.size() < 2) return 0.0;

    double mean = computeMean(series);
    double sq_sum = 0.0;
    for (double x : series) {
        double diff = x - mean;
        sq_sum += diff * diff;
    }

    return std::sqrt(sq_sum / static_cast<double>(series.size() - 1));
}

double TrajectoryComparator::computeMax(const std::vector<double>& series) {
    if (series.empty()) return 0.0;
    return *std::max_element(series.begin(), series.end());
}

double TrajectoryComparator::computeCorrelation(
    const std::vector<double>& series1,
    const std::vector<double>& series2)
{
    size_t n = std::min(series1.size(), series2.size());
    if (n < 2) return 0.0;

    double mean1 = 0.0, mean2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        mean1 += series1[i];
        mean2 += series2[i];
    }
    mean1 /= static_cast<double>(n);
    mean2 /= static_cast<double>(n);

    double cov = 0.0, var1 = 0.0, var2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff1 = series1[i] - mean1;
        double diff2 = series2[i] - mean2;
        cov += diff1 * diff2;
        var1 += diff1 * diff1;
        var2 += diff2 * diff2;
    }

    if (var1 < 1e-12 || var2 < 1e-12) return 0.0;

    return cov / std::sqrt(var1 * var2);
}

std::vector<double> TrajectoryComparator::extractComponent(
    const std::vector<Eigen::Vector2d>& series,
    int component)
{
    std::vector<double> result;
    result.reserve(series.size());

    for (const auto& vec : series) {
        result.push_back(vec(component));
    }

    return result;
}

// ============================================================================
// Output Methods
// ============================================================================

void TrajectoryComparator::writeToCSV(
    const TrajectoryData& traj,
    const std::string& output_path)
{
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output CSV: " + output_path);
    }

    // Write header
    file << "time,x,y,vx,vy,ax_obs,ay_obs,ax_geo,ay_geo,error_percent\n";

    // Write data
    size_t n = traj.time.size();
    for (size_t i = 0; i < n; ++i) {
        file << std::scientific << std::setprecision(10);
        file << traj.time[i] << ",";
        file << traj.position[i](0) << "," << traj.position[i](1) << ",";

        if (i < traj.velocity.size()) {
            file << traj.velocity[i](0) << "," << traj.velocity[i](1) << ",";
        } else {
            file << "0,0,";
        }

        if (i < traj.acceleration_observed.size()) {
            file << traj.acceleration_observed[i](0) << ","
                 << traj.acceleration_observed[i](1) << ",";
        } else {
            file << "0,0,";
        }

        if (i < traj.acceleration_geodesic.size()) {
            file << traj.acceleration_geodesic[i](0) << ","
                 << traj.acceleration_geodesic[i](1) << ",";
        } else {
            file << "0,0,";
        }

        if (i < traj.acceleration_error.size()) {
            file << traj.acceleration_error[i] * 100.0;  // Percent
        } else {
            file << "0";
        }

        file << "\n";
    }

    file.close();
}

std::string TrajectoryComparator::generateReport(const TrajectoryData& traj) {
    std::ostringstream report;

    report << "========================================\n";
    report << "Geodesic Deviation Analysis Report\n";
    report << "========================================\n\n";

    report << "Trajectory Statistics:\n";
    report << "  Number of timesteps: " << traj.time.size() << "\n";
    if (!traj.time.empty()) {
        report << "  Time range: [" << traj.time.front() << ", "
               << traj.time.back() << "]\n";
    }
    report << "\n";

    report << "Acceleration Error Metrics:\n";
    report << std::fixed << std::setprecision(4);
    report << "  Mean error:       " << (traj.mean_error * 100.0) << "%\n";
    report << "  Std deviation:    " << (traj.std_error * 100.0) << "%\n";
    report << "  Maximum error:    " << (traj.max_error * 100.0) << "%\n";
    report << "\n";

    report << "Correlation Coefficients:\n";
    report << "  ρ(a_obs_x, a_geo_x): " << traj.correlation_x << "\n";
    report << "  ρ(a_obs_y, a_geo_y): " << traj.correlation_y << "\n";
    report << "\n";

    report << "Validation Status:\n";
    report << "  Threshold: " << (ERROR_THRESHOLD * 100.0) << "%\n";
    report << "  Status: " << (traj.passes_validation ? "PASS" : "FAIL") << "\n";
    report << "\n";

    if (traj.passes_validation) {
        report << "✓ Geodesic hypothesis validated!\n";
        report << "  Observed particle motion agrees with geodesic prediction\n";
        report << "  within " << (ERROR_THRESHOLD * 100.0) << "% error tolerance.\n";
    } else {
        report << "✗ Geodesic hypothesis not validated.\n";
        report << "  Mean error (" << (traj.mean_error * 100.0) << "%) exceeds\n";
        report << "  threshold (" << (ERROR_THRESHOLD * 100.0) << "%).\n";
    }

    report << "========================================\n";

    return report.str();
}
