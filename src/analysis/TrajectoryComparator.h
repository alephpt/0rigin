#ifndef TRAJECTORY_COMPARATOR_H
#define TRAJECTORY_COMPARATOR_H

#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>

/**
 * TrajectoryComparator - Compare observed vs geodesic particle trajectories
 *
 * Physics Framework:
 * ------------------
 * Hypothesis: Particle motion follows geodesics in (t,x,y,R) spacetime
 *
 * Observed trajectory:
 * - Position: r_obs(t) from quantum expectation values ⟨x⟩, ⟨y⟩
 * - Velocity: v_obs = dr_obs/dt
 * - Acceleration: a_obs = dv_obs/dt
 *
 * Geodesic prediction:
 * - Geodesic equation: d²x^i/dτ² + Γ^i_μν (dx^μ/dτ)(dx^ν/dτ) = 0
 * - Predicted acceleration: a_geo = -Γ^i_μν v^μ v^ν
 *
 * Validation Metric:
 * - Acceleration error: ε = |a_obs - a_geo| / |a_geo|
 * - Pass criterion: ⟨ε⟩ < 10% (average over trajectory)
 *
 * Statistical Analysis:
 * - Mean error: μ_ε = mean(ε(t))
 * - Standard deviation: σ_ε = std(ε(t))
 * - Maximum error: max(ε(t))
 * - Correlation coefficient: ρ(a_obs, a_geo)
 */
class TrajectoryComparator {
public:
    /**
     * Trajectory data structure
     * Contains time series of position, velocity, acceleration
     */
    struct TrajectoryData {
        std::vector<double> time;
        std::vector<Eigen::Vector2d> position;
        std::vector<Eigen::Vector2d> velocity;
        std::vector<Eigen::Vector2d> acceleration_observed;
        std::vector<Eigen::Vector2d> acceleration_geodesic;

        // Validation metrics
        std::vector<double> acceleration_error;       // |a_obs - a_geo| / |a_geo|
        std::vector<double> acceleration_magnitude_obs;
        std::vector<double> acceleration_magnitude_geo;

        // Statistical summary
        double mean_error;
        double std_error;
        double max_error;
        double correlation_x;  // ρ(a_obs_x, a_geo_x)
        double correlation_y;  // ρ(a_obs_y, a_geo_y)
        bool passes_validation;  // mean_error < 10%
    };

    /**
     * Build trajectory from observables CSV data
     *
     * Reads time series of position expectation values and computes
     * velocity and acceleration via finite differences.
     *
     * @param csv_path Path to observables.csv file
     * @return TrajectoryData with position, velocity, acceleration_observed
     */
    static TrajectoryData buildFromObservables(const std::string& csv_path);

    /**
     * Compute geodesic accelerations for entire trajectory
     *
     * For each (position, velocity) pair in trajectory, compute
     * predicted geodesic acceleration using GeometryAnalyzer.
     *
     * @param traj Trajectory data (position, velocity filled)
     * @param R_field_timeseries Time series of R-fields (vector of vectors)
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param L_domain Domain size (Planck units)
     */
    static void computeGeodesicAccelerations(
        TrajectoryData& traj,
        const std::vector<std::vector<double>>& R_field_timeseries,
        int Nx, int Ny,
        double L_domain);

    /**
     * Compute acceleration error metrics
     *
     * Computes:
     * - ε(t) = |a_obs - a_geo| / |a_geo| for each timestep
     * - Statistical summary: mean, std, max, correlation
     *
     * @param traj Trajectory with both observed and geodesic accelerations
     */
    static void computeErrorMetrics(TrajectoryData& traj);

    /**
     * Compute percent deviation between observed and geodesic acceleration
     *
     * Returns time series of relative errors:
     * ε(t) = |a_obs(t) - a_geo(t)| / |a_geo(t)| × 100%
     *
     * @param traj Trajectory data
     * @return Vector of percent deviations
     */
    static std::vector<double> computePercentDeviation(const TrajectoryData& traj);

    /**
     * Compute average acceleration error over trajectory
     *
     * Returns: ⟨ε⟩ = (1/N) Σ |a_obs - a_geo| / |a_geo|
     *
     * @param traj Trajectory data
     * @return Average percent error
     */
    static double computeAverageError(const TrajectoryData& traj);

    /**
     * Compute Pearson correlation coefficient between two time series
     *
     * ρ = cov(X,Y) / (σ_X σ_Y)
     *
     * @param series1 First time series
     * @param series2 Second time series
     * @return Correlation coefficient ρ ∈ [-1, 1]
     */
    static double computeCorrelation(
        const std::vector<double>& series1,
        const std::vector<double>& series2);

    /**
     * Write trajectory comparison to CSV
     *
     * Format: time,x,y,vx,vy,ax_obs,ay_obs,ax_geo,ay_geo,error
     *
     * @param traj Trajectory data
     * @param output_path Output CSV file path
     */
    static void writeToCSV(const TrajectoryData& traj, const std::string& output_path);

    /**
     * Generate validation report
     *
     * Creates human-readable summary of trajectory comparison:
     * - Mean/std/max acceleration error
     * - Correlation coefficients
     * - Pass/fail status
     *
     * @param traj Trajectory data
     * @return Report string
     */
    static std::string generateReport(const TrajectoryData& traj);

private:
    /**
     * Compute velocity from position time series
     * Uses centered finite differences: v(t) ≈ (r(t+dt) - r(t-dt)) / (2dt)
     *
     * @param position Time series of positions
     * @param time Time vector
     * @return Velocity time series
     */
    static std::vector<Eigen::Vector2d> computeVelocity(
        const std::vector<Eigen::Vector2d>& position,
        const std::vector<double>& time);

    /**
     * Compute acceleration from velocity time series
     * Uses centered finite differences: a(t) ≈ (v(t+dt) - v(t-dt)) / (2dt)
     *
     * @param velocity Time series of velocities
     * @param time Time vector
     * @return Acceleration time series
     */
    static std::vector<Eigen::Vector2d> computeAcceleration(
        const std::vector<Eigen::Vector2d>& velocity,
        const std::vector<double>& time);

    /**
     * Compute mean of time series
     */
    static double computeMean(const std::vector<double>& series);

    /**
     * Compute standard deviation of time series
     */
    static double computeStdDev(const std::vector<double>& series);

    /**
     * Compute maximum of time series
     */
    static double computeMax(const std::vector<double>& series);

    /**
     * Extract vector component from 2D vector time series
     * component: 0 for x, 1 for y
     */
    static std::vector<double> extractComponent(
        const std::vector<Eigen::Vector2d>& series,
        int component);

    // Validation threshold
    static constexpr double ERROR_THRESHOLD = 0.10;  // 10% mean error
};

#endif // TRAJECTORY_COMPARATOR_H
