/**
 * EMValidator.h
 *
 * Electromagnetic field validation for SMFT framework.
 *
 * Purpose: Validate that EM fields extracted from Kuramoto phase
 * satisfy Maxwell equations and exhibit flux quantization.
 *
 * Key Validations:
 *   1. Maxwell equation residuals
 *   2. Flux quantization around vortices
 *   3. Energy-momentum conservation
 *   4. Gauge invariance checks
 */

#pragma once

#include "../physics/EMFieldComputer.h"
#include "../simulations/EMObservables.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include <complex>

class EMValidator {
public:
    /**
     * Maxwell equation validation results
     */
    struct MaxwellValidation {
        // Residuals for each equation
        double gauss_law_residual;       // |∇·E - 4πρ|
        double ampere_law_residual;      // |∇×B - 4πJ - ∂_tE|
        double faraday_law_residual;     // |∇×E + ∂_tB|
        double no_monopole_residual;     // |∇·B|

        // Relative errors (normalized by field magnitudes)
        double gauss_relative_error;
        double ampere_relative_error;
        double faraday_relative_error;

        // Spatial distribution metrics
        double max_violation;            // Maximum violation across grid
        double rms_violation;            // RMS of all violations
        int violation_count;             // Number of points above threshold

        // Violation maps
        Eigen::MatrixXd gauss_violation_map;
        Eigen::MatrixXd ampere_violation_map;
        Eigen::MatrixXd faraday_violation_map;
        Eigen::MatrixXd monopole_violation_map;

        // Constructor
        MaxwellValidation(int Nx, int Ny)
            : gauss_law_residual(0), ampere_law_residual(0),
              faraday_law_residual(0), no_monopole_residual(0),
              gauss_relative_error(0), ampere_relative_error(0),
              faraday_relative_error(0),
              max_violation(0), rms_violation(0), violation_count(0),
              gauss_violation_map(Nx, Ny),
              ampere_violation_map(Nx, Ny),
              faraday_violation_map(Nx, Ny),
              monopole_violation_map(Nx, Ny) {
            gauss_violation_map.setZero();
            ampere_violation_map.setZero();
            faraday_violation_map.setZero();
            monopole_violation_map.setZero();
        }
    };

    /**
     * Flux quantization results
     */
    struct FluxQuantization {
        double flux;                // Measured flux ∮A·dl
        double expected_flux;        // Expected flux (h/q)·W
        int winding_number;          // Extracted winding W
        double quantization_error;   // |flux - expected|/expected
        bool is_quantized;          // Whether flux is quantized

        // Radius dependence
        std::vector<double> radii;
        std::vector<double> fluxes;

        FluxQuantization()
            : flux(0), expected_flux(0), winding_number(0),
              quantization_error(0), is_quantized(false) {}
    };

private:
    // Grid parameters
    int Nx_, Ny_;
    double dx_, dy_, dt_;

    // Validation thresholds
    double maxwell_tolerance_;       // Tolerance for Maxwell equations
    double flux_tolerance_;          // Tolerance for flux quantization
    double violation_threshold_;     // Threshold for violation counting

    // Helper: Compute spatial derivatives with periodic BC
    Eigen::MatrixXd computeGradientX(const Eigen::MatrixXd& field) const;
    Eigen::MatrixXd computeGradientY(const Eigen::MatrixXd& field) const;
    Eigen::MatrixXd computeDivergence(const Eigen::MatrixXd& field_x,
                                      const Eigen::MatrixXd& field_y) const;
    double computeCurlZ(const Eigen::MatrixXd& field_x,
                        const Eigen::MatrixXd& field_y,
                        int ix, int iy) const;

    // Helper: Create contour for line integral
    std::vector<Eigen::Vector2d> createCircularContour(
        double center_x, double center_y,
        double radius, int n_points) const;

public:
    /**
     * Constructor
     * @param Nx, Ny: Grid dimensions
     * @param dx, dy: Grid spacing
     * @param dt: Timestep
     * @param maxwell_tol: Tolerance for Maxwell validation (default 0.01)
     * @param flux_tol: Tolerance for flux quantization (default 0.1)
     */
    EMValidator(int Nx, int Ny,
                double dx, double dy, double dt,
                double maxwell_tol = 0.01,
                double flux_tol = 0.1);

    /**
     * Verify Maxwell equations numerically
     * @param fields: Current EM fields
     * @param fields_prev: Previous EM fields (for time derivatives)
     * @param charge_density: Charge density ρ
     * @param current_x, current_y: Current density components
     * @return MaxwellValidation struct with residuals and maps
     */
    MaxwellValidation verifyMaxwellEquations(
        const EMFieldComputer::EMFields& fields,
        const EMFieldComputer::EMFields& fields_prev,
        const Eigen::MatrixXd& charge_density,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y);

    /**
     * Compute magnetic flux through closed loop
     * @param fields: EM fields
     * @param center_x, center_y: Center of loop [Planck lengths]
     * @param radius: Loop radius [Planck lengths]
     * @param n_points: Number of points on contour (default 100)
     * @return Flux ∮A·dl
     */
    double computeFlux(const EMFieldComputer::EMFields& fields,
                       double center_x, double center_y,
                       double radius, int n_points = 100);

    /**
     * Verify flux quantization around vortex
     * @param fields: EM fields
     * @param vortex_x, vortex_y: Vortex center [grid indices]
     * @param expected_winding: Expected winding number
     * @param charge: Effective charge (default 1.0)
     * @return FluxQuantization struct with results
     */
    FluxQuantization computeFluxQuantization(
        const EMFieldComputer::EMFields& fields,
        int vortex_x, int vortex_y,
        int expected_winding,
        double charge = 1.0);

    /**
     * Scan flux vs radius to verify quantization stability
     * @param fields: EM fields
     * @param center_x, center_y: Center point [grid indices]
     * @param r_min, r_max: Radius range [Planck lengths]
     * @param n_radii: Number of radii to test
     * @return Vector of (radius, flux) pairs
     */
    std::vector<std::pair<double, double>> scanFluxVsRadius(
        const EMFieldComputer::EMFields& fields,
        int center_x, int center_y,
        double r_min, double r_max,
        int n_radii = 20);

    /**
     * Compute energy-momentum conservation
     * @param fields: Current fields
     * @param fields_prev: Previous fields
     * @param current_x, current_y: Current density
     * @return Energy change rate dU/dt + ∇·S + J·E
     */
    double computeEnergyConservation(
        const EMFieldComputer::EMFields& fields,
        const EMFieldComputer::EMFields& fields_prev,
        const Eigen::MatrixXd& current_x,
        const Eigen::MatrixXd& current_y);

    /**
     * Check gauge invariance under θ → θ + const
     * @param fields1: Fields from phase θ
     * @param fields2: Fields from phase θ + const
     * @return Maximum difference in field strengths
     */
    double checkGaugeInvariance(
        const EMFieldComputer::EMFields& fields1,
        const EMFieldComputer::EMFields& fields2);

    /**
     * Write validation report to file
     * @param output_dir: Output directory
     * @param maxwell: Maxwell validation results
     * @param flux: Flux quantization results
     * @param test_name: Name of test for report
     */
    void writeValidationReport(
        const std::string& output_dir,
        const MaxwellValidation& maxwell,
        const FluxQuantization& flux,
        const std::string& test_name = "EM_Validation");

    /**
     * Write spatial violation maps to CSV files
     * @param output_dir: Output directory
     * @param maxwell: Maxwell validation with violation maps
     */
    void writeViolationMaps(
        const std::string& output_dir,
        const MaxwellValidation& maxwell);

    /**
     * Analyze grid convergence of Maxwell violations
     * @param violations_by_grid: Vector of (grid_size, violation) pairs
     * @return Convergence order (should be ~2 for O(h²))
     */
    static double analyzeGridConvergence(
        const std::vector<std::pair<int, double>>& violations_by_grid);

    /**
     * Extract vortex positions from B field
     * Finds peaks in |B_z| above threshold
     * @param B_z: Magnetic field z-component
     * @param threshold: Detection threshold
     * @return Vector of (x, y) vortex positions [grid indices]
     */
    static std::vector<std::pair<int, int>> findVortices(
        const Eigen::MatrixXd& B_z,
        double threshold);

    /**
     * Compute effective fine structure constant
     * @param measured_force: Measured Lorentz force
     * @param theoretical_force: Expected force with α = 1/137
     * @return α_effective
     */
    static double computeEffectiveAlpha(
        double measured_force,
        double theoretical_force);
};