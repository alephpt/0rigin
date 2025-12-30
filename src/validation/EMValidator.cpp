/**
 * EMValidator.cpp
 *
 * Implementation of electromagnetic field validation.
 */

#include "EMValidator.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>

EMValidator::EMValidator(int Nx, int Ny,
                         double dx, double dy, double dt,
                         double maxwell_tol, double flux_tol)
    : Nx_(Nx), Ny_(Ny),
      dx_(dx), dy_(dy), dt_(dt),
      maxwell_tolerance_(maxwell_tol),
      flux_tolerance_(flux_tol),
      violation_threshold_(0.01) {}

Eigen::MatrixXd EMValidator::computeGradientX(const Eigen::MatrixXd& field) const {
    Eigen::MatrixXd grad(Nx_, Ny_);

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            int ip = (i + 1) % Nx_;
            int im = (i - 1 + Nx_) % Nx_;

            grad(i, j) = (field(ip, j) - field(im, j)) / (2.0 * dx_);
        }
    }

    return grad;
}

Eigen::MatrixXd EMValidator::computeGradientY(const Eigen::MatrixXd& field) const {
    Eigen::MatrixXd grad(Nx_, Ny_);

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            int jp = (j + 1) % Ny_;
            int jm = (j - 1 + Ny_) % Ny_;

            grad(i, j) = (field(i, jp) - field(i, jm)) / (2.0 * dy_);
        }
    }

    return grad;
}

Eigen::MatrixXd EMValidator::computeDivergence(const Eigen::MatrixXd& field_x,
                                               const Eigen::MatrixXd& field_y) const {
    Eigen::MatrixXd div(Nx_, Ny_);

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            int ip = (i + 1) % Nx_;
            int im = (i - 1 + Nx_) % Nx_;
            int jp = (j + 1) % Ny_;
            int jm = (j - 1 + Ny_) % Ny_;

            double dx_component = (field_x(ip, j) - field_x(im, j)) / (2.0 * dx_);
            double dy_component = (field_y(i, jp) - field_y(i, jm)) / (2.0 * dy_);

            div(i, j) = dx_component + dy_component;
        }
    }

    return div;
}

double EMValidator::computeCurlZ(const Eigen::MatrixXd& field_x,
                                 const Eigen::MatrixXd& field_y,
                                 int ix, int iy) const {
    int ixp = (ix + 1) % Nx_;
    int ixm = (ix - 1 + Nx_) % Nx_;
    int iyp = (iy + 1) % Ny_;
    int iym = (iy - 1 + Ny_) % Ny_;

    double dFy_dx = (field_y(ixp, iy) - field_y(ixm, iy)) / (2.0 * dx_);
    double dFx_dy = (field_x(ix, iyp) - field_x(ix, iym)) / (2.0 * dy_);

    return dFy_dx - dFx_dy;  // (∇×F)_z
}

std::vector<Eigen::Vector2d> EMValidator::createCircularContour(
    double center_x, double center_y,
    double radius, int n_points) const {

    std::vector<Eigen::Vector2d> contour;
    contour.reserve(n_points);

    for (int i = 0; i < n_points; ++i) {
        double theta = 2.0 * M_PI * i / n_points;
        double x = center_x + radius * std::cos(theta);
        double y = center_y + radius * std::sin(theta);
        contour.push_back(Eigen::Vector2d(x, y));
    }

    return contour;
}

EMValidator::MaxwellValidation EMValidator::verifyMaxwellEquations(
    const EMFieldComputer::EMFields& fields,
    const EMFieldComputer::EMFields& fields_prev,
    const Eigen::MatrixXd& charge_density,
    const Eigen::MatrixXd& current_x,
    const Eigen::MatrixXd& current_y) {

    MaxwellValidation result(Nx_, Ny_);

    // Compute divergence of E
    Eigen::MatrixXd div_E = computeDivergence(fields.E_x, fields.E_y);

    // Compute divergence of B (should be zero)
    Eigen::MatrixXd div_B(Nx_, Ny_);
    div_B.setZero();  // In 2D, B only has z-component, so ∇·B = 0 automatically

    // Compute curl of B and E
    Eigen::MatrixXd curl_B_z(Nx_, Ny_);
    Eigen::MatrixXd curl_E_z(Nx_, Ny_);

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            // For B field (only has z-component in 2D), curl has x,y components
            // But we need ∇×B in z-direction which involves B_x, B_y (which are 0)
            // So we compute numerical curl differently

            // Actually in 2D, we need to be careful:
            // B = (0, 0, B_z), so ∇×B = (∂B_z/∂y, -∂B_z/∂x, 0)
            // We need the in-plane components for Ampere's law

            // Curl of E: (∇×E)_z = ∂E_y/∂x - ∂E_x/∂y
            curl_E_z(i, j) = computeCurlZ(fields.E_x, fields.E_y, i, j);
        }
    }

    // Time derivatives
    Eigen::MatrixXd dE_x_dt = (fields.E_x - fields_prev.E_x) / dt_;
    Eigen::MatrixXd dE_y_dt = (fields.E_y - fields_prev.E_y) / dt_;
    Eigen::MatrixXd dB_z_dt = (fields.B_z - fields_prev.B_z) / dt_;

    // Validate each Maxwell equation
    double gauss_sum = 0, ampere_sum = 0, faraday_sum = 0, monopole_sum = 0;
    double field_norm = 0;
    int violations = 0;

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            // Gauss's law: ∇·E = 4πρ
            double gauss_violation = std::abs(div_E(i, j) - 4.0 * M_PI * charge_density(i, j));
            result.gauss_violation_map(i, j) = gauss_violation;
            gauss_sum += gauss_violation * gauss_violation;

            // Ampere's law (2D version): ∇×B needs careful treatment
            // In 2D with B_z only: ∂B_z/∂y = 4πJ_x + ∂E_x/∂t
            //                      -∂B_z/∂x = 4πJ_y + ∂E_y/∂t
            int ip = (i + 1) % Nx_;
            int im = (i - 1 + Nx_) % Nx_;
            int jp = (j + 1) % Ny_;
            int jm = (j - 1 + Ny_) % Ny_;

            double dBz_dx = (fields.B_z(ip, j) - fields.B_z(im, j)) / (2.0 * dx_);
            double dBz_dy = (fields.B_z(i, jp) - fields.B_z(i, jm)) / (2.0 * dy_);

            double ampere_x_violation = std::abs(dBz_dy - 4.0 * M_PI * current_x(i, j) - dE_x_dt(i, j));
            double ampere_y_violation = std::abs(-dBz_dx - 4.0 * M_PI * current_y(i, j) - dE_y_dt(i, j));
            double ampere_violation = std::sqrt(ampere_x_violation * ampere_x_violation +
                                               ampere_y_violation * ampere_y_violation);

            result.ampere_violation_map(i, j) = ampere_violation;
            ampere_sum += ampere_violation * ampere_violation;

            // Faraday's law: ∇×E = -∂B/∂t
            double faraday_violation = std::abs(curl_E_z(i, j) + dB_z_dt(i, j));
            result.faraday_violation_map(i, j) = faraday_violation;
            faraday_sum += faraday_violation * faraday_violation;

            // No monopoles: ∇·B = 0 (automatically satisfied in 2D)
            result.monopole_violation_map(i, j) = 0.0;

            // Count violations above threshold
            double total_violation = gauss_violation + ampere_violation + faraday_violation;
            if (total_violation > violation_threshold_) {
                violations++;
            }

            // Track maximum
            result.max_violation = std::max(result.max_violation, total_violation);

            // Field normalization
            field_norm += fields.E_x(i, j) * fields.E_x(i, j) +
                         fields.E_y(i, j) * fields.E_y(i, j) +
                         fields.B_z(i, j) * fields.B_z(i, j);
        }
    }

    // Compute RMS residuals
    int N_total = Nx_ * Ny_;
    result.gauss_law_residual = std::sqrt(gauss_sum / N_total);
    result.ampere_law_residual = std::sqrt(ampere_sum / N_total);
    result.faraday_law_residual = std::sqrt(faraday_sum / N_total);
    result.no_monopole_residual = 0.0;  // Always zero in 2D

    // Compute relative errors
    double field_rms = std::sqrt(field_norm / N_total);
    if (field_rms > 1e-12) {
        result.gauss_relative_error = result.gauss_law_residual / field_rms;
        result.ampere_relative_error = result.ampere_law_residual / field_rms;
        result.faraday_relative_error = result.faraday_law_residual / field_rms;
    }

    // Overall RMS violation
    result.rms_violation = std::sqrt((gauss_sum + ampere_sum + faraday_sum) / (3.0 * N_total));
    result.violation_count = violations;

    return result;
}

double EMValidator::computeFlux(const EMFieldComputer::EMFields& fields,
                                double center_x, double center_y,
                                double radius, int n_points) {

    // Create contour
    auto contour = createCircularContour(center_x, center_y, radius, n_points);

    // Line integral ∮A·dl
    double flux = 0.0;

    for (int i = 0; i < n_points; ++i) {
        int i_next = (i + 1) % n_points;

        // Current and next points
        double x1 = contour[i](0);
        double y1 = contour[i](1);
        double x2 = contour[i_next](0);
        double y2 = contour[i_next](1);

        // Midpoint for field evaluation
        double x_mid = 0.5 * (x1 + x2);
        double y_mid = 0.5 * (y1 + y2);

        // Convert to grid indices
        int ix = static_cast<int>(x_mid / dx_) % Nx_;
        int iy = static_cast<int>(y_mid / dy_) % Ny_;
        if (ix < 0) ix += Nx_;
        if (iy < 0) iy += Ny_;

        // Get vector potential at midpoint
        double A_x = fields.A_x(ix, iy);
        double A_y = fields.A_y(ix, iy);

        // Line element
        double dl_x = x2 - x1;
        double dl_y = y2 - y1;

        // Contribution to integral
        flux += A_x * dl_x + A_y * dl_y;
    }

    return flux;
}

EMValidator::FluxQuantization EMValidator::computeFluxQuantization(
    const EMFieldComputer::EMFields& fields,
    int vortex_x, int vortex_y,
    int expected_winding,
    double charge) {

    FluxQuantization result;

    // Convert grid indices to physical coordinates
    double center_x = vortex_x * dx_;
    double center_y = vortex_y * dy_;

    // Test multiple radii to ensure quantization
    std::vector<double> test_radii = {3.0 * dx_, 5.0 * dx_, 10.0 * dx_, 15.0 * dx_};

    for (double radius : test_radii) {
        double flux = computeFlux(fields, center_x, center_y, radius);
        result.radii.push_back(radius);
        result.fluxes.push_back(flux);
    }

    // Use largest radius for final measurement
    result.flux = result.fluxes.back();

    // Expected flux: Φ = (h/q)·W = 2π·W (in natural units with h=2π, q=1)
    result.expected_flux = 2.0 * M_PI * expected_winding / charge;

    // Extract winding number
    result.winding_number = static_cast<int>(std::round(result.flux * charge / (2.0 * M_PI)));

    // Quantization error
    if (std::abs(result.expected_flux) > 1e-12) {
        result.quantization_error = std::abs(result.flux - result.expected_flux) /
                                    std::abs(result.expected_flux);
    } else {
        result.quantization_error = std::abs(result.flux);
    }

    // Check if quantized within tolerance
    result.is_quantized = (result.quantization_error < flux_tolerance_) &&
                          (result.winding_number == expected_winding);

    return result;
}

std::vector<std::pair<double, double>> EMValidator::scanFluxVsRadius(
    const EMFieldComputer::EMFields& fields,
    int center_x, int center_y,
    double r_min, double r_max,
    int n_radii) {

    std::vector<std::pair<double, double>> flux_vs_radius;

    double center_x_phys = center_x * dx_;
    double center_y_phys = center_y * dy_;

    for (int i = 0; i < n_radii; ++i) {
        double radius = r_min + (r_max - r_min) * i / (n_radii - 1);
        double flux = computeFlux(fields, center_x_phys, center_y_phys, radius);
        flux_vs_radius.push_back({radius, flux});
    }

    return flux_vs_radius;
}

double EMValidator::computeEnergyConservation(
    const EMFieldComputer::EMFields& fields,
    const EMFieldComputer::EMFields& fields_prev,
    const Eigen::MatrixXd& current_x,
    const Eigen::MatrixXd& current_y) {

    // Compute field energy change
    double U_curr = EMFieldComputer::computeFieldEnergy(fields, dx_, dy_);
    double U_prev = EMFieldComputer::computeFieldEnergy(fields_prev, dx_, dy_);
    double dU_dt = (U_curr - U_prev) / dt_;

    // Compute Poynting vector divergence
    auto poynting = EMFieldComputer::computePoyntingVector(fields, dx_, dy_);
    Eigen::MatrixXd div_S = computeDivergence(poynting.S_x, poynting.S_y);

    // Compute J·E (work done by fields)
    double work = 0.0;
    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            work += current_x(i, j) * fields.E_x(i, j) +
                   current_y(i, j) * fields.E_y(i, j);
        }
    }
    work *= dx_ * dy_;

    // Energy conservation: dU/dt + ∇·S = -J·E
    double div_S_total = div_S.sum() * dx_ * dy_;
    double conservation_residual = dU_dt + div_S_total + work;

    return conservation_residual;
}

double EMValidator::checkGaugeInvariance(
    const EMFieldComputer::EMFields& fields1,
    const EMFieldComputer::EMFields& fields2) {

    double max_diff = 0.0;

    for (int i = 0; i < Nx_; ++i) {
        for (int j = 0; j < Ny_; ++j) {
            // Field strengths should be identical
            double E_diff = std::abs(fields1.E_x(i, j) - fields2.E_x(i, j)) +
                           std::abs(fields1.E_y(i, j) - fields2.E_y(i, j));
            double B_diff = std::abs(fields1.B_z(i, j) - fields2.B_z(i, j));

            max_diff = std::max(max_diff, E_diff + B_diff);
        }
    }

    return max_diff;
}

void EMValidator::writeValidationReport(
    const std::string& output_dir,
    const MaxwellValidation& maxwell,
    const FluxQuantization& flux,
    const std::string& test_name) {

    std::string filename = output_dir + "/em_validation_report.txt";
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing\n";
        return;
    }

    file << "=====================================\n";
    file << "   EM VALIDATION REPORT\n";
    file << "   Test: " << test_name << "\n";
    file << "=====================================\n\n";

    file << "Grid Parameters:\n";
    file << "  Nx × Ny: " << Nx_ << " × " << Ny_ << "\n";
    file << "  dx, dy: " << dx_ << ", " << dy_ << " [Planck lengths]\n";
    file << "  dt: " << dt_ << " [Planck time]\n\n";

    file << "Maxwell Equation Validation:\n";
    file << "  Gauss Law Residual: " << std::scientific << std::setprecision(3)
         << maxwell.gauss_law_residual << "\n";
    file << "  Ampere Law Residual: " << maxwell.ampere_law_residual << "\n";
    file << "  Faraday Law Residual: " << maxwell.faraday_law_residual << "\n";
    file << "  No Monopole Residual: " << maxwell.no_monopole_residual << "\n\n";

    file << "  Relative Errors:\n";
    file << "    Gauss: " << maxwell.gauss_relative_error << "\n";
    file << "    Ampere: " << maxwell.ampere_relative_error << "\n";
    file << "    Faraday: " << maxwell.faraday_relative_error << "\n\n";

    file << "  Violation Statistics:\n";
    file << "    Max Violation: " << maxwell.max_violation << "\n";
    file << "    RMS Violation: " << maxwell.rms_violation << "\n";
    file << "    Points Above Threshold: " << maxwell.violation_count
         << " / " << Nx_ * Ny_ << "\n\n";

    // Maxwell validation pass/fail
    bool maxwell_pass = (maxwell.gauss_law_residual < maxwell_tolerance_) &&
                        (maxwell.ampere_law_residual < maxwell_tolerance_) &&
                        (maxwell.faraday_law_residual < maxwell_tolerance_);

    file << "  Maxwell Validation: " << (maxwell_pass ? "PASS" : "FAIL") << "\n\n";

    file << "Flux Quantization:\n";
    file << "  Measured Flux: " << flux.flux << "\n";
    file << "  Expected Flux: " << flux.expected_flux << "\n";
    file << "  Winding Number: " << flux.winding_number << "\n";
    file << "  Quantization Error: " << flux.quantization_error << "\n";
    file << "  Is Quantized: " << (flux.is_quantized ? "YES" : "NO") << "\n\n";

    if (!flux.radii.empty()) {
        file << "  Flux vs Radius:\n";
        for (size_t i = 0; i < flux.radii.size(); ++i) {
            file << "    r = " << std::fixed << std::setprecision(2)
                 << flux.radii[i] << " : Φ = " << std::scientific
                 << std::setprecision(4) << flux.fluxes[i] << "\n";
        }
        file << "\n";
    }

    file << "Overall Result: "
         << (maxwell_pass && flux.is_quantized ? "PASS" : "FAIL") << "\n";

    file.close();
    std::cout << "Validation report written to " << filename << "\n";
}

void EMValidator::writeViolationMaps(
    const std::string& output_dir,
    const MaxwellValidation& maxwell) {

    // Write Gauss law violations
    std::string gauss_file = output_dir + "/gauss_violations.csv";
    std::ofstream file_gauss(gauss_file);
    if (file_gauss.is_open()) {
        for (int j = 0; j < Ny_; ++j) {
            for (int i = 0; i < Nx_; ++i) {
                file_gauss << maxwell.gauss_violation_map(i, j);
                if (i < Nx_ - 1) file_gauss << ",";
            }
            file_gauss << "\n";
        }
        file_gauss.close();
    }

    // Write Ampere law violations
    std::string ampere_file = output_dir + "/ampere_violations.csv";
    std::ofstream file_ampere(ampere_file);
    if (file_ampere.is_open()) {
        for (int j = 0; j < Ny_; ++j) {
            for (int i = 0; i < Nx_; ++i) {
                file_ampere << maxwell.ampere_violation_map(i, j);
                if (i < Nx_ - 1) file_ampere << ",";
            }
            file_ampere << "\n";
        }
        file_ampere.close();
    }

    // Write Faraday law violations
    std::string faraday_file = output_dir + "/faraday_violations.csv";
    std::ofstream file_faraday(faraday_file);
    if (file_faraday.is_open()) {
        for (int j = 0; j < Ny_; ++j) {
            for (int i = 0; i < Nx_; ++i) {
                file_faraday << maxwell.faraday_violation_map(i, j);
                if (i < Nx_ - 1) file_faraday << ",";
            }
            file_faraday << "\n";
        }
        file_faraday.close();
    }

    std::cout << "Violation maps written to " << output_dir << "/\n";
}

double EMValidator::analyzeGridConvergence(
    const std::vector<std::pair<int, double>>& violations_by_grid) {

    if (violations_by_grid.size() < 2) {
        return 0.0;
    }

    // Compute log-log slope to determine convergence order
    std::vector<double> log_h, log_err;

    for (const auto& [grid_size, violation] : violations_by_grid) {
        double h = 1.0 / grid_size;  // Assume domain size = 1
        log_h.push_back(std::log(h));
        log_err.push_back(std::log(violation));
    }

    // Linear regression on log-log data
    int n = log_h.size();
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;

    for (int i = 0; i < n; ++i) {
        sum_x += log_h[i];
        sum_y += log_err[i];
        sum_xy += log_h[i] * log_err[i];
        sum_xx += log_h[i] * log_h[i];
    }

    double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);

    return slope;  // Should be ~2 for O(h²) convergence
}

std::vector<std::pair<int, int>> EMValidator::findVortices(
    const Eigen::MatrixXd& B_z,
    double threshold) {

    std::vector<std::pair<int, int>> vortices;
    int Nx = B_z.rows();
    int Ny = B_z.cols();

    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            double B_center = std::abs(B_z(i, j));

            // Check if local maximum
            bool is_max = (B_center > threshold) &&
                         (B_center > std::abs(B_z(i-1, j))) &&
                         (B_center > std::abs(B_z(i+1, j))) &&
                         (B_center > std::abs(B_z(i, j-1))) &&
                         (B_center > std::abs(B_z(i, j+1)));

            if (is_max) {
                vortices.push_back({i, j});
            }
        }
    }

    return vortices;
}

double EMValidator::computeEffectiveAlpha(
    double measured_force,
    double theoretical_force) {

    if (std::abs(theoretical_force) < 1e-12) {
        return 0.0;
    }

    // α_eff = (measured/theoretical) * (1/137)
    return (measured_force / theoretical_force) / 137.0;
}