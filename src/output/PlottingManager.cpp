// matplotlib-cpp must be included first before other headers
#include "matplotlibcpp.h"

#include "PlottingManager.h"
#include "simulations/ObservableComputer.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>

namespace plt = matplotlibcpp;

// Static member initialization
const std::map<std::string, std::string> PlottingManager::PlotStyle::VELOCITY_COLORS = {
    {"0.0", "black"},
    {"0.3", "blue"},
    {"0.5", "green"},
    {"0.7", "red"},
    {"0.8", "purple"},
    {"0.9", "orange"}
};

const std::map<std::string, std::string> PlottingManager::PlotStyle::GRID_COLORS = {
    {"16", "dodgerblue"},
    {"32", "mediumseagreen"},
    {"64", "forestgreen"},
    {"128", "darkorange"},
    {"256", "crimson"},
    {"512", "darkviolet"}
};

const std::string PlottingManager::PlotStyle::COLORMAP_THETA = "twilight";
const std::string PlottingManager::PlotStyle::COLORMAP_R = "viridis";
const std::string PlottingManager::PlotStyle::COLORMAP_ERROR = "RdYlGn_r";

PlottingManager::PlottingManager(const std::string& outputDir)
    : outputDir_(outputDir), initialized_(false) {}

PlottingManager::~PlottingManager() {
    // Clean up matplotlib if needed
    if (initialized_) {
        try {
            plt::close();
        } catch (...) {
            // Ignore errors during cleanup
        }
    }
}

bool PlottingManager::initialize() {
    if (initialized_) {
        return true;
    }

    try {
        // Set non-interactive backend for headless operation
        plt::backend("Agg");

        // Test if matplotlib works
        plt::figure();
        plt::close();

        initialized_ = true;
        std::cout << "[PlottingManager] Initialized matplotlib backend (Agg)" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[PlottingManager] Failed to initialize: " << e.what() << std::endl;
        return false;
    }
}

void PlottingManager::plotObservables6Panel(
    const std::vector<ObservableComputer::Observables>& obsTimeSeries,
    double velocity,
    int gridSize,
    int N
) {
    if (!initialize()) {
        std::cerr << "[PlottingManager] Cannot plot - initialization failed" << std::endl;
        return;
    }

    // Extract time series data
    std::vector<double> time, norm, E_total, E_kin, E_pot;
    std::vector<double> p_mag, pos_x, pos_y;
    std::vector<double> R_avg, R_min, R_max, gamma_measured;

    for (const auto& obs : obsTimeSeries) {
        time.push_back(obs.time);
        norm.push_back(obs.norm);
        E_total.push_back(obs.energy_total);
        E_kin.push_back(obs.energy_kinetic);
        E_pot.push_back(obs.energy_potential);

        double px = obs.momentum_x.real();
        double py = obs.momentum_y.real();
        p_mag.push_back(std::sqrt(px*px + py*py));

        pos_x.push_back(obs.position_x.real());
        pos_y.push_back(obs.position_y.real());

        R_avg.push_back(obs.R_avg);
        R_min.push_back(obs.R_min);
        R_max.push_back(obs.R_max);

        // Compute gamma from observables
        double E = obs.energy_total;
        double p = std::sqrt(px*px + py*py);
        double m_eff_sq = std::max(E*E - p*p, 0.0);
        double m_eff = std::sqrt(m_eff_sq);
        double gamma = (obs.R_avg > 1e-10) ? m_eff / (1.0 * obs.R_avg) : 1.0;  // DELTA = 1.0
        gamma_measured.push_back(gamma);
    }

    // Theoretical gamma
    double gamma_theory = (std::abs(velocity) < 1e-10) ? 1.0 :
                          1.0 / std::sqrt(1.0 - velocity*velocity);

    try {
        // NOTE: subplot() doesn't work reliably with matplotlib-cpp
        // Create individual plots instead

        // Plot 1: Norm
        plt::figure();
        plt::figure_size(PlotStyle::FIGURE_WIDTH_SINGLE * 100,
                         PlotStyle::FIGURE_HEIGHT_SINGLE * 100);
        std::map<std::string, std::string> norm_style = {{"color", "blue"}, {"linewidth", "2"}};
        plt::plot(time, norm, norm_style);
        plt::axhline(1.0);
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("Norm");
        plt::title("Wavefunction Norm");
        plt::grid(true);

        // Subplot 2: Energy
        plt::subplot(2, 3, 2);
        plt::named_plot("Total", time, E_total, "k-");
        plt::named_plot("Kinetic", time, E_kin, "b--");
        plt::named_plot("Potential", time, E_pot, "r--");
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("Energy (m_P c^2)");
        plt::title("Energy Conservation");
        plt::legend();
        plt::grid(true);

        // Subplot 3: Momentum
        plt::subplot(2, 3, 3);
        plt::plot(time, p_mag, "g-");
        double p_expected = gamma_theory * 1.0 * velocity;  // DELTA = 1.0
        if (velocity > 0) {
            plt::axhline(p_expected);
        }
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("|p| (m_P c)");
        plt::title("Momentum Magnitude");
        plt::grid(true);

        // Subplot 4: Trajectory
        plt::subplot(2, 3, 4);
        plt::plot(pos_x, pos_y, "b-");
        if (!pos_x.empty() && !pos_y.empty()) {
            plt::plot({pos_x[0]}, {pos_y[0]}, "go");
            plt::plot({pos_x.back()}, {pos_y.back()}, "ro");
        }
        plt::xlabel("<x> (grid units)");
        plt::ylabel("<y> (grid units)");
        plt::title("Wavepacket Trajectory");
        plt::grid(true);

        // Subplot 5: R-field
        plt::subplot(2, 3, 5);
        plt::named_plot("R_avg", time, R_avg, "k-");
        plt::named_plot("R_max", time, R_max, "b--");
        plt::named_plot("R_min", time, R_min, "r--");
        plt::axhline(0.5);
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("R");
        plt::title("Kuramoto R-field");
        plt::legend();
        plt::grid(true);

        // Subplot 6: Gamma
        plt::subplot(2, 3, 6);
        plt::named_plot("Measured", time, gamma_measured, "b-");
        plt::axhline(gamma_theory);
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("gamma");
        plt::title("Gamma Factor");
        plt::legend();
        plt::grid(true);

        // Overall title
        std::ostringstream title;
        title << "Grid " << gridSize << "x" << gridSize
              << ", v=" << std::fixed << std::setprecision(1) << velocity << "c"
              << ", N=" << N;
        plt::suptitle(title.str());

        // Save
        std::string filename = makeFilename("observables", "6panel", gridSize, velocity, N);
        std::string fullpath = outputDir_ + "/" + filename;
        plt::save(fullpath, PlotStyle::DPI);
        plt::close();

        std::cout << "[PlottingManager] Saved: " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[PlottingManager] Error creating 6-panel plot: " << e.what() << std::endl;
    }
}

void PlottingManager::plotConservationLaws(
    const std::vector<ObservableComputer::Observables>& obsTimeSeries,
    double velocity,
    int gridSize,
    int N
) {
    if (!initialize()) {
        std::cerr << "[PlottingManager] Cannot plot - initialization failed" << std::endl;
        return;
    }

    // Extract conservation data
    std::vector<double> time, norm_error, energy_rel_error;
    double E0 = 0.0;

    for (size_t i = 0; i < obsTimeSeries.size(); ++i) {
        const auto& obs = obsTimeSeries[i];
        if (i == 0) {
            E0 = obs.energy_total;
        }

        time.push_back(obs.time);
        norm_error.push_back(obs.norm_error);

        double E_rel = (E0 != 0.0) ? (obs.energy_total - E0) / E0 : 0.0;
        energy_rel_error.push_back(E_rel);
    }

    try {
        plt::figure_size(PlotStyle::FIGURE_WIDTH_SINGLE * 100,
                         PlotStyle::FIGURE_HEIGHT_SINGLE * 100);

        // Subplot 1: Norm error
        plt::subplot(2, 1, 1);
        plt::plot(time, norm_error, "b-");
        plt::axhline(0.0);
        plt::axhline(0.005);
        plt::axhline(-0.005);
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("Norm Error");
        plt::title("Norm Conservation: ||Psi||^2 - 1");
        plt::grid(true);

        // Subplot 2: Energy error
        plt::subplot(2, 1, 2);
        plt::plot(time, energy_rel_error, "r-");
        plt::axhline(0.0);
        plt::axhline(0.01);
        plt::axhline(-0.01);
        plt::xlabel("Time (l_P/c)");
        plt::ylabel("Relative Energy Error");
        plt::title("Energy Conservation: (E - E_0) / E_0");
        plt::grid(true);

        // Overall title
        std::ostringstream title;
        title << "Conservation Laws: Grid " << gridSize << "x" << gridSize
              << ", v=" << std::fixed << std::setprecision(1) << velocity << "c"
              << ", N=" << N;
        plt::suptitle(title.str());

        // Save
        std::string filename = makeFilename("conservation", "laws", gridSize, velocity, N);
        std::string fullpath = outputDir_ + "/" + filename;
        plt::save(fullpath, PlotStyle::DPI);
        plt::close();

        std::cout << "[PlottingManager] Saved: " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[PlottingManager] Error creating conservation plot: " << e.what() << std::endl;
    }
}

void PlottingManager::plotSpatialFields(
    const std::vector<std::complex<float>>& theta_field,
    const std::vector<float>& R_field,
    int Nx,
    int Ny,
    double velocity,
    int gridSize,
    int N,
    double time
) {
    if (!initialize()) {
        std::cerr << "[PlottingManager] Cannot plot - initialization failed" << std::endl;
        return;
    }

    // Extract real theta values (phase) and R values
    std::vector<std::vector<double>> theta_grid(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> R_grid(Ny, std::vector<double>(Nx));

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int idx = i * Ny + j;
            theta_grid[j][i] = std::arg(theta_field[idx]);
            R_grid[j][i] = static_cast<double>(R_field[idx]);
        }
    }

    try {
        // NOTE: matplotlib-cpp imshow requires NumPy, which we disabled
        // For now, skip spatial field plotting or implement alternative visualization
        std::cout << "[PlottingManager] Spatial field plotting requires NumPy (currently disabled)" << std::endl;
        std::cout << "[PlottingManager] Skipping spatial field plot. Use Python scripts for 2D visualization." << std::endl;
        return;

        /* TODO: Enable NumPy support or implement alternative 2D plotting
        plt::figure_size(1200, 500);

        // Subplot 1: Theta field
        plt::subplot(1, 2, 1);
        std::map<std::string, std::string> theta_map = {
            {"cmap", PlotStyle::COLORMAP_THETA},
            {"origin", "lower"},
            {"aspect", "auto"}
        };
        plt::imshow(theta_grid, theta_map);
        plt::xlabel("x (grid units)");
        plt::ylabel("y (grid units)");
        std::ostringstream theta_title;
        theta_title << "theta(x,y,t=" << std::fixed << std::setprecision(2) << time << ") - Phase Field";
        plt::title(theta_title.str());
        plt::colorbar();

        // Subplot 2: R field
        plt::subplot(1, 2, 2);
        std::map<std::string, std::string> R_map = {
            {"cmap", PlotStyle::COLORMAP_R},
            {"origin", "lower"},
            {"aspect", "auto"}
        };
        plt::imshow(R_grid, R_map);
        plt::xlabel("x (grid units)");
        plt::ylabel("y (grid units)");
        std::ostringstream R_title;
        R_title << "R(x,y,t=" << std::fixed << std::setprecision(2) << time << ") - Order Parameter";
        plt::title(R_title.str());
        plt::colorbar();
        */

        // Overall title
        std::ostringstream title;
        title << "Spatial Fields: Grid " << gridSize << "x" << gridSize
              << ", v=" << std::fixed << std::setprecision(1) << velocity << "c"
              << ", N=" << N;
        plt::suptitle(title.str());

        // Save
        std::string filename = makeFilename("spatial_fields", "2panel", gridSize, velocity, N);
        std::string fullpath = outputDir_ + "/" + filename;
        plt::save(fullpath, PlotStyle::DPI);
        plt::close();

        std::cout << "[PlottingManager] Saved: " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[PlottingManager] Error creating spatial fields plot: " << e.what() << std::endl;
    }
}

void PlottingManager::plotGammaValidation(
    const std::vector<double>& time,
    const std::vector<double>& gamma_measured,
    double gamma_theory,
    double velocity,
    int gridSize,
    int N
) {
    if (!initialize()) {
        std::cerr << "[PlottingManager] Cannot plot - initialization failed" << std::endl;
        return;
    }

    try {
        plt::figure_size(PlotStyle::FIGURE_WIDTH_SINGLE * 100,
                         PlotStyle::FIGURE_HEIGHT_SINGLE * 100);

        plt::named_plot("Measured", time, gamma_measured, "b-");
        plt::axhline(gamma_theory);

        plt::xlabel("Time (l_P/c)");
        plt::ylabel("gamma");
        plt::title("Relativistic Gamma Factor Validation");
        plt::legend();
        plt::grid(true);

        // Add text annotation
        std::ostringstream annotation;
        annotation << "v = " << std::fixed << std::setprecision(2) << velocity << "c\n"
                   << "gamma_theory = " << std::setprecision(4) << gamma_theory;
        plt::text(time[time.size()/10], gamma_theory * 1.1, annotation.str());

        // Overall title
        std::ostringstream title;
        title << "Grid " << gridSize << "x" << gridSize << ", N=" << N;
        plt::suptitle(title.str());

        // Save
        std::string filename = makeFilename("gamma", "validation", gridSize, velocity, N);
        std::string fullpath = outputDir_ + "/" + filename;
        plt::save(fullpath, PlotStyle::DPI);
        plt::close();

        std::cout << "[PlottingManager] Saved: " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[PlottingManager] Error creating gamma validation plot: " << e.what() << std::endl;
    }
}

std::string PlottingManager::makeFilename(
    const std::string& scenario,
    const std::string& plotType,
    int gridSize,
    double velocity,
    int N,
    const std::string& extension
) {
    std::ostringstream ss;
    ss << scenario << "_" << plotType << "_"
       << gridSize << "x" << gridSize << "_"
       << "v" << std::fixed << std::setprecision(1) << velocity << "_"
       << "N" << N << "." << extension;
    return ss.str();
}

void PlottingManager::configurePlotStyle() {
    // Reserved for future style configuration
}

std::string PlottingManager::getVelocityColor(double velocity) const {
    std::ostringstream key;
    key << std::fixed << std::setprecision(1) << velocity;

    auto it = PlotStyle::VELOCITY_COLORS.find(key.str());
    if (it != PlotStyle::VELOCITY_COLORS.end()) {
        return it->second;
    }
    return "gray";  // Default color
}

std::string PlottingManager::getGridColor(int gridSize) const {
    auto it = PlotStyle::GRID_COLORS.find(std::to_string(gridSize));
    if (it != PlotStyle::GRID_COLORS.end()) {
        return it->second;
    }
    return "gray";  // Default color
}
