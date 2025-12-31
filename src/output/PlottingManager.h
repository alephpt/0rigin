#ifndef PLOTTING_MANAGER_H
#define PLOTTING_MANAGER_H

#include <string>
#include <vector>
#include <map>
#include <complex>
#include "simulations/ObservableComputer.h"  // Need full definition

// Forward declare to avoid including matplotlib-cpp in header
namespace matplotlibcpp {}

// Forward declarations
struct ValidationResult;

class PlottingManager {
public:
    struct PlotStyle {
        // Standard colors for velocities
        static const std::map<std::string, std::string> VELOCITY_COLORS;
        static const std::map<std::string, std::string> GRID_COLORS;

        // Standard sizes
        static const int DPI = 300;
        static const int FIGURE_WIDTH_SINGLE = 10;
        static const int FIGURE_HEIGHT_SINGLE = 8;
        static const int FIGURE_WIDTH_MULTI = 15;
        static const int FIGURE_HEIGHT_MULTI = 10;

        // Colormaps
        static const std::string COLORMAP_THETA;
        static const std::string COLORMAP_R;
        static const std::string COLORMAP_ERROR;
    };

    PlottingManager(const std::string& outputDir);
    ~PlottingManager();

    // Initialize Python and matplotlib
    bool initialize();

    // Standard plot sets
    void plotObservables6Panel(
        const std::vector<ObservableComputer::Observables>& obsTimeSeries,
        double velocity,
        int gridSize,
        int N
    );

    void plotConservationLaws(
        const std::vector<ObservableComputer::Observables>& obsTimeSeries,
        double velocity,
        int gridSize,
        int N
    );

    void plotSpatialFields(
        const std::vector<std::complex<float>>& theta_field,
        const std::vector<float>& R_field,
        int Nx,
        int Ny,
        double velocity,
        int gridSize,
        int N,
        double time = 0.0
    );

    void plotGammaValidation(
        const std::vector<double>& time,
        const std::vector<double>& gamma_measured,
        double gamma_theory,
        double velocity,
        int gridSize,
        int N
    );

private:
    std::string outputDir_;
    bool initialized_;

    // Utility functions
    std::string makeFilename(
        const std::string& scenario,
        const std::string& plotType,
        int gridSize,
        double velocity,
        int N,
        const std::string& extension = "png"
    );

    void configurePlotStyle();

    // Get color for velocity
    std::string getVelocityColor(double velocity) const;
    std::string getGridColor(int gridSize) const;
};

#endif // PLOTTING_MANAGER_H
