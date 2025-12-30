#pragma once

#include <vector>
#include <map>
#include <string>
#include <memory>
#include <complex>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>

namespace SMFT {

/**
 * UniversalityClassifier: Comprehensive finite-size scaling analysis
 *
 * This class implements the complete FSS strategy for determining the universality
 * class of the SMFT phase transition through extraction of critical exponents
 * and data collapse analysis.
 */
class UniversalityClassifier {
public:
    /**
     * Structure to hold critical exponents with uncertainties
     */
    struct CriticalExponents {
        // Primary exponents
        double beta;            // Order parameter exponent
        double beta_err;        // Uncertainty in beta
        double nu;              // Correlation length exponent
        double nu_err;          // Uncertainty in nu
        double gamma;           // Susceptibility exponent
        double gamma_err;       // Uncertainty in gamma
        double eta;             // Anomalous dimension
        double eta_err;         // Uncertainty in eta

        // Derived exponents
        double alpha;           // Specific heat exponent
        double alpha_err;       // Uncertainty in alpha
        double delta;           // Critical isotherm exponent
        double delta_err;       // Uncertainty in delta

        // Critical point
        double sigma_c;         // Critical noise strength
        double sigma_c_err;     // Uncertainty in sigma_c

        // Quality metrics
        double chi2_collapse;   // Chi-squared for data collapse
        double collapse_quality; // Overall collapse quality (0-1)

        CriticalExponents() :
            beta(0), beta_err(0), nu(0), nu_err(0),
            gamma(0), gamma_err(0), eta(0), eta_err(0),
            alpha(0), alpha_err(0), delta(0), delta_err(0),
            sigma_c(0), sigma_c_err(0),
            chi2_collapse(0), collapse_quality(0) {}
    };

    /**
     * Data point for a single grid size and noise value
     */
    struct DataPoint {
        int L;                  // Grid size
        double sigma;           // Noise strength
        double R_mean;          // Mean order parameter
        double R_err;           // Error in R
        double R2_mean;         // Mean of R^2
        double R4_mean;         // Mean of R^4
        double chi;             // Susceptibility
        double chi_err;         // Error in susceptibility
        double U_L;             // Binder cumulant
        double U_L_err;         // Error in Binder cumulant
        int samples;            // Number of samples

        DataPoint() : L(0), sigma(0), R_mean(0), R_err(0),
                     R2_mean(0), R4_mean(0), chi(0), chi_err(0),
                     U_L(0), U_L_err(0), samples(0) {}
    };

    /**
     * Structure for correlation function data
     */
    struct CorrelationData {
        std::vector<double> r;          // Distance values
        std::vector<double> G_r;        // Correlation function G(r)
        std::vector<double> G_r_err;    // Error in G(r)
        double xi;                      // Correlation length
        double xi_err;                  // Error in correlation length
    };

    // Constructor/Destructor
    UniversalityClassifier();
    ~UniversalityClassifier() = default;

    // Data input methods
    void addDataPoint(const DataPoint& point);
    void addGridSizeData(int L, const std::vector<double>& sigma,
                        const std::vector<double>& R,
                        const std::vector<double>& R2,
                        const std::vector<double>& R4);
    void addCorrelationData(int L, double sigma, const CorrelationData& corr_data);

    // Critical point determination
    double findCriticalPoint();
    double findCriticalPointFromBinder();
    double findCriticalPointFromSusceptibility();

    // Exponent extraction methods
    CriticalExponents extractExponents();
    double extractBeta(double sigma_c);
    double extractNu(double sigma_c);
    double extractGamma(double sigma_c);
    double extractEta(double sigma_c);

    // Data collapse analysis
    bool performDataCollapse(double sigma_c, double beta_nu, double nu,
                            double& chi2, double& quality);
    void optimizeCollapse(double& sigma_c, double& beta, double& nu);

    // Scaling relation checks
    bool verifyScalingRelations(const CriticalExponents& exp);
    double checkFisherRelation(const CriticalExponents& exp);
    double checkRushbrookeRelation(const CriticalExponents& exp);
    double checkJosephsonRelation(const CriticalExponents& exp);

    // Universality classification
    std::string classifyUniversality();
    std::string classifyFromExponents(const CriticalExponents& exp);
    double compareToKnownClass(const CriticalExponents& exp, const std::string& class_name);

    // Reporting and output
    void writeDataCollapseReport(const std::string& output_dir);
    void writeCriticalExponentReport(const std::string& output_dir);
    void writeClassificationReport(const std::string& output_dir);
    void generatePythonScripts(const std::string& output_dir);

    // Data export for external analysis
    void exportDataForPython(const std::string& filename);
    void exportCollapseData(const std::string& filename);

    // Known universality classes (2D)
    static CriticalExponents getIsingExponents();
    static CriticalExponents getXYExponents();
    static CriticalExponents getPotts3Exponents();
    static CriticalExponents getMeanFieldExponents();

private:
    // Data storage
    std::map<int, std::vector<DataPoint>> data_by_size;  // L -> data points
    std::map<std::pair<int, double>, CorrelationData> correlation_data;
    std::vector<int> grid_sizes;

    // Analysis parameters
    double sigma_min, sigma_max;
    int n_bootstrap;
    double convergence_tol;

    // Helper methods for fitting
    void linearFit(const std::vector<double>& x, const std::vector<double>& y,
                  double& slope, double& intercept, double& error);
    void powerLawFit(const std::vector<double>& x, const std::vector<double>& y,
                     double& exponent, double& amplitude, double& error);
    double computeChi2(const std::vector<double>& observed,
                       const std::vector<double>& expected);

    // Bootstrap error estimation
    void bootstrapError(const std::vector<double>& data, double& mean, double& error);
    std::vector<double> resampleData(const std::vector<double>& data);

    // Finite-size scaling functions
    double scalingFunction(double x, const std::vector<double>& params);
    void fitScalingFunction(const std::vector<double>& x, const std::vector<double>& y,
                           std::vector<double>& params, double& chi2);

    // Binder cumulant analysis
    std::pair<double, double> findBinderCrossing(int L1, int L2);
    void analyzeBinderCrossings(std::vector<std::pair<double, double>>& crossings);

    // Susceptibility peak analysis
    std::pair<double, double> findSusceptibilityPeak(int L);
    void analyzeSusceptibilityPeaks(std::vector<std::pair<int, double>>& peaks);

    // Correlation length extraction
    double extractCorrelationLength(const CorrelationData& data);
    double fitExponentialCorrelation(const std::vector<double>& r,
                                    const std::vector<double>& G);

    // Utility methods
    void sortDataBySigma(int L);
    std::vector<double> interpolateData(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& x_new);
    double findIntersection(const std::vector<double>& x1, const std::vector<double>& y1,
                           const std::vector<double>& x2, const std::vector<double>& y2);
};

// Helper class for managing test configurations
class UniversalityTestConfig {
public:
    struct TestParameters {
        std::vector<int> grid_sizes;      // [32, 64, 128, 256, 512]
        double sigma_min;                  // 0.0
        double sigma_max;                  // 2.0
        int n_sigma_points;                // 41
        int equilibration_steps;          // 5000
        int measurement_steps;             // 10000
        double dt;                         // 0.01
        bool measure_correlations;         // true
        int correlation_max_r;             // L/4
        std::string output_base;          // "output/universality_scan"
    };

    UniversalityTestConfig();

    void generateConfigFiles(const std::string& config_dir);
    void generateGridConfig(int L, const std::string& filename);
    TestParameters getDefaultParameters();

private:
    TestParameters params;

    void writeYAMLConfig(int L, double sigma, const std::string& filename);
    std::string generateConfigName(int L, double sigma);
};

} // namespace SMFT