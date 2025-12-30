#ifndef ENERGY_BUDGET_H
#define ENERGY_BUDGET_H

#include "../SMFTEngine.h"
#include "../DiracEvolution.h"
#include "../simulations/ObservableComputer.h"
#include <vector>
#include <string>
#include <map>
#include <fstream>

/**
 * EnergyBudget - Comprehensive energy budget tracking for SMFT simulations
 *
 * Following Phase 1 Task 2 strategy, this class tracks ALL energy components to
 * distinguish numerical error from physical damping. The complete energy budget:
 *
 * E_total = E_Dirac + E_Kuramoto + E_Coupling + E_EM
 *
 * Where:
 * - E_Dirac = T_Dirac + V_Dirac (kinetic + mass potential)
 * - E_Kuramoto = E_gradient + E_sync (phase gradients + synchronization)
 * - E_Coupling = λ∫ψ†ψ·R (Kuramoto-Dirac interaction)
 * - E_EM = ∫(E² + B²)/(8π) (electromagnetic field energy if enabled)
 *
 * Energy conservation equation:
 * dE_total/dt = -P_dissipated
 *
 * Where P_dissipated = γ·∫|∇θ|² dx (power lost to Kuramoto damping)
 *
 * For undamped systems (γ=0), energy should be conserved to machine precision
 * with Strang splitting (already implemented).
 */
class EnergyBudget {
public:
    /**
     * Complete energy component breakdown
     */
    struct EnergyComponents {
        // Dirac field components
        double T_dirac_kinetic;      // ∫ψ†(-iα·∇)ψ dx
        double V_dirac_mass;          // ∫ψ†(β·m)ψ dx where m = Δ·R

        // Kuramoto field components
        double E_kuramoto_gradient;   // (1/2)∫|∇θ|² dx
        double E_kuramoto_sync;       // -K·Σcos(θᵢ - θⱼ)

        // Coupling energy
        double E_coupling;            // λ∫|ψ|²·R dx

        // EM field energy (if enabled)
        double E_em_field;            // ∫(E² + B²)/(8π) dx

        // Total energy
        double E_total;               // Sum of all components

        // Dissipation tracking
        double P_dissipated;          // γ·∫|∇θ|² dx (instantaneous power loss)
        double E_dissipated_cumulative; // ∫P_dissipated dt (total energy lost)

        // Conservation metrics
        double dE_dt;                 // Time derivative of total energy
        double conservation_error;    // |dE/dt + P_dissipated| (should be ~0)
    };

    /**
     * Time series data for energy evolution
     */
    struct TimeSeries {
        std::vector<double> time;
        std::vector<EnergyComponents> components;

        // Derived metrics
        std::vector<double> relative_error;    // |E(t) - E₀|/E₀
        std::vector<double> energy_drift;      // E(t) - E₀ + ∫P dt

        void clear() {
            time.clear();
            components.clear();
            relative_error.clear();
            energy_drift.clear();
        }
    };

    /**
     * Constructor
     * @param engine SMFTEngine instance
     * @param output_dir Directory for energy analysis outputs
     */
    EnergyBudget(SMFTEngine* engine, const std::string& output_dir);

    /**
     * Set physics parameters from TestConfig
     * @param delta Mass gap parameter
     * @param coupling Kuramoto-Dirac coupling strength
     * @param K Kuramoto coupling strength
     * @param damping Kuramoto damping coefficient
     * @param L_domain Domain size in Planck lengths
     */
    void setPhysicsParameters(float delta, float coupling, float K, float damping, float L_domain);

    /**
     * Set grid dimensions
     * @param Nx Grid width
     * @param Ny Grid height
     */
    void setGridDimensions(int Nx, int Ny) { _Nx = Nx; _Ny = Ny; }

    /**
     * Compute all energy components for current state
     *
     * @param dirac_evolution DiracEvolution instance (can be null for Kuramoto-only)
     * @param time Current simulation time
     * @param compute_em Whether to compute EM field energy
     * @return Complete energy component breakdown
     */
    EnergyComponents computeComponents(
        DiracEvolution* dirac_evolution,
        double time,
        bool compute_em = false);

    /**
     * Track energy components over time
     * Adds current state to time series
     *
     * @param t Current time
     * @param components Energy components at time t
     */
    void trackTimeSeries(double t, const EnergyComponents& components);

    /**
     * Analyze energy conservation over full simulation
     * Generates comprehensive energy budget report
     *
     * @param output_file Path to write analysis results
     * @return True if energy is conserved within tolerance
     */
    bool analyzeConservation(const std::string& output_file);

    /**
     * Write energy time series to CSV file
     *
     * @param csv_file Path to output CSV file
     */
    void writeTimeSeriesCSV(const std::string& csv_file);

    /**
     * Compute Kuramoto field energy components
     *
     * @param theta Phase field θ(x,y)
     * @param R_field Synchronization field R(x,y)
     * @param K Kuramoto coupling strength
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return {E_gradient, E_sync} tuple
     */
    static std::pair<double, double> computeKuramotoEnergy(
        const std::vector<float>& theta,
        const std::vector<float>& R_field,
        float K,
        float dx, float dy);

    /**
     * Compute coupling energy between Dirac and Kuramoto fields
     *
     * @param psi Dirac spinor (4-component)
     * @param R_field Synchronization field
     * @param coupling Coupling strength λ
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return E_coupling = λ∫|ψ|²·R dx
     */
    static double computeCouplingEnergy(
        const std::vector<std::complex<double>>& psi,
        const std::vector<float>& R_field,
        float coupling,
        float dx, float dy);

    /**
     * Compute dissipation rate from Kuramoto damping
     *
     * @param theta Phase field θ(x,y)
     * @param damping Damping coefficient γ
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return P_dissipated = γ·∫|∇θ|² dx
     */
    static double computeDissipationRate(
        const std::vector<float>& theta,
        float damping,
        float dx, float dy);

    /**
     * Perform Richardson extrapolation for dt→0 limit
     * Tests true energy conservation by running with dt, dt/2, dt/4
     *
     * @param dt_values Vector of timestep values to test
     * @param output_file Path to write extrapolation results
     * @return Extrapolated energy conservation error at dt=0
     */
    double richardsonExtrapolation(
        const std::vector<double>& dt_values,
        const std::string& output_file);

    /**
     * Compare Strang splitting vs first-order Euler
     * Validates that Strang provides superior energy conservation
     *
     * @param output_file Path to write comparison results
     * @return Ratio of Euler error to Strang error
     */
    double compareIntegrators(const std::string& output_file);

    /**
     * Get conservation metrics
     *
     * @return Map of metric name to value
     */
    std::map<std::string, double> getConservationMetrics() const;

    /**
     * Reset time series data
     */
    void reset() { _time_series.clear(); }

private:
    SMFTEngine* _engine;
    std::string _output_dir;
    TimeSeries _time_series;

    // Initial energy for conservation tracking
    double _E0;
    bool _E0_set;

    // Grid parameters cached for efficiency
    int _Nx, _Ny;
    float _dx, _dy;

    // Physics parameters
    float _delta;
    float _coupling;
    float _K;
    float _damping;

    /**
     * Helper: Compute phase gradient energy using finite differences
     */
    double computePhaseGradientEnergy(
        const std::vector<float>& theta,
        float dx, float dy);

    /**
     * Helper: Compute synchronization energy from phase differences
     */
    double computeSynchronizationEnergy(
        const std::vector<float>& theta,
        float K);

    /**
     * Helper: Validate energy conservation criterion
     * Returns true if |dE/dt + P_dissipated| < tolerance
     */
    bool validateConservation(double tolerance = 1e-6);

    /**
     * Helper: Fit exponential decay for damped systems
     * Returns decay rate γ_measured
     */
    double fitExponentialDecay();
};

#endif // ENERGY_BUDGET_H