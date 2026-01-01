/**
 * GeodesicIntegrator.h
 *
 * Geodesic equation verification in SMFT curved spacetime
 *
 * Theory:
 * - SMFT metric: g_μν = R²(x,y) × diag(-(1-v²), 1, 1, 0)
 * - Christoffel symbols: Γ^μ_νλ computed from ∂g_μν
 * - Geodesic equation: d²x^μ/dτ² + Γ^μ_νλ(dx^ν/dτ)(dx^λ/dτ) = 0
 * - Particle trajectory from Dirac wavepacket center should follow geodesics
 *
 * Implementation:
 * - Christoffel computation from g_μν derivatives
 * - RK4 numerical geodesic integrator
 * - Trajectory comparison: Dirac evolution vs. geodesic equation
 * - Quality gate: <1% deviation from analytical geodesic
 */

#ifndef GEODESIC_INTEGRATOR_H
#define GEODESIC_INTEGRATOR_H

#include <vector>
#include <array>
#include <cmath>
#include <complex>

/**
 * GeodesicIntegrator - Christoffel symbol computation and geodesic integration
 *
 * Computes curved spacetime Christoffel symbols from SMFT metric and integrates
 * particle geodesics for comparison with Dirac wavepacket trajectories.
 */
class GeodesicIntegrator {
public:
    /**
     * Metric tensor components at a point
     * g_μν = R² × diag(-(1-v²), 1, 1, 0)
     */
    struct MetricTensor {
        double g00;  // -(1-v²)R²
        double g11;  // R²
        double g22;  // R²
        double g33;  // 0 (no time dimension in 2D spatial metric)

        // Off-diagonal (all zero for diagonal metric)
        double g01, g02, g03;
        double g12, g13, g23;
    };

    /**
     * Christoffel symbols Γ^μ_νλ
     * Stored as: christoffel[mu][nu][lambda]
     */
    struct ChristoffelSymbols {
        double value[4][4][4];  // All 64 components

        ChristoffelSymbols() {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        value[i][j][k] = 0.0;
        }
    };

    /**
     * Geodesic trajectory point
     */
    struct GeodesicPoint {
        double t;          // parameter along geodesic
        double x, y;       // spatial position
        double vx, vy;     // spatial velocity (dx/dt, dy/dt)
        double acceleration_magnitude;  // |d²x/dt²|
    };

    /**
     * Constructor
     * @param Nx Grid size x
     * @param Ny Grid size y
     * @param delta Vacuum potential Δ
     */
    GeodesicIntegrator(int Nx, int Ny, double delta = 2.5);

    /**
     * Destructor
     */
    ~GeodesicIntegrator();

    /**
     * Compute metric tensor at point (x, y)
     * @param x X coordinate
     * @param y Y coordinate
     * @param R_field Sync field R(x,y) at each grid point
     * @param v_field Velocity field |v| = sqrt(1 - 1/R²) for g_00 component
     * @return Metric tensor g_μν
     */
    MetricTensor computeMetric(double x, double y,
                               const std::vector<double>& R_field,
                               const std::vector<double>& v_field) const;

    /**
     * Compute Christoffel symbols at point (x, y) from metric derivatives
     * Uses finite differences: Γ^μ_νλ = (1/2)g^μρ(∂_ν g_ρλ + ∂_λ g_νρ - ∂_ρ g_νλ)
     *
     * @param x X coordinate
     * @param y Y coordinate
     * @param R_field Sync field R(x,y)
     * @param v_field Velocity field for g_00
     * @param h Finite difference step (default 0.1)
     * @return Christoffel symbols Γ^μ_νλ
     */
    ChristoffelSymbols computeChristoffel(double x, double y,
                                          const std::vector<double>& R_field,
                                          const std::vector<double>& v_field,
                                          double h = 0.1) const;

    /**
     * Geodesic acceleration from geodesic equation
     * a^μ = d²x^μ/dτ² = -Γ^μ_νλ(dx^ν/dτ)(dx^λ/dτ)
     *
     * @param pos Position {x, y}
     * @param vel Velocity {vx, vy}
     * @param christoffel Christoffel symbols at position
     * @param ax Output acceleration x
     * @param ay Output acceleration y
     */
    void computeGeodesicAcceleration(const std::array<double, 2>& pos,
                                     const std::array<double, 2>& vel,
                                     const ChristoffelSymbols& christoffel,
                                     double& ax, double& ay) const;

    /**
     * Integrate geodesic using RK4 method
     * Solves: d²x^μ/dτ² + Γ^μ_νλ(dx^ν/dτ)(dx^λ/dτ) = 0
     *
     * @param initial_pos Initial position {x0, y0}
     * @param initial_vel Initial velocity {vx0, vy0}
     * @param R_field Sync field
     * @param v_field Velocity field
     * @param dt Time step
     * @param num_steps Number of integration steps
     * @return Vector of geodesic trajectory points
     */
    std::vector<GeodesicPoint> integrateGeodesic(
        const std::array<double, 2>& initial_pos,
        const std::array<double, 2>& initial_vel,
        const std::vector<double>& R_field,
        const std::vector<double>& v_field,
        double dt, int num_steps) const;

    /**
     * Compare Dirac wavepacket trajectory to geodesic prediction
     * Computes deviation as: |trajectory - geodesic| / |geodesic|
     *
     * @param dirac_trajectory Positions from Dirac evolution {t, x, y}
     * @param geodesic_trajectory Positions from geodesic equation {t, x, y}
     * @return Deviation for each time step
     */
    std::vector<double> compareTrajectories(
        const std::vector<std::array<double, 3>>& dirac_trajectory,
        const std::vector<std::array<double, 3>>& geodesic_trajectory) const;

    /**
     * Get grid size
     */
    int getNx() const { return _Nx; }
    int getNy() const { return _Ny; }

private:
    int _Nx, _Ny;
    double _delta;

    /**
     * Bilinear interpolation for field values
     */
    double interpolateField(double x, double y,
                           const std::vector<double>& field) const;

    /**
     * Finite difference gradient of field
     */
    void fieldGradient(double x, double y,
                      const std::vector<double>& field,
                      double& grad_x, double& grad_y,
                      double h = 0.1) const;

    /**
     * Finite difference second gradient (Laplacian diagonal)
     */
    void fieldSecondGradient(double x, double y,
                            const std::vector<double>& field,
                            double& d2_xx, double& d2_yy,
                            double h = 0.1) const;

    /**
     * Check if point is within grid bounds
     */
    bool isInBounds(double x, double y) const {
        return x >= 0 && x <= _Nx && y >= 0 && y <= _Ny;
    }

    /**
     * Clamp coordinate to valid range
     */
    void clampCoordinate(double& x, double& y) const {
        x = std::max(0.0, std::min(static_cast<double>(_Nx - 1), x));
        y = std::max(0.0, std::min(static_cast<double>(_Ny - 1), y));
    }
};

#endif // GEODESIC_INTEGRATOR_H
