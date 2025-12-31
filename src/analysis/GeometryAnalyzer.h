#ifndef GEOMETRY_ANALYZER_H
#define GEOMETRY_ANALYZER_H

#include <eigen3/Eigen/Dense>
#include <vector>

/**
 * GeometryAnalyzer - Compute spacetime geometry and geodesics from R-field
 *
 * Physics Framework:
 * ------------------
 * Effective metric: ds² = -R²dt² + dx² + dy²
 * - R(x,y) is the Kuramoto sync field (order parameter)
 * - Metric signature: (-,+,+) in (t,x,y) coordinates
 * - g_μν = diag(-R², 1, 1) (diagonal metric)
 *
 * Christoffel Symbols:
 * -------------------
 * Γ^λ_μν = (1/2) g^λσ (∂_μ g_νσ + ∂_ν g_μσ - ∂_σ g_μν)
 *
 * For diagonal metric g_μν = diag(-R², 1, 1):
 * - Γ^x_00 = (1/2) g^xx ∂_x g_00 = (1/2)(1)(∂_x(-R²)) = -R ∂R/∂x
 * - Γ^y_00 = (1/2) g^yy ∂_y g_00 = (1/2)(1)(∂_y(-R²)) = -R ∂R/∂y
 *
 * Geodesic Equation:
 * -----------------
 * d²x^μ/dτ² + Γ^μ_αβ (dx^α/dτ)(dx^β/dτ) = 0
 *
 * For slowly varying R (|∂R/∂x| << 1):
 * - Dominant acceleration: a^i ≈ -Γ^i_00 = R ∂R/∂x^i
 * - Geodesic "force" from spacetime curvature
 *
 * References:
 * ----------
 * Carroll, "Spacetime and Geometry" (2004)
 * - Eq. 3.31: Christoffel symbols
 * - Eq. 3.44: Geodesic equation
 */
class GeometryAnalyzer {
public:
    /**
     * Metric tensor components g_μν
     * For diagonal metric: g_tt = -R², g_xx = g_yy = 1
     */
    struct Metric {
        double g_tt;  // -R²
        double g_xx;  // +1
        double g_yy;  // +1
        double g_tx;  // 0 (diagonal)
        double g_ty;  // 0 (diagonal)

        // Inverse metric g^μν (for raising indices)
        double g_inv_tt;  // -1/R²
        double g_inv_xx;  // +1
        double g_inv_yy;  // +1
    };

    /**
     * Christoffel symbols Γ^λ_μν
     * Only non-zero components stored for efficiency
     */
    struct ChristoffelSymbols {
        // Time-time components (dominant for geodesic acceleration)
        double Gamma_x_tt;  // Γ^x_00 = -R ∂R/∂x
        double Gamma_y_tt;  // Γ^y_00 = -R ∂R/∂y

        // Mixed time-space components
        double Gamma_t_tx;  // Γ^t_0x = (1/R) ∂R/∂x
        double Gamma_t_ty;  // Γ^t_0y = (1/R) ∂R/∂y

        // Space-space components (second order in ∂R)
        double Gamma_x_xx;  // Γ^x_xx = R ∂R/∂x
        double Gamma_y_yy;  // Γ^y_yy = R ∂R/∂y
        double Gamma_x_xy;  // Γ^x_xy (off-diagonal)
        double Gamma_y_xy;  // Γ^y_xy (off-diagonal)
    };

    /**
     * Compute metric tensor at given R value
     *
     * @param R Sync field value (order parameter)
     * @return Metric struct with g_μν and g^μν components
     */
    static Metric computeMetric(double R);

    /**
     * Get 3x3 metric matrix in (t,x,y) coordinates
     *
     * @param R Sync field value
     * @return 3x3 Eigen matrix: g_μν = diag(-R², 1, 1)
     */
    static Eigen::Matrix3d getMetricMatrix(double R);

    /**
     * Compute Christoffel symbols at grid point (ix, iy)
     *
     * Uses centered finite differences for derivatives:
     * ∂R/∂x ≈ (R[i+1,j] - R[i-1,j]) / (2Δx)
     *
     * @param R_field Sync field (size Nx*Ny, row-major)
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param ix Grid index in x (1 ≤ ix < Nx-1 for centered stencil)
     * @param iy Grid index in y (1 ≤ iy < Ny-1 for centered stencil)
     * @param dx Grid spacing in x (Planck units)
     * @param dy Grid spacing in y (Planck units)
     * @return ChristoffelSymbols struct with non-zero components
     */
    static ChristoffelSymbols computeChristoffel(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        int ix, int iy,
        double dx, double dy);

    /**
     * Compute Christoffel symbols at arbitrary position (x, y)
     * using bilinear interpolation of R-field
     *
     * @param R_field Sync field (size Nx*Ny, row-major)
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param x Physical position in x (Planck units)
     * @param y Physical position in y (Planck units)
     * @param L_domain Domain size (Planck units)
     * @return ChristoffelSymbols at interpolated position
     */
    static ChristoffelSymbols computeChristoffelInterpolated(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double x, double y,
        double L_domain);

    /**
     * Compute geodesic acceleration a^i = -Γ^i_μν v^μ v^ν
     *
     * For 4-velocity v^μ = (1, vx, vy):
     * - a^x ≈ -Γ^x_00 - 2 Γ^x_0x v^x - 2 Γ^x_0y v^y - Γ^x_xx v^x² - Γ^x_yy v^y²
     * - a^y ≈ -Γ^y_00 - 2 Γ^y_0x v^x - 2 Γ^y_0y v^y - Γ^y_xx v^x² - Γ^y_yy v^y²
     *
     * Dominant term: a^i ≈ -Γ^i_00 for non-relativistic velocities
     *
     * @param R_field Sync field (size Nx*Ny)
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param x Particle position in x (Planck units)
     * @param y Particle position in y (Planck units)
     * @param vx Particle velocity in x (c = 1)
     * @param vy Particle velocity in y (c = 1)
     * @param L_domain Domain size (Planck units)
     * @return (a_x, a_y) geodesic acceleration
     */
    static Eigen::Vector2d computeGeodesicAcceleration(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double x, double y,
        double vx, double vy,
        double L_domain);

    /**
     * Compute Riemann curvature scalar R at grid point
     *
     * For 2+1D spacetime with diagonal metric:
     * R = g^μν R_μν (scalar curvature)
     *
     * Approximate for weak curvature:
     * R ≈ -2∇²(ln R) / R² (in (x,y) space)
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param ix Grid index in x
     * @param iy Grid index in y
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return Riemann scalar curvature
     */
    static double computeRiemannScalar(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        int ix, int iy,
        double dx, double dy);

    /**
     * Compute Gaussian curvature K for 2D spatial slice
     *
     * For diagonal spatial metric h_ij = diag(1, 1):
     * K = 0 (flat Euclidean space)
     *
     * For effective 2+1D metric, spatial curvature induced by
     * time component: K ≈ (∇²R) / R (leading order)
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param ix Grid index in x
     * @param iy Grid index in y
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return Gaussian curvature
     */
    static double computeGaussianCurvature(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        int ix, int iy,
        double dx, double dy);

private:
    /**
     * Compute R-field gradient using centered finite differences
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param ix Grid index in x
     * @param iy Grid index in y
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return (∂R/∂x, ∂R/∂y)
     */
    static Eigen::Vector2d computeRGradient(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        int ix, int iy,
        double dx, double dy);

    /**
     * Compute R-field Laplacian using 5-point stencil
     *
     * ∇²R = ∂²R/∂x² + ∂²R/∂y²
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param ix Grid index in x
     * @param iy Grid index in y
     * @param dx Grid spacing in x
     * @param dy Grid spacing in y
     * @return ∇²R
     */
    static double computeLaplacian(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        int ix, int iy,
        double dx, double dy);

    /**
     * Bilinear interpolation of R-field at position (x, y)
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param x Physical position in x
     * @param y Physical position in y
     * @param L_domain Domain size
     * @return Interpolated R value
     */
    static double interpolateR(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double x, double y,
        double L_domain);

    /**
     * Bilinear interpolation of R-field gradient at position (x, y)
     *
     * @param R_field Sync field
     * @param Nx Grid size in x
     * @param Ny Grid size in y
     * @param x Physical position in x
     * @param y Physical position in y
     * @param L_domain Domain size
     * @return (∂R/∂x, ∂R/∂y) interpolated
     */
    static Eigen::Vector2d interpolateRGradient(
        const std::vector<double>& R_field,
        int Nx, int Ny,
        double x, double y,
        double L_domain);

    // Singularity guard: minimum R value to avoid metric singularity
    static constexpr double R_MIN = 0.1;

    // Grid boundary guard: minimum distance from edge for centered stencil
    static constexpr int BOUNDARY_GUARD = 1;
};

#endif // GEOMETRY_ANALYZER_H
