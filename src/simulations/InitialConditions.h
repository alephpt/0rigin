#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <vector>
#include <cmath>

/**
 * InitialConditions - Utilities for creating specialized initial conditions
 *
 * Purpose: Generate domain-specific and defect-based initial conditions
 * for Phase 3 vacuum structure tests (Tests 3.1 and 3.2).
 *
 * All coordinates are in grid units unless specified as physical.
 * Physical coordinates can be converted using: grid_coord = physical / lattice_spacing
 */
class InitialConditions {
public:
    /**
     * Domain Split Initialization (Test 3.2: Vacuum Energy Density)
     *
     * Creates two regions with different order parameters:
     * - Left domain (x < x_boundary): R = R_left (synchronized)
     * - Right domain (x >= x_boundary): R = R_right (desynchronized)
     *
     * Physics: Tests energy density scaling ρ ∝ R^n
     *
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param x_boundary Domain boundary position (grid units)
     * @param R_left Order parameter for left domain (typically 1.0)
     * @param R_right Order parameter for right domain (scan: 0.0 to 1.0)
     * @param transition_width Smooth transition width in grid units (default: 1.0)
     * @return R-field vector (size Nx*Ny) with domain split
     */
    static std::vector<float> domainSplit(
        int Nx, int Ny,
        float x_boundary,
        float R_left,
        float R_right,
        float transition_width = 1.0f);

    /**
     * Linear Defect Sheets (Test 3.1: Casimir-like Force)
     *
     * Creates parallel linear desynchronization regions (defect sheets):
     * - R(x) = R_min at x = x1 and x = x2 (defect lines)
     * - R(x) = R_max elsewhere (background)
     * - Smooth transition via tanh profile
     *
     * Physics: Measures attractive force F(d) ∝ 1/d^α between sheets
     *
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param x1 Position of first defect line (grid units)
     * @param x2 Position of second defect line (grid units)
     * @param R_min Order parameter at defect core (typically 0.0)
     * @param R_max Order parameter in background (typically 1.0)
     * @param defect_width Width of defect transition region (grid units)
     * @return R-field vector (size Nx*Ny) with linear defects
     */
    static std::vector<float> linearDefects(
        int Nx, int Ny,
        float x1, float x2,
        float R_min,
        float R_max,
        float defect_width);

    /**
     * Helper: Smooth step function using tanh
     * Interpolates from 0 to 1 over distance 2*width
     *
     * @param x Position
     * @param center Center of transition
     * @param width Half-width of transition
     * @return Value in [0, 1]
     */
    static float smoothStep(float x, float center, float width);

    /**
     * Vortex Core Initialization (Phase 2 Scenario 2.3: Relativistic Mass)
     *
     * Creates a radial R-field profile with core depression:
     * - R(r) = R_min + (R_max - R_min) * tanh(r / r_core)
     * - At r=0 (vortex center): R → R_min (incoherent)
     * - At r>>r_core: R → R_max (synchronized)
     *
     * Physics: Vortex topological singularity forces R→0 at core
     *
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param center_x Vortex center x-coordinate (grid units)
     * @param center_y Vortex center y-coordinate (grid units)
     * @param core_radius Vortex core radius (grid units)
     * @param R_min Order parameter at vortex center (typically 0.0)
     * @param R_max Order parameter far from core (typically 1.0)
     * @return R-field vector (size Nx*Ny) with vortex core profile
     */
    static std::vector<float> vortexCore(
        int Nx, int Ny,
        float center_x, float center_y,
        float core_radius,
        float R_min = 0.0f,
        float R_max = 1.0f);

    /**
     * Vortex Pair Initialization (Sprint 2: Multi-Defect Interactions)
     *
     * Creates vortex-antivortex pair with opposite winding numbers:
     * - R-field: Overlapping tanh profiles at each core
     * - Phase field: θ(r) = W₁·arg(z-z₁) + W₂·arg(z-z₂)
     *
     * Physics: Tests defect interaction, annihilation dynamics, EM topology
     *
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param x1 First vortex center x (grid units)
     * @param y1 First vortex center y (grid units)
     * @param W1 First vortex winding number (typically +1)
     * @param x2 Second vortex center x (grid units)
     * @param y2 Second vortex center y (grid units)
     * @param W2 Second vortex winding number (typically -1)
     * @param core_radius Vortex core radius (grid units)
     * @param R_min Order parameter at vortex centers (typically 0.0)
     * @param R_max Order parameter far from cores (typically 1.0)
     * @return R-field vector with overlapping vortex cores
     */
    static std::vector<float> vortexPair(
        int Nx, int Ny,
        float x1, float y1, int W1,
        float x2, float y2, int W2,
        float core_radius,
        float R_min = 0.0f,
        float R_max = 1.0f);

    /**
     * Helper: Convert physical coordinates to grid coordinates
     *
     * @param physical_coord Coordinate in Planck lengths
     * @param L_domain Domain size in Planck lengths
     * @param N_grid Grid size
     * @return Grid coordinate
     */
    static float physicalToGrid(float physical_coord, float L_domain, int N_grid);

    /**
     * Helper: Convert grid coordinates to physical coordinates
     *
     * @param grid_coord Grid coordinate
     * @param L_domain Domain size in Planck lengths
     * @param N_grid Grid size
     * @return Physical coordinate in Planck lengths
     */
    static float gridToPhysical(float grid_coord, float L_domain, int N_grid);
};

#endif // INITIAL_CONDITIONS_H
