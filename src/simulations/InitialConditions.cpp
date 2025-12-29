#include "simulations/InitialConditions.h"
#include <algorithm>

std::vector<float> InitialConditions::domainSplit(
    int Nx, int Ny,
    float x_boundary,
    float R_left,
    float R_right,
    float transition_width) {

    std::vector<float> R_field(Nx * Ny);

    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;
            float x = static_cast<float>(ix);

            // Smooth interpolation from R_left to R_right at boundary
            // Using smoothStep: 0 (far left) → 1 (far right)
            float t = smoothStep(x, x_boundary, transition_width);

            // Interpolate: R = R_left * (1-t) + R_right * t
            R_field[idx] = R_left * (1.0f - t) + R_right * t;
        }
    }

    return R_field;
}

std::vector<float> InitialConditions::linearDefects(
    int Nx, int Ny,
    float x1, float x2,
    float R_min,
    float R_max,
    float defect_width) {

    std::vector<float> R_field(Nx * Ny);

    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;
            float x = static_cast<float>(ix);

            // Distance from each defect line
            float d1 = std::abs(x - x1);
            float d2 = std::abs(x - x2);

            // Defect profile: R → R_min at defect, R → R_max far away
            // Use minimum distance to nearest defect
            float d_min = std::min(d1, d2);

            // Smooth profile: R(d) = R_min + (R_max - R_min) * tanh(d / defect_width)
            // At d=0: R ≈ R_min
            // At d >> defect_width: R → R_max
            float profile = std::tanh(d_min / defect_width);
            R_field[idx] = R_min + (R_max - R_min) * profile;
        }
    }

    return R_field;
}

std::vector<float> InitialConditions::vortexCore(
    int Nx, int Ny,
    float center_x, float center_y,
    float core_radius,
    float R_min,
    float R_max) {

    std::vector<float> R_field(Nx * Ny);

    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            int idx = iy * Nx + ix;

            // Distance from vortex center
            float dx = static_cast<float>(ix) - center_x;
            float dy = static_cast<float>(iy) - center_y;
            float r = std::sqrt(dx * dx + dy * dy);

            // Vortex core profile: R(r) = R_min + (R_max - R_min) * tanh(r / r_core)
            // At r=0: R → R_min (incoherent at vortex singularity)
            // At r >> r_core: R → R_max (synchronized background)
            float profile = std::tanh(r / core_radius);
            R_field[idx] = R_min + (R_max - R_min) * profile;
        }
    }

    return R_field;
}

float InitialConditions::smoothStep(float x, float center, float width) {
    if (width < 1e-6f) {
        // Sharp transition
        return (x >= center) ? 1.0f : 0.0f;
    }

    // Smooth transition using tanh
    // tanh((x - center) / width) ranges from -1 to +1
    // Map to [0, 1]: (tanh(...) + 1) / 2
    float arg = (x - center) / width;
    return 0.5f * (std::tanh(arg) + 1.0f);
}

float InitialConditions::physicalToGrid(float physical_coord, float L_domain, int N_grid) {
    // Lattice spacing: a = L_domain / N_grid
    float a = L_domain / static_cast<float>(N_grid);
    return physical_coord / a;
}

float InitialConditions::gridToPhysical(float grid_coord, float L_domain, int N_grid) {
    // Lattice spacing: a = L_domain / N_grid
    float a = L_domain / static_cast<float>(N_grid);
    return grid_coord * a;
}
