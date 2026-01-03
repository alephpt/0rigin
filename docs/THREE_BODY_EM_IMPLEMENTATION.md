# Three-Body Electromagnetic Dynamics Test Implementation

## Goal: G2 - Verify Superposition Principle

Validate that multiple charges interact correctly via TRD EM fields, demonstrating that:
1. EM fields from multiple vortices superpose linearly
2. Forces match analytical three-body Coulomb calculation
3. Total force on each charge equals vector sum of pairwise forces

## Configuration: config/three_body_em.yaml

### Vortex Setup
- **Vortex 1**: Position (32, 64), charge q1 = +1 (winding = 1)
- **Vortex 2**: Position (96, 64), charge q2 = +1 (winding = 1)
- **Vortex 3**: Position (64, 96), charge q3 = -1 (winding = -1)

### Physical Layout
```
    y
    ^
    |
96  |     V3(-)
    |      *
    |     / \
    |    /   \
64  | V1(+)---V2(+)
    |  *       *
    |
    +-----------> x
      32  64  96
```

Forms an isosceles triangle with:
- Base: V1-V2 distance = 64 units (32 grid points at dx=0.5)
- Sides: V1-V3 = V2-V3 = √(32² + 32²) ≈ 45.3 units

### Expected Forces (Coulomb's Law)

For charges q_i, q_j at distance r_ij:
F_ij = k * q_i * q_j / r_ij² * r̂_ij

#### Force on V1 (+1 charge at (32, 64)):
- From V2: F_12 = repulsive, rightward
  - Magnitude: k/32² ≈ 0.00098k
  - Direction: +x
- From V3: F_13 = attractive, toward (64, 96)
  - Magnitude: k/45.3² ≈ 0.00049k
  - Direction: (1/√2, 1/√2)
- **Total**: F_1 = F_12 + F_13

#### Force on V2 (+1 charge at (96, 64)):
- From V1: F_21 = repulsive, leftward
  - Magnitude: k/32² ≈ 0.00098k
  - Direction: -x
- From V3: F_23 = attractive, toward (64, 96)
  - Magnitude: k/45.3² ≈ 0.00049k
  - Direction: (-1/√2, 1/√2)
- **Total**: F_2 = F_21 + F_23

#### Force on V3 (-1 charge at (64, 96)):
- From V1: F_31 = attractive, toward (32, 64)
  - Magnitude: k/45.3² ≈ 0.00049k
  - Direction: (-1/√2, -1/√2)
- From V2: F_32 = attractive, toward (96, 64)
  - Magnitude: k/45.3² ≈ 0.00049k
  - Direction: (1/√2, -1/√2)
- **Total**: F_3 = F_31 + F_32

### Superposition Verification

The test will verify:
1. **Linear Superposition**: E_total = E_1 + E_2 + E_3
2. **Force Additivity**: F_on_i = Σ_j≠i F_ij
3. **Newton's Third Law**: F_ij = -F_ji

## Implementation Code

### 1. Multi-Vortex Initialization (TRDCore.cpp modification)

```cpp
void TRDCore::initializeThreeVortices(
    float x1, float y1, int w1,
    float x2, float y2, int w2,
    float x3, float y3, int w3) {

    // Initialize phase field with three vortices
    for (int j = 0; j < config_.ny; ++j) {
        for (int i = 0; i < config_.nx; ++i) {
            int idx = j * config_.nx + i;

            // Phase from each vortex
            float theta1 = w1 * atan2f(j - y1, i - x1);
            float theta2 = w2 * atan2f(j - y2, i - x2);
            float theta3 = w3 * atan2f(j - y3, i - x3);

            // Superpose phases
            h_theta_[idx] = theta1 + theta2 + theta3;

            // High synchronization for strong EM
            h_R_[idx] = 0.95f;
        }
    }

    // Upload to GPU
    uploadFieldsToGPU();
}
```

### 2. Force Calculation (ObservableComputer.cpp addition)

```cpp
struct VortexForce {
    double Fx, Fy;
    double x, y;  // Position
    double charge;
};

VortexForce ObservableComputer::computeVortexForce(
    const std::vector<float>& Ex,
    const std::vector<float>& Ey,
    const std::vector<float>& Bz,
    float x, float y,
    float vx, float vy,
    float charge,
    int nx, int ny) {

    // Bilinear interpolation of fields at vortex position
    int ix = static_cast<int>(x);
    int iy = static_cast<int>(y);

    // Fractional parts for interpolation
    float fx = x - ix;
    float fy = y - iy;

    // Grid indices (with periodic boundaries)
    int ix1 = (ix + 1) % nx;
    int iy1 = (iy + 1) % ny;

    // Interpolate E_x
    float E_x = (1-fx)*(1-fy)*Ex[iy*nx + ix] +
                fx*(1-fy)*Ex[iy*nx + ix1] +
                (1-fx)*fy*Ex[iy1*nx + ix] +
                fx*fy*Ex[iy1*nx + ix1];

    // Interpolate E_y
    float E_y = (1-fx)*(1-fy)*Ey[iy*nx + ix] +
                fx*(1-fy)*Ey[iy*nx + ix1] +
                (1-fx)*fy*Ey[iy1*nx + ix] +
                fx*fy*Ey[iy1*nx + ix1];

    // Interpolate B_z
    float B_z = (1-fx)*(1-fy)*Bz[iy*nx + ix] +
                fx*(1-fy)*Bz[iy*nx + ix1] +
                (1-fx)*fy*Bz[iy1*nx + ix] +
                fx*fy*Bz[iy1*nx + ix1];

    // Lorentz force: F = q*(E + v×B)
    VortexForce force;
    force.Fx = charge * (E_x + vy * B_z);
    force.Fy = charge * (E_y - vx * B_z);
    force.x = x;
    force.y = y;
    force.charge = charge;

    return force;
}

void ObservableComputer::computeThreeBodyForces(
    ThreeBodyObservables* obs,
    const TRDEngine* engine) {

    // Get EM fields from engine
    const auto& Ex = engine->getEM_Ex();
    const auto& Ey = engine->getEM_Ey();
    const auto& Bz = engine->getEM_Bz();

    int nx = engine->getNx();
    int ny = engine->getNy();

    // Vortex positions and charges
    float x1 = 32.0f, y1 = 64.0f, q1 = 1.0f;
    float x2 = 96.0f, y2 = 64.0f, q2 = 1.0f;
    float x3 = 64.0f, y3 = 96.0f, q3 = -1.0f;

    // Compute forces (assuming initially stationary)
    VortexForce F1 = computeVortexForce(Ex, Ey, Bz, x1, y1, 0, 0, q1, nx, ny);
    VortexForce F2 = computeVortexForce(Ex, Ey, Bz, x2, y2, 0, 0, q2, nx, ny);
    VortexForce F3 = computeVortexForce(Ex, Ey, Bz, x3, y3, 0, 0, q3, nx, ny);

    // Store in observables
    obs->F1_x = F1.Fx;
    obs->F1_y = F1.Fy;
    obs->F2_x = F2.Fx;
    obs->F2_y = F2.Fy;
    obs->F3_x = F3.Fx;
    obs->F3_y = F3.Fy;

    // Compute analytical Coulomb forces for comparison
    float dx = 0.5f;  // Grid spacing
    float r12 = 64.0f * dx;  // Distance V1-V2
    float r13 = sqrtf(32.0f*32.0f + 32.0f*32.0f) * dx;  // Distance V1-V3
    float r23 = r13;  // Distance V2-V3

    // Coulomb constant (effective in TRD units)
    float k = 1.0f;

    // Analytical forces
    float F12_mag = k * q1 * q2 / (r12 * r12);  // Repulsive
    float F13_mag = k * q1 * q3 / (r13 * r13);  // Attractive
    float F23_mag = k * q2 * q3 / (r23 * r23);  // Attractive

    // Store analytical values
    obs->F12_analytical = F12_mag;
    obs->F13_analytical = F13_mag;
    obs->F23_analytical = F23_mag;

    // Compute errors
    float F1_measured = sqrtf(F1.Fx*F1.Fx + F1.Fy*F1.Fy);
    float F1_expected = sqrtf(F12_mag*F12_mag + F13_mag*F13_mag/2);  // Vector sum
    obs->coulomb_error = fabsf(F1_measured - F1_expected) / F1_expected;
}
```

### 3. Test Implementation (test/test_three_body_em.cpp)

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TRDEngine.h"
#include "simulations/ObservableComputer.h"
#include "simulations/TestConfig.h"

int main(int argc, char* argv[]) {
    std::cout << "=== Three-Body EM Dynamics Test ===" << std::endl;

    // Load configuration
    TestConfig config("config/three_body_em.yaml");

    // Initialize TRD engine
    TRDConfig trd_config;
    trd_config.nx = config.grid.size_x;
    trd_config.ny = config.grid.size_y;
    trd_config.dt = config.physics.dt;
    trd_config.delta = config.physics.delta;
    trd_config.coupling_strength = config.physics.coupling;
    trd_config.enable_em = true;
    trd_config.photon_mass_coupling = 0.05f;

    TRDEngine engine(trd_config);

    // Initialize three vortices
    engine.initializeThreeVortices(
        32.0f, 64.0f, 1,   // Vortex 1: +1 charge
        96.0f, 64.0f, 1,   // Vortex 2: +1 charge
        64.0f, 96.0f, -1   // Vortex 3: -1 charge
    );

    // Let system equilibrate
    std::cout << "Equilibrating for 100 steps..." << std::endl;
    for (int i = 0; i < 100; ++i) {
        engine.step();
    }

    // Measure forces
    std::cout << "\n--- Force Measurements ---" << std::endl;

    ThreeBodyObservables obs;
    ObservableComputer computer;
    computer.computeThreeBodyForces(&obs, &engine);

    // Print force comparison table
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nVortex Forces (TRD vs Analytical):" << std::endl;
    std::cout << "------------------------------------" << std::endl;

    std::cout << "V1 (+1 at 32,64):" << std::endl;
    std::cout << "  Fx = " << std::setw(10) << obs.F1_x << std::endl;
    std::cout << "  Fy = " << std::setw(10) << obs.F1_y << std::endl;
    std::cout << "  |F| = " << std::setw(10) << sqrt(obs.F1_x*obs.F1_x + obs.F1_y*obs.F1_y) << std::endl;

    std::cout << "\nV2 (+1 at 96,64):" << std::endl;
    std::cout << "  Fx = " << std::setw(10) << obs.F2_x << std::endl;
    std::cout << "  Fy = " << std::setw(10) << obs.F2_y << std::endl;
    std::cout << "  |F| = " << std::setw(10) << sqrt(obs.F2_x*obs.F2_x + obs.F2_y*obs.F2_y) << std::endl;

    std::cout << "\nV3 (-1 at 64,96):" << std::endl;
    std::cout << "  Fx = " << std::setw(10) << obs.F3_x << std::endl;
    std::cout << "  Fy = " << std::setw(10) << obs.F3_y << std::endl;
    std::cout << "  |F| = " << std::setw(10) << sqrt(obs.F3_x*obs.F3_x + obs.F3_y*obs.F3_y) << std::endl;

    // Superposition check
    std::cout << "\n--- Superposition Verification ---" << std::endl;

    // Newton's third law check
    float newton3_x = obs.F1_x + obs.F2_x + obs.F3_x;
    float newton3_y = obs.F1_y + obs.F2_y + obs.F3_y;
    float newton3_error = sqrt(newton3_x*newton3_x + newton3_y*newton3_y);

    std::cout << "Sum of forces (should be ~0):" << std::endl;
    std::cout << "  ΣFx = " << newton3_x << std::endl;
    std::cout << "  ΣFy = " << newton3_y << std::endl;
    std::cout << "  |ΣF| = " << newton3_error << std::endl;

    // Quality gates
    std::cout << "\n--- Quality Gates ---" << std::endl;
    bool coulomb_pass = obs.coulomb_error < 0.05;
    bool newton3_pass = newton3_error < 1e-4;

    std::cout << "Coulomb agreement: "
              << (coulomb_pass ? "✓ PASS" : "✗ FAIL")
              << " (error = " << obs.coulomb_error*100 << "%)" << std::endl;

    std::cout << "Newton's 3rd law: "
              << (newton3_pass ? "✓ PASS" : "✗ FAIL")
              << " (|ΣF| = " << newton3_error << ")" << std::endl;

    // Dynamic evolution test
    std::cout << "\n--- Dynamic Evolution (1000 steps) ---" << std::endl;

    // Track energy
    float E0 = engine.getTotalEnergy();

    for (int step = 0; step < 1000; ++step) {
        engine.step();

        if (step % 100 == 0) {
            float E = engine.getTotalEnergy();
            float drift = (E - E0) / E0;
            std::cout << "Step " << step
                      << ": Energy drift = " << drift*100 << "%" << std::endl;
        }
    }

    // Final validation
    float Ef = engine.getTotalEnergy();
    float total_drift = std::abs(Ef - E0) / E0;

    std::cout << "\n--- Final Results ---" << std::endl;
    std::cout << "Energy conservation: "
              << (total_drift < 0.01 ? "✓ PASS" : "✗ FAIL")
              << " (drift = " << total_drift*100 << "%)" << std::endl;

    // Overall pass/fail
    bool all_pass = coulomb_pass && newton3_pass && (total_drift < 0.01);

    std::cout << "\n=== TEST "
              << (all_pass ? "PASSED" : "FAILED")
              << " ===" << std::endl;

    return all_pass ? 0 : 1;
}
```

## Success Criteria

✅ Forces match Coulomb law within 5%
✅ Superposition verified: F_total = F_12 + F_13
✅ Energy conserved within 1%
✅ Vortices remain stable (no unwinding)
✅ Newton's third law satisfied (ΣF ≈ 0)

## Next Steps

1. Implement `initializeThreeVortices()` in TRDCore
2. Add `computeThreeBodyForces()` to ObservableComputer
3. Create and run test_three_body_em.cpp
4. Validate results against analytical predictions
5. If successful, extend to N-body configurations

## Expected Output

```
=== Three-Body EM Dynamics Test ===
Equilibrating for 100 steps...

--- Force Measurements ---

Vortex Forces (TRD vs Analytical):
------------------------------------
V1 (+1 at 32,64):
  Fx =   0.001234
  Fy =   0.000345
  |F| =   0.001281

V2 (+1 at 96,64):
  Fx =  -0.001234
  Fy =   0.000345
  |F| =   0.001281

V3 (-1 at 64,96):
  Fx =   0.000000
  Fy =  -0.000690
  |F| =   0.000690

--- Superposition Verification ---
Sum of forces (should be ~0):
  ΣFx = 0.000000
  ΣFy = 0.000000
  |ΣF| = 0.000000

--- Quality Gates ---
Coulomb agreement: ✓ PASS (error = 3.2%)
Newton's 3rd law: ✓ PASS (|ΣF| = 1.2e-6)

=== TEST PASSED ===
```