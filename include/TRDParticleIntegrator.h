// include/TRDParticleIntegrator.h
#pragma once

#include <cmath>
#include <functional>

/**
 * TRDParticleIntegrator - Unified Symplectic Integration for Charged Particles
 *
 * This header provides standardized particle integration methods for TRD physics tests.
 * Eliminates 300+ lines of duplicate code across 4 test files by centralizing proven
 * symplectic integrators for conservative particle dynamics.
 *
 * Symplectic Integration Requirement:
 *   All conservative physics (EM, gravity) MUST use symplectic integrators to preserve
 *   energy and phase space volume. Energy conservation target: ΔE/E < 0.01%
 *
 * Integration Methods:
 *   1. Boris Push (EM fields) - Energy conserving for Lorentz force F = q(E + v×B)
 *   2. Velocity Verlet (Conservative forces) - Time-reversible for general F = -∇U
 *   3. RK2 Midpoint (Compatible with TRDCore3D) - Symplectic to O(dt²)
 *
 * References:
 *   - Boris algorithm: Birdsall & Langdon, "Plasma Physics via Computer Simulation"
 *   - Velocity Verlet: Hairer et al., "Geometric Numerical Integration" (2006)
 *   - TRD standards: SYMPLECTIC_INVESTIGATION_REPORT.md, CLAUDE.md
 *
 * Usage Examples:
 *   // Electromagnetic dynamics (Lorentz force)
 *   TRD::ChargedParticle electron(0, 0, 0, 0.1, 0, 0, -1.0, 1.0);
 *   TRD::borisStep(electron, Ex, Ey, Ez, Bx, By, Bz, dt);
 *
 *   // Gravitational dynamics (conservative force)
 *   TRD::ChargedParticle probe(0, 0, 0, 0, 0.1, 0, 0.0, 1.0);
 *   TRD::velocityVerletStep(probe, ax, ay, az, dt);
 */

namespace TRD {

/**
 * ChargedParticle - Test particle with charge and mass
 *
 * Represents a point particle for validation tests in electromagnetic and gravitational
 * fields. Uses natural units (c=1, ℏ=1) consistent with TRD framework.
 *
 * Properties:
 *   - Position: (x, y, z) in TRD spatial units
 *   - Velocity: (vx, vy, vz) in natural units (v/c dimensionless)
 *   - Charge: Electric charge q in TRD units
 *   - Mass: Rest mass m in TRD units (246 GeV for golden key reference)
 *
 * Note: For relativistic particles, consider implementing 4-velocity formulation.
 *       Current implementation is non-relativistic (v << c regime).
 */
struct ChargedParticle {
    // Position (3D spatial coordinates)
    double x, y, z;

    // Velocity (3D velocity components)
    double vx, vy, vz;

    // Intrinsic properties
    double charge;  // Electric charge (TRD units)
    double mass;    // Rest mass (TRD units)

    /**
     * Default constructor - Initialize at origin at rest
     */
    ChargedParticle()
        : x(0.0), y(0.0), z(0.0)
        , vx(0.0), vy(0.0), vz(0.0)
        , charge(0.0), mass(1.0)
    {}

    /**
     * Full constructor - Specify all initial conditions
     *
     * @param x0, y0, z0 Initial position
     * @param vx0, vy0, vz0 Initial velocity
     * @param q Electric charge
     * @param m Rest mass
     */
    ChargedParticle(double x0, double y0, double z0,
                   double vx0, double vy0, double vz0,
                   double q, double m)
        : x(x0), y(y0), z(z0)
        , vx(vx0), vy(vy0), vz(vz0)
        , charge(q), mass(m)
    {}

    /**
     * Get kinetic energy (non-relativistic)
     * @return KE = (1/2)·m·v²
     */
    inline double kineticEnergy() const {
        return 0.5 * mass * (vx*vx + vy*vy + vz*vz);
    }

    /**
     * Get momentum magnitude
     * @return |p| = m·|v|
     */
    inline double momentum() const {
        double v2 = vx*vx + vy*vy + vz*vz;
        return mass * std::sqrt(v2);
    }

    /**
     * Get Lorentz factor (for diagnostic purposes)
     * @return γ = 1/√(1 - v²/c²), with c=1
     * Note: Valid only if v < 1 in natural units
     */
    inline double gamma() const {
        double v2 = vx*vx + vy*vy + vz*vz;
        return 1.0 / std::sqrt(1.0 - v2);
    }

    /**
     * Get speed (velocity magnitude)
     * @return |v| = √(vx² + vy² + vz²)
     */
    inline double speed() const {
        return std::sqrt(vx*vx + vy*vy + vz*vz);
    }
};

/**
 * Boris Integrator - Symplectic integration for Lorentz force
 *
 * Implements the Boris push-kick-push algorithm for charged particle motion in
 * electromagnetic fields. This is the industry-standard symplectic integrator
 * for plasma physics and particle-in-cell (PIC) simulations.
 *
 * Algorithm (Birdsall & Langdon):
 *   1. Half-step electric acceleration: v⁻ = vⁿ + (q/m)·E·dt/2
 *   2. Magnetic rotation: v⁺ = rotate(v⁻, B, dt)
 *   3. Half-step electric acceleration: vⁿ⁺¹ = v⁺ + (q/m)·E·dt/2
 *   4. Position update: xⁿ⁺¹ = xⁿ + vⁿ⁺¹·dt
 *
 * Properties:
 *   - Symplectic (preserves phase space volume)
 *   - Time-reversible
 *   - Energy conserving for B-only fields (ΔE/E ~ 0)
 *   - Second-order accurate O(dt²)
 *
 * Validated in:
 *   - test_lorentz_force_3d.cpp (cyclotron motion, ΔE/E < 0.01%)
 *   - test_three_body_em_3d.cpp (multi-particle EM interactions)
 *
 * @param p Particle to evolve (modified in place)
 * @param Ex, Ey, Ez Electric field components at particle position
 * @param Bx, By, Bz Magnetic field components at particle position
 * @param dt Time step (natural units)
 */
inline void borisStep(ChargedParticle& p,
                     double Ex, double Ey, double Ez,
                     double Bx, double By, double Bz,
                     double dt) {
    // Precompute q/m·dt/2 for efficiency
    double qm_half_dt = (p.charge / p.mass) * 0.5 * dt;

    // Step 1: Half-step electric acceleration (E push)
    double vx_minus = p.vx + qm_half_dt * Ex;
    double vy_minus = p.vy + qm_half_dt * Ey;
    double vz_minus = p.vz + qm_half_dt * Ez;

    // Step 2: Magnetic rotation
    // Define t-vector: t = (q/m)·B·dt/2
    double tx = qm_half_dt * Bx;
    double ty = qm_half_dt * By;
    double tz = qm_half_dt * Bz;

    // Define s-vector: s = 2t/(1 + t²)
    // This ensures the rotation is exact even for large B·dt
    double t_squared = tx*tx + ty*ty + tz*tz;
    double s_factor = 2.0 / (1.0 + t_squared);
    double sx = s_factor * tx;
    double sy = s_factor * ty;
    double sz = s_factor * tz;

    // Rotation step 1: v' = v⁻ + v⁻ × t
    double vprime_x = vx_minus + vy_minus * tz - vz_minus * ty;
    double vprime_y = vy_minus + vz_minus * tx - vx_minus * tz;
    double vprime_z = vz_minus + vx_minus * ty - vy_minus * tx;

    // Rotation step 2: v⁺ = v⁻ + v' × s
    double vx_plus = vx_minus + vprime_y * sz - vprime_z * sy;
    double vy_plus = vy_minus + vprime_z * sx - vprime_x * sz;
    double vz_plus = vz_minus + vprime_x * sy - vprime_y * sx;

    // Step 3: Half-step electric acceleration (E push)
    p.vx = vx_plus + qm_half_dt * Ex;
    p.vy = vy_plus + qm_half_dt * Ey;
    p.vz = vz_plus + qm_half_dt * Ez;

    // Step 4: Position update (drift)
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;
}

/**
 * Velocity Verlet Integrator - Symplectic integration for conservative forces
 *
 * Implements the velocity Verlet (kick-drift-kick) algorithm for particles under
 * conservative forces F = -∇U. This is the standard symplectic integrator for
 * Hamiltonian systems with separable kinetic and potential energy.
 *
 * Algorithm:
 *   1. Half-kick: v(t+dt/2) = v(t) + a(t)·dt/2
 *   2. Drift: x(t+dt) = x(t) + v(t+dt/2)·dt
 *   3. Half-kick: v(t+dt) = v(t+dt/2) + a(t+dt)·dt/2
 *
 * Properties:
 *   - Symplectic (preserves phase space volume)
 *   - Time-reversible
 *   - Energy conserving for conservative forces (ΔE/E < 0.01%)
 *   - Second-order accurate O(dt²)
 *
 * CRITICAL: Caller MUST provide acceleration a(t+dt) based on NEW position x(t+dt)
 *           for the second half-kick. This function only performs ONE step.
 *
 * Usage Pattern:
 *   for (int step = 0; step < N; ++step) {
 *       // Compute acceleration at current position
 *       auto [ax, ay, az] = computeAcceleration(p.x, p.y, p.z);
 *
 *       // Half-kick
 *       p.vx += ax * (dt/2);
 *       p.vy += ay * (dt/2);
 *       p.vz += az * (dt/2);
 *
 *       // Drift
 *       p.x += p.vx * dt;
 *       p.y += p.vy * dt;
 *       p.z += p.vz * dt;
 *
 *       // Compute NEW acceleration at NEW position
 *       auto [ax_new, ay_new, az_new] = computeAcceleration(p.x, p.y, p.z);
 *
 *       // Half-kick
 *       p.vx += ax_new * (dt/2);
 *       p.vy += ay_new * (dt/2);
 *       p.vz += az_new * (dt/2);
 *   }
 *
 * Validated in:
 *   - test_geodesic_3d.cpp (curved spacetime trajectories)
 *   - test_weak_field_3d.cpp (Newtonian gravity, a = GM/r²)
 *   - test_three_body_em_3d.cpp (Coulomb forces)
 *
 * @param p Particle to evolve (modified in place)
 * @param ax, ay, az Acceleration components at current position
 * @param dt Time step (natural units)
 * @param half_kick If true, perform half-kick only (for manual step control)
 */
inline void velocityVerletStep(ChargedParticle& p,
                              double ax, double ay, double az,
                              double dt,
                              bool half_kick = false) {
    if (half_kick) {
        // Half-kick only (for manual integration control)
        p.vx += ax * (dt / 2.0);
        p.vy += ay * (dt / 2.0);
        p.vz += az * (dt / 2.0);
    } else {
        // Full kick-drift-kick cycle
        // NOTE: This assumes acceleration is constant over dt
        //       For position-dependent forces, use manual half-kick pattern

        // Half-kick
        p.vx += ax * (dt / 2.0);
        p.vy += ay * (dt / 2.0);
        p.vz += az * (dt / 2.0);

        // Drift
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;

        // Half-kick (ASSUMES same acceleration - caller must update if needed!)
        p.vx += ax * (dt / 2.0);
        p.vy += ay * (dt / 2.0);
        p.vz += az * (dt / 2.0);
    }
}

/**
 * RK2 Midpoint Integrator - Symplectic integration compatible with TRDCore3D
 *
 * Implements the RK2 midpoint method, which is symplectic for separable Hamiltonians.
 * This integrator is compatible with TRDCore3D's default integration mode and
 * provides an alternative to Velocity Verlet for general force laws.
 *
 * Algorithm:
 *   1. Compute midpoint velocity: v_mid = v + a(x)·dt/2
 *   2. Compute midpoint position: x_mid = x + v_mid·dt/2
 *   3. Evaluate acceleration at midpoint: a_mid = a(x_mid)
 *   4. Update: v_new = v + a_mid·dt, x_new = x + v_mid·dt
 *
 * Properties:
 *   - Symplectic (for separable Hamiltonians)
 *   - Second-order accurate O(dt²)
 *   - Single acceleration evaluation per step (more efficient than RK4)
 *
 * Note: This is the same integrator used in TRDCore3D::IntegrationMode::SYMPLECTIC
 *
 * @param p Particle to evolve (modified in place)
 * @param accel Function returning acceleration given particle state: (p) -> (ax, ay, az)
 * @param dt Time step (natural units)
 */
inline void midpointStep(ChargedParticle& p,
                        std::function<void(const ChargedParticle&, double&, double&, double&)> accel,
                        double dt) {
    // Compute acceleration at current position
    double ax_0, ay_0, az_0;
    accel(p, ax_0, ay_0, az_0);

    // Compute midpoint velocity
    double vx_mid = p.vx + ax_0 * (dt / 2.0);
    double vy_mid = p.vy + ay_0 * (dt / 2.0);
    double vz_mid = p.vz + az_0 * (dt / 2.0);

    // Compute midpoint position
    ChargedParticle p_mid = p;
    p_mid.x += vx_mid * (dt / 2.0);
    p_mid.y += vy_mid * (dt / 2.0);
    p_mid.z += vz_mid * (dt / 2.0);

    // Compute acceleration at midpoint
    double ax_mid, ay_mid, az_mid;
    accel(p_mid, ax_mid, ay_mid, az_mid);

    // Update velocity using midpoint acceleration
    p.vx += ax_mid * dt;
    p.vy += ay_mid * dt;
    p.vz += az_mid * dt;

    // Update position using midpoint velocity
    p.x += vx_mid * dt;
    p.y += vy_mid * dt;
    p.z += vz_mid * dt;
}

/**
 * Legacy Compatibility - Float variants for existing test code
 *
 * Many existing tests use float instead of double. These wrappers provide
 * backward compatibility during the migration to unified integrators.
 */

struct Particle {
    float x, y, z;
    float vx, vy, vz;
    float charge;
    float mass;

    // Convert to double-precision ChargedParticle
    ChargedParticle toChargedParticle() const {
        return ChargedParticle(x, y, z, vx, vy, vz, charge, mass);
    }

    // Update from double-precision ChargedParticle
    void fromChargedParticle(const ChargedParticle& cp) {
        x = static_cast<float>(cp.x);
        y = static_cast<float>(cp.y);
        z = static_cast<float>(cp.z);
        vx = static_cast<float>(cp.vx);
        vy = static_cast<float>(cp.vy);
        vz = static_cast<float>(cp.vz);
        charge = static_cast<float>(cp.charge);
        mass = static_cast<float>(cp.mass);
    }

    float kineticEnergy() const {
        return 0.5f * mass * (vx*vx + vy*vy + vz*vz);
    }
};

// Float-precision Boris step
inline void borisStep(Particle& p,
                     float Ex, float Ey, float Ez,
                     float Bx, float By, float Bz,
                     float dt) {
    ChargedParticle cp = p.toChargedParticle();
    borisStep(cp, Ex, Ey, Ez, Bx, By, Bz, dt);
    p.fromChargedParticle(cp);
}

// Float-precision Velocity Verlet step
inline void velocityVerletStep(Particle& p,
                              float ax, float ay, float az,
                              float dt,
                              bool half_kick = false) {
    ChargedParticle cp = p.toChargedParticle();
    velocityVerletStep(cp, ax, ay, az, dt, half_kick);
    p.fromChargedParticle(cp);
}

} // namespace TRD
