/**
 * @file precision.glsl
 * @brief Floating point precision management for MSFT Vulkan compute shaders
 *
 * HAZARD B MITIGATION: Float32 vs Float64 precision
 * ===================================================
 *
 * PROBLEM: Dirac equation evolution must preserve unitarity
 * - Unitary evolution: d|ψ|²/dt = 0
 * - FP32 accumulation error → norm drift → physics violation
 * - Long time evolution (t >> 1000 steps) amplifies error
 *
 * STRATEGY (priority order):
 * 1. Check for FP64 support (GL_ARB_gpu_shader_fp64)
 * 2. If available: use double for critical accumulators
 * 3. If not: mitigate FP32 error via:
 *    - Kahan summation for accumulation
 *    - RK4 with adaptive step size
 *    - Periodic norm conservation projection
 *    - Compensated arithmetic for critical operations
 *
 * CRITICAL OPERATIONS REQUIRING HIGH PRECISION:
 * - Spinor norm: Σᵢ |ψᵢ|² (accumulation of 4 terms)
 * - RK4 k-vector accumulation: ψ + dt/6·(k₁+2k₂+2k₃+k₄)
 * - Hamiltonian action: H·ψ (matrix-vector products)
 * - Phase evolution: e^(-iHt)ψ (exponential error growth)
 */

#ifndef MSFT_PRECISION_GLSL
#define MSFT_PRECISION_GLSL

// ============================================================================
// PRECISION DETECTION
// ============================================================================

// Check for FP64 support at compile time
#ifdef GL_ARB_gpu_shader_fp64
    #define HAS_FP64 1
    #extension GL_ARB_gpu_shader_fp64 : enable
#else
    #define HAS_FP64 0
#endif

// Precision configuration (set via push constants or specialization constants)
#ifndef MSFT_USE_FP64
    #define MSFT_USE_FP64 HAS_FP64
#endif

// ============================================================================
// PRECISION TYPES
// ============================================================================

#if MSFT_USE_FP64 && HAS_FP64
    // Use FP64 for critical operations
    #define highp_float double
    #define highp_vec2  dvec2
    #define highp_vec3  dvec3
    #define highp_vec4  dvec4

    #define FLOAT_EPS 1e-15
#else
    // Fallback to FP32 with compensated arithmetic
    #define highp_float float
    #define highp_vec2  vec2
    #define highp_vec3  vec3
    #define highp_vec4  vec4

    #define FLOAT_EPS 1e-7
#endif

// ============================================================================
// KAHAN SUMMATION (COMPENSATED SUMMATION)
// ============================================================================

/**
 * @brief Kahan summation state
 *
 * Reduces FP32 accumulation error from O(nε) to O(ε)
 * where n = number of terms, ε = machine epsilon
 *
 * Critical for:
 * - Spinor norm calculation
 * - RK4 k-vector accumulation
 */
struct KahanSum {
    highp_float sum;           // Running sum
    highp_float compensation;  // Lost low-order bits
};

/**
 * @brief Initialize Kahan sum
 */
KahanSum kahan_init() {
    KahanSum ks;
    ks.sum = 0.0;
    ks.compensation = 0.0;
    return ks;
}

/**
 * @brief Add value to Kahan sum
 *
 * Formula:
 *   y = value - compensation
 *   t = sum + y
 *   compensation = (t - sum) - y  // Capture lost bits
 *   sum = t
 */
void kahan_add(inout KahanSum ks, highp_float value) {
    highp_float y = value - ks.compensation;
    highp_float t = ks.sum + y;
    ks.compensation = (t - ks.sum) - y;
    ks.sum = t;
}

/**
 * @brief Get final Kahan sum value
 */
highp_float kahan_result(KahanSum ks) {
    return ks.sum;
}

// ============================================================================
// NORM CONSERVATION
// ============================================================================

/**
 * @brief Conservative spinor norm computation using Kahan summation
 *
 * Computes: ||ψ||² = Σᵢ |ψᵢ|² with minimal error
 */
highp_float compute_spinor_norm2_conservative(vec2[4] spinor) {
    KahanSum sum = kahan_init();

    for (int i = 0; i < 4; i++) {
        highp_float term = highp_float(dot(spinor[i], spinor[i]));
        kahan_add(sum, term);
    }

    return kahan_result(sum);
}

/**
 * @brief Project spinor to unit norm (norm conservation)
 *
 * ψ → ψ / ||ψ||
 *
 * CRITICAL: Must be called periodically during evolution
 * to prevent FP32 drift from violating probability conservation
 *
 * Recommended: Every 10-100 RK4 steps
 */
void normalize_spinor_conservative(inout vec2[4] spinor) {
    highp_float norm2 = compute_spinor_norm2_conservative(spinor);

    // Avoid division by zero
    if (norm2 < FLOAT_EPS) {
        return;
    }

    float norm = float(sqrt(norm2));
    for (int i = 0; i < 4; i++) {
        spinor[i] /= norm;
    }
}

// ============================================================================
// RK4 COMPENSATED ARITHMETIC
// ============================================================================

/**
 * @brief RK4 final step with Kahan summation
 *
 * ψ(t+dt) = ψ(t) + dt/6 · (k₁ + 2k₂ + 2k₃ + k₄)
 *
 * Uses Kahan summation for weighted k-vector accumulation
 * to minimize FP32 error
 */
vec2 rk4_combine_conservative(
    vec2 psi_component,
    vec2 k1, vec2 k2, vec2 k3, vec2 k4,
    highp_float dt
) {
    // Real part accumulation
    KahanSum sum_real = kahan_init();
    kahan_add(sum_real, highp_float(k1.x));
    kahan_add(sum_real, 2.0 * highp_float(k2.x));
    kahan_add(sum_real, 2.0 * highp_float(k3.x));
    kahan_add(sum_real, highp_float(k4.x));

    // Imaginary part accumulation
    KahanSum sum_imag = kahan_init();
    kahan_add(sum_imag, highp_float(k1.y));
    kahan_add(sum_imag, 2.0 * highp_float(k2.y));
    kahan_add(sum_imag, 2.0 * highp_float(k3.y));
    kahan_add(sum_imag, highp_float(k4.y));

    // Apply time step
    vec2 increment = vec2(
        float(kahan_result(sum_real) * dt / 6.0),
        float(kahan_result(sum_imag) * dt / 6.0)
    );

    return psi_component + increment;
}

// ============================================================================
// ADAPTIVE TIME STEPPING
// ============================================================================

/**
 * @brief Estimate local truncation error for adaptive RK4
 *
 * Uses embedded Runge-Kutta method (RK4 vs RK5)
 * to estimate error and adapt time step
 *
 * Returns: Estimated error magnitude
 */
float estimate_rk4_error(
    vec2 psi_component,
    vec2 k1, vec2 k2, vec2 k3, vec2 k4,
    vec2 k5  // 5th-order estimate
) {
    // RK4 solution
    vec2 y_rk4 = psi_component + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;

    // RK5 solution (if available)
    vec2 y_rk5 = psi_component + (k1 + 4.0*k2 + k3 + 4.0*k4 + k5) / 10.0;

    // Error estimate
    return length(y_rk5 - y_rk4);
}

/**
 * @brief Adapt time step based on error tolerance
 *
 * @param dt_current Current time step
 * @param error Estimated truncation error
 * @param tolerance Desired error tolerance (e.g., 1e-6)
 * @return Adjusted time step
 */
highp_float adapt_time_step(
    highp_float dt_current,
    float error,
    float tolerance
) {
    const float safety_factor = 0.9;  // Conservative scaling
    const float max_scale = 2.0;      // Don't increase too fast
    const float min_scale = 0.3;      // Don't decrease too fast

    if (error < FLOAT_EPS) {
        // Error negligible, allow growth
        return dt_current * max_scale;
    }

    // PI controller for smooth adaptation
    float scale = safety_factor * pow(tolerance / error, 0.2);
    scale = clamp(scale, min_scale, max_scale);

    return dt_current * highp_float(scale);
}

// ============================================================================
// ENERGY CONSERVATION DIAGNOSTICS
// ============================================================================

/**
 * @brief Compute Dirac Hamiltonian energy: E = ⟨ψ|H|ψ⟩
 *
 * For free Dirac: E = ⟨ψ|(-i·α·∇ + β·m)|ψ⟩
 * With MSFT mass: m = m(x,t) = Δ·R(x,t)·e^(iθγ⁵)
 *
 * Should be conserved in unitary evolution.
 * Drift indicates precision issues.
 */
highp_float compute_dirac_energy_conservative(
    vec2[4] psi,
    vec2[4] H_psi  // Hamiltonian action on psi
) {
    KahanSum energy = kahan_init();

    for (int i = 0; i < 4; i++) {
        // ⟨ψᵢ|H|ψᵢ⟩ = Re(ψᵢ* · (Hψ)ᵢ)
        highp_float real_part = highp_float(
            psi[i].x * H_psi[i].x + psi[i].y * H_psi[i].y
        );
        kahan_add(energy, real_part);
    }

    return kahan_result(energy);
}

/**
 * @brief Monitor unitarity violation
 *
 * Computes: ||H†|| - ||H||
 * Should be zero for Hermitian Hamiltonian.
 *
 * Non-zero value indicates numerical error accumulation.
 */
float check_hermiticity_violation(
    vec2[4] psi,
    vec2[4] H_psi,
    vec2[4] H_dag_psi  // H† acting on psi
) {
    highp_float norm_H = compute_spinor_norm2_conservative(H_psi);
    highp_float norm_Hdag = compute_spinor_norm2_conservative(H_dag_psi);

    return float(abs(norm_Hdag - norm_H) / (norm_H + FLOAT_EPS));
}

// ============================================================================
// MSFT-SPECIFIC PRECISION MANAGEMENT
// ============================================================================

/**
 * @brief Precision budget for MSFT simulation
 *
 * Determines when to apply conservation corrections
 * based on accumulated error
 */
struct PrecisionBudget {
    uint step_count;              // Total steps elapsed
    float max_norm_drift;         // Max |1 - ||ψ||²|
    float max_energy_drift;       // Max |ΔE/E₀|
    bool needs_normalization;     // Trigger norm projection
    bool needs_step_reduction;    // Trigger dt reduction
};

/**
 * @brief Update precision budget after RK4 step
 */
void update_precision_budget(
    inout PrecisionBudget budget,
    float norm2,
    highp_float energy,
    highp_float energy_initial
) {
    const float norm_tolerance = 1e-4;      // 0.01% norm drift
    const float energy_tolerance = 1e-3;    // 0.1% energy drift

    budget.step_count++;

    // Check norm conservation
    float norm_drift = abs(1.0 - norm2);
    budget.max_norm_drift = max(budget.max_norm_drift, norm_drift);

    if (norm_drift > norm_tolerance) {
        budget.needs_normalization = true;
    }

    // Check energy conservation
    float energy_drift = float(abs(energy - energy_initial) / (abs(energy_initial) + FLOAT_EPS));
    budget.max_energy_drift = max(budget.max_energy_drift, energy_drift);

    if (energy_drift > energy_tolerance) {
        budget.needs_step_reduction = true;
    }
}

#endif // MSFT_PRECISION_GLSL
