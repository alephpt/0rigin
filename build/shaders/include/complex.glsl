/**
 * @file complex.glsl
 * @brief Complex number arithmetic for MSFT Vulkan compute shaders
 *
 * HAZARD A MITIGATION: Proper complex arithmetic in GLSL
 * ========================================================
 *
 * GLSL does NOT have native complex number support.
 * We use vec2(real, imag) representation with explicit operations.
 *
 * CRITICAL for MSFT/Dirac:
 * - Dirac spinors are complex: Ψ = (ψ₁, ψ₂, ψ₃, ψ₄) where ψᵢ ∈ ℂ
 * - Gamma matrices mix components: γ^μ : ℂ⁴ → ℂ⁴
 * - Mass operator includes complex phase: e^(iθγ⁵)
 *
 * Physics Requirements:
 * - Preserve unitarity: |ψ|² conserved (mitigated in Hazard B)
 * - Lorentz invariance: γ^μ anticommutation preserved
 * - Chiral symmetry: γ⁵ = iγ⁰γ¹γ²γ³ structure
 */

#ifndef MSFT_COMPLEX_GLSL
#define MSFT_COMPLEX_GLSL

// Type aliases for clarity
#define complex vec2
#define Complex vec2

// Constants
const complex C_ZERO = vec2(0.0, 0.0);
const complex C_ONE  = vec2(1.0, 0.0);
const complex C_I    = vec2(0.0, 1.0);

/**
 * @brief Complex conjugate: z* = a - ib
 */
complex conj(complex z) {
    return vec2(z.x, -z.y);
}

/**
 * @brief Complex norm squared: |z|² = a² + b²
 */
float norm2(complex z) {
    return dot(z, z);  // Optimized using dot product
}

/**
 * @brief Complex modulus: |z| = √(a² + b²)
 */
float cabs(complex z) {
    return length(z);  // Optimized using GLSL builtin
}

/**
 * @brief Complex addition: (a + ib) + (c + id) = (a+c) + i(b+d)
 */
complex cadd(complex z1, complex z2) {
    return z1 + z2;  // Component-wise
}

/**
 * @brief Complex subtraction: (a + ib) - (c + id) = (a-c) + i(b-d)
 */
complex csub(complex z1, complex z2) {
    return z1 - z2;  // Component-wise
}

/**
 * @brief Complex multiplication: (a+ib)(c+id) = (ac-bd) + i(ad+bc)
 *
 * CRITICAL: This is where GLSL naive approaches fail!
 * Cannot use: z1 * z2 (gives component-wise, NOT complex multiplication)
 */
complex cmul(complex z1, complex z2) {
    float a = z1.x, b = z1.y;
    float c = z2.x, d = z2.y;
    return vec2(a*c - b*d, a*d + b*c);
}

/**
 * @brief Scalar multiplication: r * (a + ib) = ra + i(rb)
 */
complex smul(float r, complex z) {
    return r * z;  // Component-wise OK for scalar
}

/**
 * @brief Complex division: z1/z2 = (a+ib)/(c+id)
 *
 * Formula: (a+ib)/(c+id) = [(ac+bd) + i(bc-ad)]/(c²+d²)
 */
complex cdiv(complex z1, complex z2) {
    float a = z1.x, b = z1.y;
    float c = z2.x, d = z2.y;
    float denom = c*c + d*d;

    // Avoid division by zero
    if (denom < 1e-20) {
        return C_ZERO;
    }

    return vec2((a*c + b*d) / denom, (b*c - a*d) / denom);
}

/**
 * @brief Complex exponential: e^(iθ) = cos(θ) + i·sin(θ)
 *
 * CRITICAL for MSFT mass operator: e^(iθγ⁵)
 */
complex cexp_i(float theta) {
    return vec2(cos(theta), sin(theta));
}

/**
 * @brief General complex exponential: e^z = e^a(cos b + i sin b)
 */
complex cexp(complex z) {
    float exp_a = exp(z.x);
    return vec2(exp_a * cos(z.y), exp_a * sin(z.y));
}

/**
 * @brief Complex phase (argument): arg(a + ib) = atan2(b, a)
 */
float carg(complex z) {
    return atan(z.y, z.x);
}

/**
 * @brief Complex power (integer exponent): z^n
 *
 * Uses binary exponentiation for efficiency
 * Non-recursive implementation to avoid GLSL recursion issues
 */
complex cpow_int(complex z, int n) {
    if (n == 0) return C_ONE;
    if (n == 1) return z;

    // Handle negative exponents by computing positive power then inverting
    bool invert = false;
    if (n < 0) {
        invert = true;
        n = -n;
    }

    complex result = C_ONE;
    complex base = z;

    while (n > 0) {
        if ((n & 1) == 1) {
            result = cmul(result, base);
        }
        base = cmul(base, base);
        n >>= 1;
    }

    if (invert) {
        result = cdiv(C_ONE, result);
    }

    return result;
}

/**
 * @brief Complex square root: √z
 *
 * Formula: √(a+ib) = ±(√r·cos(θ/2) + i·√r·sin(θ/2))
 * where r = |z|, θ = arg(z)
 */
complex csqrt(complex z) {
    float r = cabs(z);
    float theta = carg(z);
    float sqrt_r = sqrt(r);
    return vec2(sqrt_r * cos(theta/2.0), sqrt_r * sin(theta/2.0));
}

/**
 * @brief Natural logarithm: ln(z) = ln|z| + i·arg(z)
 */
complex clog(complex z) {
    return vec2(log(cabs(z)), carg(z));
}

/**
 * @brief Complex power: z1^z2
 *
 * Formula: z1^z2 = exp(z2 · ln(z1))
 */
complex cpow(complex z1, complex z2) {
    return cexp(cmul(z2, clog(z1)));
}

/**
 * @brief Normalize complex number: z/|z|
 */
complex cnormalize(complex z) {
    float mag = cabs(z);
    if (mag < 1e-20) return C_ZERO;
    return z / mag;
}

// ============================================================================
// MSFT-SPECIFIC OPERATIONS
// ============================================================================

/**
 * @brief MSFT mass operator phase factor: e^(iθγ⁵)
 *
 * For a single spinor component:
 * e^(iθγ⁵) = cos(θ)·I + i·sin(θ)·γ⁵
 *
 * When γ⁵ eigenvalue is ±1:
 * e^(iθγ⁵)|±⟩ = e^(±iθ)|±⟩
 */
complex MSFT_chiral_phase(float theta, float gamma5_eigenvalue) {
    // gamma5_eigenvalue = ±1 for chiral eigenstates
    float effective_theta = theta * gamma5_eigenvalue;
    return cexp_i(effective_theta);
}

/**
 * @brief Dirac gamma matrix element multiplication
 *
 * γ^μ has complex entries in general representation.
 * This computes: (γ^μ)_ij · ψ_j where both are complex.
 */
complex gamma_mul_spinor(complex gamma_element, complex spinor_component) {
    return cmul(gamma_element, spinor_component);
}

/**
 * @brief Inner product of two complex 4-spinors
 *
 * ⟨ψ|φ⟩ = Σᵢ ψᵢ* · φᵢ
 */
complex spinor_inner_product(complex[4] psi, complex[4] phi) {
    complex result = C_ZERO;
    for (int i = 0; i < 4; i++) {
        result = cadd(result, cmul(conj(psi[i]), phi[i]));
    }
    return result;
}

/**
 * @brief Spinor density: ρ = |ψ|² = Σᵢ |ψᵢ|²
 *
 * Returns REAL number (imaginary part must be zero by construction)
 */
float spinor_density(complex[4] psi) {
    float rho = 0.0;
    for (int i = 0; i < 4; i++) {
        rho += norm2(psi[i]);
    }
    return rho;
}

#endif // MSFT_COMPLEX_GLSL
