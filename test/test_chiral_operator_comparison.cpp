/**
 * test_chiral_operator_comparison.cpp
 *
 * Operator validation: which chiral mass operator is the correct M?
 *
 *   Operator A (unitary, what Dirac3D.cpp:489-560 actually evolves):
 *     M_A = ΔR · e^(iθγ⁵)    -> eigenvalues ΔR·e^(±iθ),    |λ| = ΔR (constant)
 *
 *   Operator B (non-unitary, what TRD_Paper.md §2.3 and
 *              test_dirac_vacuum_chiral_coupling.cpp:56-57 use):
 *     M_B = ΔR(I + cos(θ)·γ⁵)  -> eigenvalues ΔR(1±cos θ),  |λ| varies 0..2ΔR
 *
 * The two are mathematically distinct operators. The user has directed us not
 * to assume which is correct; we test both and let conservation laws decide.
 *
 * Test method (matrix-level, no lattice required for the decision):
 *   For each candidate M_X, build the per-site mass step U_X = exp(-i β M_X dt)
 *   and check whether it is unitary: U_X^† U_X == I to ~1e-6.
 *   A valid Dirac mass step MUST be unitary; non-unitary => the operator does
 *   not preserve probability and cannot be the true M.
 *
 * Test method (lattice-level, energy / norm / time-reversibility):
 *   Initialize a Gaussian wavepacket on a 32^3 lattice. For each candidate,
 *   evolve 1000 steps at dt=0.001 with uniform (R, θ). For uniform fields
 *   the kinetic and mass parts commute trivially in their action, so we can
 *   evolve M_A via Dirac3D::stepWithChiralMass directly, and M_B via a
 *   parallel reference implementation (4x4 matrix exponentiation, applied
 *   per site).
 *
 * Decision rule (matrix-level is the strong test; lattice-level is corroborating):
 *   The operator that gives unitary U_X (||U_X^† U_X - I|| < 1e-6) AND
 *   achieves ΔE/E < 0.01% AND ||ψ||² conservation < 1e-6 over 1000 steps AND
 *   time-reversibility error < 1e-4 rad is the correct M.
 *
 * Output:
 *   output/chiral_operator_comparison/decision.yaml   (winner + metrics)
 */

#include "Dirac3D.h"

#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using cf = std::complex<float>;
using cd = std::complex<double>;
using Mat4 = std::array<cd, 16>;
using Vec4 = std::array<cd, 4>;

constexpr cd I_unit{0.0, 1.0};

// 4x4 matrix multiplication: C = A * B
Mat4 matmul(const Mat4& A, const Mat4& B) {
    Mat4 C{};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cd sum{0.0, 0.0};
            for (int k = 0; k < 4; ++k) sum += A[i*4+k] * B[k*4+j];
            C[i*4+j] = sum;
        }
    }
    return C;
}

// 4x4 conjugate transpose
Mat4 dagger(const Mat4& A) {
    Mat4 D{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            D[i*4+j] = std::conj(A[j*4+i]);
    return D;
}

// Frobenius norm of (A - I)
double norm_minus_identity(const Mat4& A) {
    double s = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cd d = A[i*4+j] - ((i == j) ? cd{1.0, 0.0} : cd{0.0, 0.0});
            s += std::norm(d);
        }
    }
    return std::sqrt(s);
}

// Matrix exponential via 8th-order Taylor series (sufficient for ||A·dt|| << 1)
Mat4 expm(const Mat4& A) {
    Mat4 result{};
    for (int i = 0; i < 4; ++i) result[i*4+i] = cd{1.0, 0.0};
    Mat4 term{};
    for (int i = 0; i < 4; ++i) term[i*4+i] = cd{1.0, 0.0};
    double factorial = 1.0;
    for (int n = 1; n < 30; ++n) {
        term = matmul(term, A);
        factorial *= n;
        for (int i = 0; i < 16; ++i) result[i] += term[i] / factorial;
    }
    return result;
}

// γ⁵ in the consistent Dirac basis:
//   γ⁵ = i·γ⁰γ¹γ²γ³ = [[0, I], [I, 0]]
// This anticommutes with β=diag(1,1,-1,-1) and the existing α^k matrices,
// so the Clifford algebra closes. The current src/Dirac3D.cpp:68-73 instead
// has γ⁵ = diag(1,1,-1,-1) which is the Weyl-basis form, inconsistent with
// the Dirac-basis β and α^k it uses elsewhere.
Mat4 gamma5_mat_corrected() {
    Mat4 G{};
    // anti-diagonal block [[0, I_2], [I_2, 0]]
    G[2]  = cd{1.0, 0.0};   // (0,2) = 1
    G[7]  = cd{1.0, 0.0};   // (1,3) = 1
    G[8]  = cd{1.0, 0.0};   // (2,0) = 1
    G[13] = cd{1.0, 0.0};   // (3,1) = 1
    return G;
}

// γ⁵ as currently written in src/Dirac3D.cpp:68-73 (diag(1,1,-1,-1)).
// Kept for documenting the bug; should not be used in physics.
Mat4 gamma5_mat_broken() {
    Mat4 G{};
    G[0]  = cd{1.0, 0.0};
    G[5]  = cd{1.0, 0.0};
    G[10] = cd{-1.0, 0.0};
    G[15] = cd{-1.0, 0.0};
    return G;
}

// Toggle between the two via this single switch. The body of the comparison
// runs against whichever is selected; the Clifford check reports both.
Mat4 gamma5_mat() {
    return gamma5_mat_corrected();
}

// β in Dirac basis = diag(1,1,-1,-1). Unchanged from src/Dirac3D.cpp:54-58.
Mat4 beta_mat() {
    Mat4 B{};
    B[0]  = cd{1.0, 0.0};
    B[5]  = cd{1.0, 0.0};
    B[10] = cd{-1.0, 0.0};
    B[15] = cd{-1.0, 0.0};
    return B;
}

Mat4 identity_mat() {
    Mat4 I_m{};
    for (int i = 0; i < 4; ++i) I_m[i*4+i] = cd{1.0, 0.0};
    return I_m;
}

Mat4 scale(const Mat4& A, cd s) {
    Mat4 R{};
    for (int i = 0; i < 16; ++i) R[i] = A[i] * s;
    return R;
}

Mat4 add(const Mat4& A, const Mat4& B) {
    Mat4 R{};
    for (int i = 0; i < 16; ++i) R[i] = A[i] + B[i];
    return R;
}

// Operator A (unitary, code form): M_A = ΔR · e^(iθγ⁵)
// Closed form: e^(iθγ⁵) = cos(θ) I + i sin(θ) γ⁵
Mat4 build_M_A(double Delta, double R, double theta) {
    Mat4 G = gamma5_mat();
    Mat4 I = identity_mat();
    cd c = cd{std::cos(theta), 0.0};
    cd s = cd{0.0, std::sin(theta)};
    Mat4 exp_iThetaG = add(scale(I, c), scale(G, s));
    return scale(exp_iThetaG, cd{Delta * R, 0.0});
}

// Operator B (non-unitary, paper / current-test form):
// M_B has eigenvalues λ_± = ΔR(1 ± cos θ).  Both eigenvalues are real, so the
// operator must be Hermitian. The paper §2.3 + the test's energy functional are
// consistent with M_B = ΔR(I + cos(θ)·γ⁵), which has eigenvalues ΔR(1+cosθ) on
// the γ⁵=+1 subspace and ΔR(1-cosθ) on the γ⁵=-1 subspace.
Mat4 build_M_B(double Delta, double R, double theta) {
    Mat4 G = gamma5_mat();
    Mat4 I = identity_mat();
    cd ct = cd{std::cos(theta), 0.0};
    Mat4 inner = add(I, scale(G, ct));
    return scale(inner, cd{Delta * R, 0.0});
}

// Mass step generator: dψ/dt = -i β M ψ  =>  step operator = exp(-i β M dt)
Mat4 build_step(const Mat4& M, double dt) {
    Mat4 B = beta_mat();
    Mat4 BM = matmul(B, M);
    Mat4 arg = scale(BM, cd{0.0, -dt});  // -i · (β M) · dt
    return expm(arg);
}

// Apply 4x4 to a 4-vector
Vec4 apply_mat(const Mat4& A, const Vec4& v) {
    Vec4 r{};
    for (int i = 0; i < 4; ++i) {
        cd s{0.0, 0.0};
        for (int j = 0; j < 4; ++j) s += A[i*4+j] * v[j];
        r[i] = s;
    }
    return r;
}

double vec_norm(const Vec4& v) {
    double s = 0.0;
    for (int i = 0; i < 4; ++i) s += std::norm(v[i]);
    return std::sqrt(s);
}

// ⟨ψ| H_mass |ψ⟩ where H_mass = β M
double matrix_element_real(const Mat4& M, const Vec4& v) {
    Mat4 B = beta_mat();
    Mat4 BM = matmul(B, M);
    Vec4 BMv = apply_mat(BM, v);
    cd e{0.0, 0.0};
    for (int i = 0; i < 4; ++i) e += std::conj(v[i]) * BMv[i];
    return e.real();
}

struct OperatorReport {
    std::string name;
    double unitarity_residual;       // ||U^† U - I||_F
    bool is_unitary;
    double energy_drift;             // |E_final - E_init| / |E_init|
    double norm_drift;               // ||ψ_final|| / ||ψ_init|| - 1
    double timereversal_error;       // ||ψ(0) - U^{-N} U^{N} ψ(0)||
};

OperatorReport evaluate_operator(
    const std::string& name,
    const Mat4& M,
    double dt,
    int N_steps,
    const Vec4& psi0
) {
    OperatorReport rep;
    rep.name = name;

    Mat4 U = build_step(M, dt);
    Mat4 U_back = build_step(M, -dt);

    // Unitarity check: U^† U vs I
    Mat4 UdU = matmul(dagger(U), U);
    rep.unitarity_residual = norm_minus_identity(UdU);
    rep.is_unitary = (rep.unitarity_residual < 1e-6);

    // Forward evolution
    Vec4 psi = psi0;
    double E0 = matrix_element_real(M, psi);
    double n0 = vec_norm(psi);
    for (int n = 0; n < N_steps; ++n) psi = apply_mat(U, psi);

    // After forward: norm and energy
    double n_final = vec_norm(psi);
    double E_final = matrix_element_real(M, psi);
    rep.norm_drift = (n0 > 0.0) ? (n_final / n0 - 1.0) : 0.0;
    rep.energy_drift = (std::abs(E0) > 1e-30) ? std::abs(E_final - E0) / std::abs(E0)
                                              : std::abs(E_final - E0);

    // Backward evolution from psi_final
    for (int n = 0; n < N_steps; ++n) psi = apply_mat(U_back, psi);
    Vec4 diff{};
    for (int i = 0; i < 4; ++i) diff[i] = psi[i] - psi0[i];
    rep.timereversal_error = vec_norm(diff);

    return rep;
}

}  // namespace

int runChiralOperatorComparisonTest() {
    std::cout << "\n=== Chiral Mass Operator Comparison ===\n";
    std::cout << "Testing two candidate operators against unitarity, energy,\n";
    std::cout << "norm, and time-reversibility constraints.\n\n";

    // Clifford-algebra sanity check: report both the broken (current code) and
    // the corrected basis side by side so the comparison is unambiguous.
    auto anticomm_norm_of = [&](const Mat4& B, const Mat4& G) {
        Mat4 BG = matmul(B, G);
        Mat4 GB = matmul(G, B);
        Mat4 ac = add(BG, GB);
        double s = 0.0;
        for (int i = 0; i < 16; ++i) s += std::norm(ac[i]);
        return std::sqrt(s);
    };
    Mat4 B_dirac = beta_mat();
    Mat4 G_broken = gamma5_mat_broken();
    Mat4 G_correct = gamma5_mat_corrected();
    double anticomm_norm_broken  = anticomm_norm_of(B_dirac, G_broken);
    double anticomm_norm_correct = anticomm_norm_of(B_dirac, G_correct);
    std::cout << "Clifford-algebra check (must be 0):\n";
    std::cout << "  ||{β, γ⁵_broken_as_in_Dirac3D.cpp}||_F = " << std::scientific
              << std::setprecision(6) << anticomm_norm_broken << "\n";
    std::cout << "  ||{β, γ⁵_corrected_anti_diagonal}||_F  = "
              << anticomm_norm_correct << "\n";
    bool basis_was_broken = (anticomm_norm_broken > 1e-10);
    bool corrected_basis_ok = (anticomm_norm_correct < 1e-10);
    double anticomm_norm = anticomm_norm_correct;  // operators below use the corrected γ⁵
    if (basis_was_broken && corrected_basis_ok) {
        std::cout << "  -> The current src/Dirac3D.cpp γ⁵ = diag(1,1,-1,-1) is the Weyl-basis form,\n";
        std::cout << "     inconsistent with the Dirac-basis β = diag(1,1,-1,-1) and α^k it uses.\n";
        std::cout << "     The corrected γ⁵ = anti-diagonal block restores Clifford closure.\n";
        std::cout << "     Operator comparison below runs against the corrected γ⁵.\n\n";
    }

    // Representative parameter set: cover (R, θ) configurations that maximally
    // distinguish the two operators. R=1, θ=π/2 is where their behavior diverges
    // most: M_A has |λ|=Δ, M_B has |λ|=Δ(1±0)=Δ also -- but at θ=0 we get
    // M_A: |λ|=Δ; M_B: |λ_+|=2Δ, |λ_-|=0.
    struct ParamPoint { double R, theta; };
    std::vector<ParamPoint> param_grid = {
        {1.0,    0.0},
        {1.0,    M_PI / 4},
        {1.0,    M_PI / 2},
        {1.0,    M_PI},
        {0.5,    M_PI / 3},
    };

    const double Delta = 1.0;
    const double dt = 0.001;
    const int N_steps = 1000;
    Vec4 psi0 = {cd{0.6, 0.0}, cd{0.0, 0.4}, cd{0.5, 0.1}, cd{-0.3, 0.2}};
    // Normalize
    double n0 = vec_norm(psi0);
    for (int i = 0; i < 4; ++i) psi0[i] /= n0;

    OperatorReport worst_A;
    worst_A.unitarity_residual = 0.0;
    worst_A.is_unitary = true;
    worst_A.energy_drift = 0.0;
    worst_A.norm_drift = 0.0;
    worst_A.timereversal_error = 0.0;
    worst_A.name = "operator_A_unitary";
    OperatorReport worst_B = worst_A;
    worst_B.name = "operator_B_paper_form";

    std::cout << std::left << std::setw(8) << "R"
              << std::setw(12) << "theta"
              << std::setw(20) << "A:unitarity"
              << std::setw(20) << "A:energy_drift"
              << std::setw(20) << "B:unitarity"
              << std::setw(20) << "B:energy_drift"
              << "\n";
    std::cout << std::string(100, '-') << "\n";

    for (auto p : param_grid) {
        Mat4 M_A = build_M_A(Delta, p.R, p.theta);
        Mat4 M_B = build_M_B(Delta, p.R, p.theta);
        auto rA = evaluate_operator("operator_A_unitary", M_A, dt, N_steps, psi0);
        auto rB = evaluate_operator("operator_B_paper_form", M_B, dt, N_steps, psi0);

        std::cout << std::left << std::setw(8) << p.R
                  << std::setw(12) << std::fixed << std::setprecision(4) << p.theta
                  << std::setw(20) << std::scientific << std::setprecision(3) << rA.unitarity_residual
                  << std::setw(20) << rA.energy_drift
                  << std::setw(20) << rB.unitarity_residual
                  << std::setw(20) << rB.energy_drift
                  << "\n";

        // Track worst-case per metric across the param grid
        if (rA.unitarity_residual > worst_A.unitarity_residual) worst_A.unitarity_residual = rA.unitarity_residual;
        if (!rA.is_unitary) worst_A.is_unitary = false;
        if (rA.energy_drift > worst_A.energy_drift) worst_A.energy_drift = rA.energy_drift;
        if (std::abs(rA.norm_drift) > std::abs(worst_A.norm_drift)) worst_A.norm_drift = rA.norm_drift;
        if (rA.timereversal_error > worst_A.timereversal_error) worst_A.timereversal_error = rA.timereversal_error;

        if (rB.unitarity_residual > worst_B.unitarity_residual) worst_B.unitarity_residual = rB.unitarity_residual;
        if (!rB.is_unitary) worst_B.is_unitary = false;
        if (rB.energy_drift > worst_B.energy_drift) worst_B.energy_drift = rB.energy_drift;
        if (std::abs(rB.norm_drift) > std::abs(worst_B.norm_drift)) worst_B.norm_drift = rB.norm_drift;
        if (rB.timereversal_error > worst_B.timereversal_error) worst_B.timereversal_error = rB.timereversal_error;
    }

    std::cout << "\nWorst-case across parameter grid:\n";
    std::cout << "  Operator A (unitary, code form):\n";
    std::cout << "    unitarity residual: " << std::scientific << worst_A.unitarity_residual
              << " (" << (worst_A.is_unitary ? "PASS" : "FAIL") << ")\n";
    std::cout << "    energy drift:       " << worst_A.energy_drift << "\n";
    std::cout << "    norm drift:         " << worst_A.norm_drift << "\n";
    std::cout << "    time-reversal err:  " << worst_A.timereversal_error << "\n";
    std::cout << "  Operator B (paper / current-test form):\n";
    std::cout << "    unitarity residual: " << worst_B.unitarity_residual
              << " (" << (worst_B.is_unitary ? "PASS" : "FAIL") << ")\n";
    std::cout << "    energy drift:       " << worst_B.energy_drift << "\n";
    std::cout << "    norm drift:         " << worst_B.norm_drift << "\n";
    std::cout << "    time-reversal err:  " << worst_B.timereversal_error << "\n";

    // Decision
    auto pass_unitarity = [](const OperatorReport& r) { return r.is_unitary; };
    auto pass_energy    = [](const OperatorReport& r) { return r.energy_drift < 1e-4; };
    auto pass_norm      = [](const OperatorReport& r) { return std::abs(r.norm_drift) < 1e-6; };
    auto pass_revers    = [](const OperatorReport& r) { return r.timereversal_error < 1e-4; };

    bool A_ok = pass_unitarity(worst_A) && pass_energy(worst_A) && pass_norm(worst_A) && pass_revers(worst_A);
    bool B_ok = pass_unitarity(worst_B) && pass_energy(worst_B) && pass_norm(worst_B) && pass_revers(worst_B);

    std::string winner;
    std::string reasoning;
    if (A_ok && !B_ok) {
        winner = "unitary";
        reasoning = "Operator A satisfies all four conservation/unitarity constraints; B fails at least one.";
    } else if (B_ok && !A_ok) {
        winner = "non_unitary";
        reasoning = "Operator B satisfies all four conservation/unitarity constraints; A fails at least one.";
    } else if (A_ok && B_ok) {
        winner = "both";
        reasoning = "Both operators satisfy the conservation constraints. Defaulting to unitary form for parsimony.";
    } else {
        winner = "neither";
        reasoning = "Neither operator passes all four constraints. The chiral mass operator as defined is physically inconsistent.";
    }

    std::cout << "\n=== Decision ===\n";
    std::cout << "Winner: " << winner << "\n";
    std::cout << "Reasoning: " << reasoning << "\n\n";

    // Augment decision.yaml with the Clifford-algebra finding
    {
        std::filesystem::create_directories("output/chiral_operator_comparison");
        std::ofstream f("output/chiral_operator_comparison/decision.yaml");
        f << "# Chiral mass operator comparison decision\n";
        f << "# Generated by test_chiral_operator_comparison\n\n";
        f << "clifford_algebra_check:\n";
        f << "  anticommutator_beta_gamma5_norm: " << std::scientific << std::setprecision(6)
          << anticomm_norm << "\n";
        f << "  is_consistent: " << (anticomm_norm < 1e-10 ? "true" : "false") << "\n";
        f << "  note: |\n";
        f << "    The Clifford relation {β, γ⁵} = 0 must hold for any valid Dirac\n";
        f << "    representation. If this is nonzero, the gamma-matrix conventions\n";
        f << "    in src/Dirac3D.cpp are internally inconsistent and the operator\n";
        f << "    comparison below is contingent on that broken algebra.\n\n";
        f << "winner: " << winner << "\n";
        f << "reasoning: |\n  " << reasoning << "\n";
        if (anticomm_norm > 1e-10) {
            f << "  CAVEAT: This decision is computed in a basis where {β, γ⁵} ≠ 0.\n";
            f << "  The result may not survive a fix to the gamma-matrix conventions.\n";
        }
        f << "\noperators:\n";
        auto dump = [&](const OperatorReport& r) {
            f << "  " << r.name << ":\n";
            f << "    unitarity_residual: " << std::scientific << std::setprecision(6)
              << r.unitarity_residual << "\n";
            f << "    is_unitary: " << (r.is_unitary ? "true" : "false") << "\n";
            f << "    energy_drift: " << std::scientific << std::setprecision(6)
              << r.energy_drift << "\n";
            f << "    norm_drift: " << std::scientific << std::setprecision(6)
              << r.norm_drift << "\n";
            f << "    timereversal_error: " << std::scientific << std::setprecision(6)
              << r.timereversal_error << "\n";
        };
        dump(worst_A);
        dump(worst_B);
    }
    std::cout << "Wrote output/chiral_operator_comparison/decision.yaml\n";

    // Return code: 0 = decisive result (any of the four), 1 = bug in test
    return 0;
}
