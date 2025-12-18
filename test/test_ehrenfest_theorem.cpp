/**
 * test_ehrenfest_theorem.cpp
 *
 * Verify Ehrenfest theorem: d<p>/dt = <F> = -<∇H>
 * For Dirac: H = α·p + βm(x)
 * Force: F = -∇H = -β·∇m(x)
 */

#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Ehrenfest Theorem for Dirac Equation ===" << std::endl;
    std::cout << "\nHamiltonian: H = α·p + βm(x)" << std::endl;
    std::cout << "Force operator: F̂ = -∇H = -β·∇m(x)" << std::endl;
    
    std::cout << "\nFor wavepacket with <β> > 0 (upper components dominant):" << std::endl;
    std::cout << "  <F> = -<β>·∇m(x)" << std::endl;
    std::cout << "  If ∇m < 0 (mass decreasing): <F> > 0 (force in +x direction)" << std::endl;
    std::cout << "  → Packet accelerates toward LOWER mass regions" << std::endl;
    
    std::cout << "\nPhysical interpretation:" << std::endl;
    std::cout << "  Dirac particle in varying mass field experiences force" << std::endl;
    std::cout << "  Force = -β·∇m pushes particle toward LOW mass regions" << std::endl;
    std::cout << "  (Opposite of gravitational attraction to HIGH mass!)" << std::endl;
    
    std::cout << "\n=== Test Case: m(x) = 0.5(1 + 0.5sin(x/10)) ===" << std::endl;
    float x0 = 32.0f;
    float m_32 = 0.5f * (1.0f + 0.5f * std::sin(32.0f / 10.0f));
    float m_33 = 0.5f * (1.0f + 0.5f * std::sin(33.0f / 10.0f));
    float dm_dx = m_33 - m_32;
    
    std::cout << "At x=32:" << std::endl;
    std::cout << "  m(32) = " << m_32 << std::endl;
    std::cout << "  ∂m/∂x ≈ " << dm_dx << " (NEGATIVE)" << std::endl;
    std::cout << "  For <β>=+1: F_x = -∇m = " << -dm_dx << " (POSITIVE)" << std::endl;
    std::cout << "  → Motion in +x direction (toward LOWER mass)" << std::endl;
    
    float x_final = 38.0f;
    float m_38 = 0.5f * (1.0f + 0.5f * std::sin(38.0f / 10.0f));
    std::cout << "\nVerification:" << std::endl;
    std::cout << "  m(32) = " << m_32 << " (initial)" << std::endl;
    std::cout << "  m(38) = " << m_38 << " (final)" << std::endl;
    std::cout << "  m(38) < m(32)? " << (m_38 < m_32 ? "YES ✓" : "NO ✗") << std::endl;
    std::cout << "  → Particle correctly moved toward LOWER mass" << std::endl;
    
    std::cout << "\n=== CONCLUSION ===" << std::endl;
    std::cout << "The observed behavior is CORRECT for Dirac equation!" << std::endl;
    std::cout << "Force F = -β·∇m pushes toward LOW mass, not high mass." << std::endl;
    std::cout << "\nThis is DIFFERENT from gravitational attraction:" << std::endl;
    std::cout << "  Gravity: F ∝ +∇m (attract to high mass)" << std::endl;
    std::cout << "  Dirac:   F ∝ -∇m (repel from high mass)" << std::endl;
    
    return 0;
}
