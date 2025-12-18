
===== PHASE 2.2 VALIDATION REPORT =====
Scenario: Traveling Wave Surfing
Date: 2025-12-18 14:46:49

[Test Configuration]
Grid: 128 x 64
Wave vector: kx = 0.300, ky = 0.000
Kuramoto coupling: K = 3.00
Damping: γ = 0.020
Particle coupling: λ = 1.00
Timestep: dt = 0.0200
Total time: T = 9.98

[Wave Properties]
Analytical wave velocity: v_w = 0.825688 grid units/time
Wavelength: λ = 20.94 grid units

[Particle Dynamics]
Initial position: x_0 = 32.02
Final position: x_f = 37.41
Net displacement: Δx = 5.39
Average velocity in locked region: <v_x> = 0.848605
Expected wave velocity: v_w = 0.825688
Velocity match: 2.78%

[Validation Results]

1. Velocity Locking: ✗ FAIL
   Average |v_p - v_w|/v_w = 0.4156 (threshold: 0.20)
   Maximum |v_p - v_w|/v_w = 0.6822

2. Average Velocity Match: ✓ PASS
   |<v_p> - v_w|/v_w = 0.0278 (threshold: 0.10)

3. Lock Duration: ✗ FAIL
   Maximum sustained lock = 2.28 time units (threshold: 3.0)

[Conservation Laws]
Norm error: max |Ψ|² - 1 = 1.077315e-04
Energy drift: ΔE/E₀ = 1.779229e-02

[Overall Result]
Phase 2.2 Traveling Wave Surfing: ✗✗✗ FAIL ✗✗✗

[Interpretation]

The particle does NOT exhibit sustained wave surfing behavior.

Issues identified:
- Velocity locking too weak (avg 0.4156 > 0.20)
- Lock duration too short (2.28 < 3.0)

Recommendations:
- Increase Kuramoto coupling K for stronger wave stability
- Increase particle coupling λ for stronger wave-particle interaction
- Reduce timestep dt for better numerical accuracy
- Verify phase gradient creates stable traveling wave

===== END OF REPORT =====
