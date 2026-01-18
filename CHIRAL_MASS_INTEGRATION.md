# Chiral Mass Step Integration

## Implementation Complete

Successfully extended the existing `Dirac3D` class with chiral mass support.

### Files Modified
1. **include/Dirac3D.h** - Added method declaration (line 134)
2. **src/Dirac3D.cpp** - Added implementation (lines 318-365)

### Implementation Details

The new `applyChiralMassStep()` method:
- **Signature**: `void applyChiralMassStep(const std::vector<float>& R_field, const std::vector<float>& theta_field, float Delta, float dt)`
- **Location**: Immediately after the existing `applyMassStep()` method
- **Pattern**: Follows identical structure to existing mass step implementation

### Key Features
1. **Chiral Decomposition**: M = Δ·R·e^(iθγ^5) = m_S·I + i·m_P·γ^5
   - m_S = Δ·R·cos(θ) (scalar mass)
   - m_P = Δ·R·sin(θ) (pseudoscalar mass)

2. **Exact Match to GPU Shader**: Implementation matches dirac_rk4.comp lines 209-223

3. **Eigenvalue Decomposition**:
   - Upper components (0,1): γ^5 eigenvalue = +1
   - Lower components (2,3): γ^5 eigenvalue = -1

### Integration with Existing Step Method

The existing `step()` method uses Strang splitting:
```cpp
void Dirac3D::step(const std::vector<float>& mass_field, float dt) {
    applyKineticHalfStep(dt / 2.0f);
    applyMassStep(mass_field, dt);     // <-- Can be replaced with applyChiralMassStep
    applyKineticHalfStep(dt / 2.0f);
}
```

To use chiral mass instead of scalar mass, create a new overloaded `step()` method:
```cpp
void Dirac3D::stepChiral(const std::vector<float>& R_field,
                         const std::vector<float>& theta_field,
                         float Delta, float dt) {
    applyKineticHalfStep(dt / 2.0f);
    applyChiralMassStep(R_field, theta_field, Delta, dt);
    applyKineticHalfStep(dt / 2.0f);
}
```

### Validation Requirements
- Energy conservation: ΔE/E < 0.01%
- Time reversibility: < 1e-4 rad phase error
- Chiral symmetry preservation

### No Duplication
- Reused existing `applyMassStep()` pattern
- No new files created
- No duplicate matrix definitions
- Leverages existing FFT and kinetic step infrastructure