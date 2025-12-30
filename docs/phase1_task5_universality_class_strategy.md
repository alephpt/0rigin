# Finite-Size Scaling Strategy for SMFT Universality Class Identification

## Executive Summary

This document outlines a comprehensive finite-size scaling (FSS) strategy to definitively identify the universality class of the SMFT phase transition. Current measurements show β = 0.099 ± 0.004, statistically different from 2D Ising (β = 0.125), necessitating a complete FSS analysis to extract multiple critical exponents and classify the system.

## 1. Theoretical Framework

### 1.1 Finite-Size Scaling Formalism

Near critical points, observables exhibit universal scaling behavior controlled by the correlation length ξ:

**Core FSS Ansatz:**
```
A_L(σ) = L^(ζ/ν) f[(σ - σ_c)L^(1/ν)]
```

where:
- `A_L`: Observable at system size L
- `σ`: Noise strength (control parameter)
- `σ_c`: Critical noise strength
- `ζ`: Observable-specific critical exponent
- `ν`: Correlation length exponent
- `f`: Universal scaling function

### 1.2 Critical Exponent Relations

**Independent exponents:** Only two exponents are independent (typically ν and η)

**Scaling relations:**
- Fisher: `γ = ν(2 - η)`
- Rushbrooke: `α + 2β + γ = 2`
- Josephson: `νd = 2 - α` (hyperscaling)
- Widom: `γ = β(δ - 1)`

### 1.3 Key Observables

1. **Order Parameter:** `⟨R⟩_L ~ L^(-β/ν) g[(σ - σ_c)L^(1/ν)]`
2. **Susceptibility:** `χ_L = L^2[⟨R²⟩ - ⟨R⟩²] ~ L^(γ/ν) h[(σ - σ_c)L^(1/ν)]`
3. **Binder Cumulant:** `U_L = 1 - ⟨R⁴⟩/(3⟨R²⟩²)`
4. **Correlation Function:** `G(r) ~ r^(-(d-2+η))` at criticality
5. **Correlation Length:** `ξ_L ~ L` at σ_c

## 2. Test Matrix Specification

### 2.1 Grid Sizes

**Primary set (5 sizes):**
- L = 32 (1024 sites)
- L = 64 (4096 sites)
- L = 128 (16384 sites)
- L = 256 (65536 sites)
- L = 512 (262144 sites)

**Extended set (optional, for verification):**
- L = 16 (quick tests)
- L = 1024 (if computationally feasible)

### 2.2 Noise Scan Parameters

**Coarse scan (Phase 1):**
- Range: σ ∈ [0.5, 1.0]
- Points: 21 equally spaced
- Purpose: Locate critical region

**Fine scan (Phase 2):**
- Range: σ ∈ [σ_c - 0.1, σ_c + 0.1]
- Points: 41 equally spaced
- Purpose: Precise exponent extraction

**Total simulations:** 5 sizes × 41 points = 205 production runs

### 2.3 Computational Cost Estimate

**Per simulation:**
- Equilibration: 5000 steps
- Measurement: 10000 steps
- Time: ~5 min (L=32) to ~2 hours (L=512)

**Total runtime:**
- Coarse scan: ~50 CPU hours
- Fine scan: ~200 CPU hours
- **Total: ~250 CPU hours** (parallelizable to ~10 wall hours on cluster)

## 3. Analysis Methods

### 3.1 Extract β (Order Parameter Exponent)

**Method 1: Power-law fit at σ_c**
```
⟨R⟩_L(σ_c) ~ L^(-β/ν)
```
- Plot log(⟨R⟩) vs log(L) at critical point
- Slope gives -β/ν

**Method 2: Scaling collapse**
```
⟨R⟩_L * L^(β/ν) vs (σ - σ_c)L^(1/ν)
```
- Optimize β/ν for best data collapse
- Use χ² minimization

### 3.2 Extract ν (Correlation Length Exponent)

**Method 1: Binder cumulant crossing**
- Plot U_L(σ) for different L
- Curves cross at σ_c (size-independent)
- Width of crossing region scales as L^(-1/ν)

**Method 2: Susceptibility peak scaling**
```
χ_max ~ L^(γ/ν)
σ(χ_max) - σ_c ~ L^(-1/ν)
```

**Method 3: Derivative peak**
```
|d⟨R⟩/dσ|_max ~ L^((1+β)/ν)
```

### 3.3 Extract η (Anomalous Dimension)

**Method 1: Correlation function at criticality**
```
G(r, σ_c) ~ r^(-(d-2+η))
```
- Measure spatial correlations ⟨R(0)R(r)⟩
- Power-law fit gives η

**Method 2: Susceptibility scaling**
```
χ(σ_c) ~ L^(2-η)
```
- Use Fisher relation: γ = ν(2-η)

### 3.4 Extract γ (Susceptibility Exponent)

**Direct measurement:**
```
χ_max ~ L^(γ/ν)
```
- Plot log(χ_max) vs log(L)
- Slope gives γ/ν

### 3.5 Data Collapse Validation

**Procedure:**
1. Fix σ_c from Binder cumulant crossing
2. Optimize ν and β simultaneously:
   - Plot ⟨R⟩_L * L^(β/ν) vs (σ - σ_c)L^(1/ν)
   - Minimize scatter between curves
3. Validate with susceptibility:
   - Plot χ_L * L^(-γ/ν) vs (σ - σ_c)L^(1/ν)
   - Should collapse with γ = ν(2-η)

**Quality metrics:**
- χ² per degree of freedom < 2
- Residual scatter < 5%
- Visual inspection of collapse quality

## 4. Universality Class Decision Tree

### 4.1 Known 2D Universality Classes

| Class | β | ν | η | γ | Notes |
|-------|---|---|---|---|-------|
| **2D Ising** | 0.125 (exact) | 1.0 (exact) | 0.25 (exact) | 1.75 | Z₂ symmetry |
| **2D XY** | - | 0.67 | ~0.04 | ~1.3 | BKT transition |
| **2D 3-state Potts** | 1/9 | 5/6 | 4/15 | 13/9 | Z₃ symmetry |
| **2D 4-state Potts** | 1/12 | 2/3 | 1/4 | 7/6 | First-order |
| **Mean Field** | 0.5 | 0.5 | 0 | 1.0 | d ≥ 4 or long-range |

### 4.2 Classification Algorithm

```
1. Measure all exponents (β, ν, η, γ)
2. Compare with known classes:

   IF |β - 0.125| < 0.01 AND |ν - 1.0| < 0.05:
      → 2D Ising (confirm with η = 0.25)

   ELIF |ν - 0.67| < 0.05 AND η < 0.1:
      → 2D XY/BKT (check for logarithmic corrections)

   ELIF |β - 1/9| < 0.01 AND |ν - 5/6| < 0.05:
      → 2D 3-state Potts

   ELIF β > 0.3:
      → Check for mean-field behavior

   ELSE:
      → Novel universality class
      → Report all exponents with uncertainties
      → Search for theoretical models with similar exponents
```

### 4.3 Novel Class Indicators

**If SMFT shows novel universality:**
- β ≈ 0.099 (already measured, differs from Ising)
- Expected deviations:
  - ν might differ due to non-local interactions
  - η could indicate modified scaling dimensions
  - γ/ν ratio reveals susceptibility anomalies

**Possible explanations:**
1. **Crossover behavior:** System interpolates between classes
2. **Logarithmic corrections:** Near upper critical dimension
3. **Long-range interactions:** Modified scaling from coupling
4. **Non-equilibrium effects:** Dynamic universality class

## 5. Success Criteria

### 5.1 Minimal Requirements
- [ ] σ_c determined to ±0.001 precision
- [ ] β measured to ±0.005 uncertainty
- [ ] ν measured to ±0.01 uncertainty
- [ ] Data collapse quality: χ²/dof < 2
- [ ] Consistent exponents via multiple methods

### 5.2 Optimal Goals
- [ ] All exponents (β, ν, η, γ, α, δ) determined
- [ ] Scaling relations verified to < 5% accuracy
- [ ] Correction-to-scaling exponents extracted
- [ ] Universal amplitude ratios computed
- [ ] Definitive classification or proof of novel class

## 6. Publication-Quality Figure Specifications

### Figure 1: Order Parameter Scaling
- **Panel A:** ⟨R⟩ vs σ for all L (raw data)
- **Panel B:** Log-log plot of ⟨R⟩(σ_c) vs L
- **Panel C:** Data collapse of ⟨R⟩L^(β/ν) vs (σ-σ_c)L^(1/ν)
- **Insets:** Residuals from fits

### Figure 2: Binder Cumulant Analysis
- **Main:** U_L vs σ showing crossing point
- **Inset 1:** Zoom of crossing region
- **Inset 2:** σ_cross vs 1/L extrapolation

### Figure 3: Critical Exponent Summary
- **Panel A:** Susceptibility scaling χ_max vs L
- **Panel B:** Correlation function G(r) at criticality
- **Panel C:** FSS collapse for χ
- **Panel D:** Exponent comparison table/chart

### Figure 4: Universality Classification
- **Main:** Radar plot comparing exponents to known classes
- **Table:** Measured vs theoretical exponents
- **Error bars:** 95% confidence intervals

### Figure 5: Finite-Size Effects
- **Panel A:** Effective exponents β_eff(L)
- **Panel B:** Correction-to-scaling analysis
- **Panel C:** Systematic error estimates

## 7. Implementation Roadmap

### Phase 1: Infrastructure (Week 1)
1. Extend PhaseTransitionAnalyzer.cpp:
   - Add susceptibility calculation
   - Add Binder cumulant
   - Add correlation function measurement
2. Create FSS analysis scripts:
   - Data collapse optimization
   - Multi-exponent fitting
   - Bootstrap error analysis

### Phase 2: Coarse Scan (Week 2)
1. Run L = 32, 64, 128 with 21 noise points
2. Preliminary σ_c from Binder crossing
3. Quick exponent estimates

### Phase 3: Production Runs (Weeks 3-4)
1. Fine scan around σ_c for all sizes
2. Include L = 256, 512
3. Collect correlation functions

### Phase 4: Analysis (Week 5)
1. Extract all exponents
2. Perform data collapses
3. Verify scaling relations
4. Classification decision

### Phase 5: Documentation (Week 6)
1. Generate publication figures
2. Write detailed report
3. Theoretical interpretation
4. Prepare for publication

## 8. Risk Mitigation

### Computational Challenges
- **Mitigation:** Start with smaller grids, use parallelization
- **Fallback:** Reduce to 3 grid sizes if needed

### Equilibration Issues
- **Mitigation:** Extend equilibration time for large L
- **Check:** Monitor autocorrelation times

### Ambiguous Classification
- **Mitigation:** Measure additional observables
- **Consult:** Theoretical physicists for interpretation

## 9. Expected Outcomes

### Scenario 1: Confirms 2D Ising
- All exponents match within error bars
- Publish confirmation of universality
- Focus on physical mechanism

### Scenario 2: Identifies Different Known Class
- Clear match to XY or Potts
- Investigate symmetry breaking pattern
- Theoretical model revision

### Scenario 3: Novel Universality Class
- **High Impact Discovery**
- Comprehensive exponent table
- Search for theoretical framework
- Multiple follow-up studies needed

## 10. Conclusion

This finite-size scaling strategy provides a systematic approach to definitively classify the SMFT universality class. The measured β = 0.099 ± 0.004 already suggests deviation from 2D Ising, making this investigation critical for understanding the fundamental physics. With ~250 CPU hours of computation and systematic analysis, we can either confirm membership in a known universality class or establish SMFT as exhibiting novel critical behavior - either outcome would be scientifically significant.

The comprehensive test matrix (5 sizes × 41 points) combined with multiple analysis methods ensures robust exponent extraction with quantified uncertainties. The decision tree provides clear classification criteria, while the publication-quality figures will effectively communicate results to the scientific community.

## References

1. Cardy, J. "Scaling and Renormalization in Statistical Physics" (1996)
2. Fisher, M.E. "The theory of equilibrium critical phenomena" Rep. Prog. Phys. 30, 615 (1967)
3. Binder, K. "Finite size scaling analysis of Ising model block distribution functions" Z. Phys. B 43, 119 (1981)
4. Pelissetto, A. & Vicari, E. "Critical phenomena and renormalization-group theory" Phys. Rep. 368, 549 (2002)
5. Hasenbusch, M. "Monte Carlo studies of the Ising model" series (2001-2020)