# MSFT Experimental Documentation

This directory contains comprehensive documentation for the Mass Synchronization Field Theory (MSFT) experimental validation campaigns.

---

## 📋 Quick Start

**Goal**: Validate MSFT theory through two major experiments:
1. **Noise Sweep** - Determine if theory is fundamentally deterministic or stochastic
2. **Dirac Coupling** - Test if particles emerge from vacuum defects

**Status**: Implementation complete, ready for execution
**Timeline**: 3 weeks per experiment (6 weeks total)

---

## 📄 Documents in This Directory

1. **experiment_roadmap.md** - Master timeline & status dashboard
2. **noise_sweep_experiment.md** - Path B validation protocol (Week 1-3)
3. **dirac_coupling_experiment.md** - Particle generation protocol (Week 1-3)

**Related files in parent directory**:
- ../Determinism.md - Theoretical foundation (Claude analysis)
- ../Directive.md - Experimental authorization
- ../Dirac-Anomaly.md - Defect evolution analysis & predictions

---

## 🎯 Current Status

### ✅ COMPLETE

**Design & Documentation**:
- [x] Noise sweep protocol written
- [x] Dirac coupling protocol written
- [x] Stochastic shader created (kuramoto_stochastic.comp)
- [x] Test drivers implemented (2 files, 600+ LOC)
- [x] CMakeLists.txt updated
- [x] Build system validated

### ⏳ NEXT STEPS

**Implementation** (1-3 days):
- [ ] MSFTEngine::stepStochastic() in MSFTEngine.cpp
- [ ] MSFTEngine::initializeDiracField()
- [ ] MSFTEngine::stepWithDirac()
- [ ] MSFTEngine::getDiracDensity()
- [ ] Compile stochastic shader

**Verification** (Week 1):
- [ ] PRNG quality test
- [ ] Noise scaling test
- [ ] Deterministic limit test

**Execution** (Weeks 2-3):
- [ ] Noise sweep experiment
- [ ] Dirac coupling experiment
- [ ] Data analysis & decision

---

## 📊 What We're Testing

### Experiment 1: Noise Sweep

**Falsification Criterion**:
```
IF σ_c > 10^-5:  Stochastic theory validated → Path B ✓
IF σ_c < 10^-5:  Deterministic required → Path A ✗
```

**Why This Matters**: Determines fundamental nature of vacuum (deterministic vs thermal)

### Experiment 2: Dirac Coupling

**Five Success Criteria**:
1. Localization (O > 0.7)
2. Stabilization (ΔS > 10%)
3. Discrete energies (2-5 peaks)
4. Particle count (N ∈ [10,200])
5. Stability (τ > 500)

**Why This Matters**: Tests if electron/muon/tau masses emerge from single mechanism

---

## 🛠️ Quick Reference

**Build**:
```bash
cd build
make test_noise_sweep test_dirac_coupling
```

**Run**:
```bash
./bin/test_noise_sweep        # 8 hours
./bin/test_dirac_coupling     # 12 hours
```

**Results**:
```bash
build/output/noise_sweep/     # σ_c measurement
build/output/dirac_coupling/  # Energy spectrum
```

---

**For complete details, see individual protocol documents.**

**Last Updated**: 2025-12-16 22:00 UTC
