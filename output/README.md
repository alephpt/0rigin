# Output Directory - Simulation Results

This directory contains results from various SMFT (Synchronization Mass Field Theory) validation experiments.

All experiment directories follow the naming convention: `YYYYMMDD_HHMMSS_experiment_name`

## Directory Structure

### `20251218_123158_vortex_validation/` ✅ **LATEST & VALIDATED**
**Vortex operator splitting validation with non-uniform R field (Dec 18, 2024)**

Complete validation using phase defect initial conditions to test coupled dynamics.

**Contents:**
- Analysis markdown and visualizations
- Comparison of N=1 vs N=100 operator splitting
- Spatial R(x,y) field snapshots showing vortex structure
- Demonstrates N-dependent coupling effects

**Key Findings:**
- Vortex creates non-uniform R field (σ_R ≈ 0.118, 1000× larger than uniform)
- N=100 prevents vortex core from deepening (ΔR_core = -0.053)
- Energy difference: 0.062% between N=1 and N=100
- Properly tests coupled Dirac-Kuramoto dynamics

### `20251218_122006_operator_splitting_10k_validation/` ⚠️ **NULL TEST**
**Gaussian wavepacket with nearly uniform R field (Dec 18, 2024)**

Initial 10k timestep validation - identified as null test due to uniform R field.

**Key Results:**
- Nearly uniform R field (σ_R ≈ 0.0001)
- Shows DIVERGENCE not convergence (Error increased 1.6×)
- ∇R ≈ 0 → No coupling force → Free particle evolution
- Cannot validate operator splitting for coupled dynamics

---

### `phase1_full_validation/`
**Phase 1 validation results (earlier work)**

Initial validation of the coupled Dirac-Kuramoto system.

---

### `coupled_dynamics_validation/`
**Coupled dynamics validation (earlier work)**

Testing of the synchronization field and mass field coupling.

---

### `coupled_64x64_N1/`, `coupled_64x64_N10/`, `coupled_64x64_N100/`
**Old operator splitting tests (DEPRECATED - had recursion bug)**

These directories contain results from tests that crashed due to the recursion bug. Superseded by `operator_splitting_10k_validation/`.

---

### `primitive/`
**Historical experimental data**

Early experimental results and prototypes.

---

## Quick Start

To visualize the latest 10k validation results:

```bash
cd operator_splitting_10k_validation
python3 visualize_operator_splitting_10k.py
```

This will display statistics and regenerate the plots.

---

## Analysis Scripts

- `visualize_operator_splitting_10k.py` - Latest 10k validation (in `operator_splitting_10k_validation/`)
- `visualize_phase1.py` - Phase 1 visualization

---

## What to Use

**For operator splitting research**: Use `operator_splitting_10k_validation/` - this is the validated, bug-free implementation.

**For historical reference**: Other directories preserved for comparison and debugging history.
