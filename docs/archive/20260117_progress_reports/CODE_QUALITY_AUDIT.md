# Code Quality Audit Report
**Date**: January 17, 2026
**Purpose**: Identify and eliminate duplicates, fragmented code, and unnecessary files

## 1. DUPLICATE IMPLEMENTATIONS FOUND

### Dirac Solvers (2 separate implementations)
1. **Dirac3D** (`include/Dirac3D.h`, `src/Dirac3D.cpp`) - 3D+1 dimensional solver
2. **DiracEvolution** (`src/DiracEvolution.h`, `src/DiracEvolution.cpp`) - 2D solver

**Action Required**: Consolidate into single unified Dirac solver supporting 2D/3D

### Conservative Solver (1 implementation - OK)
- **ConservativeSolver** (`include/ConservativeSolver.h`, `src/ConservativeSolver.cpp`)
- Single implementation, no duplicates

### TRD Core Classes (2 versions)
1. **TRDCore3D** (`include/TRDCore3D.h`) - 3D implementation
2. **TRDCore** (`include/TRDCore.h`) - Original 2D implementation

**Action Required**: Verify if both needed or consolidate

## 2. TEST FILE CHAOS

### Statistics
- **Total test files found**: 114 test_*.cpp files
- **In test/ directory**: 88 files (proper location)
- **In root directory**: 13 files (WRONG LOCATION)
- **Test executables in build/**: 20 executables

### Misplaced Test Files in Root
```
./test_dirac_debug.cpp
./test_msft_validation.cpp
./test_em_wave_simple.cpp
./test_msft_gpu.cpp
./test_grid_resolution.cpp
./test_msft_phase2.cpp
./test_strang_validation.cpp
./test_velocity_verlet_chiral.cpp
./test_descriptor_bindings.cpp
./test_laplacian_accuracy.cpp
./test_spatial_order_comparison.cpp
./test_integrator_comparison.cpp
./test_spatial_order_timestep.cpp
```

### Duplicate Test Executables (violates single executable rule)
Found 20 test executables in build/ when there should be NONE (all tests should use `./trd --test`)

## 3. PARTIAL IMPLEMENTATIONS & TODOs

### TODOs Found (9 in TRDEngine3D.cpp)
All in `src/TRDEngine3D.cpp`:
- GPU dispatch for kuramoto3d.comp
- GPU dispatch for sync_field3d.comp
- Allocate Ex, Ey, Ez, Bx, By, Bz buffers
- Dispatch maxwell_evolve_E.comp and maxwell_evolve_B.comp
- Allocate phi, Ax, Ay, Az, A0 buffers
- Implement gauge transformation
- Initialize phase field with vortex ring
- Allocate 4-component spinor buffer
- Dispatch dirac3d.comp

**Action Required**: These are incomplete GPU implementations - either complete or remove

## 4. UNNECESSARY FILES

### Backup Files
- `./main.cpp.backup` - DELETE

### Build Files in Wrong Location
- None found in src/ or include/ (GOOD)

### Dead Shader Files
- No dead shaders found (already cleaned)

## 5. CONSOLIDATION PLAN

### Phase 1: Clean Test Files
1. Move 13 test files from root to test/
2. Remove test executables from CMakeLists.txt
3. Delete test executables from build/
4. Integrate all tests into `./trd --test` framework

### Phase 2: Merge Duplicate Implementations
1. Consolidate Dirac3D and DiracEvolution into unified solver
2. Merge TRDCore and TRDCore3D if possible
3. Update all references

### Phase 3: Complete or Remove Stubs
1. Review 9 TODOs in TRDEngine3D.cpp
2. Either implement GPU features or remove stubs
3. Clean up incomplete code paths

### Phase 4: Delete Unnecessary Files
1. Delete main.cpp.backup
2. Clean any build artifacts
3. Remove deprecated code

## 6. VERIFICATION METRICS

### Before
- Test files in root: 13
- Test executables: 20
- Duplicate Dirac implementations: 2
- TODOs: 9
- Backup files: 1

### Target After
- Test files in root: 0
- Test executables: 0 (only ./trd)
- Duplicate implementations: 0
- TODOs: 0
- Backup files: 0

## 7. IMMEDIATE ACTIONS

### Critical (Block Further Development)
1. Move test files from root to test/
2. Remove test executables from CMakeLists.txt
3. Delete main.cpp.backup

### High Priority
1. Consolidate Dirac implementations
2. Address TODOs in TRDEngine3D.cpp

### Medium Priority
1. Review TRDCore vs TRDCore3D
2. Clean up worktree duplicates

## 8. RISKS & CONSIDERATIONS

- **Test Migration Risk**: Moving to single executable may break existing workflows
- **GPU TODOs**: Incomplete GPU code may be needed for future features
- **Worktree Files**: 26 test files in worktrees/feature-wave1d-gauge-covariant-em

## 9. RECOMMENDATION

Start with Phase 1 (test cleanup) as it has lowest risk and highest impact. Then proceed with duplicate consolidation only after verifying nothing depends on the separate implementations.