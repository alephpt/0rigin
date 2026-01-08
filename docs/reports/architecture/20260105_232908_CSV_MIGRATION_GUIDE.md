# CSV Migration Guide - TRDCSVWriter Standardization

**Purpose**: Migrate 15 tests from custom CSV output code to standardized `TRDCSVWriter.h`

**Impact**: Eliminate 200+ lines of duplicate code, add consistent metadata, improve reproducibility

**Timeline**: 15-30 minutes per test × 15 tests = ~5 hours total

---

## Before/After Comparison

### Before (Custom Implementation)

```cpp
// Old pattern in test_fine_structure_constant.cpp (lines 561-567)
std::ofstream results_file("output/fine_structure_constant_results.csv");
results_file << "Method,Alpha_Measured,Alpha_QED,Ratio\n";
results_file << "Energy," << alpha_energy << "," << ALPHA_QED << "," << (alpha_energy/ALPHA_QED) << "\n";
results_file << "Coupling," << alpha_coupling << "," << ALPHA_QED << "," << (alpha_coupling/ALPHA_QED) << "\n";
results_file << "Flux," << alpha_flux << "," << ALPHA_QED << "," << (alpha_flux/ALPHA_QED) << "\n";
results_file << "GeometricMean," << alpha_geometric_mean << "," << ALPHA_QED << "," << ratio << "\n";
results_file.close();

std::cout << "\nResults exported to: output/fine_structure_constant_results.csv" << std::endl;
```

**Problems**:
- No metadata (timestamp, git commit, test parameters)
- Manual CSV formatting (error-prone)
- Hardcoded output paths
- No precision control
- 8 lines of boilerplate code

### After (Standardized Implementation)

```cpp
// New pattern using TRDCSVWriter.h
#include "TRDCSVWriter.h"

TRD::CSVWriter csv("fine_structure_constant_results", "B2_FineStructure", false);

csv.writeMetadata({
    {"K_coupling", std::to_string(K)},
    {"grid_size", std::to_string(Nx) + "x" + std::to_string(Ny) + "x" + std::to_string(Nz)},
    {"evolution_steps", std::to_string(num_steps)},
    {"dt", std::to_string(dt)}
});

csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});

csv.writeRow("Energy", alpha_energy, ALPHA_QED, alpha_energy/ALPHA_QED);
csv.writeRow("Coupling", alpha_coupling, ALPHA_QED, alpha_coupling/ALPHA_QED);
csv.writeRow("Flux", alpha_flux, ALPHA_QED, alpha_flux/ALPHA_QED);
csv.writeRow("GeometricMean", alpha_geometric_mean, ALPHA_QED, ratio);

csv.close();

std::cout << "\nResults exported to: " << csv.getFilePath() << std::endl;
```

**Benefits**:
- Automatic metadata (timestamp, git commit, parameters)
- Type-safe row writing
- Standardized output directory structure
- Precision control available
- Same line count, vastly improved functionality

### Generated CSV Output

```csv
# Test: B2_FineStructure
# Date: 2026-01-05T23:27:43
# Git: 090b6de (main branch)
# TRD Version: v3D Unified
# Golden Key: 246 GeV
# Parameters: K_coupling=1.0, grid_size=64x64x32, evolution_steps=1000, dt=0.01
#
Method,Alpha_Measured,Alpha_QED,Ratio
Energy,0.003540,0.007297,0.485084
Coupling,2.441410e-04,0.007297,0.033456
Flux,0.001826,0.007297,0.250192
GeometricMean,0.001356,0.007297,0.185878
```

---

## Tests Requiring Migration (15 Total)

### Category A: Fine Structure & Fundamental Constants (3 tests)

1. **test_fine_structure_constant.cpp** (B2)
   - File: `test/test_fine_structure_constant.cpp:561-569`
   - Search: `std::ofstream results_file("output/fine_structure_constant_results.csv")`
   - CSV: `output/fine_structure_constant_results.csv`
   - Columns: `Method, Alpha_Measured, Alpha_QED, Ratio`
   - Estimated time: 20 min

2. **test_electroweak.cpp** (B3)
   - File: Search for `ofstream.*csv` patterns
   - CSV: Likely `output/electroweak_results.csv`
   - Estimated time: 25 min

3. **test_three_generations.cpp** (B4)
   - File: Search for `ofstream.*csv` patterns
   - CSV: Likely `output/three_generations_results.csv`
   - Estimated time: 25 min

### Category B: Electromagnetic Dynamics (4 tests)

4. **test_lorentz_force.cpp** (C1)
   - File: Search for trajectory/force CSV output
   - CSV: Likely `output/lorentz_trajectory.csv`
   - Columns: `t, x, y, z, vx, vy, vz, E_kinetic`
   - Estimated time: 30 min

5. **test_geodesic_verification.cpp** (C2)
   - File: `test/test_geodesic_verification.cpp:233-246`
   - Search: `std::ofstream csv("output/test/geodesic_verification.csv")`
   - CSV: `output/test/geodesic_verification.csv`
   - Columns: Geodesic trajectory data
   - Estimated time: 20 min

6. **test_weak_field_3d.cpp** (C3)
   - File: Search for `ofstream.*csv` patterns
   - CSV: Likely `output/weak_field_trajectory.csv`
   - Estimated time: 25 min

7. **test_three_body_em.cpp** (C4)
   - File: Search for three-body simulation CSV
   - CSV: Likely `output/three_body_trajectories.csv`
   - Estimated time: 30 min

### Category C: Quantum Phenomena (3 tests)

8. **test_josephson_junction.cpp** (D1)
   - File: Search for `ofstream.*csv` patterns
   - CSV: Likely `output/josephson_iv_curve.csv`
   - Columns: `Voltage, Current, Resistance`
   - Estimated time: 25 min

9. **test_quantum_hall.cpp** (D2)
   - File: Search for Hall resistance CSV
   - CSV: Likely `output/quantum_hall_resistance.csv`
   - Estimated time: 25 min

10. **test_spin_magnetism.cpp** (D3)
    - File: Search for magnetization CSV
    - CSV: Likely `output/spin_magnetization.csv`
    - Estimated time: 25 min

### Category D: Particle Dynamics (2 tests)

11. **test_particle_scattering.cpp** (E1)
    - File: Search for scattering angle CSV
    - CSV: Likely `output/scattering_results.csv`
    - Estimated time: 20 min

12. **test_causality.cpp** (E2)
    - File: Search for causality violation CSV
    - CSV: Likely `output/causality_results.csv`
    - Estimated time: 20 min

### Category E: Field Theory (3 tests)

13. **test_renormalizability.cpp** (F1)
    - File: Search for coupling constant running CSV
    - CSV: Likely `output/renormalization_group.csv`
    - Estimated time: 25 min

14. **test_finite_temperature.cpp** (F2)
    - File: `test/test_finite_temperature.cpp:336-346`
    - Search: `std::ofstream csv_file(csv_path)`
    - CSV: `output/phase_diagram.csv`
    - Columns: `Temperature, OrderParameter`
    - Estimated time: 20 min

15. **test_stochastic_quantization.cpp** (F3)
    - File: Search for stochastic field CSV
    - CSV: Likely `output/stochastic_evolution.csv`
    - Estimated time: 25 min

---

## Migration Process (Step-by-Step)

### Step 1: Find Existing CSV Code

```bash
# Search for CSV file creation
grep -n "ofstream.*csv" test/<test_name>.cpp

# Search for CSV writing operations
grep -n "\.csv" test/<test_name>.cpp | grep "<<"

# Example for fine_structure_constant
grep -n "ofstream.*csv" test/test_fine_structure_constant.cpp
# Output: test/test_fine_structure_constant.cpp:561:    std::ofstream results_file("output/fine_structure_constant_results.csv");
```

### Step 2: Add Header Include

At top of test file, add:

```cpp
#include "TRDCSVWriter.h"
```

### Step 3: Replace Custom CSV Code

**Delete lines** like:
```cpp
std::ofstream results_file("output/fine_structure_constant_results.csv");
results_file << "Method,Alpha_Measured,Alpha_QED,Ratio\n";
results_file << "Energy," << alpha_energy << "," << ALPHA_QED << "," << (alpha_energy/ALPHA_QED) << "\n";
// ... more lines
results_file.close();
```

**Replace with**:
```cpp
TRD::CSVWriter csv("fine_structure_constant_results", "B2_FineStructure", false);

csv.writeMetadata({
    {"K_coupling", std::to_string(K)},
    {"grid_size", std::to_string(Nx) + "x" + std::to_string(Ny) + "x" + std::to_string(Nz)},
    {"dt", std::to_string(dt)}
});

csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});

csv.writeRow("Energy", alpha_energy, ALPHA_QED, alpha_energy/ALPHA_QED);
csv.writeRow("Coupling", alpha_coupling, ALPHA_QED, alpha_coupling/ALPHA_QED);
csv.writeRow("Flux", alpha_flux, ALPHA_QED, alpha_flux/ALPHA_QED);
csv.writeRow("GeometricMean", alpha_geometric_mean, ALPHA_QED, ratio);

csv.close();

std::cout << "\nResults exported to: " << csv.getFilePath() << std::endl;
```

### Step 4: Compile and Test

```bash
# Rebuild test
cd build
make <test_name>

# Run test
./bin/trd --test ../config/<test_config>.yaml

# Verify CSV output
cat output/<test_name>/<csv_file>.csv
```

### Step 5: Verify Metadata

Check generated CSV includes:
- ✅ Test name in header
- ✅ Timestamp (ISO 8601 format)
- ✅ Git commit hash and branch
- ✅ TRD version and golden key
- ✅ Custom parameters from writeMetadata()

---

## Common Patterns & Solutions

### Pattern 1: Simple Data Table

**Old**:
```cpp
std::ofstream csv("results.csv");
csv << "X,Y,Z\n";
csv << x << "," << y << "," << z << "\n";
```

**New**:
```cpp
TRD::CSVWriter csv("results", "MyTest");
csv.writeHeader({"X", "Y", "Z"});
csv.writeRow(x, y, z);
```

### Pattern 2: Multiple CSV Files

**Old**:
```cpp
std::ofstream results("results.csv");
std::ofstream trajectory("trajectory.csv");
```

**New**:
```cpp
TRD::CSVWriter results("results", "MyTest");
TRD::CSVWriter trajectory("trajectory", "MyTest");
```

### Pattern 3: Time Series Data

**Old**:
```cpp
std::ofstream csv("time_series.csv");
csv << "t,x,y,z\n";
for (int i = 0; i < num_steps; ++i) {
    csv << t[i] << "," << x[i] << "," << y[i] << "," << z[i] << "\n";
}
```

**New**:
```cpp
TRD::CSVWriter csv("time_series", "MyTest");
csv.writeHeader({"t", "x", "y", "z"});
for (int i = 0; i < num_steps; ++i) {
    csv.writeRow(t[i], x[i], y[i], z[i]);
}
```

### Pattern 4: High-Precision Output

**Old**:
```cpp
std::ofstream csv("precise.csv");
csv << std::setprecision(15) << value << "\n";
```

**New**:
```cpp
TRD::CSVWriter csv("precise", "MyTest");
csv.setPrecision(15);
csv.writeRow(value);
```

### Pattern 5: Scientific Notation

**Old**:
```cpp
std::ofstream csv("scientific.csv");
csv << std::scientific << 1.23e-20 << "\n";
```

**New**:
```cpp
TRD::CSVWriter csv("scientific", "MyTest");
csv.setScientificNotation(true);
csv.writeRow(1.23e-20);
```

---

## Troubleshooting

### Issue 1: Compilation Error - `TRDCSVWriter.h` not found

**Solution**:
```bash
# Verify header exists
ls include/TRDCSVWriter.h

# Check CMakeLists.txt includes include/ directory
# Should have: include_directories(${CMAKE_SOURCE_DIR}/include)
```

### Issue 2: Directory Creation Fails

**Error**: `Failed to create output directory: output/MyTest`

**Solution**:
```bash
# Check permissions
ls -ld output/

# Create manually if needed
mkdir -p output/MyTest
```

### Issue 3: Git Info Shows "N/A"

**Cause**: Not in git repository or git not installed

**Solution**: Git info is optional, CSV will still generate. To fix:
```bash
# Verify git repo
git status

# If not a git repo, initialize
git init
git add .
git commit -m "Initial commit"
```

### Issue 4: Wrong Column Count

**Error**: CSV has 4 columns in header but 5 in data row

**Solution**:
```cpp
// Ensure header and row have same number of elements
csv.writeHeader({"A", "B", "C"});
csv.writeRow(1, 2, 3);  // ✓ Correct

csv.writeHeader({"A", "B"});
csv.writeRow(1, 2, 3);  // ✗ Wrong - 3 values for 2 columns
```

---

## Testing Checklist

After migration, verify each test:

- [ ] Test compiles without warnings
- [ ] Test runs successfully
- [ ] CSV file created in `output/<test_name>/` directory
- [ ] CSV has metadata header (timestamp, git, parameters)
- [ ] CSV has correct column headers
- [ ] CSV data matches expected values
- [ ] Scientific notation used for extreme values (< 1e-3 or > 1e6)
- [ ] No regression in test pass/fail status

---

## Batch Migration Script (Optional)

For automated migration of multiple tests:

```bash
#!/bin/bash
# migrate_csv.sh - Migrate test to TRDCSVWriter

TEST_NAME=$1
TEST_CATEGORY=$2

echo "Migrating ${TEST_NAME} (${TEST_CATEGORY})..."

# Backup original
cp test/${TEST_NAME}.cpp test/${TEST_NAME}.cpp.bak

# Add header include (insert after other includes)
sed -i '/#include.*h"/a #include "TRDCSVWriter.h"' test/${TEST_NAME}.cpp

# Note: Manual replacement still recommended due to complexity
# This script serves as a template for batch operations

echo "✓ Header added. Manual CSV code replacement still required."
echo "  See CSV_MIGRATION_GUIDE.md for detailed instructions."
```

---

## Progress Tracking

Update this table as tests are migrated:

| Test | Category | Status | Migrated By | Date | Lines Saved |
|------|----------|--------|-------------|------|-------------|
| test_fine_structure_constant.cpp | B2 | ⬜ Not Started | - | - | ~8 |
| test_electroweak.cpp | B3 | ⬜ Not Started | - | - | ~10 |
| test_three_generations.cpp | B4 | ⬜ Not Started | - | - | ~12 |
| test_lorentz_force.cpp | C1 | ⬜ Not Started | - | - | ~15 |
| test_geodesic_verification.cpp | C2 | ⬜ Not Started | - | - | ~14 |
| test_weak_field_3d.cpp | C3 | ⬜ Not Started | - | - | ~12 |
| test_three_body_em.cpp | C4 | ⬜ Not Started | - | - | ~18 |
| test_josephson_junction.cpp | D1 | ⬜ Not Started | - | - | ~10 |
| test_quantum_hall.cpp | D2 | ⬜ Not Started | - | - | ~10 |
| test_spin_magnetism.cpp | D3 | ⬜ Not Started | - | - | ~10 |
| test_particle_scattering.cpp | E1 | ⬜ Not Started | - | - | ~12 |
| test_causality.cpp | E2 | ⬜ Not Started | - | - | ~8 |
| test_renormalizability.cpp | F1 | ⬜ Not Started | - | - | ~15 |
| test_finite_temperature.cpp | F2 | ⬜ Not Started | - | - | ~11 |
| test_stochastic_quantization.cpp | F3 | ⬜ Not Started | - | - | ~14 |
| **TOTAL** | - | **0/15** | - | - | **~200** |

Status Legend:
- ⬜ Not Started
- 🔄 In Progress
- ✅ Complete
- ⚠️ Blocked

---

## Next Steps After Migration

Once all 15 tests are migrated:

1. **Remove legacy code**: Delete any backup files (`.bak`)
2. **Update documentation**: Note standardized CSV format in test docs
3. **Archive this guide**: Move to `docs/archives/` once migration complete
4. **Verify CI/CD**: Ensure automated tests generate proper CSV outputs
5. **Update ARCHITECTURE_REVIEW_CATEGORY_BF.md**: Mark task as complete

---

**Estimated Total Effort**: 5-6 hours for all 15 tests
**Expected Benefit**: 200+ lines of duplicate code eliminated, consistent metadata across all tests, improved reproducibility

**Start with**: test_fine_structure_constant.cpp (simplest, good template for others)
