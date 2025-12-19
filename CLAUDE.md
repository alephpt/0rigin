1. keep the code clean, modular, reusable - whenever possible
2. make sure tests are distinct injections of our SMFT engine
3. make sure documents and outputs are done in the correct directories with the correct subdirectory labeling
4. we want expectations, methodologies, results and summaries documented and anything that can be visualized should be visualized.
5. unless you have confirmed that the GPU code is working, use the CPU code
6. don't make redundant files, fix what we have. this should be a clean and professional academic environment
7. don't overwrite things that will be crucial for research or evidentiary in proving methodologies and results

---

## SMFT Test Runner Usage

### Basic Usage

```bash
# Run a test from config file
./build/bin/smft --test config/your_test.yaml

# Output will be automatically logged to:
# output/YYYYMMDD_HHMMSS_testname_gridsize/console.log
# (No need to manually pipe output to a log file)
```

### Grid Convergence Testing

To test multiple grid sizes from a single config file, specify `grid_sizes` array in your YAML:

```yaml
grid:
  size_x: 128          # Default size (used if grid_sizes is empty)
  size_y: 128
  L_domain: 100.0      # Domain size in Planck lengths
  grid_sizes: [16, 32, 64]  # Test all three grid sizes automatically
```

This will create separate output directories for each grid size:
- `output/YYYYMMDD_HHMMSS_testname_16x16/`
- `output/YYYYMMDD_HHMMSS_testname_32x32/`
- `output/YYYYMMDD_HHMMSS_testname_64x64/`

### Output Structure

Each test run creates a timestamped directory containing:

```
output/YYYYMMDD_HHMMSS_testname_gridsize/
├── console.log              # Complete console output (automatically generated)
├── test_report.txt          # Validation results summary
├── N_1/
│   ├── observables.csv      # Time series data (norm, energy, etc.)
│   ├── theta_field_t0.csv   # Spatial snapshots of phase field
│   ├── R_field_t0.csv       # Spatial snapshots of order parameter
│   └── ...
├── N_10/
│   └── ...
└── N_100/
    └── ...
```

### Console Logging

The SMFT Test Runner **automatically logs all console output** to `console.log` in each output directory.

**You no longer need to manually redirect output** with `> logfile.log` or `| tee logfile.log`.

Just run:
```bash
./build/bin/smft --test config/your_test.yaml
```

The console output will:
- Display on your terminal in real-time
- Simultaneously write to `output/YYYYMMDD_HHMMSS_testname/console.log`

For grid convergence tests with multiple grid sizes, each grid gets its own console.log:
- `output/YYYYMMDD_HHMMSS_testname_16x16/console.log`
- `output/YYYYMMDD_HHMMSS_testname_32x32/console.log`
- `output/YYYYMMDD_HHMMSS_testname_64x64/console.log`

### Validation and Reporting

After test completion, check:

1. **test_report.txt** - Contains:
   - Test configuration summary
   - Validation results (norm conservation, energy conservation)
   - Convergence analysis between N values
   - Overall PASS/FAIL status

2. **console.log** - Contains:
   - Complete execution log with all diagnostic messages
   - Physical parameter initialization details
   - Step-by-step progress output
   - Validation checkpoints

### Example Workflow

```bash
# 1. Run grid convergence test (16x16, 32x32, 64x64)
./build/bin/smft --test config/defect_localization_validation.yaml

# 2. Check latest output directory
ls -lt output/ | head -5

# 3. View validation report
cat output/20251219_112601_defect_localization_64x64/test_report.txt

# 4. Review full console log
less output/20251219_112601_defect_localization_64x64/console.log

# 5. Analyze time series data
python3 visualize.py output/20251219_112601_defect_localization_64x64/
```

### Grid-Independent Initialization

All physical parameters in config files use **Planck units** (ℏ = c = G = Δ = 1):

```yaml
initial_conditions:
  dirac:
    type: "gaussian"
    x0_physical: 60.0        # Position in Planck lengths (grid-independent)
    y0_physical: 50.0
    sigma_physical: 3.0      # Width in Planck lengths

  kuramoto:
    vortex_core_radius: 3.0  # Radius in Planck lengths
    vortex_center_x: 50.0    # Center in Planck lengths
    vortex_center_y: 50.0
```

The test runner automatically converts these physical coordinates to grid units based on:
- `L_domain`: Physical domain size (default: 100 ℓ_P)
- `grid_sizes`: Grid resolution array

This ensures **grid-independent results** - the same physical system is simulated at different resolutions.
