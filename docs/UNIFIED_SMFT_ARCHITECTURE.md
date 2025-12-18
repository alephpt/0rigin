# Unified SMFT Architecture

## Overview

The SMFT application now has a unified architecture with a single `./smft` executable that supports both interactive visualization and automated testing modes.

## Usage

```bash
# Interactive visualization mode (default)
./smft

# Test mode with YAML configuration
./smft --test config/timesync_validation.yaml
./smft --test config/quick_validation.yaml

# Help
./smft --help
```

## Architecture

### Single Binary: `./smft`

**Location**: `build/bin/smft`

**Modes**:
1. **Interactive Mode** (no arguments)
   - Full GPU-accelerated visualization
   - Real-time parameter tuning
   - Visual analysis of SMFT dynamics

2. **Test Mode** (`--test <config.yaml>`)
   - Automated quantitative validation
   - Configuration-driven testing
   - CSV output and comprehensive reports

### Directory Structure

```
src/
├── SMFT.{h,cpp}                    # Interactive visualization main class
├── SMFTEngine.{h,cpp}              # Core physics engine
├── SMFTPipelineFactory.{h,cpp}     # GPU pipeline management
├── SMFTBufferManager.{h,cpp}       # Vulkan buffer operations
├── SMFTCompute.{h,cpp}             # GPU compute dispatch
├── SMFTDescriptorManager.{h,cpp}   # Vulkan descriptor management
├── DiracEvolution.{h,cpp}          # Dirac field evolution (CPU)
├── output/
│   └── OutputManager.{h,cpp}       # Systematic output generation
└── simulations/                    # Test framework
    ├── TestConfig.{h,cpp}          # YAML configuration parser
    ├── SMFTTestRunner.{h,cpp}      # Test execution and validation
    └── ObservableComputer.{h,cpp}  # Observable computation

config/
├── timesync_validation.yaml        # Full validation (N=1,10,100)
└── quick_validation.yaml           # Quick test (N=1,10)

main.cpp                             # Unified entry point with mode dispatch
```

## Test Framework Features

### Quantitative Validation

1. **Norm Conservation**: ||Ψ||² - 1 < 10⁻⁴
2. **Energy Conservation**: |ΔE/E₀| < 1%
3. **Convergence**: Compare observables across N ratios < 5% error

### Observables Computed

- Norm: ||Ψ||² (should be ≈ 1.0)
- Energy: Total, kinetic, potential
- Position expectation: <Ψ|x|Ψ>, <Ψ|y|Ψ>
- Momentum expectation: <Ψ|p_x|Ψ>, <Ψ|p_y|Ψ>
- Sync field stats: R_avg, R_max, R_min, R_variance
- Validation flags: norm_valid, energy_valid

### Configuration Format (YAML)

```yaml
test_name: "timesync_validation"
description: "Validate operator splitting with N=1, 10, 100"

grid:
  size_x: 64
  size_y: 64

physics:
  delta: 2.5
  coupling: 0.1
  dt: 0.01
  total_steps: 100

operator_splitting:
  enabled: true
  substep_ratios: [1, 10, 100]

validation:
  norm_tolerance: 1.0e-4
  energy_tolerance: 1.0e-2
  convergence_tolerance: 0.05

output:
  directory: "output/timesync_validation"
  save_every: 10
  formats: ["csv"]
```

### Output Structure

```
output/
└── timesync_validation/
    ├── N_1/
    │   └── observables.csv
    ├── N_10/
    │   └── observables.csv
    ├── N_100/
    │   └── observables.csv
    └── test_report.txt
```

## Implementation Details

### Command-Line Parsing (main.cpp)

```cpp
int main(int argc, char* argv[]) {
    if (argc == 1) {
        return runInteractiveMode();  // Default: visualization
    }
    
    if (argv[1] == "--test") {
        return runTestMode(argv[2]);  // Test mode with config
    }
    
    if (argv[1] == "--help") {
        printUsage(argv[0]);
        return 0;
    }
}
```

### Test Execution Flow

1. **Load Configuration**: Parse YAML with TestConfig
2. **Initialize**: Create Nova, SMFTEngine, set parameters
3. **Run Tests**: For each N ratio:
   - Initialize fields (Dirac, Kuramoto)
   - Run evolution for total_steps
   - Compute observables every save_every steps
   - Validate norm, energy conservation
4. **Convergence Check**: Compare observables across N ratios
5. **Generate Report**: Summary with pass/fail status

### Observable Computation

ObservableComputer provides centralized physics calculations:

```cpp
auto obs = ObservableComputer::compute(
    dirac,           // DiracEvolution state
    R_field,         // Sync field
    delta,           // Mass gap
    time,            // Current time
    E0,              // Initial energy
    norm_tolerance,  // 1e-4
    energy_tolerance // 1e-2
);
```

## Benefits of Unified Architecture

1. **Single Binary**: No confusion about which executable to run
2. **Consistent Interface**: Same engine for visualization and testing
3. **Modular Design**: Test framework in `src/simulations/` is self-contained
4. **Flexible Testing**: Easy to add new test configurations
5. **Professional**: Standard CLI with `--help`, `--test` flags

## Future Extensions

Potential additions to `./smft`:

```bash
./smft --benchmark <config>         # Performance benchmarking
./smft --profile <config>           # GPU profiling mode
./smft --export <config>            # Batch export mode
./smft --analyze <output_dir>       # Post-processing analysis
```

All can be added as new modes dispatched from main.cpp.
