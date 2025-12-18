# Empty Plot Bug Investigation Report

## Overview
Two demo scripts are generating images with empty subplots:
1. **local_coupling_demo.png** - Oscillator Phases subplot is empty
2. **SMFT_basic_evolution.png** - Bottom 2 panels (Final Sync Field & Final Mediator Field) are empty

## Bug #1: Empty Oscillator Phases Plot in local_coupling_demo.png

### Location
- **File**: `/home/persist/neotec/0rigin/examples/field_theory/SMFT_full_demo.py`
- **Lines**: 186-201 (Oscillator Phases plotting section)

### Root Cause
The demo_local_coupling() function is trying to access `result['phases']` (line 187), but the LocalFieldCoupling.evolve_coupled_system() method returns phases as a 2D array with shape (n_snapshots, N_oscillators). The code extracts `final_phases = result['phases'][-1]` correctly, but the scatter plot is not displaying any points.

### Investigation Findings
Looking at the scatter plot code (lines 188-201):
- The scatter plot uses `positions[:, 0]` and `positions[:, 1]` for x,y coordinates
- The `positions` variable is created locally in demo_local_coupling() at line 132-134
- BUT: The positions are NOT passed back from evolve_coupled_system()
- The scatter plot tries to use local positions with phases from evolution, but there's likely a mismatch

### Specific Issue
The positions used for plotting (lines 188-189) are the INITIAL positions created before evolution, but they need to match what was actually used during evolution. The evolve_coupled_system doesn't return positions, so we can't verify they match.

## Bug #2: Empty Field Plots in SMFT_basic_evolution.png

### Location
- **File**: `/home/persist/neotec/0rigin/examples/field_theory/SMFT_demo.py`
- **Lines**: 70-92 (Final field plotting section)

### Root Cause
The SMFTSystem.evolve() method is experiencing numerical instability leading to NaN values:
- Console output shows: "Final R: nan" and "Energy conservation: nan"
- This indicates the simulation diverged, producing NaN values in the fields
- When imshow() tries to plot NaN arrays, it shows empty plots

### Evidence
From the console output when running SMFT_demo.py:
```
RuntimeWarning: overflow encountered in square
RuntimeWarning: overflow encountered in add
RuntimeWarning: invalid value encountered in reduce
RuntimeWarning: overflow encountered in divide
Final R: nan
Energy conservation: nan
```

### Specific Issues
1. **Line 72**: `solution['sync_field'][-1]` contains NaN values
2. **Line 84**: `solution['mediator_field'][-1]` contains NaN values
3. The simulation is numerically unstable with current parameters (dt=0.01, diffusion_coeff=0.02, mediator_mass=10.0)

## Additional Findings

### SMFT_demo.py Issues
- Uses relative paths for saving images (line 95: 'SMFT_basic_evolution.png')
- Should use absolute paths like SMFT_full_demo.py does

### SMFT_full_demo.py Issues
- Missing data validation before plotting
- No checks for NaN/infinity values
- No error handling for failed simulations

## Suggested Fixes

### Fix for Bug #1 (Empty Oscillator Phases)
1. **Option A**: Return positions from evolve_coupled_system()
   - Add 'positions': positions to the return dict at line 208-215 in local_coupling.py
2. **Option B**: Use the positions that were actually stored
   - The positions variable in demo needs to be consistent with what was used

### Fix for Bug #2 (Empty Field Plots)
1. **Reduce time step**: Change dt from 0.01 to 0.001 or 0.005
2. **Adjust parameters**: Reduce diffusion_coeff or adjust mediator_mass
3. **Add stability checks**: Detect NaN/overflow during evolution and adapt
4. **Add data validation**: Check for NaN before plotting:
   ```python
   if np.isnan(solution['sync_field'][-1]).any():
       print("Warning: NaN values in sync field")
       # Use earlier snapshot or handle gracefully
   ```

### General Improvements
1. Add bounds checking in evolution methods
2. Implement adaptive time stepping when instability detected
3. Add validation before plotting to detect empty/NaN data
4. Use consistent absolute paths for all file saves
5. Add error recovery mechanisms in demos

## Code Snippets Showing Bugs

### Bug #1 - Missing positions data
```python
# Line 187-189 in SMFT_full_demo.py
final_phases = result['phases'][-1]  # Gets phases correctly
scatter = ax_phase.scatter(
    positions[:, 0],  # Uses LOCAL positions, not from result
    positions[:, 1],  # These may not match evolution positions
    c=final_phases,
    ...
)
```

### Bug #2 - NaN field data
```python
# Lines 71-72 in SMFT_demo.py
im1 = axes[1, 0].imshow(
    solution['sync_field'][-1].T,  # This contains NaN values
    ...
)
# Lines 83-84
im2 = axes[1, 1].imshow(
    solution['mediator_field'][-1].T,  # This also contains NaN
    ...
)
```

## Testing Recommendations
1. Run with smaller time steps to verify stability
2. Add print statements to check array contents before plotting
3. Verify positions consistency between creation and usage
4. Test with simpler parameters first (fewer oscillators, smaller grid)
5. Add progress monitoring during evolution to catch divergence early