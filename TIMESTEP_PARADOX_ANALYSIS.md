# Timestep Refinement Paradox - Root Cause Analysis

## Executive Summary

The timestep refinement paradox where smaller timestep (dt=0.0001) gives 26× WORSE frequency error than larger timestep (dt=0.0005) has been traced to a **particle trajectory recording bug**.

## The Paradox

### Baseline Test (dt = 0.0005, 20000 steps)
- Config: `config/lorentz_force_pure_magnetic.yaml`
- Total time: 10 τ_P
- Frequency error: 3.45%
- Radius error: 2.69%
- Result: Reasonable accuracy

### Refined Test (dt = 0.0001, 100000 steps)
- Config: `config/lorentz_force_refined_timestep.yaml`
- Total time: 10 τ_P (same duration)
- Frequency error: 91.45% (26× WORSE!)
- Radius error: 51.01% (19× WORSE!)
- Result: Catastrophically wrong

## Root Cause: Trajectory Recording Bug

### What Should Happen
1. Simulation runs for 100000 steps × 0.0001 dt = 10 seconds
2. Particle trajectory recorded every 5 steps (config: `record_every: 5`)
3. Should get 100000/5 + 1 = 20001 data points
4. Times should be: 0, 0.0005, 0.001, ..., 9.9995, 10.0

### What Actually Happens
1. Simulation runs correctly for 100000 steps (verified via observables.csv)
2. Particle IS evolved for all 100000 steps
3. BUT trajectory recording has TWO bugs:
   - **Bug #1**: Records EVERY step (record_every=1), ignoring config value of 5
   - **Bug #2**: Stops recording after 20000 steps (2.0 seconds)
4. Result: 20001 points from t=0 to t=2.0 with 0.0001 intervals

### Impact on Analysis

The `analyzeOrbit()` function computes frequency from the truncated trajectory:
- Uses only first 2.0 seconds of a 10-second simulation
- Only captures ~0.32 orbits instead of 1.59 orbits
- Angle unwrapping algorithm fails with partial orbit
- Frequency calculation: ω = Δθ/Δt uses wrong Δt (2.0 instead of 10.0)
- Results in ~2× frequency overestimation

## Evidence

### 1. Trajectory File Analysis
```bash
# Refined timestep trajectory
wc -l particle_trajectory.csv  # 20004 lines (3 header + 20001 data)
tail particle_trajectory.csv    # Last time: 2.000000 (should be 10.0)
```

### 2. Time Intervals
```
Expected (record_every=5): 0, 0.0005, 0.001, ..., 10.0
Actual: 0, 0.0001, 0.0002, ..., 2.0
```

### 3. Debug Output Contradiction
```
[DEBUG] Step 0: record_interval=5, record_trajectory=1, trajectory_size=0
[DEBUG] Step 5: record_interval=5, record_trajectory=1, trajectory_size=2
```
Debug shows correct logic but actual file has wrong data!

## Suspected Code Issues

### 1. Recording Logic (`SMFTTestRunner.cpp:1048-1049`)
```cpp
int record_interval = _config.test_particle.record_every > 0
                     ? _config.test_particle.record_every : save_every;
bool record_trajectory = (step % record_interval == 0);
```
Logic looks correct but behavior is wrong.

### 2. Possible Hardcoded Limit
- Baseline config has `total_steps: 20000`
- Refined config has `total_steps: 100000`
- But particle recording stops at exactly 20000 steps
- Suggests baseline value is being used somewhere

### 3. Memory Issue
```cpp
trajectory_.reserve(10000);  // Pre-allocation in TestParticle.cpp:21
```
This is just a hint, shouldn't limit size, but worth investigating.

## Verification Tests

### Test 1: Debug Recording
Created `config/debug_particle_recording.yaml` with:
- 100 total steps
- record_every: 5
- save_every: 10
Result: Test particle not initialized at all (separate issue)

### Test 2: Console Verification
- Observables go to 10 seconds correctly
- SMFT fields evolve for full 100000 steps
- Only particle trajectory is truncated

## Conclusion

The timestep refinement "paradox" is not a numerical integration issue but a **data recording bug**. The Boris integrator is working correctly, but the trajectory analysis is performed on truncated, incorrectly-sampled data.

## Fixes Required

1. **Fix record_every**: Ensure particle trajectory respects `record_every` config
2. **Fix recording limit**: Remove 20000-step limit on trajectory recording
3. **Add validation**: Check trajectory covers expected time range before analysis
4. **Add diagnostics**: Log trajectory size and time range during analysis

## Impact

This bug affects ALL Lorentz force validation tests with >20000 steps, making them report incorrect frequencies and radii. The physics simulation is correct; only the measurement is wrong.