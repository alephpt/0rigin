# Phase 2 Implementation Coordination Plan
## Preventing Agent Conflicts and System Crashes

**Date**: December 29, 2025
**Purpose**: Ensure 6 parallel implementations don't interfere with each other

---

## File System Partitioning

Each task gets isolated workspace to prevent conflicts:

```
Task 1 (Phase Transition):
  - Code: src/validation/FiniteSizeScaling.{h,cpp}
  - Config: config/fss_analysis/
  - Output: output/phase2_task1_fss/
  - Build: Sequential after all code ready

Task 2 (Energy Conservation):
  - Code: src/validation/EnergyBudget.{h,cpp}
  - Config: config/energy_verification/
  - Output: output/phase2_task2_energy/
  - Build: Sequential after all code ready

Task 3 (EM Forces):
  - Code: src/validation/TestParticle.{h,cpp}, src/validation/EMValidator.{h,cpp}
  - Config: config/em_verification/
  - Output: output/phase2_task3_em/
  - Build: Sequential after all code ready

Task 4 (Metric Derivation):
  - Code: docs/metric_derivation/ (theoretical work - no code conflicts)
  - Scripts: analysis/metric_validation/
  - Output: output/phase2_task4_metric/
  - Build: N/A (analysis only)

Task 5 (Universality Class):
  - Code: src/validation/UniversalityClassifier.{h,cpp}
  - Config: config/universality_scan/
  - Output: output/phase2_task5_universality/
  - Build: Sequential after all code ready

Task 6 (Vacuum Energy):
  - Code: docs/vacuum_energy_resolution/ (theoretical work - no code conflicts)
  - Scripts: analysis/vacuum_energy/
  - Output: output/phase2_task6_vacuum/
  - Build: N/A (theory + analysis)
```

---

## Build Coordination

**CRITICAL: Only ONE build at a time**

### Strategy:
1. All agents implement their code changes to isolated files
2. Agents create `.ready` marker when code complete
3. **Main Claude** coordinates sequential builds:
   ```bash
   # Wait for all .ready markers
   while [ ! -f /tmp/task1.ready ] || [ ! -f /tmp/task2.ready ] || ...; do
     sleep 5
   done

   # Sequential builds
   echo "Building Task 1..." && cmake --build build && touch /tmp/task1.built
   echo "Building Task 2..." && cmake --build build && touch /tmp/task2.built
   # etc.
   ```

### Build Order:
1. Task 2 (Energy - core SMFTEngine changes)
2. Task 3 (EM - new validation classes)
3. Task 1 (FSS - analysis only, minimal changes)
4. Task 5 (Universality - analysis only)
5. Tasks 4 & 6 don't require builds (theory/analysis)

---

## Test Execution Coordination

**CRITICAL: Only ONE test running at a time to prevent GPU/memory crashes**

### Strategy: Test Queue System

Create `/tmp/test_queue.sh`:
```bash
#!/bin/bash
# Test coordination queue - ensures serial execution

LOCKFILE="/tmp/smft_test.lock"

run_test() {
    local task=$1
    local config=$2
    local output=$3

    # Wait for lock
    while [ -f "$LOCKFILE" ]; do
        echo "[$task] Waiting for test lock..."
        sleep 10
    done

    # Acquire lock
    touch "$LOCKFILE"
    echo "[$task] Lock acquired, running test..."

    # Run test
    ./build/bin/smft --test "$config" 2>&1 | tee "$output"

    # Release lock
    rm "$LOCKFILE"
    echo "[$task] Test complete, lock released"
}

export -f run_test
```

### Test Schedule:
```bash
# Task 2: Energy conservation tests (6 configs × ~5 min = 30 min)
run_test "Task2_T1" config/energy_verification/undamped.yaml output/phase2_task2_energy/T1.log
run_test "Task2_T2" config/energy_verification/damped.yaml output/phase2_task2_energy/T2.log
# ... (sequential)

# Task 3: EM verification (3 tests × ~10 min = 30 min)
run_test "Task3_A" config/em_verification/lorentz_force.yaml output/phase2_task3_em/lorentz.log
run_test "Task3_B" config/em_verification/maxwell_check.yaml output/phase2_task3_em/maxwell.log
run_test "Task3_C" config/em_verification/flux_quantization.yaml output/phase2_task3_em/flux.log

# Task 1: FSS analysis (5 grids × 41 points × 2 min = ~7 hours)
for L in 32 64 128 256 512; do
    run_test "Task1_L${L}" config/fss_analysis/grid_${L}.yaml output/phase2_task1_fss/L${L}.log
done

# Task 5: Universality (similar to Task 1, ~7 hours)
# ... (sequential)
```

**Total Test Time**: ~15 hours (serial execution)

---

## Agent Instructions for Phase 2

### General Rules:
1. ✅ **DO**: Work in your assigned file namespace
2. ✅ **DO**: Create `.ready` marker when code complete
3. ✅ **DO**: Use test queue system for any test runs
4. ❌ **DON'T**: Modify files outside your namespace
5. ❌ **DON'T**: Run `cmake --build` yourself (Main Claude handles)
6. ❌ **DON'T**: Run tests in parallel

### Task-Specific Namespaces:

**Task 1 (Data Analyst - FSS):**
- Files: `src/validation/FiniteSizeScaling.{h,cpp}`
- Configs: `config/fss_analysis/*.yaml`
- Marker: `/tmp/task1.ready`

**Task 2 (Developer - Energy):**
- Files: `src/validation/EnergyBudget.{h,cpp}`
- Configs: `config/energy_verification/*.yaml`
- Marker: `/tmp/task2.ready`

**Task 3 (Integration - EM):**
- Files: `src/validation/TestParticle.{h,cpp}`, `src/validation/EMValidator.{h,cpp}`
- Configs: `config/em_verification/*.yaml`
- Marker: `/tmp/task3.ready`

**Task 4 (Developer - Metric):**
- Files: `docs/metric_derivation/*.md`, `analysis/metric_validation/*.py`
- No build needed
- Marker: `/tmp/task4.ready`

**Task 5 (Data Analyst - Universality):**
- Files: `src/validation/UniversalityClassifier.{h,cpp}`
- Configs: `config/universality_scan/*.yaml`
- Marker: `/tmp/task5.ready`

**Task 6 (Developer - Vacuum):**
- Files: `docs/vacuum_energy_resolution/*.md`, `analysis/vacuum_energy/*.py`
- No build needed
- Marker: `/tmp/task6.ready`

---

## Phase 2 Execution Timeline

### Parallel Implementation (Agents work simultaneously):
**Hours 0-4**: All 6 agents implement their solutions
- No conflicts (isolated namespaces)
- No builds yet
- Create `.ready` markers when done

### Sequential Build Phase (Main Claude coordinates):
**Hours 4-5**: Sequential builds
```bash
# Clean build to catch any conflicts
cd /home/persist/neotec/0rigin
rm -rf build/*
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build each task sequentially
cmake --build build  # Task 2 (Energy)
cmake --build build  # Task 3 (EM)
cmake --build build  # Task 1 (FSS)
cmake --build build  # Task 5 (Universality)
```

### Sequential Test Phase (Queue system):
**Hours 5-20**: Tests run one at a time
- Energy tests: ~30 min
- EM tests: ~30 min
- FSS tests: ~7 hours
- Universality tests: ~7 hours
- Metric/Vacuum: Analysis only (parallel OK)

---

## Conflict Detection

**Main Claude monitors**:
```bash
# Check for file conflicts before Phase 2 starts
check_conflicts() {
    echo "Checking for namespace conflicts..."

    # Each task should only touch their files
    # If git diff shows changes outside namespace -> ABORT

    git diff --name-only | while read file; do
        case "$file" in
            src/validation/FiniteSizeScaling.*)
                [ "$CURRENT_TASK" != "1" ] && echo "ERROR: Task $CURRENT_TASK touched Task 1 files"
                ;;
            src/validation/EnergyBudget.*)
                [ "$CURRENT_TASK" != "2" ] && echo "ERROR: Task $CURRENT_TASK touched Task 2 files"
                ;;
            # ... etc
        esac
    done
}
```

---

## Emergency Procedures

**If agent violates namespace**:
1. Stop agent immediately
2. `git checkout` to revert changes
3. Reassign task with clearer instructions

**If test causes system crash**:
1. Kill all background processes: `pkill -f smft`
2. Clear GPU state: `nvidia-smi -r` (if NVIDIA) or reboot
3. Resume test queue from last successful test

**If build fails**:
1. Identify which task's code caused failure
2. Agent fixes in isolated namespace
3. Retry build for that task only

---

## Success Criteria

Phase 2 completes successfully when:
- ✅ All 6 tasks create `.ready` markers
- ✅ Sequential builds complete without errors
- ✅ All tests execute serially without crashes
- ✅ No file conflicts detected
- ✅ Output directories contain expected results

---

## Next Steps

After Phase 2 completes:
1. Main Claude reviews all outputs
2. Deploy Phase 3 QA agents (one per task)
3. QA agents verify against Phase 1 expectations
4. Final assessment and recommendations

**Coordination is Key**: Serial execution prevents crashes, isolated namespaces prevent conflicts.
