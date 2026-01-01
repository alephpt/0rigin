#!/bin/bash
# Real-time Batch B Validation Monitor
# Usage: watch -n 10 ./monitor_batch_b.sh

clear
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║        SMFT Batch B Comprehensive Validation - Live Monitor          ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""
date
echo ""

# Find latest batch_b validation directory
LATEST_DIR=$(ls -td output/batch_b_validation_* 2>/dev/null | head -1)

if [ -z "$LATEST_DIR" ]; then
    echo "⚠️  No Batch B validation directory found yet..."
    echo ""
    echo "Waiting for first test to start..."
    exit 0
fi

echo "📂 Output Directory: $LATEST_DIR"
echo ""

# Detect which test is currently running
CURRENT_TEST="Unknown"
if pgrep -f "scenario_2.4A" > /dev/null; then
    CURRENT_TEST="[1/7] Phase 2.4A - Velocity Threshold"
elif pgrep -f "scenario_2.5A" > /dev/null; then
    CURRENT_TEST="[2/7] Phase 2.5A - Klein-Gordon"
elif pgrep -f "scenario_2.6A" > /dev/null; then
    CURRENT_TEST="[3/7] Phase 2.6A - Energy-Momentum"
elif pgrep -f "vacuum_energy" > /dev/null; then
    CURRENT_TEST="[4/7] Phase 3.2 - Vacuum Energy"
elif pgrep -f "time_dilation" > /dev/null; then
    CURRENT_TEST="[5/7] Phase 4.1 - Time Dilation"
elif pgrep -f "scenario_2.6B" > /dev/null; then
    CURRENT_TEST="[6/7] Phase 5 - EM Coupling"
elif pgrep -f "phase_2.3" > /dev/null; then
    CURRENT_TEST="[7/7] Regression - Phase 2.3"
else
    CURRENT_TEST="🏁 All tests complete or between tests"
fi

echo "🔄 Current Test: $CURRENT_TEST"
echo ""

# Test completion status
echo "═══════════════════════════════════════════════════════════════════════"
echo "TEST COMPLETION STATUS"
echo "═══════════════════════════════════════════════════════════════════════"

declare -A tests=(
    ["2.4A"]="Phase 2.4A - Velocity Threshold"
    ["2.5A"]="Phase 2.5A - Klein-Gordon"
    ["2.6A"]="Phase 2.6A - Energy-Momentum"
    ["3.2"]="Phase 3.2 - Vacuum Energy"
    ["4.1"]="Phase 4.1 - Time Dilation"
    ["5"]="Phase 5 - EM Coupling"
    ["reg"]="Regression - Phase 2.3"
)

for key in "2.4A" "2.5A" "2.6A" "3.2" "4.1" "5" "reg"; do
    name="${tests[$key]}"
    
    # Check for log file
    case $key in
        "2.4A") logfile="${LATEST_DIR}/2.4A_velocity.log" ;;
        "2.5A") logfile="${LATEST_DIR}/2.5A_klein_gordon.log" ;;
        "2.6A") logfile="${LATEST_DIR}/2.6A_energy.log" ;;
        "3.2") logfile="${LATEST_DIR}/3.2_vacuum_energy.log" ;;
        "4.1") logfile="${LATEST_DIR}/4.1_time_dilation.log" ;;
        "5") logfile="${LATEST_DIR}/5_em_coupling.log" ;;
        "reg") logfile="${LATEST_DIR}/regression_2.3.log" ;;
    esac
    
    if [ ! -f "$logfile" ]; then
        echo "⏳ $name - Pending"
    elif tail -50 "$logfile" 2>/dev/null | grep -q "TESTS PASSED"; then
        echo "✅ $name - PASSED"
    elif tail -50 "$logfile" 2>/dev/null | grep -q "TESTS FAILED"; then
        echo "❌ $name - FAILED"
    elif tail -50 "$logfile" 2>/dev/null | grep -q "SOME TESTS FAILED"; then
        echo "⚠️  $name - PARTIAL"
    elif tail -5 "$logfile" 2>/dev/null | grep -q "Step.*norm"; then
        # Extract current step from running test
        step=$(tail -50 "$logfile" | grep "Step [0-9]" | tail -1 | sed -n 's/.*Step \([0-9]*\)\/\([0-9]*\).*/\1\/\2/p')
        echo "🔄 $name - Running ($step)"
    else
        echo "🔄 $name - In Progress"
    fi
done

echo ""
echo "═══════════════════════════════════════════════════════════════════════"
echo "LATEST OUTPUT (Last 15 lines)"
echo "═══════════════════════════════════════════════════════════════════════"

# Show latest output from current test log
for logfile in "${LATEST_DIR}"/*.log; do
    if [ -f "$logfile" ]; then
        # Check if file was modified in last 60 seconds (active test)
        if [ $(( $(date +%s) - $(stat -c %Y "$logfile" 2>/dev/null || stat -f %m "$logfile" 2>/dev/null) )) -lt 60 ]; then
            echo "From: $(basename $logfile)"
            tail -15 "$logfile" 2>/dev/null | sed 's/^/  /'
            break
        fi
    fi
done

echo ""
echo "═══════════════════════════════════════════════════════════════════════"
echo "RUNTIME ESTIMATE"
echo "═══════════════════════════════════════════════════════════════════════"

# Count completed tests
COMPLETED=$(ls -1 "${LATEST_DIR}"/*.log 2>/dev/null | wc -l)
TOTAL=7

# Estimate based on test times (2+2+4+2+3+4+1 = 18h)
declare -A test_times=(["2.4A"]=2 ["2.5A"]=2 ["2.6A"]=4 ["3.2"]=2 ["4.1"]=3 ["5"]=4 ["reg"]=1)
ELAPSED_EST=0

for key in "2.4A" "2.5A" "2.6A" "3.2" "4.1" "5" "reg"; do
    case $key in
        "2.4A") logfile="${LATEST_DIR}/2.4A_velocity.log" ;;
        "2.5A") logfile="${LATEST_DIR}/2.5A_klein_gordon.log" ;;
        "2.6A") logfile="${LATEST_DIR}/2.6A_energy.log" ;;
        "3.2") logfile="${LATEST_DIR}/3.2_vacuum_energy.log" ;;
        "4.1") logfile="${LATEST_DIR}/4.1_time_dilation.log" ;;
        "5") logfile="${LATEST_DIR}/5_em_coupling.log" ;;
        "reg") logfile="${LATEST_DIR}/regression_2.3.log" ;;
    esac
    
    if [ -f "$logfile" ] && tail -50 "$logfile" 2>/dev/null | grep -qE "(PASSED|FAILED|PARTIAL)"; then
        ELAPSED_EST=$((ELAPSED_EST + ${test_times[$key]}))
    fi
done

REMAINING_EST=$((18 - ELAPSED_EST))

echo "Progress: ${COMPLETED}/7 tests logged"
echo "Estimated elapsed: ~${ELAPSED_EST}h"
echo "Estimated remaining: ~${REMAINING_EST}h"
echo ""
echo "Expected completion: $(date -d "+${REMAINING_EST} hours" 2>/dev/null || date -v+${REMAINING_EST}H 2>/dev/null || echo "N/A")"

echo ""
echo "═══════════════════════════════════════════════════════════════════════"
echo "Refresh: watch -n 10 ./monitor_batch_b.sh"
echo "═══════════════════════════════════════════════════════════════════════"
