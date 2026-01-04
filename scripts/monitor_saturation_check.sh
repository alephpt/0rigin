#!/bin/bash
# Monitor B1 Saturation Check Progress

LOG_FILE="/tmp/b1_saturation_256.log"
CSV_FILE="/home/persist/neotec/0rigin/analysis/b1_saturation_check_256.csv"

echo "==========================================="
echo "B1 SATURATION CHECK MONITOR"
echo "==========================================="
echo ""

# Check if process is running
PID=$(pgrep -f "trd --test.*saturation_check" 2>/dev/null)
if [ -n "$PID" ]; then
    echo "Status: RUNNING (PID $PID)"

    # Extract latest separation being tested
    CURRENT=$(tail -100 "$LOG_FILE" | grep "Testing separation" | tail -1)
    echo "Current: $CURRENT"

    # Count completed separations
    COMPLETED=$(grep -c "Testing separation" "$LOG_FILE")
    echo "Progress: $COMPLETED / 6 separations"
    echo ""

    # Show latest results
    echo "Latest Results:"
    tail -50 "$LOG_FILE" | grep -E "(separation d=|m1=|m2_m1=)" | tail -5

    echo ""
    echo "Estimated time remaining: $(( (6 - COMPLETED) * 25 )) minutes"

else
    echo "Status: COMPLETE or NOT RUNNING"

    # Check for CSV output
    if [ -f "$CSV_FILE" ]; then
        echo ""
        echo "Results available in: $CSV_FILE"
        echo ""
        echo "Run analysis:"
        echo "  python3 scripts/analyze_b1_saturation.py"
    else
        echo ""
        echo "No CSV output found. Check log for errors:"
        echo "  tail -100 $LOG_FILE"
    fi
fi

echo ""
echo "==========================================="
echo "Monitor commands:"
echo "  Watch log:     tail -f $LOG_FILE"
echo "  Check status:  ps aux | grep saturation_check"
echo "  Kill process:  pkill -f saturation_check"
echo "==========================================="
