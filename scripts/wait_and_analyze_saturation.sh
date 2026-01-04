#!/bin/bash
# Wait for saturation check completion and auto-analyze

LOG_FILE="/tmp/b1_saturation_256.log"
CSV_FILE="/home/persist/neotec/0rigin/analysis/b1_saturation_check_256.csv"

echo "Waiting for B1 Saturation Check completion..."
echo "CSV output: $CSV_FILE"
echo ""

# Wait for CSV file to appear
while [ ! -f "$CSV_FILE" ]; do
    # Check if process is still running
    if ! pgrep -f "trd --test.*saturation_check" > /dev/null; then
        echo "WARNING: Process not running but no CSV found!"
        echo "Check log for errors: tail -50 $LOG_FILE"
        exit 1
    fi

    # Show progress
    COMPLETED=$(grep -c "Testing separation" "$LOG_FILE" 2>/dev/null || echo 0)
    echo -ne "\rProgress: $COMPLETED / 6 separations completed..."
    sleep 10
done

echo ""
echo "✓ Saturation check complete! CSV generated."
echo ""

# Wait a moment for file to be fully written
sleep 2

# Display CSV contents
echo "=========================================="
echo "RAW RESULTS:"
echo "=========================================="
cat "$CSV_FILE"
echo ""

# Run analysis
echo "=========================================="
echo "RUNNING SATURATION ANALYSIS..."
echo "=========================================="
python3 scripts/analyze_b1_saturation.py "$CSV_FILE"

echo ""
echo "=========================================="
echo "DELIVERABLES:"
echo "=========================================="
echo "1. Raw data:        $CSV_FILE"
echo "2. Analysis report: (printed above)"
echo "3. Visualization:   analysis/b1_saturation_check_plot.png"
echo "4. Documentation:   docs/B1_SATURATION_ANALYSIS.md"
echo ""
echo "Next steps: Review verdict and update docs/B1_SATURATION_ANALYSIS.md with results"
