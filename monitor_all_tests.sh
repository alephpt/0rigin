#!/bin/bash
# Monitor grid convergence test (single test runs all grid sizes)

echo "========================================="
echo "Grid Convergence Test Monitor"
echo "========================================="
echo ""

log_file="grid_convergence_16_32_64.log"

if [ ! -f "$log_file" ]; then
    echo "  Status: NOT STARTED (log file missing: $log_file)"
    echo ""
    echo "To start the test:"
    echo "  build/bin/smft --test config/defect_localization_validation.yaml > grid_convergence_16_32_64.log 2>&1 &"
    exit 1
fi

echo "Log file: $log_file"
echo ""

# Check if completed
if grep -q "Grid Convergence Tests Complete" "$log_file" 2>/dev/null; then
    echo "  Status: ✓ COMPLETED"
    echo ""

    # Show which grid sizes were tested
    if grep -q "Testing grid sizes:" "$log_file"; then
        grid_sizes=$(grep "Testing grid sizes:" "$log_file" | head -1 | sed 's/.*: //')
        echo "  Grid sizes tested: $grid_sizes"
    fi

    # Show summary for each grid size
    for grid in 16 32 64; do
        if grep -q "Grid Size: ${grid}×${grid}" "$log_file"; then
            echo ""
            echo "  === ${grid}×${grid} Grid ==="

            # Check if this grid size completed
            # Look for completion between this grid header and the next
            grid_section=$(sed -n "/Grid Size: ${grid}×${grid}/,/Grid Size: [0-9]/p" "$log_file")

            if echo "$grid_section" | grep -q "Convergence validation PASSED"; then
                echo "    Convergence: ✓ PASS"
            elif echo "$grid_section" | grep -q "Convergence validation FAILED"; then
                echo "    Convergence: ✗ FAIL"
            fi

            # Show grid spacing info
            if echo "$grid_section" | grep -q "Grid spacing:"; then
                spacing=$(echo "$grid_section" | grep "Grid spacing:" | head -1 | awk '{print $3" "$4}')
                core=$(echo "$grid_section" | grep "Core radius:" | head -1 | grep -oP '\([0-9.]+ grid points\)')
                echo "    Grid spacing: $spacing"
                echo "    Core radius: $core"
            fi
        fi
    done

else
    echo "  Status: RUNNING"
    echo ""

    # Determine which grid size is currently running
    last_grid=$(grep -oP "Grid Size: \K[0-9]+×[0-9]+" "$log_file" | tail -1)
    if [ -n "$last_grid" ]; then
        echo "  Current grid: $last_grid"
    fi

    # Get latest progress (Step number)
    latest_step=$(grep -oP "Step \K[0-9]+" "$log_file" 2>/dev/null | tail -1)
    if [ -n "$latest_step" ]; then
        total_steps=10000
        progress=$((latest_step * 100 / total_steps))
        echo "  Progress: $latest_step / $total_steps ($progress%)"

        # Estimate which N value is running (each grid runs N=1,10,100)
        steps_per_n=10000
        n_index=$((latest_step / steps_per_n))
        case $n_index in
            0) echo "  Running: N=1" ;;
            1) echo "  Running: N=10" ;;
            2) echo "  Running: N=100" ;;
        esac
    fi
fi

echo ""
echo "========================================="
echo ""

if grep -q "Grid Convergence Tests Complete" "$log_file" 2>/dev/null; then
    echo "✓ Grid convergence test complete!"
    echo ""
    echo "Output directories created:"
    ls -dt output/20*_defect_localization_*x* 2>/dev/null | head -10 || echo "  (check output/ directory)"
    echo ""
    echo "Next steps:"
    echo "  1. Analyze results across grid sizes"
    echo "  2. Check if R_var(t=0) converged"
    echo "  3. Check if dynamics (ΔR_var) converged"
else
    echo "⏳ Test still running. Re-run this script to check progress."
    echo ""
    echo "Monitor log:"
    echo "  tail -f $log_file"
fi

echo ""
