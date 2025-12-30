#!/bin/bash
# Universality Class Finite-Size Scaling Analysis
# Comprehensive scan across multiple grid sizes to extract critical exponents

echo "=========================================="
echo "SMFT UNIVERSALITY CLASS ANALYSIS"
echo "=========================================="
echo ""
echo "This will run FSS analysis with:"
echo "  - Grid sizes: L = 32, 64, 128, 256, 512"
echo "  - Noise points: 41 values from σ = 0.0 to 2.0"
echo "  - Total simulations: 5 × 41 = 205"
echo ""
echo "Estimated runtime: ~10 hours (parallelizable)"
echo ""

# Check if SMFT executable exists
if [ ! -f "../../build/bin/smft" ]; then
    echo "ERROR: SMFT executable not found!"
    echo "Please build the project first: cd ../.. && cmake --build build"
    exit 1
fi

# Create output directory
OUTPUT_BASE="../../output/universality_scan"
mkdir -p $OUTPUT_BASE

# Function to run a single grid size
run_grid_size() {
    local L=$1
    echo ""
    echo "=========================================="
    echo "Running L = $L scan..."
    echo "=========================================="

    CONFIG="universality_L${L}.yaml"

    if [ ! -f "$CONFIG" ]; then
        echo "ERROR: Config file $CONFIG not found!"
        return 1
    fi

    # Run the simulation
    ../../build/bin/smft --test $CONFIG

    if [ $? -eq 0 ]; then
        echo "✓ L = $L completed successfully"
    else
        echo "✗ L = $L failed!"
        return 1
    fi
}

# Option to run in parallel
if [ "$1" == "--parallel" ]; then
    echo "Running in PARALLEL mode..."
    echo ""

    # Run all grid sizes in parallel
    for L in 32 64 128 256 512; do
        run_grid_size $L &
    done

    # Wait for all background jobs
    wait

else
    echo "Running in SEQUENTIAL mode..."
    echo "(Use --parallel flag for faster execution)"
    echo ""

    # Run each grid size sequentially
    for L in 32 64 128 256 512; do
        run_grid_size $L
        if [ $? -ne 0 ]; then
            echo "Stopping due to error in L = $L"
            exit 1
        fi
    done
fi

echo ""
echo "=========================================="
echo "ALL SIMULATIONS COMPLETE!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Analyze data:"
echo "   python3 ../../analysis/universality_analysis/data_collapse.py $OUTPUT_BASE"
echo ""
echo "2. Extract exponents:"
echo "   python3 ../../analysis/universality_analysis/exponent_extraction.py $OUTPUT_BASE"
echo ""
echo "3. Classify universality:"
echo "   python3 ../../analysis/universality_analysis/classification.py $OUTPUT_BASE/critical_exponents.txt"
echo ""