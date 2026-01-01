#!/bin/bash
# K-Parameter Systematic Scan: Measure B_max(K) from K=0.0 to K=2.0
# Critical research to resolve synchronization-EM paradox

set -e

OUTPUT_DIR="output/k_scan_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# K values to test
K_VALUES=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0)

echo "=========================================="
echo "K-Parameter Scan: Synchronization vs. EM"
echo "=========================================="
echo "Output: $OUTPUT_DIR"
echo ""

for K in "${K_VALUES[@]}"; do
    echo "Testing K = $K"

    # Create temporary config with current K value
    CONFIG="/tmp/k_scan_${K}.yaml"
    cp config/stuckelberg_integration_test.yaml "$CONFIG"

    # Update K parameter in config
    sed -i "s/K: .*/K: $K/" "$CONFIG"

    # Run test and extract final observables
    LOG="$OUTPUT_DIR/k_${K}.log"
    ./build/bin/smft --test "$CONFIG" > "$LOG" 2>&1

    # Extract final B_max and R_avg
    B_max=$(grep "EM_B_max" "$LOG" | tail -1 | awk '{print $NF}')
    R_avg=$(grep "R_avg" "$LOG" | tail -1 | awk '{print $NF}')

    echo "  K=$K → B_max=$B_max, R_avg=$R_avg"
    echo "$K,$B_max,$R_avg" >> "$OUTPUT_DIR/k_scan_results.csv"
done

echo ""
echo "Scan complete. Results: $OUTPUT_DIR/k_scan_results.csv"
echo ""
echo "Analysis:"
python3 - <<'EOF'
import pandas as pd
import sys

data = pd.read_csv("$OUTPUT_DIR/k_scan_results.csv", names=['K', 'B_max', 'R_avg'])

print("\nK-Parameter Scan Results:")
print(data.to_string(index=False))

print("\n" + "="*60)
print("Critical Questions:")
print("="*60)

# Check if any K>0 preserves EM
strong_em = data[data['B_max'] > 1.0]
if len(strong_em) > 0:
    print(f"✓ Strong EM (B_max>1.0) found at K={strong_em['K'].tolist()}")
else:
    print("✗ No K>0 preserves strong EM fields")

# Check synchronization threshold
synced = data[data['R_avg'] > 0.9]
if len(synced) > 0:
    print(f"✓ Synchronization (R>0.9) found at K={synced['K'].tolist()}")

# Check coexistence
coexist = data[(data['B_max'] > 1.0) & (data['R_avg'] > 0.9)]
if len(coexist) > 0:
    print(f"✅ COEXISTENCE at K={coexist['K'].tolist()}")
    print("   → Option 1 (abandon sync) NOT necessary")
else:
    print("❌ NO COEXISTENCE found")
    print("   → Proceed to Option 2 (dynamic K) or Option 3 (multi-component)")
EOF
