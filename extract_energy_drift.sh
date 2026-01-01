#!/bin/bash
# Extract energy drift from test console.log files

if [ -z "$1" ]; then
    echo "Usage: $1 <output_directory>"
    exit 1
fi

OUTPUT_DIR="$1"
CONSOLE_LOG="$OUTPUT_DIR/console.log"

if [ ! -f "$CONSOLE_LOG" ]; then
    echo "ERROR: Console log not found at $CONSOLE_LOG"
    exit 1
fi

# Extract energy drift percentage from validation report
grep -i "energy drift" "$CONSOLE_LOG" | tail -1 || echo "Energy drift: NOT FOUND"

# Extract final energy conservation metric
grep -i "energy conservation" "$CONSOLE_LOG" | tail -1 || echo "Energy conservation: NOT FOUND"

# Extract validation status
grep -i "VALIDATION:" "$CONSOLE_LOG" | tail -1 || echo "Validation status: NOT FOUND"
