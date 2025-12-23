#!/bin/bash
# Archive Legacy Outputs Script
# Archives 1.7GB of pre-automation outputs before clean rerun
# Created: 2025-12-22
# Part of: Clean Validation Campaign (Phase: phase_1766476820454_rllxro102)

set -e  # Exit on error

ARCHIVE_DATE=$(date +%Y%m%d_%H%M%S)
ARCHIVE_DIR="archive/pre_automation_${ARCHIVE_DATE}"
COMPRESS_NAME="pre_automation_${ARCHIVE_DATE}.tar.gz"

echo "========================================="
echo "SMFT Legacy Output Archival"
echo "========================================="
echo ""
echo "Archive directory: ${ARCHIVE_DIR}"
echo "Compressed archive: archive/${COMPRESS_NAME}"
echo ""

# Create archive directory structure
echo "[1/6] Creating archive directory..."
mkdir -p "archive/${ARCHIVE_DIR}/outputs"
mkdir -p "archive/${ARCHIVE_DIR}/scripts"
mkdir -p "archive/${ARCHIVE_DIR}/reports"

# Archive outputs (1.7GB, 155 directories)
echo "[2/6] Archiving output directories (1.7GB, this may take a few minutes)..."
if [ -d "output" ]; then
    mv output/ "archive/${ARCHIVE_DIR}/outputs/"
    echo "  ✓ Moved output/ (155 directories)"
else
    echo "  ⚠ Warning: output/ directory not found, skipping"
fi

# Archive temporary Python analysis scripts
echo "[3/6] Archiving Python analysis scripts..."
if ls analyze_*.py 1> /dev/null 2>&1; then
    mv analyze_*.py "archive/${ARCHIVE_DIR}/scripts/"
    echo "  ✓ Moved analyze_*.py scripts"
else
    echo "  ⚠ No analyze_*.py scripts found"
fi

if ls visualize_*.py 1> /dev/null 2>&1; then
    mv visualize_*.py "archive/${ARCHIVE_DIR}/scripts/"
    echo "  ✓ Moved visualize_*.py scripts"
else
    echo "  ⚠ No visualize_*.py scripts found"
fi

# Archive old phase reports (mark as pre-automation)
echo "[4/6] Archiving pre-automation reports..."
if ls docs/PHASE_2.3_*.md 1> /dev/null 2>&1; then
    cp docs/PHASE_2.3_*.md "archive/${ARCHIVE_DIR}/reports/" 2>/dev/null || true
    echo "  ✓ Copied Phase 2.3 reports"
fi
if ls docs/PHASE_2.4_*.md 1> /dev/null 2>&1; then
    cp docs/PHASE_2.4_*.md "archive/${ARCHIVE_DIR}/reports/" 2>/dev/null || true
    echo "  ✓ Copied Phase 2.4 reports"
fi

# Create archive README
echo "[5/6] Creating archive documentation..."
cat > "archive/${ARCHIVE_DIR}/README.md" << EOF
# SMFT Legacy Outputs Archive

**Archived**: $(date)
**Reason**: Clean rerun with automated validation + C++ plotting system

## Contents

### outputs/
Legacy output directories (155 total, 1.7GB):
- Timestamped directories (20251218_*, 20251219_*, etc.)
- Mixed manual validation results
- Inconsistent naming conventions
- Python-script dependent analysis

### scripts/
Temporary Python analysis/visualization scripts:
- analyze_*.py - Manual validation scripts
- visualize_*.py - Manual plotting scripts
- **Note**: Replaced by C++ automated validation + plotting

### reports/
Pre-automation phase reports:
- PHASE_2.3_*.md - Manual analysis summaries
- PHASE_2.4_*.md - Investigation reports
- **Note**: New reports generated with automated system

## Why Archived

The project has transitioned to:
1. **Automated validation** - ScenarioValidator integrated into ./smft
2. **C++ plotting** - matplotlib-cpp for publication-quality plots
3. **Semantic organization** - Phase-based structure, not timestamps
4. **Consistent workflow** - Single executable, no manual scripts

## Restoration (if needed)

```bash
cd archive
tar -xzf $(basename "$COMPRESS_NAME")
# Outputs will be in pre_automation_*/outputs/
```
EOF

# Compress archive
echo "[6/6] Compressing archive..."
cd archive
if [ -d "pre_automation_${ARCHIVE_DATE}" ]; then
    tar -czf "${COMPRESS_NAME}" "pre_automation_${ARCHIVE_DATE}/"
    COMPRESSED_SIZE=$(du -h "${COMPRESS_NAME}" | cut -f1)
    echo "  ✓ Created ${COMPRESS_NAME} (${COMPRESSED_SIZE})"

    # Remove uncompressed archive directory
    rm -rf "pre_automation_${ARCHIVE_DATE}/"
    echo "  ✓ Removed uncompressed directory"
else
    echo "  ⚠ Warning: Archive directory not found for compression"
fi
cd ..

# Recreate clean output directory
echo ""
echo "[Cleanup] Creating clean output directory structure..."
mkdir -p output

echo ""
echo "========================================="
echo "Archive Complete!"
echo "========================================="
echo ""
echo "Summary:"
echo "  Archive location: archive/${COMPRESS_NAME}"
echo "  Compressed size: ${COMPRESSED_SIZE}"
echo "  Original size: ~1.7GB"
echo "  Compression ratio: ~3-4x"
echo ""
echo "✓ Ready for clean validation campaign"
echo ""
