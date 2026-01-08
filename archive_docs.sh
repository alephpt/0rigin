#!/bin/bash
# Documentation Archive Script v2
# Created: 2026-01-08
# Purpose: Move 122 scattered .md files to organized structure with timestamps

set -e

REPO_ROOT="/home/persist/neotec/0rigin"
cd "$REPO_ROOT"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "================================================================"
echo "TRD Documentation Archive - Timestamp Preservation v2"
echo "================================================================"
echo ""

# Create archive structure
echo -e "${YELLOW}Creating archive directory structure...${NC}"
mkdir -p docs/reports/validation/wave-{A,B,C,D,E,F,G,H}
mkdir -p docs/reports/implementation
mkdir -p docs/reports/architecture
mkdir -p docs/analysis
mkdir -p docs/archive
mkdir -p docs/archive/duplicates

echo -e "${GREEN}✓ Directory structure created${NC}"
echo ""

# Initialize archive map
ARCHIVE_MAP="docs/ARCHIVE_MAP.md"
cat > "$ARCHIVE_MAP" << 'EOF'
# Documentation Archive Map
Generated: 2026-01-08

This file maps original root-level .md files to their archived locations with timestamps.

## Archive Structure

- `docs/reports/validation/wave-X/` - Wave validation reports (A1-H3)
- `docs/reports/implementation/` - Implementation summaries
- `docs/reports/architecture/` - Architecture and cleanup reports
- `docs/analysis/` - Agent-generated analysis reports
- `docs/archive/` - Miscellaneous documentation
- `docs/archive/duplicates/` - Duplicate file versions

## Timestamp Format

Files archived with creation timestamp: `YYYYMMDD_HHMMSS_originalname.md`

## File Mappings

EOF

# Function to get timestamp and format
get_timestamp() {
    local file="$1"
    local epoch=$(stat -c '%Y' "$file" 2>/dev/null || echo "0")
    date -d "@$epoch" '+%Y%m%d_%H%M%S' 2>/dev/null || echo "19700101_000000"
}

# Function to archive file (use mv not git mv to avoid staging issues)
archive_file() {
    local file="$1"
    local dest_dir="$2"
    local prefix="${3:-}"

    if [[ ! -f "$file" ]]; then
        return
    fi

    local timestamp=$(get_timestamp "$file")
    local basename=$(basename "$file")
    local new_name="${timestamp}_${basename}"

    if [[ -n "$prefix" ]]; then
        new_name="${timestamp}_${prefix}_${basename}"
    fi

    local dest_path="${dest_dir}/${new_name}"

    echo "  $file → $dest_path"
    mv "$file" "$dest_path"

    # Add to archive map
    echo "- \`$file\` → \`$dest_path\`" >> "$ARCHIVE_MAP"
}

# Counters
count_wave_a=0
count_wave_b=0
count_wave_c=0
count_wave_d=0
count_wave_e=0
count_wave_f=0
count_wave_g=0
count_wave_h=0
count_impl=0
count_arch=0
count_analysis=0
count_archive=0
count_duplicates=0

# Keep list of files to preserve in root
KEEP_FILES=(
    "README.md"
    "CLAUDE.md"
    "TODO.md"
    "CONTRIBUTING.md"
)

# Function to check if file should be kept
should_keep() {
    local file="$1"
    for keep in "${KEEP_FILES[@]}"; do
        if [[ "$file" == "$keep" ]]; then
            return 0
        fi
    done
    return 1
}

echo -e "${YELLOW}Archiving Wave A validation reports...${NC}"
echo "### Wave A Validation Reports" >> "$ARCHIVE_MAP"
for file in A1_*.md A2_*.md A3_*.md A4_*.md A5_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-A"
        ((count_wave_a++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave B validation reports...${NC}"
echo "### Wave B Validation Reports" >> "$ARCHIVE_MAP"
for file in B1_*.md B2_*.md B3_*.md B4_*.md B5_*.md B6_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-B"
        ((count_wave_b++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave C validation reports...${NC}"
echo "### Wave C Validation Reports" >> "$ARCHIVE_MAP"
for file in C1_*.md C2_*.md C3_*.md C4_*.md C5_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-C"
        ((count_wave_c++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave D validation reports...${NC}"
echo "### Wave D Validation Reports" >> "$ARCHIVE_MAP"
for file in D1_*.md D2_*.md D3_*.md D4_*.md D5_*.md; do
    if [[ -f "$file" ]]; then
        # Handle D1 duplicates - keep SUMMARY, archive others
        if [[ "$file" == "D1_EXPERIMENTAL_PREDICTIONS_SUMMARY.md" ]]; then
            archive_file "$file" "docs/reports/validation/wave-D"
        elif [[ "$file" == D1_EXPERIMENTAL_PREDICTIONS*.md ]]; then
            archive_file "$file" "docs/archive/duplicates" "dup"
            ((count_duplicates++))
        else
            archive_file "$file" "docs/reports/validation/wave-D"
        fi
        ((count_wave_d++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave E validation reports...${NC}"
echo "### Wave E Validation Reports" >> "$ARCHIVE_MAP"
for file in E1_*.md E2_*.md E3_*.md E4_*.md E5_*.md; do
    if [[ -f "$file" ]]; then
        # Handle E1 duplicates - keep COMPREHENSIVE
        if [[ "$file" == "E1_RENORMALIZABILITY_COMPREHENSIVE.md" ]]; then
            archive_file "$file" "docs/reports/validation/wave-E"
        elif [[ "$file" == E1_RENORMALIZABILITY*.md ]]; then
            archive_file "$file" "docs/archive/duplicates" "dup"
            ((count_duplicates++))
        # Handle E3 duplicates - keep COMPREHENSIVE
        elif [[ "$file" == "E3_CAUSALITY_COMPREHENSIVE.md" ]]; then
            archive_file "$file" "docs/reports/validation/wave-E"
        elif [[ "$file" == E3_CAUSALITY*.md ]]; then
            archive_file "$file" "docs/archive/duplicates" "dup"
            ((count_duplicates++))
        else
            archive_file "$file" "docs/reports/validation/wave-E"
        fi
        ((count_wave_e++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave F validation reports...${NC}"
echo "### Wave F Validation Reports" >> "$ARCHIVE_MAP"
for file in F2_*.md F3_*.md F4_*.md F5_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-F"
        ((count_wave_f++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave G validation reports...${NC}"
echo "### Wave G Validation Reports" >> "$ARCHIVE_MAP"
for file in G4_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-G"
        ((count_wave_g++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving Wave H validation reports...${NC}"
echo "### Wave H Validation Reports" >> "$ARCHIVE_MAP"
for file in H1_*.md H2_*.md H3_*.md; do
    if [[ -f "$file" ]]; then
        archive_file "$file" "docs/reports/validation/wave-H"
        ((count_wave_h++))
    fi
done
echo ""

echo -e "${YELLOW}Archiving implementation summaries...${NC}"
echo "### Implementation Summaries" >> "$ARCHIVE_MAP"
for pattern in "*IMPLEMENTATION*.md" "*COMPLETE.md" "*SUMMARY*.md" "PHASE*.md"; do
    for file in $pattern; do
        if [[ -f "$file" ]] && ! should_keep "$file"; then
            archive_file "$file" "docs/reports/implementation"
            ((count_impl++))
        fi
    done
done
echo ""

echo -e "${YELLOW}Archiving architecture reports...${NC}"
echo "### Architecture Reports" >> "$ARCHIVE_MAP"
for pattern in "*CLEANUP*.md" "*CONSOLIDATION*.md" "SYMPLECTIC*.md" "*ARCHITECTURE*.md" "BUILD_SYSTEM*.md" "MIGRATION*.md"; do
    for file in $pattern; do
        if [[ -f "$file" ]]; then
            archive_file "$file" "docs/reports/architecture"
            ((count_arch++))
        fi
    done
done
echo ""

echo -e "${YELLOW}Archiving analysis reports...${NC}"
echo "### Analysis Reports" >> "$ARCHIVE_MAP"
for pattern in "TRD_CODEBASE_ANALYSIS*.md" "TRD_VALIDATION_REPORT*.md" "*ANALYSIS*.md" "*DIAGNOSTIC*.md" "*INVESTIGATION*.md" "*REPORT*.md"; do
    for file in $pattern; do
        if [[ -f "$file" ]]; then
            archive_file "$file" "docs/analysis"
            ((count_analysis++))
        fi
    done
done
echo ""

echo -e "${YELLOW}Archiving miscellaneous documentation...${NC}"
echo "### Miscellaneous Documentation" >> "$ARCHIVE_MAP"
# Archive remaining .md files (excluding formal docs)
for file in *.md; do
    if [[ -f "$file" ]] && ! should_keep "$file"; then
        archive_file "$file" "docs/archive"
        ((count_archive++))
    fi
done
echo ""

# Summary
total_moved=$((count_wave_a + count_wave_b + count_wave_c + count_wave_d + count_wave_e + count_wave_f + count_wave_g + count_wave_h + count_impl + count_arch + count_analysis + count_archive))

echo "================================================================"
echo -e "${GREEN}Archive Complete${NC}"
echo "================================================================"
echo ""
echo "Files Moved by Category:"
echo "  Wave A:           $count_wave_a"
echo "  Wave B:           $count_wave_b"
echo "  Wave C:           $count_wave_c"
echo "  Wave D:           $count_wave_d"
echo "  Wave E:           $count_wave_e"
echo "  Wave F:           $count_wave_f"
echo "  Wave G:           $count_wave_g"
echo "  Wave H:           $count_wave_h"
echo "  Implementation:   $count_impl"
echo "  Architecture:     $count_arch"
echo "  Analysis:         $count_analysis"
echo "  Archive:          $count_archive"
echo "  Duplicates:       $count_duplicates"
echo "  ------------------------"
echo "  TOTAL:            $total_moved"
echo ""

remaining=$(ls -1 *.md 2>/dev/null | wc -l)
echo "Remaining root .md files: $remaining"
echo ""
echo -e "${YELLOW}Files kept in root:${NC}"
ls -1 *.md 2>/dev/null || echo "No markdown files remaining in root"
echo ""
echo "Archive map: docs/ARCHIVE_MAP.md"
echo ""

echo -e "${GREEN}Done!${NC}"
echo ""
echo "Next Steps:"
echo "1. Review docs/ARCHIVE_MAP.md for full mapping"
echo "2. Verify remaining root files are appropriate"
echo "3. Stage archive: git add docs/ *.md"
echo "4. Commit: git commit -m 'docs: Archive scattered documentation with timestamps'"
