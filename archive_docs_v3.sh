#!/bin/bash
# Documentation Archive Script v3 - Robust Version
# Created: 2026-01-08
# Purpose: Move 122 scattered .md files to organized structure with timestamps

REPO_ROOT="/home/persist/neotec/0rigin"
cd "$REPO_ROOT"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "================================================================"
echo "TRD Documentation Archive v3"
echo "================================================================"
echo ""

# Create structure
echo -e "${YELLOW}Creating directories...${NC}"
mkdir -p docs/reports/validation/{wave-A,wave-B,wave-C,wave-D,wave-E,wave-F,wave-G,wave-H}
mkdir -p docs/reports/{implementation,architecture}
mkdir -p docs/{analysis,archive,archive/duplicates}
echo -e "${GREEN}✓ Done${NC}"
echo ""

# Initialize map
cat > docs/ARCHIVE_MAP.md << 'EOF'
# Documentation Archive Map
Generated: 2026-01-08

Files archived with creation timestamp: `YYYYMMDD_HHMMSS_originalname.md`

## Archive Structure
- `docs/reports/validation/wave-X/` - Wave validation reports
- `docs/reports/implementation/` - Implementation summaries
- `docs/reports/architecture/` - Architecture reports
- `docs/analysis/` - Agent-generated analysis
- `docs/archive/` - Miscellaneous
- `docs/archive/duplicates/` - Duplicate versions

## Mappings

EOF

# Keep these in root
declare -A KEEP
KEEP["README.md"]=1
KEEP["CLAUDE.md"]=1
KEEP["TODO.md"]=1
KEEP["CONTRIBUTING.md"]=1

# Get timestamp helper
get_ts() {
    date -d "@$(stat -c '%Y' "$1")" '+%Y%m%d_%H%M%S' 2>/dev/null || echo "19700101_000000"
}

# Archive helper
archive() {
    local src="$1"
    local dst_dir="$2"
    local prefix="${3:-}"

    [[ ! -f "$src" ]] && return
    [[ "${KEEP[$src]}" == "1" ]] && return

    local ts=$(get_ts "$src")
    local name=$(basename "$src")
    local dst="$dst_dir/${ts}_${prefix}${name}"

    echo "  $src → $dst"
    mv "$src" "$dst"
    echo "- \`$src\` → \`$dst\`" >> docs/ARCHIVE_MAP.md
}

count=0

echo -e "${YELLOW}Wave A (A*)...${NC}"
echo "### Wave A" >> docs/ARCHIVE_MAP.md
for f in A*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ ^ARCHITECTURE ]] && { archive "$f" "docs/reports/validation/wave-A"; ((count++)); }
done

echo -e "${YELLOW}Wave B (B*)...${NC}"
echo "### Wave B" >> docs/ARCHIVE_MAP.md
for f in B*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ ^(BINARY|BUILD) ]] && { archive "$f" "docs/reports/validation/wave-B"; ((count++)); }
done

echo -e "${YELLOW}Wave C (C*)...${NC}"
echo "### Wave C" >> docs/ARCHIVE_MAP.md
for f in C*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ ^(CLAUDE|CONTRIBUTING|CSV) ]] && { archive "$f" "docs/reports/validation/wave-C"; ((count++)); }
done

echo -e "${YELLOW}Wave D (D*)...${NC}"
echo "### Wave D" >> docs/ARCHIVE_MAP.md
# Handle D1 duplicates
if [[ -f "D1_EXPERIMENTAL_PREDICTIONS_SUMMARY.md" ]]; then
    archive "D1_EXPERIMENTAL_PREDICTIONS_SUMMARY.md" "docs/reports/validation/wave-D"
    ((count++))
fi
for f in D1_EXPERIMENTAL_PREDICTIONS*.md; do
    [[ -f "$f" ]] && [[ "$f" != "D1_EXPERIMENTAL_PREDICTIONS_SUMMARY.md" ]] && { archive "$f" "docs/archive/duplicates" "dup_"; ((count++)); }
done
for f in D[2-9]*.md; do
    [[ -f "$f" ]] && { archive "$f" "docs/reports/validation/wave-D"; ((count++)); }
done

echo -e "${YELLOW}Wave E (E*)...${NC}"
echo "### Wave E" >> docs/ARCHIVE_MAP.md
# Handle E1 duplicates - keep COMPREHENSIVE
if [[ -f "E1_RENORMALIZABILITY_COMPREHENSIVE_REPORT.md" ]]; then
    archive "E1_RENORMALIZABILITY_COMPREHENSIVE_REPORT.md" "docs/reports/validation/wave-E"
    ((count++))
fi
for f in E1_RENORMALIZABILITY*.md; do
    [[ -f "$f" ]] && [[ "$f" != "E1_RENORMALIZABILITY_COMPREHENSIVE_REPORT.md" ]] && { archive "$f" "docs/archive/duplicates" "dup_"; ((count++)); }
done
# Handle E3 duplicates - keep COMPREHENSIVE
if [[ -f "E3_CAUSALITY_COMPREHENSIVE_REPORT.md" ]]; then
    archive "E3_CAUSALITY_COMPREHENSIVE_REPORT.md" "docs/reports/validation/wave-E"
    ((count++))
fi
for f in E3_CAUSALITY*.md; do
    [[ -f "$f" ]] && [[ "$f" != "E3_CAUSALITY_COMPREHENSIVE_REPORT.md" ]] && { archive "$f" "docs/archive/duplicates" "dup_"; ((count++)); }
done
# Other E files
for f in E[1-9]*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ (EXPERIMENTAL|EINSTEIN) ]] && ! [[ "$f" =~ ^E[13]_.*CAUSALITY ]] && ! [[ "$f" =~ ^E1_RENORMALIZABILITY ]] && { archive "$f" "docs/reports/validation/wave-E"; ((count++)); }
done

echo -e "${YELLOW}Wave F (F*)...${NC}"
echo "### Wave F" >> docs/ARCHIVE_MAP.md
for f in F[2-9]*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ ^FIELD ]] && ! [[ "$f" =~ ^FINAL ]] && { archive "$f" "docs/reports/validation/wave-F"; ((count++)); }
done

echo -e "${YELLOW}Wave G (G*)...${NC}"
echo "### Wave G" >> docs/ARCHIVE_MAP.md
for f in G[4-9]*.md; do
    [[ -f "$f" ]] && ! [[ "$f" =~ ^(GOLDEN|GPU) ]] && { archive "$f" "docs/reports/validation/wave-G"; ((count++)); }
done

echo -e "${YELLOW}Wave H (H*)...${NC}"
echo "### Wave H" >> docs/ARCHIVE_MAP.md
for f in H[1-9]*.md; do
    [[ -f "$f" ]] && { archive "$f" "docs/reports/validation/wave-H"; ((count++)); }
done

echo -e "${YELLOW}Implementation reports...${NC}"
echo "### Implementation" >> docs/ARCHIVE_MAP.md
for f in *IMPLEMENTATION*.md *COMPLETE*.md *SUMMARY*.md PHASE*.md; do
    [[ -f "$f" ]] && ! [[ "${KEEP[$f]}" == "1" ]] && { archive "$f" "docs/reports/implementation"; ((count++)); }
done

echo -e "${YELLOW}Architecture reports...${NC}"
echo "### Architecture" >> docs/ARCHIVE_MAP.md
for f in *CLEANUP*.md *CONSOLIDATION*.md SYMPLECTIC*.md ARCHITECTURE*.md BUILD_SYSTEM*.md *MIGRATION*.md; do
    [[ -f "$f" ]] && { archive "$f" "docs/reports/architecture"; ((count++)); }
done

echo -e "${YELLOW}Analysis reports...${NC}"
echo "### Analysis" >> docs/ARCHIVE_MAP.md
for f in TRD_*.md *ANALYSIS*.md *DIAGNOSTIC*.md *INVESTIGATION*.md *VALIDATION_REPORT*.md *VERIFICATION*.md; do
    [[ -f "$f" ]] && ! [[ "${KEEP[$f]}" == "1" ]] && { archive "$f" "docs/analysis"; ((count++)); }
done

echo -e "${YELLOW}Miscellaneous...${NC}"
echo "### Miscellaneous" >> docs/ARCHIVE_MAP.md
for f in *.md; do
    [[ -f "$f" ]] && ! [[ "${KEEP[$f]}" == "1" ]] && { archive "$f" "docs/archive"; ((count++)); }
done

echo ""
echo "================================================================"
echo -e "${GREEN}Complete!${NC}"
echo "================================================================"
echo "Files archived: $count"
echo "Remaining root .md files: $(ls -1 *.md 2>/dev/null | wc -l)"
echo ""
echo -e "${YELLOW}Root files:${NC}"
ls -1 *.md 2>/dev/null
echo ""
echo "Archive map: docs/ARCHIVE_MAP.md"
echo ""
echo "Next: git add docs/ && git commit -m 'docs: Archive scattered documentation'"
