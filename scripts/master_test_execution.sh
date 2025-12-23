#!/bin/bash
# Master Test Execution Script
# Comprehensive SMFT validation campaign with automated validation + C++ plotting
# Created: 2025-12-22
# Part of: Clean Validation Campaign (Phase: phase_1766476820454_rllxro102)

set -e  # Exit on error

# Configuration
LOG_DIR="logs/test_execution"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MASTER_LOG="${LOG_DIR}/master_execution_${TIMESTAMP}.log"
EMAIL_NOTIFICATIONS=false  # Set to true if you want email alerts

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Create log directory
mkdir -p "${LOG_DIR}"

# Logging function
log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "${MASTER_LOG}"
}

# Error handler
error_exit() {
    log "${RED}ERROR: $1${NC}"
    exit 1
}

# Runtime estimation function
estimate_time() {
    local hours=$1
    local days=$((hours / 24))
    local remaining_hours=$((hours % 24))
    if [ $days -gt 0 ]; then
        echo "${days}d ${remaining_hours}h"
    else
        echo "${hours}h"
    fi
}

# Test execution function
run_test() {
    local test_name=$1
    local config_file=$2
    local est_runtime=$3

    log "${BLUE}=========================================${NC}"
    log "${BLUE}Starting: ${test_name}${NC}"
    log "${BLUE}Config: ${config_file}${NC}"
    log "${BLUE}Estimated runtime: ${est_runtime}${NC}"
    log "${BLUE}=========================================${NC}"

    local start_time=$(date +%s)

    # Run test and capture output
    if ./build/bin/smft --test "${config_file}" 2>&1 | tee "${LOG_DIR}/${test_name}_${TIMESTAMP}.log"; then
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        local elapsed_min=$((elapsed / 60))

        log "${GREEN}✓ ${test_name} COMPLETE (${elapsed_min} minutes)${NC}"
        return 0
    else
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        local elapsed_min=$((elapsed / 60))

        log "${RED}✗ ${test_name} FAILED (${elapsed_min} minutes)${NC}"
        return 1
    fi
}

# Header
clear
echo "========================================="
echo "SMFT Master Test Execution"
echo "========================================="
echo ""
echo "Campaign: Clean Validation with Automated Framework"
echo "Start time: $(date)"
echo "Log file: ${MASTER_LOG}"
echo ""

log "========================================="
log "SMFT Master Test Execution Starting"
log "========================================="

# =========================================
# BATCH 1: Phase 2.4 - Breakdown Investigation
# =========================================
echo ""
echo "========================================="
echo "BATCH 1: Phase 2.4 Breakdown Investigation"
echo "Estimated total runtime: 30 hours"
echo "========================================="
echo ""

BATCH1_START=$(date +%s)

# Phase 2.4A: Velocity Threshold (2h)
run_test "phase_2.4A_velocity_threshold" \
         "config/scenario_2.4A_velocity_threshold.yaml" \
         "2h" || true

# Phase 2.4B: R-field Dynamics (4h)
run_test "phase_2.4B_R_field_dynamics" \
         "config/scenario_2.4B_R_field_dynamics.yaml" \
         "4h" || true

# Phase 2.4C: Ultra-Relativistic (24h)
echo ""
echo "${YELLOW}WARNING: Phase 2.4C will run for ~24 hours (512×512 grid)${NC}"
echo "Consider running overnight or over weekend"
read -p "Proceed with 2.4C? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    run_test "phase_2.4C_ultra_relativistic" \
             "config/scenario_2.4C_ultra_relativistic.yaml" \
             "24h" || true
else
    log "${YELLOW}⊘ Phase 2.4C skipped by user${NC}"
fi

BATCH1_END=$(date +%s)
BATCH1_ELAPSED=$((($BATCH1_END - $BATCH1_START) / 3600))
log "${GREEN}Batch 1 complete in ${BATCH1_ELAPSED} hours${NC}"

# =========================================
# BATCH 2: Phase 2.3 - Clean Rerun
# =========================================
echo ""
echo "========================================="
echo "BATCH 2: Phase 2.3 Clean Rerun"
echo "Estimated runtime: 12 hours"
echo "========================================="
echo ""

read -p "Continue with Batch 2 (Phase 2.3 rerun)? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    BATCH2_START=$(date +%s)

    run_test "phase_2.3_relativistic_mass_full" \
             "config/phase_2.3_full_validation.yaml" \
             "12h" || true

    BATCH2_END=$(date +%s)
    BATCH2_ELAPSED=$((($BATCH2_END - $BATCH2_START) / 3600))
    log "${GREEN}Batch 2 complete in ${BATCH2_ELAPSED} hours${NC}"
else
    log "${YELLOW}⊘ Batch 2 skipped by user${NC}"
fi

# =========================================
# BATCH 3: Phase 2.6 - Advanced Physics
# =========================================
echo ""
echo "========================================="
echo "BATCH 3: Phase 2.6 Advanced Physics"
echo "Estimated runtime: 45 hours"
echo "========================================="
echo ""

read -p "Continue with Batch 3 (Phase 2.6)? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    BATCH3_START=$(date +%s)

    # Phase 2.6A: Energy-Momentum Relation (15h)
    run_test "phase_2.6A_energy_momentum" \
             "config/scenario_2.6A_energy_momentum_relation.yaml" \
             "15h" || true

    # Phase 2.6C: Soliton Stability (30h)
    echo ""
    echo "${YELLOW}WARNING: Phase 2.6C will run for ~30 hours (long evolution)${NC}"
    read -p "Proceed with 2.6C? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        run_test "phase_2.6C_soliton_stability" \
                 "config/scenario_2.6C_soliton_stability.yaml" \
                 "30h" || true
    else
        log "${YELLOW}⊘ Phase 2.6C skipped by user${NC}"
    fi

    BATCH3_END=$(date +%s)
    BATCH3_ELAPSED=$((($BATCH3_END - $BATCH3_START) / 3600))
    log "${GREEN}Batch 3 complete in ${BATCH3_ELAPSED} hours${NC}"
else
    log "${YELLOW}⊘ Batch 3 skipped by user${NC}"
fi

# =========================================
# BATCH 4: Phase 2.5C + Verification
# =========================================
echo ""
echo "========================================="
echo "BATCH 4: Phase 2.5C + Verification Tests"
echo "Estimated runtime: 17 hours"
echo "========================================="
echo ""

read -p "Continue with Batch 4? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    BATCH4_START=$(date +%s)

    # Phase 2.5C: Ultra-Relativistic Limit (12h)
    run_test "phase_2.5C_ultra_relativistic_limit" \
             "config/scenario_2.5C_ultra_relativistic_limit.yaml" \
             "12h" || true

    # Phase 2.1: Defect Localization (2h)
    run_test "phase_2.1_defect_localization" \
             "config/defect_localization_validation.yaml" \
             "2h" || true

    # Phase 2.2: Traveling Wave (3h)
    run_test "phase_2.2_traveling_wave" \
             "config/traveling_wave_validation.yaml" \
             "3h" || true

    BATCH4_END=$(date +%s)
    BATCH4_ELAPSED=$((($BATCH4_END - $BATCH4_START) / 3600))
    log "${GREEN}Batch 4 complete in ${BATCH4_ELAPSED} hours${NC}"
else
    log "${YELLOW}⊘ Batch 4 skipped by user${NC}"
fi

# =========================================
# Summary
# =========================================
echo ""
echo "========================================="
echo "Master Test Execution Complete"
echo "========================================="
echo ""

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(((TOTAL_END - BATCH1_START) / 3600))

log "========================================="
log "EXECUTION SUMMARY"
log "========================================="
log "Total runtime: ${TOTAL_ELAPSED} hours"
log "End time: $(date)"
log ""
log "Check individual test results in output/ directory"
log "Check logs in ${LOG_DIR}/"
log ""
log "Next steps:"
log "  1. Review test_report.txt in each output directory"
log "  2. Check validation results (automated)"
log "  3. Review generated plots in output/*/plots/"
log "  4. Create phase-level summaries"
log "========================================="

echo ""
echo "${GREEN}All batches complete!${NC}"
echo ""
echo "Review logs: ${MASTER_LOG}"
echo ""
