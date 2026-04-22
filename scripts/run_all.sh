#!/usr/bin/env bash
# Full DADA2 pipeline on F7D{6,7,9,25,65} subset via dada2:last Docker image.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

START=$(date +%s)

bash "${SCRIPT_DIR}/01_filter_trim.sh"
bash "${SCRIPT_DIR}/02_learn_errors.sh"
bash "${SCRIPT_DIR}/03_denoise.sh"
bash "${SCRIPT_DIR}/04_merge_pairs.sh"
bash "${SCRIPT_DIR}/05_seq_table.sh"
bash "${SCRIPT_DIR}/06_remove_bimera.sh"

END=$(date +%s)
echo ""
echo "All steps complete in $(( END - START ))s."
echo "Final outputs in: ${SCRIPT_DIR}/docker-output/06_remove_bimera/"
echo "  asvs.fasta        — denoised, chimera-free ASV sequences"
echo "  asv_table.tsv     — ASV x sample count matrix"
echo "  seqtab_nochim.rds — R object for downstream analysis"
