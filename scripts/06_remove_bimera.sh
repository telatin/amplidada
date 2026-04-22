#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 6: removeBimeraDenovo ==="
echo "Input:  ${REPO_ROOT}/scripts/docker-output/05_seq_table/"
echo "Output: ${REPO_ROOT}/scripts/docker-output/06_remove_bimera/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/06_remove_bimera.R

echo ""
echo "=== Step 6 complete ==="
