#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 5: makeSequenceTable ==="
echo "Input:  ${REPO_ROOT}/scripts/docker-output/04_merge_pairs/"
echo "Output: ${REPO_ROOT}/scripts/docker-output/05_seq_table/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/05_seq_table.R

echo ""
echo "=== Step 5 complete ==="
