#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 4: mergePairs ==="
echo "Input:  ${REPO_ROOT}/scripts/docker-output/03_denoise/"
echo "Output: ${REPO_ROOT}/scripts/docker-output/04_merge_pairs/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/04_merge_pairs.R

echo ""
echo "=== Step 4 complete ==="
