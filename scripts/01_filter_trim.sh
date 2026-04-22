#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 1: filterAndTrim ==="
echo "Input:  ${REPO_ROOT}/data/reads/F7D{6,7,9,25,65}*"
echo "Output: ${REPO_ROOT}/scripts/docker-output/01_filter_trim/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/01_filter_trim.R

echo ""
echo "=== Step 1 complete ==="
