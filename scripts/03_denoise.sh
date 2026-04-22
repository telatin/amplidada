#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 3: dada (denoise) ==="
echo "Input:  ${REPO_ROOT}/scripts/docker-output/01_filter_trim/ (reads)"
echo "        ${REPO_ROOT}/scripts/docker-output/02_learn_errors/ (error models)"
echo "Output: ${REPO_ROOT}/scripts/docker-output/03_denoise/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/03_denoise.R

echo ""
echo "=== Step 3 complete ==="
