#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== Step 2: learnErrors ==="
echo "Input:  ${REPO_ROOT}/scripts/docker-output/01_filter_trim/"
echo "Output: ${REPO_ROOT}/scripts/docker-output/02_learn_errors/"
echo ""

docker run --rm \
  -v "${REPO_ROOT}:/data" \
  dada2:last \
  Rscript /data/scripts/R/02_learn_errors.R

echo ""
echo "=== Step 2 complete ==="
