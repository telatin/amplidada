#!/usr/bin/env bash
# DADA2 Nim pipeline — subset F7D{6,7,9,25,65}
# Parallel to run_all.sh (Docker/R).  Outputs in scripts/nim-output/{stepname}/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
BIN="${REPO_ROOT}/bin"
READS="${REPO_ROOT}/data/reads"
NIM_OUT="${SCRIPT_DIR}/nim-output"

SAMPLES=(F7D6 F7D7 F7D9 F7D25 F7D65)
START=$(date +%s)

# ── Step 1: filterAndTrim ─────────────────────────────────────────────────────
echo "=== Step 1: filterAndTrim ==="
echo "Input:  ${READS}/F7D{6,7,9,25,65}*"
echo "Output: ${NIM_OUT}/01_filter_trim/"

mkdir -p "${NIM_OUT}/01_filter_trim"

# Build samplesheet — use full sample name (e.g. F7D6_S322_L001) to match R pipeline
SHEET="${NIM_OUT}/01_filter_trim/samplesheet.csv"
echo "SampleID,reads_forward,reads_reverse" > "$SHEET"
for s in "${SAMPLES[@]}"; do
  r1=$(ls "${READS}/${s}_"*"_R1_"*.fastq.gz 2>/dev/null | head -1 || true)
  r2=$(ls "${READS}/${s}_"*"_R2_"*.fastq.gz 2>/dev/null | head -1 || true)
  if [[ -z "$r1" || -z "$r2" ]]; then
    echo "WARNING: reads not found for ${s}, skipping" >&2
    continue
  fi
  sample_name=$(basename "$r1" | sed 's/_R1_.*//')
  echo "${sample_name},${r1},${r2}" >> "$SHEET"
done

"${BIN}/filterAndTrim" \
  --samplesheet    "$SHEET" \
  --output-dir     "${NIM_OUT}/01_filter_trim" \
  --trunc-len 240  --trunc-len-r 160 \
  --max-n 0        --max-ee 2  --trunc-q 2 \
  --out-forward-suffix _R1_filt.fastq.gz \
  --out-reverse-suffix _R2_filt.fastq.gz \
  --threads 8      --verbose

echo "=== Step 1 complete ==="

# ── Step 2: learnErrors ───────────────────────────────────────────────────────
echo "=== Step 2: learnErrors ==="
echo "Input:  ${NIM_OUT}/01_filter_trim/"
echo "Output: ${NIM_OUT}/02_learn_errors/"

mkdir -p "${NIM_OUT}/02_learn_errors"

FWD_FILT=$(ls "${NIM_OUT}/01_filter_trim/"*_R1_filt.fastq.gz | paste -sd,)
REV_FILT=$(ls "${NIM_OUT}/01_filter_trim/"*_R2_filt.fastq.gz | paste -sd,)

echo "  Learning forward error model..."
"${BIN}/learnErrors" \
  --inputs "$FWD_FILT" \
  --out    "${NIM_OUT}/02_learn_errors/errF.tsv" \
  --self-consist --verbose

echo "  Learning reverse error model..."
"${BIN}/learnErrors" \
  --inputs "$REV_FILT" \
  --out    "${NIM_OUT}/02_learn_errors/errR.tsv" \
  --self-consist --verbose

echo "=== Step 2 complete ==="

# ── Step 3: dada2 — denoise + mergePairs, per sample ─────────────────────────
echo "=== Step 3: dada2 (denoise + mergePairs) per sample ==="
echo "Input:  ${NIM_OUT}/01_filter_trim/  +  ${NIM_OUT}/02_learn_errors/"
echo "Output: ${NIM_OUT}/03_dada2/"

mkdir -p "${NIM_OUT}/03_dada2"

for fwd in "${NIM_OUT}/01_filter_trim/"*_R1_filt.fastq.gz; do
  sample=$(basename "$fwd" _R1_filt.fastq.gz)
  rev="${NIM_OUT}/01_filter_trim/${sample}_R2_filt.fastq.gz"
  echo "  [dada2] ${sample}"
  "${BIN}/dada2" \
    --in-forward        "$fwd" \
    --in-reverse        "$rev" \
    --err-matrix-forward "${NIM_OUT}/02_learn_errors/errF.tsv" \
    --err-matrix-reverse "${NIM_OUT}/02_learn_errors/errR.tsv" \
    --out               "${NIM_OUT}/03_dada2/${sample}_merged.tsv" \
    --summary           "${NIM_OUT}/03_dada2/${sample}.summary" \
    --threads 8         --verbose
done

echo "=== Step 3 complete ==="

# ── Step 4: combine per-sample TSVs → long format ────────────────────────────
# dada2 output: merged_id / sequence / abundance  (cols 1 2 3)
# long format:  SampleID  / Sequence / Abundance
echo "=== Step 4: combine into long format ==="
echo "Output: ${NIM_OUT}/04_combine/all_samples.tsv"

mkdir -p "${NIM_OUT}/04_combine"
{
  printf "SampleID\tSequence\tAbundance\n"
  for tsv in "${NIM_OUT}/03_dada2/"*_merged.tsv; do
    sample=$(basename "$tsv" _merged.tsv)
    awk -v s="$sample" 'NR>1 && $3+0>0 { print s "\t" $2 "\t" $3 }' "$tsv"
  done
} > "${NIM_OUT}/04_combine/all_samples.tsv"

echo "  Rows (excl. header): $(( $(wc -l < "${NIM_OUT}/04_combine/all_samples.tsv") - 1 ))"
echo "=== Step 4 complete ==="

# ── Step 5: removeBimeraDenovo ────────────────────────────────────────────────
echo "=== Step 5: removeBimeraDenovo ==="
echo "Input:  ${NIM_OUT}/04_combine/all_samples.tsv"
echo "Output: ${NIM_OUT}/05_remove_bimera/seqtab_nochim.tsv"

mkdir -p "${NIM_OUT}/05_remove_bimera"
"${BIN}/removeBimerDenovo" \
  --input  "${NIM_OUT}/04_combine/all_samples.tsv" \
  --output "${NIM_OUT}/05_remove_bimera/seqtab_nochim.tsv" \
  --method consensus \
  --verbose

echo "=== Step 5 complete ==="

# ── Step 6: asvs.fasta + asv_table.tsv ───────────────────────────────────────
# Converts SampleID/Sequence/Abundance long TSV → FASTA + wide count matrix.
# Sequence order: first-appearance (matches R makeSequenceTable).
echo "=== Step 6: generate final outputs ==="
echo "Output: ${NIM_OUT}/06_final/"

mkdir -p "${NIM_OUT}/06_final"

python3 - \
    "${NIM_OUT}/05_remove_bimera/seqtab_nochim.tsv" \
    "${NIM_OUT}/06_final" \
<<'PYEOF'
import sys, os

inpath, outdir = sys.argv[1], sys.argv[2]

rows = []
with open(inpath) as f:
    hdr = f.readline().rstrip().split('\t')
    ci, cs, ca = hdr.index('SampleID'), hdr.index('Sequence'), hdr.index('Abundance')
    for line in f:
        p = line.rstrip().split('\t')
        rows.append((p[ci], p[cs], int(p[ca])))

seqs_order    = {}
samples_order = {}
for sid, seq, _ in rows:
    seqs_order.setdefault(seq, len(seqs_order))
    samples_order.setdefault(sid, len(samples_order))

seqs      = sorted(seqs_order,    key=seqs_order.__getitem__)
samples   = sorted(samples_order, key=samples_order.__getitem__)
asv_names = [f"ASV{i+1:04d}" for i in range(len(seqs))]

counts = {}
for sid, seq, abd in rows:
    counts[(sid, seq)] = abd

with open(os.path.join(outdir, "asvs.fasta"), "w") as f:
    for name, seq in zip(asv_names, seqs):
        f.write(f">{name}\n{seq}\n")

with open(os.path.join(outdir, "asv_table.tsv"), "w") as f:
    f.write("\t" + "\t".join(samples) + "\n")
    for name, seq in zip(asv_names, seqs):
        row = [str(counts.get((s, seq), 0)) for s in samples]
        f.write(name + "\t" + "\t".join(row) + "\n")

print(f"  ASVs: {len(seqs)}   Samples: {len(samples)}")
print(f"  asvs.fasta    written to: {outdir}")
print(f"  asv_table.tsv written to: {outdir}")
PYEOF

END=$(date +%s)
echo "=== Step 6 complete ==="
echo ""
echo "All steps complete in $(( END - START ))s."
echo "Final outputs: ${NIM_OUT}/06_final/"
echo "  asvs.fasta     — denoised, chimera-free ASV sequences"
echo "  asv_table.tsv  — ASV × sample count matrix"
