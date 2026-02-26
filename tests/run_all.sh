#!/usr/bin/env bash
set -euo pipefail

# ---------- Config ----------
INPUT_DIR="${1:-data/small/filter_and_trim/}"
OUT_DIR="${2:-data/small/results/}"
DADA2_BIN="${DADA2_BIN:-bin/dada2}"

# Recommended settings for filtered reads (from your sweep)
THREADS="${THREADS:-4}"
MERGE_THREADS="${MERGE_THREADS:-4}"
MIN_OVERLAP="${MIN_OVERLAP:-12}"
MAX_MISMATCH="${MAX_MISMATCH:-1}"
TRIM_OVERHANG="${TRIM_OVERHANG:-1}"   # 1=yes, 0=no

# Canonical self-consistency (optional, 1=yes, 0=no)
USE_DADA_SELF="${USE_DADA_SELF:-1}"
DADA_SELF_MIN="${DADA_SELF_MIN:-2}"
DADA_SELF_MAX="${DADA_SELF_MAX:-4}"
DADA_SELF_TOL="${DADA_SELF_TOL:-0}"

# ASV reporting filters:
# - Use ASV_MIN_LEN/ASV_MAX_LEN to force a specific merged-length window.
# - If unset, an automatic window around the dominant merged length is used.
ASV_MIN_TOTAL="${ASV_MIN_TOTAL:-10}"
ASV_MIN_LEN="${ASV_MIN_LEN:-}"
ASV_MAX_LEN="${ASV_MAX_LEN:-}"
ASV_LEN_MODE_WINDOW="${ASV_LEN_MODE_WINDOW:-2}"

# If 1, skip samples that already have outputs
RESUME="${RESUME:-1}"
# ---------- /Config ----------

mkdir -p "$OUT_DIR"/{merged,summary,logs}
[[ -x "$DADA2_BIN" ]] || { echo "ERROR: cannot execute $DADA2_BIN"; exit 1; }

shopt -s nullglob
r1_files=("$INPUT_DIR"/*_R1.filtered.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1.filtered.fastq.gz in $INPUT_DIR"; exit 1; }

echo "Found ${#r1_files[@]} R1 files in $INPUT_DIR"

ok=0
fail=0

for r1 in "${r1_files[@]}"; do
  base="$(basename "$r1")"
  sample="${base%_R1.filtered.fastq.gz}"
  r2="$INPUT_DIR/${sample}_R2.filtered.fastq.gz"

  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for $sample, skipping ($r2)"
    ((fail++)) || true
    continue
  fi

  merged_out="$OUT_DIR/merged/${sample}.merged.tsv"
  summary_out="$OUT_DIR/summary/${sample}.summary.txt"
  log_out="$OUT_DIR/logs/${sample}.log"

  if [[ "$RESUME" == "1" && -s "$merged_out" && -s "$summary_out" ]]; then
    echo "[SKIP] $sample (already done)"
    ((ok++)) || true
    continue
  fi

  cmd=(
    "$DADA2_BIN"
    --in-forward "$r1"
    --in-reverse "$r2"
    --out "$merged_out"
    --summary "$summary_out"
    --threads "$THREADS"
    --merge-threads "$MERGE_THREADS"
    --min-overlap "$MIN_OVERLAP"
    --max-mismatch "$MAX_MISMATCH"
    --max-iter 64
    --max-shuffle 12
    --quiet
  )

  if [[ "$TRIM_OVERHANG" == "1" ]]; then
    cmd+=(--trim-overhang)
  fi

  if [[ "$USE_DADA_SELF" == "1" ]]; then
    cmd+=(
      --dada-self-consist
      --dada-self-min-iter "$DADA_SELF_MIN"
      --dada-self-max-iter "$DADA_SELF_MAX"
      --dada-self-tol "$DADA_SELF_TOL"
    )
  fi

  echo "[RUN] $sample"
  if "${cmd[@]}" >"$log_out" 2>&1; then
    ((ok++)) || true
  else
    echo "[FAIL] $sample (see $log_out)"
    ((fail++)) || true
  fi
done

echo "Done. success=$ok failed=$fail"

# -------- Build long + wide ASV tables --------
raw_long_tsv="$OUT_DIR/asv_long.raw.tsv"
long_tsv="$OUT_DIR/asv_long.tsv"
seq_list="$OUT_DIR/asv_sequences.txt"
seqtab="$OUT_DIR/seqtab.tsv"
asv_fasta="$OUT_DIR/asv.fasta"
counts_tsv="$OUT_DIR/counts.tsv"
asv_map="$OUT_DIR/asv_map.tsv"
sample_list="$OUT_DIR/sample_list.txt"

echo -e "SampleID\tSequence\tAbundance" > "$raw_long_tsv"
for f in "$OUT_DIR"/merged/*.merged.tsv; do
  [[ -f "$f" ]] || continue
  sample="$(basename "$f" .merged.tsv)"
  awk -F'\t' -v s="$sample" 'NR>1{print s"\t"$2"\t"$3}' "$f" >> "$raw_long_tsv"
done

if [[ -z "$ASV_MIN_LEN" || -z "$ASV_MAX_LEN" ]]; then
  mode_len="$(awk -F'\t' '
    NR > 1 {
      l = length($2)
      w[l] += $3
    }
    END {
      best = -1
      mode = 0
      for (l in w) {
        if (w[l] > best) {
          best = w[l]
          mode = l
        }
      }
      print mode
    }
  ' "$raw_long_tsv")"

  if [[ "$mode_len" =~ ^[0-9]+$ && "$mode_len" -gt 0 ]]; then
    if [[ -z "$ASV_MIN_LEN" ]]; then
      ASV_MIN_LEN=$((mode_len - ASV_LEN_MODE_WINDOW))
    fi
    if [[ -z "$ASV_MAX_LEN" ]]; then
      ASV_MAX_LEN=$((mode_len + ASV_LEN_MODE_WINDOW))
    fi
  else
    ASV_MIN_LEN="${ASV_MIN_LEN:-0}"
    ASV_MAX_LEN="${ASV_MAX_LEN:-0}"
  fi
fi

if [[ "$ASV_MIN_LEN" =~ ^-?[0-9]+$ && "$ASV_MIN_LEN" -lt 0 ]]; then
  ASV_MIN_LEN=0
fi
if [[ "$ASV_MAX_LEN" =~ ^-?[0-9]+$ && "$ASV_MAX_LEN" -lt 0 ]]; then
  ASV_MAX_LEN=0
fi

awk -F'\t' -v OFS='\t' \
  -v min_total="$ASV_MIN_TOTAL" \
  -v min_len="$ASV_MIN_LEN" \
  -v max_len="$ASV_MAX_LEN" '
  NR == 1 {
    header = $0
    next
  }
  {
    seq = $2
    totals[seq] += $3
    row[NR] = $0
    rowSeq[NR] = seq
    rowLen[NR] = length(seq)
  }
  END {
    print header
    for (i = 2; i <= NR; i++) {
      seq = rowSeq[i]
      l = rowLen[i]
      if (totals[seq] < min_total) {
        continue
      }
      if (min_len > 0 && l < min_len) {
        continue
      }
      if (max_len > 0 && l > max_len) {
        continue
      }
      print row[i]
    }
  }
' "$raw_long_tsv" > "$long_tsv"

ls "$OUT_DIR"/merged/*.merged.tsv 2>/dev/null | xargs -n1 basename | sed 's/\.merged\.tsv$//' | sort > "$sample_list"

awk -F'\t' 'NR>1{print $2}' "$long_tsv" | sort -u > "$seq_list"

{
  printf "SampleID"
  while IFS= read -r seq; do
    printf "\t%s" "$seq"
  done < "$seq_list"
  printf "\n"
} > "$seqtab"

for f in "$OUT_DIR"/merged/*.merged.tsv; do
  [[ -f "$f" ]] || continue
  sample="$(basename "$f" .merged.tsv)"
  declare -A counts=()
  while IFS=$'\t' read -r _id seq abundance; do
    [[ "$seq" == "sequence" ]] && continue
    counts["$seq"]="$abundance"
  done < "$f"

  {
    printf "%s" "$sample"
    while IFS= read -r seq; do
      printf "\t%s" "${counts[$seq]:-0}"
    done < "$seq_list"
    printf "\n"
  } >> "$seqtab"
done

# -------- Build ASV FASTA + transposed counts table --------
# Map each unique sequence to asv1..asvN sorted by descending total abundance.
awk -F'\t' 'NR>1{sum[$2]+=$3} END{for(seq in sum) printf "%s\t%d\n", seq, sum[seq]}' "$long_tsv" \
  | sort -k2,2nr -k1,1 \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print "asv" NR, $1, $2}' > "$asv_map"

awk -F'\t' 'BEGIN{OFS="\t"} {printf ">%s totalcounts=%s\n%s\n", $1, $3, $2}' "$asv_map" > "$asv_fasta"

awk -F'\t' '
BEGIN { OFS="\t" }
FILENAME == ARGV[1] {
  asvOrder[++nAsv] = $1
  seqToAsv[$2] = $1
  next
}
FILENAME == ARGV[2] {
  samples[++nSample] = $1
  next
}
FILENAME == ARGV[3] {
  if (FNR == 1) next
  sample = $1
  seq = $2
  abundance = $3 + 0
  asv = seqToAsv[seq]
  if (asv != "") {
    counts[asv SUBSEP sample] += abundance
  }
  next
}
END {
  printf "ASV"
  for (i = 1; i <= nSample; i++) printf OFS samples[i]
  printf "\n"

  for (i = 1; i <= nAsv; i++) {
    asv = asvOrder[i]
    printf "%s", asv
    for (j = 1; j <= nSample; j++) {
      key = asv SUBSEP samples[j]
      if (key in counts) printf OFS counts[key]
      else printf OFS "0"
    }
    printf "\n"
  }
}
' "$asv_map" "$sample_list" "$long_tsv" > "$counts_tsv"

rm -f "$sample_list"

echo "ASV outputs:"
echo "  per-sample merged: $OUT_DIR/merged/"
echo "  per-sample summary: $OUT_DIR/summary/"
echo "  long table (raw): $raw_long_tsv"
echo "  long table: $long_tsv"
echo "  ASV filters: min_total=$ASV_MIN_TOTAL min_len=$ASV_MIN_LEN max_len=$ASV_MAX_LEN (auto-window=$ASV_LEN_MODE_WINDOW)"
echo "  ASV FASTA: $asv_fasta"
echo "  ASV map: $asv_map"
echo "  transposed counts: $counts_tsv"
echo "  sequence table: $seqtab"
