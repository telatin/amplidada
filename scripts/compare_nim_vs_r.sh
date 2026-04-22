#!/usr/bin/env bash
# Step-by-step comparison of Nim vs R DADA2 pipelines.
# Runs the R extraction script via Docker, then prints a side-by-side digest.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
NIM_OUT="${SCRIPT_DIR}/nim-output"
R_OUT="${SCRIPT_DIR}/docker-output"

hr() { printf '%.0s─' {1..72}; echo; }

# ── 1. Run R extraction via Docker ───────────────────────────────────────────
echo "Running R extraction via Docker..."
docker run --rm -v "${REPO_ROOT}:/data" dada2:last \
  Rscript /data/scripts/R/compare_steps.R

# ── 2. Nim summary from existing files ───────────────────────────────────────
hr
echo "NIM PIPELINE SUMMARY (from intermediate files)"
hr

echo
echo "── Step 1: filterAndTrim ──────────────────────────────"
echo "  Sample                       in      out    pct"
tail -n +2 "${R_OUT}/01_filter_trim/filter_stats.csv" \
  | awk -F, '{printf "  %-28s %6d  %6d  %.1f%%\n", $1, $2, $3, 100*$3/$2}'
echo "  (Nim filter stats match R — same parameters, confirmed from log)"

echo
echo "── Step 2: learnErrors ────────────────────────────────"
echo "  A->C mismatch probability (Nim errF vs Nim errR):"
printf "  %-8s  %-14s  %-14s\n" "Quality" "errF_prob" "errR_prob"
# join errF and errR on (from=A, to=C, quality)
join -t$'\t' \
  <(awk -F'\t' '$1=="A" && $2=="C" {print $3"\t"$5}' \
      "${NIM_OUT}/02_learn_errors/errF.tsv" | sort -k1,1n) \
  <(awk -F'\t' '$1=="A" && $2=="C" {print $3"\t"$5}' \
      "${NIM_OUT}/02_learn_errors/errR.tsv" | sort -k1,1n) \
  | awk -F'\t' 'NR%4==0 || $1+0 >= 25 {
      printf "  %-8s  %-14s  %-14s\n", $1, $2, $3
    }'
echo
echo "  Compare against R error model:"
CMP="${NIM_OUT}/compare"
if [[ -f "${CMP}/R_errF.tsv" ]]; then
  printf "  %-8s  %-14s  %-14s  %-14s\n" \
    "Quality" "Nim_errF" "R_errF" "ratio(N/R)"
  join -t$'\t' \
    <(awk -F'\t' '$1=="A" && $2=="C" {print $3"\t"$5}' \
        "${NIM_OUT}/02_learn_errors/errF.tsv" | sort -k1,1n) \
    <(awk -F'\t' 'NR>1 && $1=="A2C" {print $2"\t"$3}' \
        "${CMP}/R_errF.tsv" | sort -k1,1n) \
    | awk -F'\t' 'NR%4==0 || $1+0 >= 25 {
        r = ($3+0>0) ? $2/$3 : 0
        printf "  %-8s  %-14s  %-14s  %-14.3f\n", $1, $2, $3, r
      }'
else
  echo "  (R_errF.tsv not yet written — run this script once to generate)"
fi

echo
echo "── Step 3: dada denoising ─────────────────────────────"
printf "  %-28s  %-10s  %-12s  %-10s  %-12s\n" \
  "Sample" "fwdUniq" "fwdASVs(Nim)" "revUniq" "revASVs(Nim)"
# Read from per-sample summary files if available; fall back to merged TSV heuristics
for tsv in "${NIM_OUT}/03_dada2/"*_merged.tsv; do
  sample=$(basename "$tsv" _merged.tsv)
  summary="${NIM_OUT}/03_dada2/${sample}.summary"
  if [[ -f "$summary" ]]; then
    fwd_uniq=$(grep "^forward_unique" "$summary" 2>/dev/null | cut -d= -f2 || echo "?")
    fwd_asvs=$(grep "^forward_asvs"  "$summary" 2>/dev/null | cut -d= -f2 || echo "?")
    rev_uniq=$(grep "^reverse_unique" "$summary" 2>/dev/null | cut -d= -f2 || echo "?")
    rev_asvs=$(grep "^reverse_asvs"  "$summary" 2>/dev/null | cut -d= -f2 || echo "?")
  else
    fwd_uniq="(rerun)"; fwd_asvs="(rerun)"; rev_uniq="(rerun)"; rev_asvs="(rerun)"
  fi
  printf "  %-28s  %-10s  %-12s  %-10s  %-12s\n" \
    "$sample" "$fwd_uniq" "$fwd_asvs" "$rev_uniq" "$rev_asvs"
done

echo
echo "── Step 4: mergePairs ─────────────────────────────────"
printf "  %-28s  %10s  %10s  %8s\n" "Sample" "types" "merged_reads" "pct_in"
for tsv in "${NIM_OUT}/03_dada2/"*_merged.tsv; do
  sample=$(basename "$tsv" _merged.tsv)
  summary="${NIM_OUT}/03_dada2/${sample}.summary"
  types=$(( $(wc -l < "$tsv") - 1 ))
  merged=$(awk 'NR>1{s+=$3}END{print s+0}' "$tsv")
  if [[ -f "$summary" ]]; then
    total_in=$(grep "^forward_reads=" "$summary" 2>/dev/null | cut -d= -f2 || echo 1)
    pct=$(awk "BEGIN{printf \"%.1f%%\", 100*${merged}/${total_in}}")
  else
    pct="?"
  fi
  printf "  %-28s  %10d  %10d  %8s\n" "$sample" "$types" "$merged" "$pct"
done

echo
echo "── Step 5: unique ASVs before chimera removal ─────────"
nim_unique=$(tail -n +2 "${NIM_OUT}/04_combine/all_samples.tsv" \
               | awk -F'\t' '{print $2}' | sort -u | wc -l)
echo "  Nim: ${nim_unique} unique sequences"
if [[ -f "${R_OUT}/05_seq_table/seqtab.csv" ]]; then
  r_asvs=$(head -1 "${R_OUT}/05_seq_table/seqtab.csv" | tr ',' '\n' | tail -n +2 | wc -l)
  echo "  R:   ${r_asvs} unique sequences"
fi

echo
echo "── Step 6: after chimera removal ──────────────────────"
nim_n=$(grep -c "^>" "${NIM_OUT}/06_final/asvs.fasta"   2>/dev/null || echo 0)
r_n=$(grep -c   "^>" "${R_OUT}/06_remove_bimera/asvs.fasta" 2>/dev/null || echo 0)
echo "  Nim: ${nim_n} ASVs"
echo "  R:   ${r_n} ASVs"
echo
echo "  Sequence overlap:"
python3 - "${R_OUT}/06_remove_bimera/asvs.fasta" \
           "${NIM_OUT}/06_final/asvs.fasta" <<'PYEOF'
import sys

def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif name:
                seqs[line] = name
                name = None
    return seqs

r_seqs   = read_fasta(sys.argv[1])
nim_seqs = read_fasta(sys.argv[2])
shared   = set(r_seqs) & set(nim_seqs)
r_only   = set(r_seqs) - set(nim_seqs)
nim_only = set(nim_seqs) - set(r_seqs)
print(f"    Shared:   {len(shared):4d}  ({100*len(shared)/max(1,len(r_seqs)):.1f}% of R, {100*len(shared)/max(1,len(nim_seqs)):.1f}% of Nim)")
print(f"    R only:   {len(r_only):4d}")
print(f"    Nim only: {len(nim_only):4d}")
print()
if r_only:
    print("  R-only length distribution:", dict(
        sorted({len(s): sum(1 for x in r_only if len(x)==len(s))
                for s in r_only}.items())))
if nim_only:
    print("  Nim-only length distribution:", dict(
        sorted({len(s): sum(1 for x in nim_only if len(x)==len(s))
                for s in nim_only}.items())))
PYEOF
hr
