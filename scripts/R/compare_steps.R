#!/usr/bin/env Rscript
# Extracts per-step stats from R/DADA2 intermediate RDS files for comparison
# with Nim pipeline outputs. Run via Docker (repo mounted at /data).
suppressPackageStartupMessages(library(dada2))

cat("====================================================\n")
cat("DADA2 step-by-step comparison: R (Docker) side\n")
cat("====================================================\n\n")

r_base <- "/data/scripts/docker-output"

# ── Step 1: filterAndTrim ─────────────────────────────────────────────────────
cat("── Step 1: filterAndTrim ──────────────────────────────\n")
fs <- read.csv(file.path(r_base, "01_filter_trim", "filter_stats.csv"), row.names = 1)
for (s in rownames(fs)) {
  cat(sprintf("  %-25s in=%6d  out=%6d  (%.1f%%)\n",
    s, fs[s, 1], fs[s, 2], 100 * fs[s, 2] / fs[s, 1]))
}
cat(sprintf("  TOTAL                     in=%6d  out=%6d  (%.1f%%)\n",
  sum(fs[,1]), sum(fs[,2]), 100*sum(fs[,2])/sum(fs[,1])))
cat("\n")

# ── Step 2: learnErrors ───────────────────────────────────────────────────────
cat("── Step 2: learnErrors ────────────────────────────────\n")
errF <- readRDS(file.path(r_base, "02_learn_errors", "errF.rds"))
errR <- readRDS(file.path(r_base, "02_learn_errors", "errR.rds"))
# errF / errR are lists; the estimated error matrix is in $err_out (16 x 41)
eF <- errF$err_out
eR <- errR$err_out
qbins <- c(0, 10, 20, 25, 30, 35, 40)
cat("  Mismatch (A->C) probability at selected quality bins:\n")
cat(sprintf("  %8s  %12s  %12s\n", "Quality", "errF (R)", "errR (R)"))
for (q in qbins) {
  pF <- eF["A2C", as.character(q)]
  pR <- eR["A2C", as.character(q)]
  cat(sprintf("  %8d  %12.2e  %12.2e\n", q, pF, pR))
}
cat("\n")

# ── Step 3: dada (denoise) ────────────────────────────────────────────────────
cat("── Step 3: dada (denoise) ─────────────────────────────\n")
dadaF <- readRDS(file.path(r_base, "03_denoise", "dadaF.rds"))
dadaR <- readRDS(file.path(r_base, "03_denoise", "dadaR.rds"))
derepF <- readRDS(file.path(r_base, "03_denoise", "derepF.rds"))
derepR <- readRDS(file.path(r_base, "03_denoise", "derepR.rds"))
cat(sprintf("  %-25s  %7s  %7s  %7s  %7s\n",
  "Sample", "fwdUniq", "fwdASVs", "revUniq", "revASVs"))
for (s in names(dadaF)) {
  nFU <- length(derepF[[s]]$uniques)
  nRU <- length(derepR[[s]]$uniques)
  nFA <- length(dadaF[[s]]$denoised)
  nRA <- length(dadaR[[s]]$denoised)
  cat(sprintf("  %-25s  %7d  %7d  %7d  %7d\n", s, nFU, nFA, nRU, nRA))
}
cat("\n")

# ── Step 4: mergePairs ────────────────────────────────────────────────────────
cat("── Step 4: mergePairs ─────────────────────────────────\n")
mergers <- readRDS(file.path(r_base, "04_merge_pairs", "mergers.rds"))
cat(sprintf("  %-25s  %8s  %8s  %7s\n", "Sample", "merged", "types", "pct"))
for (s in names(mergers)) {
  m        <- mergers[[s]]
  total    <- sum(m$abundance[m$accept])
  types    <- sum(m$accept)
  total_in <- sum(derepF[[s]]$uniques)
  cat(sprintf("  %-25s  %8d  %8d  %6.1f%%\n", s, total, types,
    100 * total / total_in))
}
cat("\n")

# ── Step 5: seqtab ────────────────────────────────────────────────────────────
cat("── Step 5: makeSequenceTable ──────────────────────────\n")
seqtab <- readRDS(file.path(r_base, "05_seq_table", "seqtab.rds"))
cat(sprintf("  Samples: %d   ASVs: %d\n", nrow(seqtab), ncol(seqtab)))
cat("  Length distribution:\n")
lt <- table(nchar(colnames(seqtab)))
for (n in names(lt)) cat(sprintf("    %s bp: %d\n", n, lt[[n]]))
cat("\n")

# ── Step 6: removeBimeraDenovo ────────────────────────────────────────────────
cat("── Step 6: removeBimeraDenovo ─────────────────────────\n")
seqtab_nc <- readRDS(file.path(r_base, "06_remove_bimera", "seqtab_nochim.rds"))
cat(sprintf("  Before: %d ASVs\n", ncol(seqtab)))
cat(sprintf("  After:  %d ASVs  (%.1f%% of reads retained)\n",
  ncol(seqtab_nc), 100 * sum(seqtab_nc) / sum(seqtab)))
cat("  Per-sample reads retained:\n")
for (s in rownames(seqtab_nc)) {
  cat(sprintf("    %-25s  %d -> %d  (%.1f%%)\n",
    s, sum(seqtab[s,]), sum(seqtab_nc[s,]),
    100 * sum(seqtab_nc[s,]) / sum(seqtab[s,])))
}
cat("\n")

# ── Error model comparison: write Nim-readable TSV ───────────────────────────
cat("── Exporting R error model for direct comparison ──────\n")
out_dir <- "/data/scripts/nim-output/compare"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_err_tsv <- function(mat, path) {
  transitions <- rownames(mat)
  qualities   <- as.integer(colnames(mat))
  rows <- data.frame(
    transition = rep(transitions, each = length(qualities)),
    quality    = rep(qualities, times = length(transitions)),
    prob       = as.vector(t(mat))
  )
  write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE)
}
write_err_tsv(eF, file.path(out_dir, "R_errF.tsv"))
write_err_tsv(eR, file.path(out_dir, "R_errR.tsv"))
cat("  Written:", file.path(out_dir, "R_errF.tsv"), "\n")
cat("  Written:", file.path(out_dir, "R_errR.tsv"), "\n\n")

cat("Done.\n")
