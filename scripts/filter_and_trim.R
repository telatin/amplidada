#!/usr/bin/env Rscript

# This script will be executed in a micromamba environment with R and DADA2 installed and should
# perform DADA2 filter and Trim procedure to be able to test the results and compare with the Nim implementation

# ── Parameters ────────────────────────────────────────────────────────────────
truncLen  <- c(240, 160)   # truncation length for R1, R2
maxN      <- 0             # no ambiguous bases allowed
maxEE     <- c(2, 2)       # max expected errors for R1, R2
truncQ    <- 2             # truncate at first base with quality <= truncQ
rm.phix   <- FALSE          # remove PhiX reads
compress  <- TRUE          # gzip output files
multithread <- TRUE        # use all available cores (set FALSE on Windows)

setwd("/Users/telatina/git/dada2-nim/")
input_Dir  <- "./data/reads"
output_Dir <- "./data/filter_trim_R_out"

# ── Setup ──────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(dada2))

if (!dir.exists(output_Dir)) {
  dir.create(output_Dir, recursive = TRUE)
  cat("Created output directory:", output_Dir, "\n")
}

# ── Discover paired files ──────────────────────────────────────────────────────
fwd_reads <- sort(list.files(input_Dir, pattern = "_R1_.*\\.fastq(\\.gz)?$", full.names = TRUE))
rev_reads <- sort(list.files(input_Dir, pattern = "_R2_.*\\.fastq(\\.gz)?$", full.names = TRUE))

stopifnot("No forward reads found" = length(fwd_reads) > 0)
stopifnot("Forward/reverse count mismatch" = length(fwd_reads) == length(rev_reads))

sample_names <- sub("_R1_.*", "", basename(fwd_reads))

cat(sprintf("Found %d paired samples:\n", length(sample_names)))
for (s in sample_names) cat("  ", s, "\n")

# ── Build output file paths ────────────────────────────────────────────────────
fwd_out <- file.path(output_Dir, paste0(sample_names, "_R1_filt.fastq.gz"))
rev_out <- file.path(output_Dir, paste0(sample_names, "_R2_filt.fastq.gz"))

# ── filterAndTrim ──────────────────────────────────────────────────────────────
cat("\nRunning filterAndTrim...\n")
cat(sprintf("  truncLen  = (%d, %d)\n", truncLen[1], truncLen[2]))
cat(sprintf("  maxEE     = (%g, %g)\n", maxEE[1], maxEE[2]))
cat(sprintf("  maxN      = %d\n", maxN))
cat(sprintf("  truncQ    = %d\n", truncQ))
cat(sprintf("  rm.phix   = %s\n", rm.phix))
cat(sprintf("  multithread = %s\n\n", multithread))

out <- filterAndTrim(
  fwd        = fwd_reads,
  filt       = fwd_out,
  rev        = rev_reads,
  filt.rev   = rev_out,
  truncLen   = truncLen,
  maxN       = maxN,
  maxEE      = maxEE,
  truncQ     = truncQ,
  rm.phix    = rm.phix,
  compress   = compress,
  multithread = multithread,
  verbose    = TRUE
)

# ── Report ─────────────────────────────────────────────────────────────────────
rownames(out) <- sample_names
cat("Results (reads.in / reads.out / % retained):\n")
for (i in seq_len(nrow(out))) {
  pct <- if (out[i, "reads.in"] > 0) 100 * out[i, "reads.out"] / out[i, "reads.in"] else 0
  cat(sprintf("  %-30s %6d -> %6d  (%.1f%%)\n",
              rownames(out)[i], out[i, "reads.in"], out[i, "reads.out"], pct))
}
cat(sprintf("\nTotal: %d -> %d reads retained (%.1f%%)\n",
            sum(out[, "reads.in"]),
            sum(out[, "reads.out"]),
            100 * sum(out[, "reads.out"]) / sum(out[, "reads.in"])))

cat("\nDone. Filtered reads written to:", output_Dir, "\n")
