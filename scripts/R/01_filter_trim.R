#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

input_dir  <- "/data/data/reads"
output_dir <- "/data/scripts/docker-output/01_filter_trim"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters (match existing project settings)
truncLen    <- c(240, 160)
maxN        <- 0
maxEE       <- c(2, 2)
truncQ      <- 2
rm.phix     <- FALSE
multithread <- TRUE

# Subset: F7D6, F7D7, F7D9, F7D25, F7D65
fwd_all <- sort(list.files(input_dir, pattern = "_R1_.*\\.fastq(\\.gz)?$", full.names = TRUE))
rev_all <- sort(list.files(input_dir, pattern = "_R2_.*\\.fastq(\\.gz)?$", full.names = TRUE))

keep      <- grepl("^F7D(6|7|9|25|65)_", basename(fwd_all))
fwd_reads <- fwd_all[keep]
rev_reads <- rev_all[keep]

stopifnot("No matching samples found" = length(fwd_reads) > 0)
stopifnot("R1/R2 count mismatch"      = length(fwd_reads) == length(rev_reads))

sample_names <- sub("_R1_.*", "", basename(fwd_reads))
cat(sprintf("Samples (%d): %s\n\n", length(sample_names), paste(sample_names, collapse = ", ")))

fwd_out <- file.path(output_dir, paste0(sample_names, "_R1_filt.fastq.gz"))
rev_out <- file.path(output_dir, paste0(sample_names, "_R2_filt.fastq.gz"))
names(fwd_out) <- names(rev_out) <- sample_names

out <- filterAndTrim(
  fwd         = fwd_reads,
  filt        = fwd_out,
  rev         = rev_reads,
  filt.rev    = rev_out,
  truncLen    = truncLen,
  maxN        = maxN,
  maxEE       = maxEE,
  truncQ      = truncQ,
  rm.phix     = rm.phix,
  compress    = TRUE,
  multithread = multithread,
  verbose     = TRUE
)
rownames(out) <- sample_names

cat("\nResults (reads.in -> reads.out):\n")
for (i in seq_len(nrow(out))) {
  pct <- if (out[i, 1] > 0) 100 * out[i, 2] / out[i, 1] else 0
  cat(sprintf("  %-20s %6d -> %6d  (%.1f%%)\n",
              rownames(out)[i], out[i, 1], out[i, 2], pct))
}

saveRDS(out,         file.path(output_dir, "filter_stats.rds"))
write.csv(out,       file.path(output_dir, "filter_stats.csv"))
writeLines(sample_names, file.path(output_dir, "sample_names.txt"))
cat("\nDone:", output_dir, "\n")
