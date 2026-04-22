#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

denoise_dir <- "/data/scripts/docker-output/03_denoise"
output_dir  <- "/data/scripts/docker-output/04_merge_pairs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dadaF  <- readRDS(file.path(denoise_dir, "dadaF.rds"))
dadaR  <- readRDS(file.path(denoise_dir, "dadaR.rds"))
derepF <- readRDS(file.path(denoise_dir, "derepF.rds"))
derepR <- readRDS(file.path(denoise_dir, "derepR.rds"))

sample_names <- names(dadaF)
cat(sprintf("Merging pairs for %d samples...\n", length(sample_names)))

mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)

cat("\nMerged reads per sample:\n")
for (s in sample_names) {
  n <- sum(mergers[[s]]$accept)
  cat(sprintf("  %-20s %d\n", s, n))
}

saveRDS(mergers, file.path(output_dir, "mergers.rds"))
cat("Mergers saved to:", output_dir, "\n")
