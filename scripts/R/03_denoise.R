#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

filt_dir   <- "/data/scripts/docker-output/01_filter_trim"
err_dir    <- "/data/scripts/docker-output/02_learn_errors"
output_dir <- "/data/scripts/docker-output/03_denoise"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

errF <- readRDS(file.path(err_dir, "errF.rds"))
errR <- readRDS(file.path(err_dir, "errR.rds"))

fwd_filt <- sort(list.files(filt_dir, pattern = "_R1_filt\\.fastq\\.gz$", full.names = TRUE))
rev_filt <- sort(list.files(filt_dir, pattern = "_R2_filt\\.fastq\\.gz$", full.names = TRUE))
sample_names <- sub("_R1_filt.*", "", basename(fwd_filt))

stopifnot("No filtered reads found — run 01_filter_trim first" = length(fwd_filt) > 0)
cat(sprintf("Denoising %d samples...\n", length(fwd_filt)))

derepF <- derepFastq(fwd_filt, verbose = TRUE)
derepR <- derepFastq(rev_filt, verbose = TRUE)
names(derepF) <- names(derepR) <- sample_names

dadaF <- dada(derepF, err = errF, multithread = TRUE, verbose = TRUE)
dadaR <- dada(derepR, err = errR, multithread = TRUE, verbose = TRUE)

cat("\nASVs per sample (forward):\n")
for (s in sample_names) cat(sprintf("  %-20s %d\n", s, dadaF[[s]]$denoised))

saveRDS(dadaF,  file.path(output_dir, "dadaF.rds"))
saveRDS(dadaR,  file.path(output_dir, "dadaR.rds"))
saveRDS(derepF, file.path(output_dir, "derepF.rds"))
saveRDS(derepR, file.path(output_dir, "derepR.rds"))
cat("Denoised results saved to:", output_dir, "\n")
