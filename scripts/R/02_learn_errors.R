#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

filt_dir   <- "/data/scripts/docker-output/01_filter_trim"
output_dir <- "/data/scripts/docker-output/02_learn_errors"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fwd_filt <- sort(list.files(filt_dir, pattern = "_R1_filt\\.fastq\\.gz$", full.names = TRUE))
rev_filt <- sort(list.files(filt_dir, pattern = "_R2_filt\\.fastq\\.gz$", full.names = TRUE))
stopifnot("No filtered reads found — run 01_filter_trim first" = length(fwd_filt) > 0)

cat(sprintf("Learning errors from %d samples...\n", length(fwd_filt)))

errF <- learnErrors(fwd_filt, multithread = TRUE, verbose = TRUE)
errR <- learnErrors(rev_filt, multithread = TRUE, verbose = TRUE)

saveRDS(errF, file.path(output_dir, "errF.rds"))
saveRDS(errR, file.path(output_dir, "errR.rds"))

pdf(file.path(output_dir, "error_plots.pdf"))
print(plotErrors(errF, nominalQ = TRUE) + ggplot2::ggtitle("Forward error model"))
print(plotErrors(errR, nominalQ = TRUE) + ggplot2::ggtitle("Reverse error model"))
dev.off()

cat("Error models saved to:", output_dir, "\n")
