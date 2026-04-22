#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

merge_dir  <- "/data/scripts/docker-output/04_merge_pairs"
output_dir <- "/data/scripts/docker-output/05_seq_table"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mergers <- readRDS(file.path(merge_dir, "mergers.rds"))

cat("Building sequence table...\n")
seqtab <- makeSequenceTable(mergers)
cat(sprintf("  Samples: %d   ASVs: %d\n", nrow(seqtab), ncol(seqtab)))

cat("\nASV length distribution:\n")
print(table(nchar(getSequences(seqtab))))

saveRDS(seqtab,   file.path(output_dir, "seqtab.rds"))
write.csv(t(seqtab), file.path(output_dir, "seqtab.csv"))
cat("Sequence table saved to:", output_dir, "\n")
