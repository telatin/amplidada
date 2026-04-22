#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqtab_dir <- "/data/scripts/docker-output/05_seq_table"
output_dir <- "/data/scripts/docker-output/06_remove_bimera"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seqtab <- readRDS(file.path(seqtab_dir, "seqtab.rds"))
cat(sprintf("Before chimera removal: %d ASVs\n", ncol(seqtab)))

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)

pct_reads  <- 100 * sum(seqtab_nochim) / sum(seqtab)
pct_asvs   <- 100 * ncol(seqtab_nochim) / ncol(seqtab)
cat(sprintf("After  chimera removal: %d ASVs (%.1f%% of ASVs, %.1f%% of reads retained)\n",
            ncol(seqtab_nochim), pct_asvs, pct_reads))

# FASTA
asvs      <- getSequences(seqtab_nochim)
asv_names <- sprintf("ASV%04d", seq_along(asvs))
writeLines(as.vector(rbind(paste0(">", asv_names), asvs)),
           file.path(output_dir, "asvs.fasta"))

# TSV: rows = ASVs, cols = samples
asv_tab <- t(seqtab_nochim)
rownames(asv_tab) <- asv_names
write.table(asv_tab, file.path(output_dir, "asv_table.tsv"), sep = "\t", quote = FALSE)

saveRDS(seqtab_nochim, file.path(output_dir, "seqtab_nochim.rds"))
cat("Final output saved to:", output_dir, "\n")
cat("  asvs.fasta       —", length(asvs), "sequences\n")
cat("  asv_table.tsv    — ASVs x samples count matrix\n")
cat("  seqtab_nochim.rds\n")
