# amplidada

[![rewrites.bio - Follows best practice principles for rewriting bioinformatics tools with AI](https://rewrites.bio/badges/rewrites-bio.svg)](https://rewrites.bio)

![Amplidada](docs/amplidada.png)

Port in Nim of part of [DADA2](https://benjjneb.github.io/dada2/)
features. This project started in February 2026, inspired by principles that were later explicited by
 [rewrites.bio](https://rewrites.bio) manifesto. We now try to stick to the manifesto.

A library and a set of CLI tools:
* dada2
* derepFastq
* fastqFilter
* filterAndTrim
* learnErrors
* removeBimerDenovo


## Credits

[DADA2](https://benjjneb.github.io/dada2/) is an R package developed by **Benjamin Callahan**.

Citations:

* *DADA2: High resolution sample inference from Illumina amplicon data*. Nature Methods, 2016. [link](http://dx.doi.org/10.1038/nmeth.3869) [(Open Access link.)](http://rdcu.be/ipGh)
* *Bioconductor workflow for microbiome data analysis: from raw reads to community analyses*. F1000 Research, 2016. [link](https://f1000research.com/articles/5-1492)
* *Exact sequence variants should replace operational taxonomic units in marker-gene data analysis*. ISMEJ, 2017. [link](http://dx.doi.org/10.1038/ismej.2017.119)
* *High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution*. Nucleic Acids Research, 2019. [link](http://dx.doi.org/10.1093/nar/gkz569)

## AI Assistance

This tool was written with the assistance of AI coding agents. The initial plan was drafted using Claude Opus and the first implementation done with Claude Code.
OpenAI Codex was then used for expansion of the core functions.
Validation suite is being done against a Docker image using R 4.5.3 (2026-03-11) and *dada2 1.38.0*, with scripts to perform each step.

The AI generated the implementation; humans defined the validation criteria and verified the results.

