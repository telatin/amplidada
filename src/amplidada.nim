import amplidada/[dada, derep, fastq_filter, filter_and_trim, learn_errors, loess, loess_errfun_mod4, merge_pairs, pipeline, remove_bimera_denovo, sequence_table, types, version]

export dada, derep, fastq_filter, filter_and_trim, learn_errors, loess, loess_errfun_mod4, merge_pairs, pipeline, remove_bimera_denovo, sequence_table, types, version

proc libVersion*(): string =
  ## Returns the public AmpliDADA library version.
  AmpliDadaVersion
