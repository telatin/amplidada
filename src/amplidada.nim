import amplidada/[dada, derep, fastq_filter, filter_and_trim, learn_errors, merge_pairs, pipeline, remove_bimera_denovo, types, version]

export dada, derep, fastq_filter, filter_and_trim, learn_errors, merge_pairs, pipeline, remove_bimera_denovo, types, version

proc libVersion*(): string =
  ## Returns the public AmpliDADA library version.
  AmpliDadaVersion
