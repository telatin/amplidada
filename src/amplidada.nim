import amplidada/[dada, derep, fastq_filter, filter_and_trim, learn_errors, merge_pairs, pipeline, types, version]

export dada, derep, fastq_filter, filter_and_trim, learn_errors, merge_pairs, pipeline, types, version

proc libVersion*(): string =
  ## Returns the public AmpliDADA library version.
  AmpliDadaVersion
