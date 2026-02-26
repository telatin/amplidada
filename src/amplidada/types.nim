type
  ReadLayout* = enum
    rlSingleEnd,
    rlPairedEnd

  PipelineStep* = enum
    psFilter,
    psDereplicate,
    psDenoise,
    psMergePairs,
    psChimeraRemoval,
    psSequenceTable

