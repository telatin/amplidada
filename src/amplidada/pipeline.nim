import std/[strformat]
import amplidada/[fastq_filter, types]

type
  PipelineOptions* = object
    layout*: ReadLayout
    filter*: FastqFilterOptions

proc defaultPipelineOptions*(layout = rlSingleEnd): PipelineOptions =
  PipelineOptions(layout: layout, filter: defaultFastqFilterOptions())

proc pipelineSummary*(opts: PipelineOptions): string =
  &"layout={opts.layout} filter=({opts.filter.describe()})"

