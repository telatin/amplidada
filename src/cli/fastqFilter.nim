import std/[os, parseopt, sets, strformat, strutils]
import amplidada

proc usage() =
  echo "Usage: fastqFilter [mode options] [filter options]"
  echo ""
  echo "Direct single/paired mode:"
  echo "  -i, --in <path>              Input FASTQ/FASTQ.GZ"
  echo "  -o, --out <path>             Output FASTQ/FASTQ.GZ (use '-' for stdout)"
  echo "      --in2 <path>             Reverse input for paired-end filtering"
  echo "      --out2 <path>            Reverse output for paired-end filtering"
  echo "      --paired                 Force paired mode (requires --in2/--out2)"
  echo ""
  echo "Batch filterAndTrim mode (paired):"
  echo "      --input-dir <dir>        Discover *_R{1,2} FASTQs in a directory"
  echo "      --samplesheet <csv>      CSV with SampleID,reads_forward,reads_reverse"
  echo "      --output-dir <dir>       Output directory for filtered FASTQs"
  echo "      --r1-token <str>         Forward token for discovery (default: _R1)"
  echo "      --r2-token <str>         Reverse token for discovery (default: _R2)"
  echo "      --recursive              Recursive directory discovery"
  echo "      --out-forward-suffix <s> Output suffix (default: _R1.filtered.fastq.gz)"
  echo "      --out-reverse-suffix <s> Output suffix (default: _R2.filtered.fastq.gz)"
  echo "      --sample-threads <int>   Parallel sample workers in batch mode"
  echo "      --continue-on-error      Continue batch even if one sample fails"
  echo ""
  echo "Sample-sheet controls:"
  echo "      --sheet-delim <char>     Delimiter (default: ,)"
  echo "      --sheet-has-header       Interpret first row as header (default)"
  echo "      --sheet-no-header        Parse rows as fixed order"
  echo "      --sheet-sample-col <s>   Sample ID column name"
  echo "      --sheet-forward-col <s>  Forward reads column name"
  echo "      --sheet-reverse-col <s>  Reverse reads column name"
  echo ""
  echo "Core filters (applied to both reads unless reverse override is provided):"
  echo "      --trunc-len <int>        Truncate reads to this length (0 disables)"
  echo "      --trunc-q <int>          Truncate at first quality <= value (0 disables)"
  echo "      --max-len <int>          Drop reads longer than this (0 disables)"
  echo "      --max-ee <float>         Max expected errors"
  echo "      --max-n <int>            Max ambiguous bases (N)"
  echo "      --min-len <int>          Min read length after trimming"
  echo "      --min-q <int>            Drop if minimum Q <= value"
  echo "      --trim-left <int>        Trim fixed bases from 5' end"
  echo "      --trim-right <int>       Trim fixed bases from 3' end"
  echo ""
  echo "Reverse-only overrides (paired mode):"
  echo "      --trunc-len-r <int>"
  echo "      --trunc-q-r <int>"
  echo "      --max-len-r <int>"
  echo "      --max-ee-r <float>"
  echo "      --max-n-r <int>"
  echo "      --min-len-r <int>"
  echo "      --min-q-r <int>"
  echo "      --trim-left-r <int>"
  echo "      --trim-right-r <int>"
  echo ""
  echo "Parallelism and execution:"
  echo "      --threads <int>          Worker threads for record-level filtering"
  echo "      --chunk-size <int>       Records per processing chunk"
  echo "      --parallel-min-chunk <int>  Minimum chunk size for parallel mode"
  echo "      --check-names            Require synchronized R1/R2 read IDs"
  echo "      --no-check-names         Disable R1/R2 read ID checking"
  echo "      --summary <path>         Write machine-readable run summary"
  echo "      --quiet                  Only print fatal errors"
  echo "  -v, --version                Print version and exit"
  echo "  -h, --help                   Show this help and exit"

proc parseIntFlag(flag, value: string): int =
  try:
    result = parseInt(value)
  except ValueError:
    raise newException(ValueError, "Invalid integer for --" & flag & ": " & value)

proc parseFloatFlag(flag, value: string): float =
  try:
    result = parseFloat(value)
  except ValueError:
    raise newException(ValueError, "Invalid float for --" & flag & ": " & value)

proc parseCharFlag(flag, value: string): char =
  if value.len != 1:
    raise newException(ValueError, "Invalid char for --" & flag & ": " & value)
  value[0]

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Error: " & msg)
  quit(2)

proc normalizeArgs(args: seq[string]): seq[string] =
  let longNeedsValue = toHashSet([
    "--in", "--in2", "--out", "--out2", "--summary",
    "--input-dir", "--output-dir", "--samplesheet",
    "--r1-token", "--r2-token", "--out-forward-suffix", "--out-reverse-suffix",
    "--sheet-delim", "--sheet-sample-col", "--sheet-forward-col", "--sheet-reverse-col",
    "--sample-threads",
    "--threads", "--chunk-size", "--parallel-min-chunk",
    "--trunc-len", "--trunc-q", "--max-len", "--max-ee", "--max-n", "--min-len", "--min-q",
    "--trim-left", "--trim-right",
    "--trunc-len-r", "--trunc-q-r", "--max-len-r", "--max-ee-r", "--max-n-r",
    "--min-len-r", "--min-q-r", "--trim-left-r", "--trim-right-r"
  ])

  var i = 0
  while i < args.len:
    let a = args[i]
    if a in longNeedsValue and i + 1 < args.len and (args[i + 1].len == 0 or args[i + 1][0] != '-'):
      result.add(a & "=" & args[i + 1])
      inc i
    else:
      result.add(a)
    inc i

proc reasonSummary(stats: FastqFilterStats): string =
  var parts: seq[string] = @[]
  for reason in FilterRejectReason:
    if reason == frrNone:
      continue
    let count = stats.rejectedByReason[reason]
    if count > 0:
      parts.add(&"{reason}={count}")
  if parts.len == 0:
    return "none"
  parts.join(",")

proc printSingleSummary(stats: FastqFilterStats, opts: FastqFilterOptions) =
  echo "Mode: single-end"
  echo &"Options: {opts.describe()}"
  echo &"Reads in:  {stats.inputReads}"
  echo &"Reads out: {stats.outputReads}"
  echo &"Discarded: {stats.discardedReads} ({formatFloat(stats.discardRatio() * 100.0, ffDecimal, 2)}%)"
  echo &"Reject reasons: {reasonSummary(stats)}"

proc printPairedSummary(stats: FastqPairedFilterStats, opts: PairedFastqFilterOptions) =
  echo "Mode: paired-end"
  echo &"Options: {opts.describe()}"
  echo &"Pairs in:  {stats.inputPairs}"
  echo &"Pairs out: {stats.outputPairs}"
  echo &"Pairs discarded: {stats.discardedPairs} ({formatFloat(stats.discardRatio() * 100.0, ffDecimal, 2)}%)"
  echo &"Forward reject reasons: {reasonSummary(stats.forwardStats)}"
  echo &"Reverse reject reasons: {reasonSummary(stats.reverseStats)}"
  var pairReasonParts: seq[string] = @[]
  for reason in PairRejectReason:
    let count = stats.rejectedByPairReason[reason]
    if count > 0:
      pairReasonParts.add(&"{reason}={count}")
  if pairReasonParts.len == 0:
    echo "Pair reject reasons: none"
  else:
    echo "Pair reject reasons: " & pairReasonParts.join(",")

proc printBatchSummary(
  batch: BatchPairedFilterResult,
  samplesCount: int,
  filterOpts: PairedFastqFilterOptions,
  batchOpts: BatchFilterOptions
) =
  echo "Mode: batch paired-end (filterAndTrim style)"
  echo &"Samples: {samplesCount}"
  echo &"Failed samples: {batch.failedSamples}"
  echo &"Sample threads: {batchOpts.sampleThreads}"
  echo &"Record-level options: {filterOpts.describe()}"
  printPairedSummary(batch.totals, filterOpts)

proc writeSingleSummary(path: string, stats: FastqFilterStats, opts: FastqFilterOptions) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to write summary file: " & path)
  defer: close(f)

  f.writeLine("mode=single-end")
  f.writeLine("options=" & opts.describe())
  f.writeLine("input_reads=" & $stats.inputReads)
  f.writeLine("output_reads=" & $stats.outputReads)
  f.writeLine("discarded_reads=" & $stats.discardedReads)
  for reason in FilterRejectReason:
    if reason == frrNone:
      continue
    f.writeLine(&"reject_{reason}={stats.rejectedByReason[reason]}")

proc writePairedSummary(path: string, stats: FastqPairedFilterStats, opts: PairedFastqFilterOptions) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to write summary file: " & path)
  defer: close(f)

  f.writeLine("mode=paired-end")
  f.writeLine("options=" & opts.describe())
  f.writeLine("input_pairs=" & $stats.inputPairs)
  f.writeLine("output_pairs=" & $stats.outputPairs)
  f.writeLine("discarded_pairs=" & $stats.discardedPairs)
  for reason in FilterRejectReason:
    if reason == frrNone:
      continue
    f.writeLine(&"forward_reject_{reason}={stats.forwardStats.rejectedByReason[reason]}")
    f.writeLine(&"reverse_reject_{reason}={stats.reverseStats.rejectedByReason[reason]}")
  for reason in PairRejectReason:
    f.writeLine(&"pair_reject_{reason}={stats.rejectedByPairReason[reason]}")

proc writeBatchSummary(
  path: string,
  batch: BatchPairedFilterResult,
  samplesCount: int,
  filterOpts: PairedFastqFilterOptions,
  batchOpts: BatchFilterOptions
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to write summary file: " & path)
  defer: close(f)

  f.writeLine("mode=batch-paired")
  f.writeLine("samples=" & $samplesCount)
  f.writeLine("failed_samples=" & $batch.failedSamples)
  f.writeLine("continue_on_error=" & $batchOpts.continueOnError)
  f.writeLine("sample_threads=" & $batchOpts.sampleThreads)
  f.writeLine("filter_options=" & filterOpts.describe())
  f.writeLine("input_pairs=" & $batch.totals.inputPairs)
  f.writeLine("output_pairs=" & $batch.totals.outputPairs)
  f.writeLine("discarded_pairs=" & $batch.totals.discardedPairs)
  for i, sample in batch.samples:
    let status = if sample.success: "ok" else: "failed"
    f.writeLine(&"sample[{i}].id={sample.sampleId}")
    f.writeLine(&"sample[{i}].status={status}")
    if not sample.success:
      f.writeLine(&"sample[{i}].error={sample.errorMessage}")

proc main() =
  var input1 = ""
  var input2 = ""
  var output1 = ""
  var output2 = ""
  var inputDir = ""
  var outputDir = ""
  var sampleSheet = ""
  var summaryPath = ""
  var forcedPaired = false
  var quiet = false
  var recursive = false

  var r1Token = "_R1"
  var r2Token = "_R2"
  var outForwardSuffix = "_R1.filtered.fastq.gz"
  var outReverseSuffix = "_R2.filtered.fastq.gz"

  var pairedOpts = defaultPairedFastqFilterOptions()
  var batchOpts = defaultBatchFilterOptions()
  var sheetOpts = defaultSampleSheetOptions()

  var parser = initOptParser(normalizeArgs(commandLineParams()))
  try:
    for kind, key, value in parser.getopt():
      case kind
      of cmdLongOption, cmdShortOption:
        let flag = key.toLowerAscii()
        case flag
        of "help", "h":
          usage()
          return
        of "version", "v":
          echo AmpliDadaVersion
          return

        of "in", "i":
          input1 = value
        of "in2":
          input2 = value
        of "out", "o":
          output1 = value
        of "out2":
          output2 = value
        of "paired":
          forcedPaired = true

        of "input-dir":
          inputDir = value
        of "output-dir":
          outputDir = value
        of "samplesheet":
          sampleSheet = value
        of "r1-token":
          r1Token = value
        of "r2-token":
          r2Token = value
        of "recursive":
          recursive = true
        of "out-forward-suffix":
          outForwardSuffix = value
        of "out-reverse-suffix":
          outReverseSuffix = value
        of "sample-threads":
          batchOpts.sampleThreads = parseIntFlag(flag, value)
        of "continue-on-error":
          batchOpts.continueOnError = true

        of "sheet-delim":
          sheetOpts.delimiter = parseCharFlag(flag, value)
        of "sheet-has-header":
          sheetOpts.hasHeader = true
        of "sheet-no-header":
          sheetOpts.hasHeader = false
        of "sheet-sample-col":
          sheetOpts.sampleIdColumn = value
        of "sheet-forward-col":
          sheetOpts.forwardColumn = value
        of "sheet-reverse-col":
          sheetOpts.reverseColumn = value

        of "check-names":
          pairedOpts.checkNames = true
        of "no-check-names":
          pairedOpts.checkNames = false

        of "summary":
          summaryPath = value
        of "quiet":
          quiet = true

        of "threads":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.threads = x
          pairedOpts.reverse.threads = x
        of "chunk-size":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.chunkSize = x
          pairedOpts.reverse.chunkSize = x
        of "parallel-min-chunk":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.parallelMinChunk = x
          pairedOpts.reverse.parallelMinChunk = x

        of "trunc-len":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.truncLen = x
          pairedOpts.reverse.truncLen = x
        of "trunc-q":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.truncQ = x
          pairedOpts.reverse.truncQ = x
        of "max-len":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.maxLen = x
          pairedOpts.reverse.maxLen = x
        of "max-ee":
          let x = parseFloatFlag(flag, value)
          pairedOpts.forward.maxEE = x
          pairedOpts.reverse.maxEE = x
        of "max-n":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.maxN = x
          pairedOpts.reverse.maxN = x
        of "min-len":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.minLen = x
          pairedOpts.reverse.minLen = x
        of "min-q":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.minQ = x
          pairedOpts.reverse.minQ = x
        of "trim-left":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.trimLeft = x
          pairedOpts.reverse.trimLeft = x
        of "trim-right":
          let x = parseIntFlag(flag, value)
          pairedOpts.forward.trimRight = x
          pairedOpts.reverse.trimRight = x

        of "trunc-len-r":
          pairedOpts.reverse.truncLen = parseIntFlag(flag, value)
        of "trunc-q-r":
          pairedOpts.reverse.truncQ = parseIntFlag(flag, value)
        of "max-len-r":
          pairedOpts.reverse.maxLen = parseIntFlag(flag, value)
        of "max-ee-r":
          pairedOpts.reverse.maxEE = parseFloatFlag(flag, value)
        of "max-n-r":
          pairedOpts.reverse.maxN = parseIntFlag(flag, value)
        of "min-len-r":
          pairedOpts.reverse.minLen = parseIntFlag(flag, value)
        of "min-q-r":
          pairedOpts.reverse.minQ = parseIntFlag(flag, value)
        of "trim-left-r":
          pairedOpts.reverse.trimLeft = parseIntFlag(flag, value)
        of "trim-right-r":
          pairedOpts.reverse.trimRight = parseIntFlag(flag, value)

        else:
          fail("Unknown option: --" & key)

      of cmdArgument:
        fail("Unexpected positional argument: " & key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  let batchMode = inputDir.len > 0 or sampleSheet.len > 0

  if batchMode:
    if inputDir.len > 0 and sampleSheet.len > 0:
      fail("Batch mode accepts only one source: use either --input-dir or --samplesheet")
    if outputDir.len == 0:
      fail("Batch mode requires --output-dir")
    if input1.len > 0 or input2.len > 0 or output1.len > 0 or output2.len > 0:
      fail("Batch mode does not use --in/--in2/--out/--out2")

    try:
      validateOptions(pairedOpts)
      validateOptions(batchOpts)
      validateOptions(sheetOpts)

      var samples: seq[PairedSampleInput]
      if inputDir.len > 0:
        samples = discoverPairedSamples(inputDir, r1Token = r1Token, r2Token = r2Token, recursive = recursive)
      else:
        samples = readPairedSampleSheet(sampleSheet, sheetOpts)

      if samples.len == 0:
        fail("No paired samples found")

      let jobs = buildPairedJobs(samples, outputDir, outForwardSuffix, outReverseSuffix)
      let batch = filterAndTrimPaired(jobs, pairedOpts, batchOpts)

      if not quiet:
        printBatchSummary(batch, samples.len, pairedOpts, batchOpts)

      if summaryPath.len > 0:
        writeBatchSummary(summaryPath, batch, samples.len, pairedOpts, batchOpts)

      if batch.failedSamples > 0:
        quit(3)
      return
    except CatchableError as e:
      fail(e.msg)

  if input1.len == 0:
    fail("Missing required option: --in")
  if output1.len == 0:
    fail("Missing required option: --out")

  let pairedMode = forcedPaired or input2.len > 0 or output2.len > 0

  if pairedMode:
    if input2.len == 0:
      fail("Paired mode requires --in2")
    if output2.len == 0:
      fail("Paired mode requires --out2")
  else:
    if input2.len > 0 or output2.len > 0:
      fail("--in2/--out2 require paired mode")

  try:
    if pairedMode:
      validateOptions(pairedOpts)
      let stats = fastqPairedFilter(input1, input2, output1, output2, pairedOpts)
      if not quiet:
        printPairedSummary(stats, pairedOpts)
      if summaryPath.len > 0:
        writePairedSummary(summaryPath, stats, pairedOpts)
    else:
      validateOptions(pairedOpts.forward)
      let stats = fastqFilter(input1, output1, pairedOpts.forward)
      if not quiet:
        printSingleSummary(stats, pairedOpts.forward)
      if summaryPath.len > 0:
        writeSingleSummary(summaryPath, stats, pairedOpts.forward)
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
