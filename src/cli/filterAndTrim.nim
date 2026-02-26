import std/[os, parseopt, sets, strformat, strutils]
import amplidada

proc usage() =
  echo "Usage: filterAndTrim (--input-dir <dir> | --samplesheet <csv>) --output-dir <dir> [options]"
  echo ""
  echo "Input source (choose one):"
  echo "      --input-dir <dir>        Discover paired FASTQs in a directory"
  echo "      --samplesheet <csv>      CSV with sample/read paths"
  echo "      --output-dir <dir>       Destination directory"
  echo ""
  echo "Directory discovery options:"
  echo "      --r1-token <str>         Forward token (default: _R1)"
  echo "      --r2-token <str>         Reverse token (default: _R2)"
  echo "      --recursive              Scan subdirectories recursively"
  echo ""
  echo "Sample-sheet options:"
  echo "      --sheet-delim <char>     Delimiter (default: ,)"
  echo "      --sheet-has-header       First row is header (default)"
  echo "      --sheet-no-header        No header row"
  echo "      --sheet-sample-col <s>   Sample ID column (default: SampleID)"
  echo "      --sheet-forward-col <s>  Forward reads column (default: reads_forward)"
  echo "      --sheet-reverse-col <s>  Reverse reads column (default: reads_reverse)"
  echo ""
  echo "Output naming:"
  echo "      --out-forward-suffix <s> Default: _R1.filtered.fastq.gz"
  echo "      --out-reverse-suffix <s> Default: _R2.filtered.fastq.gz"
  echo ""
  echo "Filtering options (applied to both reads unless -r overrides are used):"
  echo "      --trunc-len <int>"
  echo "      --trunc-q <int>"
  echo "      --max-len <int>"
  echo "      --max-ee <float>"
  echo "      --max-n <int>"
  echo "      --min-len <int>"
  echo "      --min-q <int>"
  echo "      --trim-left <int>"
  echo "      --trim-right <int>"
  echo ""
  echo "Reverse-only filtering overrides:"
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
  echo "      --threads <int>          Read-level workers per sample"
  echo "      --chunk-size <int>       Reads per chunk"
  echo "      --parallel-min-chunk <int>"
  echo "      --sample-threads <int>   Parallel sample workers"
  echo "      --check-names"
  echo "      --no-check-names"
  echo "      --continue-on-error"
  echo "      --verbose                Print per-sample progress details"
  echo "      --summary <path>"
  echo "      --quiet"
  echo "  -h, --help"
  echo "  -v, --version"

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
    "--input-dir", "--samplesheet", "--output-dir", "--summary",
    "--r1-token", "--r2-token", "--out-forward-suffix", "--out-reverse-suffix",
    "--sheet-delim", "--sheet-sample-col", "--sheet-forward-col", "--sheet-reverse-col",
    "--threads", "--chunk-size", "--parallel-min-chunk", "--sample-threads",
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

proc printBatchSummary(batch: BatchPairedFilterResult, sampleCount: int,
                      filterOpts: PairedFastqFilterOptions, batchOpts: BatchFilterOptions) =
  echo "Mode: filterAndTrim batch paired-end"
  echo &"Samples: {sampleCount}"
  echo &"Failed samples: {batch.failedSamples}"
  echo &"Sample threads: {batchOpts.sampleThreads}"
  echo &"Read-level options: {filterOpts.describe()}"
  echo &"Pairs in: {batch.totals.inputPairs}"
  echo &"Pairs out: {batch.totals.outputPairs}"
  echo &"Pairs discarded: {batch.totals.discardedPairs}"
  echo &"Forward reject reasons: {reasonSummary(batch.totals.forwardStats)}"
  echo &"Reverse reject reasons: {reasonSummary(batch.totals.reverseStats)}"

proc writeBatchSummary(path: string, batch: BatchPairedFilterResult, sampleCount: int,
                      filterOpts: PairedFastqFilterOptions, batchOpts: BatchFilterOptions) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to write summary file: " & path)
  defer: close(f)

  f.writeLine("mode=batch-paired")
  f.writeLine("samples=" & $sampleCount)
  f.writeLine("failed_samples=" & $batch.failedSamples)
  f.writeLine("sample_threads=" & $batchOpts.sampleThreads)
  f.writeLine("continue_on_error=" & $batchOpts.continueOnError)
  f.writeLine("filter_options=" & filterOpts.describe())
  f.writeLine("input_pairs=" & $batch.totals.inputPairs)
  f.writeLine("output_pairs=" & $batch.totals.outputPairs)
  f.writeLine("discarded_pairs=" & $batch.totals.discardedPairs)

  for i, sample in batch.samples:
    f.writeLine(&"sample[{i}].id={sample.sampleId}")
    f.writeLine(&"sample[{i}].success={sample.success}")
    if not sample.success:
      f.writeLine(&"sample[{i}].error={sample.errorMessage}")

proc main() =
  var inputDir = ""
  var sampleSheet = ""
  var outputDir = ""
  var summaryPath = ""
  var quiet = false
  var recursive = false
  var r1Token = "_R1"
  var r2Token = "_R2"
  var outForwardSuffix = "_R1.filtered.fastq.gz"
  var outReverseSuffix = "_R2.filtered.fastq.gz"

  var filterOpts = defaultPairedFastqFilterOptions()
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
        of "input-dir":
          inputDir = value
        of "samplesheet":
          sampleSheet = value
        of "output-dir":
          outputDir = value
        of "summary":
          summaryPath = value
        of "quiet":
          quiet = true
        of "verbose":
          batchOpts.verbose = true

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

        of "threads":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.threads = x
          filterOpts.reverse.threads = x
        of "chunk-size":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.chunkSize = x
          filterOpts.reverse.chunkSize = x
        of "parallel-min-chunk":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.parallelMinChunk = x
          filterOpts.reverse.parallelMinChunk = x
        of "sample-threads":
          batchOpts.sampleThreads = parseIntFlag(flag, value)

        of "check-names":
          filterOpts.checkNames = true
        of "no-check-names":
          filterOpts.checkNames = false
        of "continue-on-error":
          batchOpts.continueOnError = true

        of "trunc-len":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.truncLen = x
          filterOpts.reverse.truncLen = x
        of "trunc-q":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.truncQ = x
          filterOpts.reverse.truncQ = x
        of "max-len":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.maxLen = x
          filterOpts.reverse.maxLen = x
        of "max-ee":
          let x = parseFloatFlag(flag, value)
          filterOpts.forward.maxEE = x
          filterOpts.reverse.maxEE = x
        of "max-n":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.maxN = x
          filterOpts.reverse.maxN = x
        of "min-len":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.minLen = x
          filterOpts.reverse.minLen = x
        of "min-q":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.minQ = x
          filterOpts.reverse.minQ = x
        of "trim-left":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.trimLeft = x
          filterOpts.reverse.trimLeft = x
        of "trim-right":
          let x = parseIntFlag(flag, value)
          filterOpts.forward.trimRight = x
          filterOpts.reverse.trimRight = x

        of "trunc-len-r":
          filterOpts.reverse.truncLen = parseIntFlag(flag, value)
        of "trunc-q-r":
          filterOpts.reverse.truncQ = parseIntFlag(flag, value)
        of "max-len-r":
          filterOpts.reverse.maxLen = parseIntFlag(flag, value)
        of "max-ee-r":
          filterOpts.reverse.maxEE = parseFloatFlag(flag, value)
        of "max-n-r":
          filterOpts.reverse.maxN = parseIntFlag(flag, value)
        of "min-len-r":
          filterOpts.reverse.minLen = parseIntFlag(flag, value)
        of "min-q-r":
          filterOpts.reverse.minQ = parseIntFlag(flag, value)
        of "trim-left-r":
          filterOpts.reverse.trimLeft = parseIntFlag(flag, value)
        of "trim-right-r":
          filterOpts.reverse.trimRight = parseIntFlag(flag, value)

        else:
          fail("Unknown option: --" & key)
      of cmdArgument:
        fail("Unexpected positional argument: " & key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  if outputDir.len == 0:
    fail("Missing required option: --output-dir")
  if quiet and batchOpts.verbose:
    fail("--quiet and --verbose are mutually exclusive")
  if (inputDir.len == 0 and sampleSheet.len == 0) or (inputDir.len > 0 and sampleSheet.len > 0):
    fail("Choose exactly one input source: --input-dir OR --samplesheet")

  try:
    validateOptions(filterOpts)
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
    let batch = filterAndTrimPaired(jobs, filterOpts, batchOpts)

    if not quiet:
      printBatchSummary(batch, samples.len, filterOpts, batchOpts)

    if summaryPath.len > 0:
      writeBatchSummary(summaryPath, batch, samples.len, filterOpts, batchOpts)

    if batch.failedSamples > 0:
      quit(3)
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
