import std/[os, parseopt, sets, strformat, strutils]
import amplidada

proc usage() =
  echo "Usage: learnErrors --in <reads.fastq[.gz]> [--in <reads2.fastq[.gz]> ...] --out <error_matrix.tsv> [options]"
  echo ""
  echo "Inputs:"
  echo "  --in <path>                   Input FASTQ/FASTQ.GZ (repeatable)"
  echo "  --inputs <p1,p2,...>          Comma-separated input FASTQ list"
  echo "  --input-list <path>           Text file with one FASTQ path per line"
  echo "  --input-dir <dir>             Include FASTQ files from a directory"
  echo "  --recursive                   Recurse when scanning --input-dir"
  echo "  --nbases <int>                Stop after this many bases (default: 100000000)"
  echo ""
  echo "Required:"
  echo "  --out <path>                  Output error matrix TSV"
  echo ""
  echo "Optional outputs:"
  echo "  --counts-out <path>           Output raw transition counts TSV"
  echo "  --summary <path>              Write run summary"
  echo ""
  echo "learnErrors options:"
  echo "  --q-min <int>                 Min quality bin (default: 0)"
  echo "  --q-max <int>                 Max quality bin (default: 40)"
  echo "  --pseudocount <float>         Additive smoothing (default: 1.0)"
  echo "  --zero-mismatch-prob <float>  Mismatch prob when no data at a bin"
  echo "  --self-consist                Enable iterative self-consistency"
  echo "  --min-iter <int>              Minimum self-consistency iterations"
  echo "  --max-iter <int>              Maximum self-consistency iterations"
  echo "  --tol <float>                 Convergence tolerance for max prob delta"
  echo "  --max-centers-per-length <int> Candidate centers per length during refinement"
  echo ""
  echo "Derep options:"
  echo "  --threads <int>               Worker threads"
  echo "  --chunk-size <int>            Reads per chunk"
  echo "  --parallel-min-chunk <int>    Min chunk for multithreading"
  echo "  --order <abundance|first>     Derep ordering"
  echo ""
  echo "Reporting:"
  echo "  --verbose                     Print progress details"
  echo "  --quiet                       Suppress summary output"
  echo ""
  echo "General:"
  echo "  -h, --help"
  echo "  -v, --version"

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Error: " & msg)
  quit(2)

proc parseIntFlag(flag, value: string): int =
  try:
    result = parseInt(value)
  except ValueError:
    raise newException(ValueError, "Invalid integer for --" & flag & ": " & value)

proc parseFloatFlag(flag, value: string): float64 =
  try:
    result = parseFloat(value)
  except ValueError:
    raise newException(ValueError, "Invalid float for --" & flag & ": " & value)

proc normalizeArgs(args: seq[string]): seq[string] =
  let longNeedsValue = toHashSet([
    "--in", "--inputs", "--input-list", "--input-dir", "--nbases",
    "--out", "--counts-out", "--summary", "--q-min", "--q-max",
    "--pseudocount", "--zero-mismatch-prob", "--threads", "--chunk-size",
    "--parallel-min-chunk", "--order", "--min-iter", "--max-iter",
    "--tol", "--max-centers-per-length"
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

proc isFastqPath(path: string): bool =
  let lower = path.toLowerAscii()
  lower.endsWith(".fastq") or lower.endsWith(".fastq.gz") or
    lower.endsWith(".fq") or lower.endsWith(".fq.gz")

proc addCsvInputs(inputPaths: var seq[string], csvValue: string) =
  for raw in csvValue.split(','):
    let p = raw.strip()
    if p.len > 0:
      inputPaths.add(p)

proc addListInputs(inputPaths: var seq[string], listPath: string) =
  if not fileExists(listPath):
    fail("Input list file does not exist: " & listPath)
  for line in lines(listPath):
    let p = line.strip()
    if p.len == 0 or p[0] == '#':
      continue
    inputPaths.add(p)

proc addDiscoveredInputs(inputPaths: var seq[string], inputDir: string, recursive: bool) =
  if not dirExists(inputDir):
    fail("Input directory does not exist: " & inputDir)

  if recursive:
    for filePath in walkDirRec(inputDir):
      if isFastqPath(filePath):
        inputPaths.add(filePath)
  else:
    for kind, filePath in walkDir(inputDir):
      if kind == pcFile and isFastqPath(filePath):
        inputPaths.add(filePath)

proc normalizeInputPaths(inputPaths: openArray[string]): seq[string] =
  var seen = initHashSet[string]()
  for raw in inputPaths:
    let stripped = raw.strip()
    if stripped.len == 0:
      continue
    let absPath = absolutePath(stripped)
    if absPath notin seen:
      seen.incl(absPath)
      result.add(absPath)

proc writeSummary(path: string, batchRes: LearnErrorsBatchResult, learnRes: LearnErrorsResult,
                 derepOpts: DerepFastqOptions, learnOpts: LearnErrorsOptions,
                 scOpts: LearnErrorsSelfConsistOptions) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open summary output file: " & path)
  defer: close(f)

  f.writeLine("input_files_requested=" & $batchRes.filesRequested)
  f.writeLine("input_files_used=" & $batchRes.filesUsed)
  f.writeLine("input_reads=" & $batchRes.readsUsed)
  f.writeLine("input_bases=" & $batchRes.basesUsed)
  f.writeLine("unique_sequences_sum=" & $batchRes.uniqueCountSum)
  f.writeLine("centers_by_length=" & $learnRes.centerCountByLength)
  f.writeLine("total_transitions=" & $learnRes.totalTransitions)
  f.writeLine("skipped_positions=" & $learnRes.skippedPositions)
  f.writeLine("self_consist_enabled=" & $batchRes.selfConsistEnabled)
  f.writeLine("self_consist_iterations=" & $batchRes.iterationsRun)
  f.writeLine("self_consist_converged=" & $batchRes.converged)
  f.writeLine("self_consist_max_abs_prob_delta=" & $batchRes.maxAbsProbDelta)
  f.writeLine("self_consist_center_changes_total=" & $batchRes.centerChangesTotal)
  f.writeLine("derep_options=" & derepOpts.describe())
  f.writeLine("learn_options=" & learnOpts.describe())
  f.writeLine("self_consist_options=" & scOpts.describe())

proc printSummary(batchRes: LearnErrorsBatchResult, learnRes: LearnErrorsResult,
                 derepOpts: DerepFastqOptions, learnOpts: LearnErrorsOptions,
                 scOpts: LearnErrorsSelfConsistOptions) =
  echo "Mode: learnErrors (phase 2)"
  echo &"Input files used: {batchRes.filesUsed}/{batchRes.filesRequested}"
  echo &"Input reads: {batchRes.readsUsed}"
  echo &"Input bases: {batchRes.basesUsed}"
  echo &"Unique sequences (sum): {batchRes.uniqueCountSum}"
  echo &"Centers by length: {learnRes.centerCountByLength}"
  echo &"Total transitions: {learnRes.totalTransitions}"
  echo &"Skipped positions: {learnRes.skippedPositions}"
  echo &"Self-consistency: enabled={batchRes.selfConsistEnabled} iterations={batchRes.iterationsRun} converged={batchRes.converged}"
  echo &"Self-consistency max prob delta: {batchRes.maxAbsProbDelta}"
  echo &"Self-consistency center changes: {batchRes.centerChangesTotal}"
  echo &"Derep options: {derepOpts.describe()}"
  echo &"Learn options: {learnOpts.describe()}"
  echo &"Self-consistency options: {scOpts.describe()}"

proc main() =
  var inputPaths: seq[string] = @[]
  var inputDir = ""
  var recursive = false
  var nbases = 100_000_000'i64
  var outputPath = ""
  var countsOutPath = ""
  var summaryPath = ""
  var quiet = false
  var verbose = false

  var derepOpts = defaultDerepFastqOptions()
  derepOpts.keepReadMap = false

  var learnOpts = defaultLearnErrorsOptions()
  var scOpts = defaultLearnErrorsSelfConsistOptions()

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
        of "in":
          inputPaths.add(value)
        of "inputs":
          addCsvInputs(inputPaths, value)
        of "input-list":
          addListInputs(inputPaths, value)
        of "input-dir":
          inputDir = value
        of "recursive":
          recursive = true
        of "nbases":
          nbases = parseIntFlag(flag, value).int64
        of "out":
          outputPath = value
        of "counts-out":
          countsOutPath = value
        of "summary":
          summaryPath = value
        of "threads":
          derepOpts.threads = parseIntFlag(flag, value)
        of "chunk-size":
          derepOpts.chunkSize = parseIntFlag(flag, value)
        of "parallel-min-chunk":
          derepOpts.parallelMinChunk = parseIntFlag(flag, value)
        of "order":
          case value.toLowerAscii()
          of "abundance":
            derepOpts.order = dsoAbundance
          of "first", "firstseen", "first-seen":
            derepOpts.order = dsoFirstSeen
          else:
            fail("Invalid value for --order: " & value)
        of "q-min":
          learnOpts.qMin = parseIntFlag(flag, value)
        of "q-max":
          learnOpts.qMax = parseIntFlag(flag, value)
        of "pseudocount":
          learnOpts.pseudocount = parseFloatFlag(flag, value)
        of "zero-mismatch-prob":
          learnOpts.zeroTotalMismatchProb = parseFloatFlag(flag, value)
        of "self-consist":
          scOpts.enabled = true
        of "min-iter":
          scOpts.minIterations = parseIntFlag(flag, value)
        of "max-iter":
          scOpts.maxIterations = parseIntFlag(flag, value)
        of "tol":
          scOpts.tol = parseFloatFlag(flag, value)
        of "max-centers-per-length":
          scOpts.maxCentersPerLength = parseIntFlag(flag, value)
        of "verbose":
          verbose = true
        of "quiet":
          quiet = true
        else:
          fail("Unknown option: --" & key)
      of cmdArgument:
        inputPaths.add(key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  if inputDir.len > 0:
    addDiscoveredInputs(inputPaths, inputDir, recursive)

  let normalizedInputs = normalizeInputPaths(inputPaths)
  if normalizedInputs.len == 0:
    fail("No input FASTQ files provided. Use --in, --inputs, --input-list, or --input-dir")
  if outputPath.len == 0:
    fail("Missing required option: --out")
  if quiet and verbose:
    fail("--quiet and --verbose are mutually exclusive")
  if nbases == 0:
    fail("--nbases must be > 0 or negative to disable cutoff")

  try:
    validateOptions(derepOpts)
    validateOptions(learnOpts)
    validateOptions(scOpts)

    if verbose:
      echo &"[learnErrors] Inputs: {normalizedInputs.len} file(s)"
      if nbases > 0:
        echo &"[learnErrors] Base target: {nbases}"
      else:
        echo "[learnErrors] Base target: disabled"
      echo &"[learnErrors] Self-consistency: {scOpts.describe()}"

    let batchRes = learnErrorsFromFastqPaths(
      inputPaths = normalizedInputs,
      derepOpts = derepOpts,
      learnOpts = learnOpts,
      nbases = nbases,
      selfConsistOpts = scOpts
    )
    let learnRes = batchRes.matrix

    if verbose:
      echo &"[learnErrors] Used {batchRes.filesUsed} file(s), reads={batchRes.readsUsed}, bases={batchRes.basesUsed}"

    if verbose:
      echo &"[learnErrors] Writing matrix: {outputPath}"
    writeLearnErrorsTsv(outputPath, learnRes)

    if countsOutPath.len > 0:
      if verbose:
        echo &"[learnErrors] Writing counts: {countsOutPath}"
      writeLearnErrorsCountsTsv(countsOutPath, learnRes)

    if summaryPath.len > 0:
      writeSummary(summaryPath, batchRes, learnRes, derepOpts, learnOpts, scOpts)

    if not quiet:
      printSummary(batchRes, learnRes, derepOpts, learnOpts, scOpts)
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
