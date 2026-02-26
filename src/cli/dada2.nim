import std/[os, parseopt, sets, strformat, strutils]
import amplidada

proc usage() =
  echo "Usage:"
  echo "  Single-end: dada2 --in <reads.fastq[.gz]> --out <asv.tsv> [options]"
  echo "  Paired-end: dada2 --in-forward <R1.fastq[.gz]> --in-reverse <R2.fastq[.gz]> --out <merged.tsv> [options]"
  echo ""
  echo "Inputs:"
  echo "  --in <path>                      Single-end input FASTQ/FASTQ.GZ"
  echo "  --in-forward <path>              Forward reads FASTQ/FASTQ.GZ"
  echo "  --in-reverse <path>              Reverse reads FASTQ/FASTQ.GZ"
  echo "  --out <path>                     Output TSV (ASVs for single-end, merged ASVs for paired-end)"
  echo ""
  echo "Error model:"
  echo "  --err-matrix <path>              Error matrix TSV for all reads (single-end, or both directions)"
  echo "  --err-matrix-forward <path>      Error matrix TSV for forward reads"
  echo "  --err-matrix-reverse <path>      Error matrix TSV for reverse reads"
  echo "  --err-out <path>                 Write matrix used for denoising (single-end)"
  echo "  --err-out-forward <path>         Write matrix used for forward denoising"
  echo "  --err-out-reverse <path>         Write matrix used for reverse denoising"
  echo "  --learn-self-consist             Learn errors with iterative self-consistency"
  echo "  --learn-min-iter <int>           Minimum learnErrors self-consistency iterations"
  echo "  --learn-max-iter <int>           Maximum learnErrors self-consistency iterations"
  echo "  --learn-tol <float>              Convergence tolerance for learnErrors self-consistency"
  echo "  --learn-max-centers-per-length <int> Candidate centers per length in learnErrors self-consistency"
  echo "  --q-min <int>                    learnErrors min quality bin"
  echo "  --q-max <int>                    learnErrors max quality bin"
  echo "  --pseudocount <float>            learnErrors additive smoothing"
  echo "  --zero-mismatch-prob <float>     learnErrors fallback mismatch probability"
  echo ""
  echo "Derep options:"
  echo "  --threads <int>                  Worker threads for derep + dada"
  echo "  --chunk-size <int>               Reads per derep chunk"
  echo "  --parallel-min-chunk <int>       Min derep chunk for multithreading"
  echo "  --order <abundance|first>        Derep ordering"
  echo ""
  echo "Denoising options:"
  echo "  --omega-a <float>                Bonferroni-corrected split threshold"
  echo "  --omega-c <float>                Post-partition correction threshold (0 = correct none)"
  echo "  --min-abundance <int>            Minimum abundance to consider split"
  echo "  --min-hamming <int>              Minimum Hamming distance from center"
  echo "  --min-fold <float>               Minimum fold over expectation to split"
  echo "  --max-clusters <int>             Max ASV clusters per length (0 = unlimited)"
  echo "  --max-iter <int>                 Max denoising outer iterations"
  echo "  --max-shuffle <int>              Max assignment shuffle rounds per iteration"
  echo "  --no-qual                        Ignore quality profiles during denoising"
  echo "  --parallel-min-groups <int>      Min length-groups for parallel denoising"
  echo "  --dada-self-consist              Run canonical dada<->error self-consistency loop"
  echo "  --dada-self-min-iter <int>       Minimum canonical self-consistency iterations"
  echo "  --dada-self-max-iter <int>       Maximum canonical self-consistency iterations"
  echo "  --dada-self-tol <float>          Convergence tolerance for canonical self-consistency"
  echo "  --dada-self-keep-omega-c         Keep configured omegaC during self-consistency iterations"
  echo "  --dada-self-no-final-denoise     Skip final denoise pass with latest learned errors"
  echo ""
  echo "Paired-end merge options:"
  echo "  --min-overlap <int>              Minimum overlap length"
  echo "  --max-mismatch <int>             Maximum mismatches in overlap"
  echo "  --trim-overhang                  Allow and trim overhangs when merging"
  echo "  --merge-threads <int>            Worker threads for mergePairs"
  echo "  --merge-parallel-min-pairs <int> Min unique ASV-pairs for parallel merge"
  echo ""
  echo "Outputs and reporting:"
  echo "  --summary <path>                 Write run summary"
  echo "  --verbose                        Print progress details"
  echo "  --quiet                          Suppress summary output"
  echo ""
  echo "General:"
  echo "  -h, --help"
  echo "  -v, --version"

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Error: " & msg)
  quit(2)

proc parseIntFlag(flag, value: string): int =
  try:
    parseInt(value)
  except ValueError:
    raise newException(ValueError, "Invalid integer for --" & flag & ": " & value)

proc parseFloatFlag(flag, value: string): float64 =
  try:
    parseFloat(value)
  except ValueError:
    raise newException(ValueError, "Invalid float for --" & flag & ": " & value)

proc normalizeArgs(args: seq[string]): seq[string] =
  let longNeedsValue = toHashSet([
    "--in", "--in-forward", "--in-reverse", "--out", "--summary",
    "--err-matrix", "--err-matrix-forward", "--err-matrix-reverse",
    "--err-out", "--err-out-forward", "--err-out-reverse",
    "--q-min", "--q-max", "--pseudocount", "--zero-mismatch-prob",
    "--learn-min-iter", "--learn-max-iter", "--learn-tol", "--learn-max-centers-per-length",
    "--dada-self-min-iter", "--dada-self-max-iter", "--dada-self-tol",
    "--threads", "--chunk-size", "--parallel-min-chunk", "--order",
    "--omega-a", "--omega-c", "--min-abundance", "--min-hamming", "--min-fold", "--max-clusters",
    "--max-iter", "--max-shuffle", "--parallel-min-groups",
    "--min-overlap", "--max-mismatch", "--merge-threads", "--merge-parallel-min-pairs"
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

proc makeLearnBatchFromLoadedMatrix(derepRes: DerepFastqResult): LearnErrorsBatchResult =
  result.filesRequested = 1
  result.filesUsed = 1
  result.readsUsed = derepRes.totalReads
  result.basesUsed = derepBaseCount(derepRes)
  result.uniqueCountSum = derepRes.uniqueCount
  result.selfConsistEnabled = false
  result.iterationsRun = 0
  result.converged = true
  result.maxAbsProbDelta = 0.0
  result.centerChangesTotal = 0

proc learnOrLoadErrorMatrix(
  label: string,
  derepRes: DerepFastqResult,
  matrixPath: string,
  outPath: string,
  learnOpts: LearnErrorsOptions,
  learnScOpts: LearnErrorsSelfConsistOptions,
  verbose: bool
): tuple[matrix: LearnErrorsResult, batch: LearnErrorsBatchResult, loadedFromFile: bool] =
  if matrixPath.len > 0:
    if verbose:
      echo &"[dada2] Loading {label} error matrix: {matrixPath}"
    result.matrix = readLearnErrorsTsv(matrixPath)
    result.batch = makeLearnBatchFromLoadedMatrix(derepRes)
    result.loadedFromFile = true
  else:
    if verbose:
      echo &"[dada2] Learning {label} error matrix"
    let scRes = learnErrorsSelfConsistentFromDereps([derepRes], learnOpts, learnScOpts)
    result.matrix = scRes.matrix
    result.batch.filesRequested = 1
    result.batch.filesUsed = 1
    result.batch.readsUsed = derepRes.totalReads
    result.batch.basesUsed = derepBaseCount(derepRes)
    result.batch.uniqueCountSum = derepRes.uniqueCount
    result.batch.selfConsistEnabled = learnScOpts.enabled
    result.batch.iterationsRun = scRes.iterationsRun
    result.batch.converged = scRes.converged
    result.batch.maxAbsProbDelta = scRes.maxAbsProbDelta
    result.batch.centerChangesTotal = scRes.centerChangesTotal
    result.loadedFromFile = false

  if outPath.len > 0:
    if verbose:
      echo &"[dada2] Writing {label} matrix: {outPath}"
    writeLearnErrorsTsv(outPath, result.matrix)

proc writeSingleSummary(
  path: string,
  inputPath: string,
  derepRes: DerepFastqResult,
  learnBatch: LearnErrorsBatchResult,
  learnOpts: LearnErrorsOptions,
  learnScOpts: LearnErrorsSelfConsistOptions,
  dadaScOpts: DadaSelfConsistOptions,
  dadaScRan: bool,
  dadaScRes: DadaSelfConsistResult,
  dadaRes: DadaResult,
  dadaOpts: DadaOptions,
  loadedFromFile: bool
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open summary output file: " & path)
  defer: close(f)

  f.writeLine("mode=single")
  f.writeLine("input=" & absolutePath(inputPath))
  f.writeLine("input_reads=" & $derepRes.totalReads)
  f.writeLine("input_unique_sequences=" & $derepRes.uniqueCount)
  f.writeLine("learn_reads_used=" & $learnBatch.readsUsed)
  f.writeLine("learn_bases_used=" & $learnBatch.basesUsed)
  f.writeLine("learn_self_consist_enabled=" & $learnBatch.selfConsistEnabled)
  f.writeLine("learn_self_consist_iterations=" & $learnBatch.iterationsRun)
  f.writeLine("learn_self_consist_converged=" & $learnBatch.converged)
  f.writeLine("dada_self_consist_enabled=" & $dadaScOpts.enabled)
  f.writeLine("dada_self_consist_ran=" & $dadaScRan)
  f.writeLine("dada_self_consist_iterations=" & $dadaScRes.iterationsRun)
  f.writeLine("dada_self_consist_converged=" & $dadaScRes.converged)
  f.writeLine("dada_self_consist_max_abs_prob_delta=" & $dadaScRes.maxAbsProbDelta)
  f.writeLine("dada_self_consist_assignment_changes=" & $dadaScRes.assignmentChanges)
  f.writeLine("err_loaded_from_file=" & $loadedFromFile)
  f.writeLine("asv_count=" & $dadaRes.asvCount)
  f.writeLine("corrected_uniques=" & $dadaRes.correctedUniques)
  f.writeLine("uncorrected_uniques=" & $dadaRes.uncorrectedUniques)
  f.writeLine("denoise_group_count=" & $dadaRes.groupCount)
  f.writeLine("denoise_iterations_total=" & $dadaRes.iterationsRun)
  f.writeLine("denoise_splits_accepted=" & $dadaRes.splitsAccepted)
  f.writeLine("learn_options=" & learnOpts.describe())
  f.writeLine("learn_self_consist_options=" & learnScOpts.describe())
  f.writeLine("dada_self_consist_options=" & dadaScOpts.describe())
  f.writeLine("dada_options=" & dadaOpts.describe())

proc writePairedSummary(
  path: string,
  inputForward: string,
  inputReverse: string,
  derepForward: DerepFastqResult,
  derepReverse: DerepFastqResult,
  learnBatchForward: LearnErrorsBatchResult,
  learnBatchReverse: LearnErrorsBatchResult,
  learnOpts: LearnErrorsOptions,
  learnScOpts: LearnErrorsSelfConsistOptions,
  dadaScOpts: DadaSelfConsistOptions,
  dadaForward: DadaResult,
  dadaReverse: DadaResult,
  forwardDadaScRan: bool,
  reverseDadaScRan: bool,
  forwardDadaScRes: DadaSelfConsistResult,
  reverseDadaScRes: DadaSelfConsistResult,
  dadaOpts: DadaOptions,
  mergeRes: MergePairsResult,
  mergeOpts: MergePairsOptions,
  forwardLoadedFromFile: bool,
  reverseLoadedFromFile: bool
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open summary output file: " & path)
  defer: close(f)

  f.writeLine("mode=paired")
  f.writeLine("input_forward=" & absolutePath(inputForward))
  f.writeLine("input_reverse=" & absolutePath(inputReverse))
  f.writeLine("forward_reads=" & $derepForward.totalReads)
  f.writeLine("reverse_reads=" & $derepReverse.totalReads)
  f.writeLine("forward_unique_sequences=" & $derepForward.uniqueCount)
  f.writeLine("reverse_unique_sequences=" & $derepReverse.uniqueCount)
  f.writeLine("forward_err_loaded_from_file=" & $forwardLoadedFromFile)
  f.writeLine("reverse_err_loaded_from_file=" & $reverseLoadedFromFile)
  f.writeLine("dada_self_consist_enabled=" & $dadaScOpts.enabled)
  f.writeLine("forward_dada_self_consist_ran=" & $forwardDadaScRan)
  f.writeLine("reverse_dada_self_consist_ran=" & $reverseDadaScRan)
  f.writeLine("forward_dada_self_consist_iterations=" & $forwardDadaScRes.iterationsRun)
  f.writeLine("reverse_dada_self_consist_iterations=" & $reverseDadaScRes.iterationsRun)
  f.writeLine("forward_dada_self_consist_converged=" & $forwardDadaScRes.converged)
  f.writeLine("reverse_dada_self_consist_converged=" & $reverseDadaScRes.converged)
  f.writeLine("forward_dada_self_consist_max_abs_prob_delta=" & $forwardDadaScRes.maxAbsProbDelta)
  f.writeLine("reverse_dada_self_consist_max_abs_prob_delta=" & $reverseDadaScRes.maxAbsProbDelta)
  f.writeLine("forward_dada_self_consist_assignment_changes=" & $forwardDadaScRes.assignmentChanges)
  f.writeLine("reverse_dada_self_consist_assignment_changes=" & $reverseDadaScRes.assignmentChanges)
  f.writeLine("forward_asv_count=" & $dadaForward.asvCount)
  f.writeLine("reverse_asv_count=" & $dadaReverse.asvCount)
  f.writeLine("forward_corrected_uniques=" & $dadaForward.correctedUniques)
  f.writeLine("reverse_corrected_uniques=" & $dadaReverse.correctedUniques)
  f.writeLine("merged_count=" & $mergeRes.merged.len)
  f.writeLine("merged_pairs=" & $mergeRes.stats.mergedPairs)
  f.writeLine("merge_candidate_pairs=" & $mergeRes.stats.candidatePairs)
  f.writeLine("merge_rejected_pairs=" & $mergeRes.stats.rejectedPairs)
  f.writeLine("merge_unsynchronized_reads=" & $mergeRes.stats.unsynchronizedReads)
  f.writeLine("merge_unique_asv_pairs=" & $mergeRes.stats.uniqueAsvPairs)
  f.writeLine("merge_unique_sequences=" & $mergeRes.stats.uniqueMergedSequences)
  for reason in MergeRejectReason:
    f.writeLine(&"merge_reject_{reason}={mergeRes.stats.rejectedByReason[reason]}")
  f.writeLine("learn_forward_reads_used=" & $learnBatchForward.readsUsed)
  f.writeLine("learn_reverse_reads_used=" & $learnBatchReverse.readsUsed)
  f.writeLine("learn_options=" & learnOpts.describe())
  f.writeLine("learn_self_consist_options=" & learnScOpts.describe())
  f.writeLine("dada_self_consist_options=" & dadaScOpts.describe())
  f.writeLine("dada_options=" & dadaOpts.describe())
  f.writeLine("merge_options=" & mergeOpts.describe())

proc printSingleSummary(
  derepRes: DerepFastqResult,
  learnBatch: LearnErrorsBatchResult,
  learnOpts: LearnErrorsOptions,
  learnScOpts: LearnErrorsSelfConsistOptions,
  dadaScOpts: DadaSelfConsistOptions,
  dadaScRan: bool,
  dadaScRes: DadaSelfConsistResult,
  dadaRes: DadaResult,
  dadaOpts: DadaOptions,
  loadedFromFile: bool
) =
  echo "Mode: single-end denoising"
  echo &"Input reads: {derepRes.totalReads}"
  echo &"Input unique sequences: {derepRes.uniqueCount}"
  echo &"Error model loaded from file: {loadedFromFile}"
  echo &"learnErrors iterations: {learnBatch.iterationsRun} converged={learnBatch.converged}"
  echo &"dada self-consistency enabled={dadaScOpts.enabled} ran={dadaScRan} iterations={dadaScRes.iterationsRun} converged={dadaScRes.converged}"
  echo &"ASV count: {dadaRes.asvCount}"
  echo &"Corrected uniques: {dadaRes.correctedUniques} uncorrected uniques: {dadaRes.uncorrectedUniques}"
  echo &"learnErrors options: {learnOpts.describe()}"
  echo &"learnErrors self-consistency options: {learnScOpts.describe()}"
  echo &"dada self-consistency options: {dadaScOpts.describe()}"
  echo &"dada options: {dadaOpts.describe()}"

proc printPairedSummary(
  derepForward: DerepFastqResult,
  derepReverse: DerepFastqResult,
  dadaForward: DadaResult,
  dadaReverse: DadaResult,
  mergeRes: MergePairsResult,
  mergeOpts: MergePairsOptions,
  forwardLoadedFromFile: bool,
  reverseLoadedFromFile: bool
) =
  echo "Mode: paired-end denoising + mergePairs"
  echo &"Forward reads: {derepForward.totalReads} unique={derepForward.uniqueCount}"
  echo &"Reverse reads: {derepReverse.totalReads} unique={derepReverse.uniqueCount}"
  echo &"Forward error matrix from file: {forwardLoadedFromFile}"
  echo &"Reverse error matrix from file: {reverseLoadedFromFile}"
  echo &"Forward ASVs: {dadaForward.asvCount} reverse ASVs: {dadaReverse.asvCount}"
  echo &"Merged sequences: {mergeRes.merged.len}"
  echo &"Merged pairs: {mergeRes.stats.mergedPairs} rejected: {mergeRes.stats.rejectedPairs} candidate: {mergeRes.stats.candidatePairs}"
  echo &"mergePairs options: {mergeOpts.describe()}"

proc defaultPairedErrOut(basePath: string; suffix: string): string =
  if basePath.len == 0:
    return ""
  let p = splitFile(basePath)
  let ext =
    if p.ext.len > 0: p.ext else: ".tsv"
  joinPath(p.dir, p.name & "." & suffix & ext)

proc main() =
  var inputPath = ""
  var inputForward = ""
  var inputReverse = ""
  var outPath = ""
  var summaryPath = ""

  var errMatrixPath = ""
  var errMatrixForwardPath = ""
  var errMatrixReversePath = ""
  var errOutPath = ""
  var errOutForwardPath = ""
  var errOutReversePath = ""

  var quiet = false
  var verbose = false

  var derepOpts = defaultDerepFastqOptions()
  var learnOpts = defaultLearnErrorsOptions()
  var learnScOpts = defaultLearnErrorsSelfConsistOptions()
  var dadaOpts = defaultDadaOptions()
  var dadaScOpts = defaultDadaSelfConsistOptions()
  var mergeOpts = defaultMergePairsOptions()

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
          inputPath = value
        of "in-forward":
          inputForward = value
        of "in-reverse":
          inputReverse = value
        of "out":
          outPath = value
        of "summary":
          summaryPath = value
        of "err-matrix":
          errMatrixPath = value
        of "err-matrix-forward":
          errMatrixForwardPath = value
        of "err-matrix-reverse":
          errMatrixReversePath = value
        of "err-out":
          errOutPath = value
        of "err-out-forward":
          errOutForwardPath = value
        of "err-out-reverse":
          errOutReversePath = value
        of "threads":
          let t = parseIntFlag(flag, value)
          derepOpts.threads = t
          dadaOpts.threads = t
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
        of "learn-self-consist":
          learnScOpts.enabled = true
        of "learn-min-iter":
          learnScOpts.minIterations = parseIntFlag(flag, value)
        of "learn-max-iter":
          learnScOpts.maxIterations = parseIntFlag(flag, value)
        of "learn-tol":
          learnScOpts.tol = parseFloatFlag(flag, value)
        of "learn-max-centers-per-length":
          learnScOpts.maxCentersPerLength = parseIntFlag(flag, value)
        of "dada-self-consist":
          dadaScOpts.enabled = true
        of "dada-self-min-iter":
          dadaScOpts.minIterations = parseIntFlag(flag, value)
        of "dada-self-max-iter":
          dadaScOpts.maxIterations = parseIntFlag(flag, value)
        of "dada-self-tol":
          dadaScOpts.tol = parseFloatFlag(flag, value)
        of "dada-self-keep-omega-c":
          dadaScOpts.forceOmegaCZero = false
        of "dada-self-no-final-denoise":
          dadaScOpts.finalDenoiseWithLatestErrors = false
        of "omega-a":
          dadaOpts.omegaA = parseFloatFlag(flag, value)
        of "omega-c":
          dadaOpts.omegaC = parseFloatFlag(flag, value)
        of "min-abundance":
          dadaOpts.minAbundance = parseIntFlag(flag, value).int64
        of "min-hamming":
          dadaOpts.minHamming = parseIntFlag(flag, value)
        of "min-fold":
          dadaOpts.minFold = parseFloatFlag(flag, value)
        of "max-clusters":
          dadaOpts.maxClusters = parseIntFlag(flag, value)
        of "max-iter":
          dadaOpts.maxIterations = parseIntFlag(flag, value)
        of "max-shuffle":
          dadaOpts.maxShuffle = parseIntFlag(flag, value)
        of "no-qual":
          dadaOpts.useQuality = false
        of "parallel-min-groups":
          dadaOpts.parallelMinGroups = parseIntFlag(flag, value)
        of "min-overlap":
          mergeOpts.minOverlap = parseIntFlag(flag, value)
        of "max-mismatch":
          mergeOpts.maxMismatch = parseIntFlag(flag, value)
        of "trim-overhang":
          mergeOpts.trimOverhang = true
        of "merge-threads":
          mergeOpts.threads = parseIntFlag(flag, value)
        of "merge-parallel-min-pairs":
          mergeOpts.parallelMinPairs = parseIntFlag(flag, value)
        of "verbose":
          verbose = true
        of "quiet":
          quiet = true
        else:
          fail("Unknown option: --" & key)
      of cmdArgument:
        fail("Unexpected positional argument: " & key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  let pairedMode = inputForward.len > 0 or inputReverse.len > 0
  if outPath.len == 0:
    fail("Missing required option: --out")
  if quiet and verbose:
    fail("--quiet and --verbose are mutually exclusive")

  if pairedMode:
    if inputForward.len == 0 or inputReverse.len == 0:
      fail("Paired-end mode requires both --in-forward and --in-reverse")
    if inputPath.len > 0:
      fail("Do not combine --in with paired-end inputs; use either single-end or paired-end mode")
  else:
    if inputPath.len == 0:
      fail("Missing required option: --in")

  try:
    validateOptions(derepOpts)
    validateOptions(learnOpts)
    validateOptions(learnScOpts)
    validateOptions(dadaOpts)
    validateOptions(dadaScOpts)
    validateOptions(mergeOpts)

    if pairedMode:
      var derepForwardOpts = derepOpts
      var derepReverseOpts = derepOpts
      derepForwardOpts.keepReadMap = true
      derepReverseOpts.keepReadMap = true

      if verbose:
        echo &"[dada2] Dereplicating forward reads: {inputForward}"
      let derepForward = derepFastq(inputForward, derepForwardOpts)
      if verbose:
        echo &"[dada2] Dereplicating reverse reads: {inputReverse}"
      let derepReverse = derepFastq(inputReverse, derepReverseOpts)

      let forwardMatrixPath =
        if errMatrixForwardPath.len > 0: errMatrixForwardPath
        else: errMatrixPath
      let reverseMatrixPath =
        if errMatrixReversePath.len > 0: errMatrixReversePath
        else: errMatrixPath

      let forwardErrOut =
        if errOutForwardPath.len > 0: errOutForwardPath
        elif errOutPath.len > 0: defaultPairedErrOut(errOutPath, "forward")
        else: ""
      let reverseErrOut =
        if errOutReversePath.len > 0: errOutReversePath
        elif errOutPath.len > 0: defaultPairedErrOut(errOutPath, "reverse")
        else: ""

      let forwardErr = learnOrLoadErrorMatrix(
        label = "forward",
        derepRes = derepForward,
        matrixPath = forwardMatrixPath,
        outPath = forwardErrOut,
        learnOpts = learnOpts,
        learnScOpts = learnScOpts,
        verbose = verbose
      )
      let reverseErr = learnOrLoadErrorMatrix(
        label = "reverse",
        derepRes = derepReverse,
        matrixPath = reverseMatrixPath,
        outPath = reverseErrOut,
        learnOpts = learnOpts,
        learnScOpts = learnScOpts,
        verbose = verbose
      )

      var forwardMatrix = forwardErr.matrix
      var reverseMatrix = reverseErr.matrix
      var dadaForward: DadaResult
      var dadaReverse: DadaResult
      var forwardDadaScRan = false
      var reverseDadaScRan = false
      var forwardDadaScRes = DadaSelfConsistResult()
      var reverseDadaScRes = DadaSelfConsistResult()

      if dadaScOpts.enabled and not forwardErr.loadedFromFile:
        if verbose:
          echo "[dada2] Running canonical self-consistency for forward reads"
        forwardDadaScRes = dadaSelfConsistentFromInitialErrors(
          derepRes = derepForward,
          initialErr = forwardMatrix,
          learnOpts = learnOpts,
          dadaOpts = dadaOpts,
          scOpts = dadaScOpts
        )
        forwardDadaScRan = true
        forwardMatrix = forwardDadaScRes.errorMatrix
        dadaForward = forwardDadaScRes.dada
        if forwardErrOut.len > 0:
          writeLearnErrorsTsv(forwardErrOut, forwardMatrix)
      else:
        if dadaScOpts.enabled and forwardErr.loadedFromFile and verbose:
          echo "[dada2] Skipping forward canonical self-consistency because --err-matrix was provided"
        if verbose:
          echo "[dada2] Denoising forward reads"
        dadaForward = dadaDenoise(derepForward, forwardMatrix, dadaOpts)

      if dadaScOpts.enabled and not reverseErr.loadedFromFile:
        if verbose:
          echo "[dada2] Running canonical self-consistency for reverse reads"
        reverseDadaScRes = dadaSelfConsistentFromInitialErrors(
          derepRes = derepReverse,
          initialErr = reverseMatrix,
          learnOpts = learnOpts,
          dadaOpts = dadaOpts,
          scOpts = dadaScOpts
        )
        reverseDadaScRan = true
        reverseMatrix = reverseDadaScRes.errorMatrix
        dadaReverse = reverseDadaScRes.dada
        if reverseErrOut.len > 0:
          writeLearnErrorsTsv(reverseErrOut, reverseMatrix)
      else:
        if dadaScOpts.enabled and reverseErr.loadedFromFile and verbose:
          echo "[dada2] Skipping reverse canonical self-consistency because --err-matrix was provided"
        if verbose:
          echo "[dada2] Denoising reverse reads"
        dadaReverse = dadaDenoise(derepReverse, reverseMatrix, dadaOpts)

      if verbose:
        echo "[dada2] Merging denoised paired reads"
      let merged = mergePairs(
        derepForward = derepForward,
        derepReverse = derepReverse,
        dadaForward = dadaForward,
        dadaReverse = dadaReverse,
        opts = mergeOpts
      )

      if verbose:
        echo &"[dada2] Writing merged output: {outPath}"
      writeMergedTsv(outPath, merged)

      if summaryPath.len > 0:
        writePairedSummary(
          path = summaryPath,
          inputForward = inputForward,
          inputReverse = inputReverse,
          derepForward = derepForward,
          derepReverse = derepReverse,
          learnBatchForward = forwardErr.batch,
          learnBatchReverse = reverseErr.batch,
          learnOpts = learnOpts,
          learnScOpts = learnScOpts,
          dadaScOpts = dadaScOpts,
          dadaForward = dadaForward,
          dadaReverse = dadaReverse,
          forwardDadaScRan = forwardDadaScRan,
          reverseDadaScRan = reverseDadaScRan,
          forwardDadaScRes = forwardDadaScRes,
          reverseDadaScRes = reverseDadaScRes,
          dadaOpts = dadaOpts,
          mergeRes = merged,
          mergeOpts = mergeOpts,
          forwardLoadedFromFile = forwardErr.loadedFromFile,
          reverseLoadedFromFile = reverseErr.loadedFromFile
        )

      if not quiet:
        printPairedSummary(
          derepForward = derepForward,
          derepReverse = derepReverse,
          dadaForward = dadaForward,
          dadaReverse = dadaReverse,
          mergeRes = merged,
          mergeOpts = mergeOpts,
          forwardLoadedFromFile = forwardErr.loadedFromFile,
          reverseLoadedFromFile = reverseErr.loadedFromFile
        )
    else:
      var singleDerepOpts = derepOpts
      singleDerepOpts.keepReadMap = false

      if verbose:
        echo &"[dada2] Dereplicating input: {inputPath}"
      let derepRes = derepFastq(inputPath, singleDerepOpts)

      let err = learnOrLoadErrorMatrix(
        label = "single-end",
        derepRes = derepRes,
        matrixPath = errMatrixPath,
        outPath = errOutPath,
        learnOpts = learnOpts,
        learnScOpts = learnScOpts,
        verbose = verbose
      )

      var finalErrMatrix = err.matrix
      var dadaRes: DadaResult
      var dadaScRan = false
      var dadaScRes = DadaSelfConsistResult()

      if dadaScOpts.enabled and not err.loadedFromFile:
        if verbose:
          echo "[dada2] Running canonical self-consistency denoising loop"
        dadaScRes = dadaSelfConsistentFromInitialErrors(
          derepRes = derepRes,
          initialErr = finalErrMatrix,
          learnOpts = learnOpts,
          dadaOpts = dadaOpts,
          scOpts = dadaScOpts
        )
        dadaScRan = true
        finalErrMatrix = dadaScRes.errorMatrix
        dadaRes = dadaScRes.dada
        if errOutPath.len > 0:
          writeLearnErrorsTsv(errOutPath, finalErrMatrix)
      else:
        if dadaScOpts.enabled and err.loadedFromFile and verbose:
          echo "[dada2] Skipping canonical self-consistency because --err-matrix was provided"
        if verbose:
          echo "[dada2] Running denoising"
        dadaRes = dadaDenoise(derepRes, finalErrMatrix, dadaOpts)

      if verbose:
        echo &"[dada2] Writing ASV table: {outPath}"
      writeAsvTsv(outPath, dadaRes)

      if summaryPath.len > 0:
        writeSingleSummary(
          path = summaryPath,
          inputPath = inputPath,
          derepRes = derepRes,
          learnBatch = err.batch,
          learnOpts = learnOpts,
          learnScOpts = learnScOpts,
          dadaScOpts = dadaScOpts,
          dadaScRan = dadaScRan,
          dadaScRes = dadaScRes,
          dadaRes = dadaRes,
          dadaOpts = dadaOpts,
          loadedFromFile = err.loadedFromFile
        )

      if not quiet:
        printSingleSummary(
          derepRes = derepRes,
          learnBatch = err.batch,
          learnOpts = learnOpts,
          learnScOpts = learnScOpts,
          dadaScOpts = dadaScOpts,
          dadaScRan = dadaScRan,
          dadaScRes = dadaScRes,
          dadaRes = dadaRes,
          dadaOpts = dadaOpts,
          loadedFromFile = err.loadedFromFile
        )
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
