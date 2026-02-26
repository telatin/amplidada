import std/[algorithm, locks, os, parsecsv, sets, strformat, strutils, tables]
import malebolgia
import amplidada/fastq_filter

type
  PairedSampleInput* = object
    sampleId*: string
    readsForward*: string
    readsReverse*: string

  PairedSampleJob* = object
    sampleId*: string
    readsForward*: string
    readsReverse*: string
    outForward*: string
    outReverse*: string

  SampleSheetOptions* = object
    delimiter*: char
    hasHeader*: bool
    sampleIdColumn*: string
    forwardColumn*: string
    reverseColumn*: string

  BatchFilterOptions* = object
    sampleThreads*: int
    parallelMinSamples*: int
    continueOnError*: bool
    verbose*: bool

  PairedSampleResult* = object
    sampleId*: string
    job*: PairedSampleJob
    success*: bool
    errorMessage*: string
    stats*: FastqPairedFilterStats

  BatchPairedFilterResult* = object
    samples*: seq[PairedSampleResult]
    totals*: FastqPairedFilterStats
    failedSamples*: int


proc defaultSampleSheetOptions*(): SampleSheetOptions =
  SampleSheetOptions(
    delimiter: ',',
    hasHeader: true,
    sampleIdColumn: "SampleID",
    forwardColumn: "reads_forward",
    reverseColumn: "reads_reverse"
  )

proc defaultBatchFilterOptions*(): BatchFilterOptions =
  BatchFilterOptions(
    sampleThreads: 1,
    parallelMinSamples: 2,
    continueOnError: false,
    verbose: false
  )

proc validateOptions*(opts: SampleSheetOptions) =
  if opts.delimiter == '\0':
    raise newException(ValueError, "Sample sheet delimiter must be a valid character")
  if opts.sampleIdColumn.strip().len == 0:
    raise newException(ValueError, "sampleIdColumn must not be empty")
  if opts.forwardColumn.strip().len == 0:
    raise newException(ValueError, "forwardColumn must not be empty")
  if opts.reverseColumn.strip().len == 0:
    raise newException(ValueError, "reverseColumn must not be empty")

proc validateOptions*(opts: BatchFilterOptions) =
  if opts.sampleThreads <= 0:
    raise newException(ValueError, "sampleThreads must be > 0")
  if opts.parallelMinSamples <= 0:
    raise newException(ValueError, "parallelMinSamples must be > 0")

proc fastqExt*(path: string): bool =
  let lower = path.toLowerAscii()
  lower.endsWith(".fastq") or lower.endsWith(".fastq.gz") or
    lower.endsWith(".fq") or lower.endsWith(".fq.gz")

proc removeFastqExt(filename: string): string =
  result = filename
  let lower = filename.toLowerAscii()
  if lower.endsWith(".fastq.gz"):
    result.setLen(filename.len - ".fastq.gz".len)
  elif lower.endsWith(".fq.gz"):
    result.setLen(filename.len - ".fq.gz".len)
  elif lower.endsWith(".fastq"):
    result.setLen(filename.len - ".fastq".len)
  elif lower.endsWith(".fq"):
    result.setLen(filename.len - ".fq".len)

proc sanitizeSampleId*(sampleId: string): string =
  var clean = newStringOfCap(sampleId.len)
  for c in sampleId:
    if c.isAlphaNumeric() or c in {'_', '-', '.'}:
      clean.add(c)
    else:
      clean.add('_')
  result = clean.strip(chars = {'_'})
  if result.len == 0:
    result = "sample"

proc deriveSampleIdFromForward*(forwardPath: string; r1Token = "_R1"): string =
  let filename = forwardPath.extractFilename()
  let idx = filename.find(r1Token)
  if idx >= 0:
    result = filename[0 ..< idx]
  else:
    result = removeFastqExt(filename)
  if result.len == 0:
    result = removeFastqExt(filename)

proc addDiscoveredCandidate(
  samples: var seq[PairedSampleInput],
  seenIds: var HashSet[string],
  filePath: string,
  r1Token: string,
  r2Token: string,
  requireMate: bool
) =
  if not fastqExt(filePath):
    return

  let fileName = filePath.extractFilename()
  let idx = fileName.find(r1Token)
  if idx < 0:
    return

  let mateName = fileName[0 ..< idx] & r2Token & fileName[idx + r1Token.len .. ^1]
  let matePath = joinPath(filePath.parentDir(), mateName)

  if not fileExists(matePath):
    if requireMate:
      raise newException(IOError, "Missing mate for " & filePath & ": expected " & matePath)
    return

  let sampleId = deriveSampleIdFromForward(filePath, r1Token)
  if sampleId in seenIds:
    raise newException(ValueError, "Duplicate discovered sample ID: " & sampleId)

  seenIds.incl(sampleId)
  samples.add(PairedSampleInput(
    sampleId: sampleId,
    readsForward: absolutePath(filePath),
    readsReverse: absolutePath(matePath)
  ))

proc discoverPairedSamples*(
  inputDir: string,
  r1Token = "_R1",
  r2Token = "_R2",
  recursive = false,
  requireMate = true
): seq[PairedSampleInput] =
  if not dirExists(inputDir):
    raise newException(IOError, "Input directory does not exist: " & inputDir)

  if r1Token.len == 0 or r2Token.len == 0:
    raise newException(ValueError, "r1Token and r2Token must not be empty")

  var seenIds = initHashSet[string]()

  if recursive:
    for filePath in walkDirRec(inputDir):
      addDiscoveredCandidate(result, seenIds, filePath, r1Token, r2Token, requireMate)
  else:
    for kind, filePath in walkDir(inputDir):
      if kind == pcFile:
        addDiscoveredCandidate(result, seenIds, filePath, r1Token, r2Token, requireMate)

  result.sort(proc(a, b: PairedSampleInput): int = cmp(a.sampleId, b.sampleId))

proc normalizeHeader(s: string): string =
  s.strip().toLowerAscii()

proc readPairedSampleSheet*(sheetPath: string, opts = defaultSampleSheetOptions()): seq[PairedSampleInput] =
  validateOptions(opts)

  if not fileExists(sheetPath):
    raise newException(IOError, "Sample sheet does not exist: " & sheetPath)

  let sheetDir = sheetPath.parentDir()
  var parser: CsvParser
  parser.open(sheetPath, separator = opts.delimiter)
  defer: parser.close()

  var colSample = 0
  var colForward = 1
  var colReverse = 2

  if opts.hasHeader:
    if not parser.readRow():
      return @[]

    var byName = initTable[string, int]()
    for i, h in parser.row:
      byName[normalizeHeader(h)] = i

    let sampleKey = normalizeHeader(opts.sampleIdColumn)
    let forwardKey = normalizeHeader(opts.forwardColumn)
    let reverseKey = normalizeHeader(opts.reverseColumn)

    if sampleKey notin byName:
      raise newException(ValueError, "Sample sheet missing column: " & opts.sampleIdColumn)
    if forwardKey notin byName:
      raise newException(ValueError, "Sample sheet missing column: " & opts.forwardColumn)
    if reverseKey notin byName:
      raise newException(ValueError, "Sample sheet missing column: " & opts.reverseColumn)

    colSample = byName[sampleKey]
    colForward = byName[forwardKey]
    colReverse = byName[reverseKey]

  var seenIds = initHashSet[string]()

  while parser.readRow():
    if parser.row.len == 0:
      continue

    if max(colSample, max(colForward, colReverse)) >= parser.row.len:
      raise newException(ValueError, "Malformed sample sheet row, expected at least 3 columns")

    var sampleId = parser.row[colSample].strip()
    var readsForward = parser.row[colForward].strip()
    var readsReverse = parser.row[colReverse].strip()

    if readsForward.len == 0 or readsReverse.len == 0:
      raise newException(ValueError, "Sample sheet contains empty read path")

    if not readsForward.isAbsolute():
      readsForward = joinPath(sheetDir, readsForward)
    if not readsReverse.isAbsolute():
      readsReverse = joinPath(sheetDir, readsReverse)

    readsForward = absolutePath(readsForward)
    readsReverse = absolutePath(readsReverse)

    if sampleId.len == 0:
      sampleId = deriveSampleIdFromForward(readsForward)

    if sampleId in seenIds:
      raise newException(ValueError, "Duplicate sample ID in sample sheet: " & sampleId)

    seenIds.incl(sampleId)
    result.add(PairedSampleInput(
      sampleId: sampleId,
      readsForward: readsForward,
      readsReverse: readsReverse
    ))

proc buildPairedJobs*(
  samples: openArray[PairedSampleInput],
  outputDir: string,
  outForwardSuffix = "_R1.filtered.fastq.gz",
  outReverseSuffix = "_R2.filtered.fastq.gz"
): seq[PairedSampleJob] =
  if outputDir.len == 0:
    raise newException(ValueError, "outputDir must not be empty")

  if outForwardSuffix.len == 0 or outReverseSuffix.len == 0:
    raise newException(ValueError, "Output suffixes must not be empty")

  for sample in samples:
    if sample.readsForward.len == 0 or sample.readsReverse.len == 0:
      raise newException(ValueError, "Sample has empty input read path")

    var sampleId = sample.sampleId
    if sampleId.len == 0:
      sampleId = deriveSampleIdFromForward(sample.readsForward)

    let cleanId = sanitizeSampleId(sampleId)
    result.add(PairedSampleJob(
      sampleId: sampleId,
      readsForward: sample.readsForward,
      readsReverse: sample.readsReverse,
      outForward: joinPath(outputDir, cleanId & outForwardSuffix),
      outReverse: joinPath(outputDir, cleanId & outReverseSuffix)
    ))

proc mergeStats(total: var FastqPairedFilterStats, x: FastqPairedFilterStats) =
  total.inputPairs += x.inputPairs
  total.outputPairs += x.outputPairs
  total.discardedPairs += x.discardedPairs

  total.forwardStats.inputReads += x.forwardStats.inputReads
  total.forwardStats.outputReads += x.forwardStats.outputReads
  total.forwardStats.discardedReads += x.forwardStats.discardedReads

  total.reverseStats.inputReads += x.reverseStats.inputReads
  total.reverseStats.outputReads += x.reverseStats.outputReads
  total.reverseStats.discardedReads += x.reverseStats.discardedReads

  for reason in FilterRejectReason:
    total.forwardStats.rejectedByReason[reason] += x.forwardStats.rejectedByReason[reason]
    total.reverseStats.rejectedByReason[reason] += x.reverseStats.rejectedByReason[reason]

  for reason in PairRejectReason:
    total.rejectedByPairReason[reason] += x.rejectedByPairReason[reason]

proc verboseLog(lock: ptr Lock; msg: string) =
  if lock != nil:
    acquire(lock[])
    defer: release(lock[])
  stdout.writeLine(msg)

proc runSingleJob(
  job: PairedSampleJob,
  filterOpts: PairedFastqFilterOptions,
  sampleIndex: int,
  totalSamples: int,
  verbose: bool,
  logLock: ptr Lock
): PairedSampleResult =
  result.sampleId = job.sampleId
  result.job = job

  if verbose:
    verboseLog(logLock, &"[filterAndTrim] Starting sample {sampleIndex + 1}/{totalSamples}: {job.sampleId}")

  try:
    if job.outForward.parentDir().len > 0:
      discard existsOrCreateDir(job.outForward.parentDir())
    if job.outReverse.parentDir().len > 0:
      discard existsOrCreateDir(job.outReverse.parentDir())

    result.stats = fastqPairedFilter(
      inputPath1 = job.readsForward,
      inputPath2 = job.readsReverse,
      outputPath1 = job.outForward,
      outputPath2 = job.outReverse,
      opts = filterOpts
    )
    result.success = true
    if verbose:
      verboseLog(logLock, &"[filterAndTrim] Completed sample {sampleIndex + 1}/{totalSamples}: {job.sampleId} pairs_in={result.stats.inputPairs} pairs_out={result.stats.outputPairs}")
  except CatchableError as e:
    result.success = false
    result.errorMessage = e.msg
    if verbose:
      verboseLog(logLock, &"[filterAndTrim] Failed sample {sampleIndex + 1}/{totalSamples}: {job.sampleId} error={e.msg}")

template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc runRange(
  jobs: ptr UncheckedArray[PairedSampleJob],
  results: ptr UncheckedArray[PairedSampleResult],
  startIndex: int,
  totalSamples: int,
  n: int,
  filterOpts: PairedFastqFilterOptions,
  verbose: bool,
  logLock: ptr Lock
) =
  for i in 0..<n:
    results[i] = runSingleJob(jobs[i], filterOpts, startIndex + i, totalSamples, verbose, logLock)

proc filterAndTrimPaired*(
  jobs: openArray[PairedSampleJob],
  filterOpts: PairedFastqFilterOptions,
  batchOpts = defaultBatchFilterOptions()
): BatchPairedFilterResult =
  validateOptions(filterOpts)
  validateOptions(batchOpts)

  var work = @jobs
  result.samples = newSeq[PairedSampleResult](work.len)
  var logLock: Lock
  var lockPtr: ptr Lock = nil

  if work.len == 0:
    return

  if batchOpts.verbose:
    initLock(logLock)
    lockPtr = addr(logLock)
    defer: deinitLock(logLock)

  if batchOpts.sampleThreads == 1 or work.len < batchOpts.parallelMinSamples:
    for i in 0..<work.len:
      result.samples[i] = runSingleJob(work[i], filterOpts, i, work.len, batchOpts.verbose, lockPtr)
  else:
    let workerCount = min(batchOpts.sampleThreads, work.len)
    let blockSize = (work.len + workerCount - 1) div workerCount

    var m = createMaster()
    m.awaitAll:
      var idx = 0
      while idx < work.len:
        let n = min(blockSize, work.len - idx)
        m.spawn runRange(atPtr(work, idx), atPtr(result.samples, idx), idx, work.len, n, filterOpts, batchOpts.verbose, lockPtr)
        idx += n

  var firstError = ""
  for sampleResult in result.samples:
    if sampleResult.success:
      mergeStats(result.totals, sampleResult.stats)
    else:
      inc result.failedSamples
      if firstError.len == 0:
        firstError = &"sample={sampleResult.sampleId}: {sampleResult.errorMessage}"

  if result.failedSamples > 0 and not batchOpts.continueOnError:
    raise newException(IOError, "Batch filtering failed: " & firstError)

proc filterAndTrimPaired*(
  samples: openArray[PairedSampleInput],
  outputDir: string,
  filterOpts: PairedFastqFilterOptions,
  batchOpts = defaultBatchFilterOptions(),
  outForwardSuffix = "_R1.filtered.fastq.gz",
  outReverseSuffix = "_R2.filtered.fastq.gz"
): BatchPairedFilterResult =
  let jobs = buildPairedJobs(samples, outputDir, outForwardSuffix, outReverseSuffix)
  filterAndTrimPaired(jobs, filterOpts, batchOpts)

proc filterAndTrimPairedDir*(
  inputDir: string,
  outputDir: string,
  filterOpts: PairedFastqFilterOptions,
  batchOpts = defaultBatchFilterOptions(),
  r1Token = "_R1",
  r2Token = "_R2",
  recursive = false,
  outForwardSuffix = "_R1.filtered.fastq.gz",
  outReverseSuffix = "_R2.filtered.fastq.gz"
): BatchPairedFilterResult =
  let samples = discoverPairedSamples(inputDir, r1Token = r1Token, r2Token = r2Token, recursive = recursive)
  filterAndTrimPaired(samples, outputDir, filterOpts, batchOpts, outForwardSuffix, outReverseSuffix)

proc filterAndTrimPairedSheet*(
  sheetPath: string,
  outputDir: string,
  filterOpts: PairedFastqFilterOptions,
  sheetOpts = defaultSampleSheetOptions(),
  batchOpts = defaultBatchFilterOptions(),
  outForwardSuffix = "_R1.filtered.fastq.gz",
  outReverseSuffix = "_R2.filtered.fastq.gz"
): BatchPairedFilterResult =
  let samples = readPairedSampleSheet(sheetPath, sheetOpts)
  filterAndTrimPaired(samples, outputDir, filterOpts, batchOpts, outForwardSuffix, outReverseSuffix)
