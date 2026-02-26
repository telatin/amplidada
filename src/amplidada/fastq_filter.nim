import std/[cpuinfo, math, strformat, strutils]
import readfx
import malebolgia

when defined(windows):
  const libz = "zlib1.dll"
elif defined(macosx):
  const libz = "libz.dylib"
else:
  const libz = "libz.so.1"

type
  GzWriteFile = pointer

proc gzopen(path: cstring, mode: cstring): GzWriteFile {.cdecl, dynlib: libz, importc: "gzopen".}
proc gzwrite(theFile: GzWriteFile, buf: pointer, length: cuint): cint {.cdecl, dynlib: libz, importc: "gzwrite".}
proc gzclose(theFile: GzWriteFile): cint {.cdecl, dynlib: libz, importc: "gzclose".}

type
  FilterRejectReason* = enum
    frrNone,
    frrTooLong,
    frrTrimmedEmpty,
    frrTruncLen,
    frrTooShort,
    frrTooManyN,
    frrLowMinQ,
    frrTooManyExpectedErrors,
    frrMalformedQuality

  PairRejectReason* = enum
    prrForwardFailed,
    prrReverseFailed,
    prrBothFailed

  FastqFilterOptions* = object
    truncLen*: int
    truncQ*: int
    maxLen*: int
    maxEE*: float
    maxN*: int
    minLen*: int
    minQ*: int
    trimLeft*: int
    trimRight*: int
    chunkSize*: int
    threads*: int
    parallelMinChunk*: int

  PairedFastqFilterOptions* = object
    forward*: FastqFilterOptions
    reverse*: FastqFilterOptions
    checkNames*: bool

  FastqFilterStats* = object
    inputReads*: int64
    outputReads*: int64
    discardedReads*: int64
    rejectedByReason*: array[FilterRejectReason, int64]

  FastqPairedFilterStats* = object
    inputPairs*: int64
    outputPairs*: int64
    discardedPairs*: int64
    forwardStats*: FastqFilterStats
    reverseStats*: FastqFilterStats
    rejectedByPairReason*: array[PairRejectReason, int64]

type
  WriterKind = enum
    wkPlain,
    wkStdout,
    wkGzip

  FastqWriter = object
    kind: WriterKind
    file: File
    gz: GzWriteFile

  FilteredRead = object
    keep: bool
    record: FQRecord
    reason: FilterRejectReason

  FilteredPair = object
    keep: bool
    read1: FQRecord
    read2: FQRecord
    r1Reason: FilterRejectReason
    r2Reason: FilterRejectReason
    pairReason: PairRejectReason


proc defaultFastqFilterOptions*(): FastqFilterOptions =
  FastqFilterOptions(
    truncLen: 0,
    truncQ: 0,
    maxLen: 0,
    maxEE: 2.0,
    maxN: 0,
    minLen: 20,
    minQ: 2,
    trimLeft: 0,
    trimRight: 0,
    chunkSize: 10000,
    threads: max(1, countProcessors()),
    parallelMinChunk: 512
  )

proc defaultPairedFastqFilterOptions*(): PairedFastqFilterOptions =
  result.forward = defaultFastqFilterOptions()
  result.reverse = defaultFastqFilterOptions()
  result.checkNames = true

proc validateOptions*(opts: FastqFilterOptions) =
  if opts.truncLen < 0:
    raise newException(ValueError, "truncLen must be >= 0")
  if opts.truncQ < 0:
    raise newException(ValueError, "truncQ must be >= 0")
  if opts.maxLen < 0:
    raise newException(ValueError, "maxLen must be >= 0")
  if opts.maxEE < 0:
    raise newException(ValueError, "maxEE must be >= 0")
  if opts.maxN < 0:
    raise newException(ValueError, "maxN must be >= 0")
  if opts.minLen < 0:
    raise newException(ValueError, "minLen must be >= 0")
  if opts.minQ < 0:
    raise newException(ValueError, "minQ must be >= 0")
  if opts.trimLeft < 0:
    raise newException(ValueError, "trimLeft must be >= 0")
  if opts.trimRight < 0:
    raise newException(ValueError, "trimRight must be >= 0")
  if opts.chunkSize <= 0:
    raise newException(ValueError, "chunkSize must be > 0")
  if opts.threads <= 0:
    raise newException(ValueError, "threads must be > 0")
  if opts.parallelMinChunk <= 0:
    raise newException(ValueError, "parallelMinChunk must be > 0")

proc validateOptions*(opts: PairedFastqFilterOptions) =
  validateOptions(opts.forward)
  validateOptions(opts.reverse)

proc describe*(opts: FastqFilterOptions): string =
  validateOptions(opts)
  &"truncLen={opts.truncLen} truncQ={opts.truncQ} maxLen={opts.maxLen} " &
    &"maxEE={opts.maxEE} maxN={opts.maxN} minLen={opts.minLen} " &
    &"minQ={opts.minQ} trimLeft={opts.trimLeft} trimRight={opts.trimRight} " &
    &"chunkSize={opts.chunkSize} threads={opts.threads} parallelMinChunk={opts.parallelMinChunk}"

proc describe*(opts: PairedFastqFilterOptions): string =
  validateOptions(opts)
  &"checkNames={opts.checkNames} forward=({opts.forward.describe()}) " &
    &"reverse=({opts.reverse.describe()})"

proc discardRatio*(stats: FastqFilterStats): float =
  if stats.inputReads == 0:
    return 0
  stats.discardedReads.float / stats.inputReads.float

proc discardRatio*(stats: FastqPairedFilterStats): float =
  if stats.inputPairs == 0:
    return 0
  stats.discardedPairs.float / stats.inputPairs.float


proc openWriter(path: string): FastqWriter =
  if path == "-":
    result.kind = wkStdout
    result.file = stdout
    return

  if path.toLowerAscii().endsWith(".gz"):
    result.kind = wkGzip
    result.gz = gzopen(path.cstring, "wb")
    if result.gz.isNil:
      raise newException(IOError, "Unable to open gzip output file: " & path)
    return

  result.kind = wkPlain
  if not open(result.file, path, fmWrite):
    raise newException(IOError, "Unable to open output file: " & path)

proc writeData(writer: var FastqWriter; text: string) =
  case writer.kind
  of wkPlain, wkStdout:
    writer.file.write(text)
  of wkGzip:
    if text.len == 0:
      return
    let written = gzwrite(writer.gz, unsafeAddr text[0], cuint(text.len))
    if written < 0:
      raise newException(IOError, "Failed writing gzip output")

proc writeRecords(writer: var FastqWriter; records: openArray[FQRecord]) =
  if records.len == 0:
    return
  var buffer = newStringOfCap(records.len * 96)
  for rec in records:
    buffer.add($rec)
    buffer.add('\n')
  writeData(writer, buffer)

proc closeWriter(writer: var FastqWriter) =
  case writer.kind
  of wkPlain:
    if writer.file != nil:
      close(writer.file)
      writer.file = nil
  of wkStdout:
    discard
  of wkGzip:
    if not writer.gz.isNil:
      discard gzclose(writer.gz)
      writer.gz = nil


proc phred(c: char): int {.inline.} =
  ord(c) - 33

proc countNs(seq: string): int =
  for c in seq:
    if c == 'N' or c == 'n':
      inc result

var eeLookupInit = false
var eeLookup: array[94, float]

proc ensureEeLookup() =
  if eeLookupInit:
    return
  for q in 0..high(eeLookup):
    eeLookup[q] = pow(10.0, -(q.float / 10.0))
  eeLookupInit = true

proc expectedErrors(quality: string): float =
  ensureEeLookup()
  for qchar in quality:
    let q = phred(qchar)
    if q < 0:
      return Inf
    if q <= high(eeLookup):
      result += eeLookup[q]
    else:
      result += pow(10.0, -(q.float / 10.0))

proc minQuality(quality: string): int =
  result = high(int)
  for qchar in quality:
    let q = phred(qchar)
    if q < result:
      result = q

proc applyFilter(rec: FQRecord; opts: FastqFilterOptions): FilteredRead =
  if rec.sequence.len != rec.quality.len:
    return FilteredRead(keep: false, reason: frrMalformedQuality)

  var seq = rec.sequence
  var qual = rec.quality

  if opts.maxLen > 0 and seq.len > opts.maxLen:
    return FilteredRead(keep: false, reason: frrTooLong)

  if opts.trimLeft > 0 or opts.trimRight > 0:
    let endPos = seq.len - opts.trimRight
    let startPos = opts.trimLeft
    if startPos >= endPos:
      return FilteredRead(keep: false, reason: frrTrimmedEmpty)
    seq = seq[startPos ..< endPos]
    qual = qual[startPos ..< endPos]

  if opts.truncQ > 0 and qual.len > 0:
    var cutPos = -1
    for i in 0..<qual.len:
      if phred(qual[i]) <= opts.truncQ:
        cutPos = i
        break
    if cutPos == 0:
      return FilteredRead(keep: false, reason: frrTrimmedEmpty)
    if cutPos > 0:
      seq.setLen(cutPos)
      qual.setLen(cutPos)

  if opts.truncLen > 0:
    if seq.len < opts.truncLen:
      return FilteredRead(keep: false, reason: frrTruncLen)
    if seq.len > opts.truncLen:
      seq.setLen(opts.truncLen)
      qual.setLen(opts.truncLen)

  if seq.len < opts.minLen:
    return FilteredRead(keep: false, reason: frrTooShort)

  if countNs(seq) > opts.maxN:
    return FilteredRead(keep: false, reason: frrTooManyN)

  let mq = minQuality(qual)
  if mq < 0:
    return FilteredRead(keep: false, reason: frrMalformedQuality)
  if mq <= opts.minQ:
    return FilteredRead(keep: false, reason: frrLowMinQ)

  if expectedErrors(qual) > opts.maxEE:
    return FilteredRead(keep: false, reason: frrTooManyExpectedErrors)

  var filtered = rec
  filtered.sequence = seq
  filtered.quality = qual
  FilteredRead(keep: true, record: filtered, reason: frrNone)

proc applyPairFilter(pair: FQPair; opts: PairedFastqFilterOptions): FilteredPair =
  let r1 = applyFilter(pair.read1, opts.forward)
  let r2 = applyFilter(pair.read2, opts.reverse)

  result.r1Reason = r1.reason
  result.r2Reason = r2.reason

  if r1.keep and r2.keep:
    result.keep = true
    result.read1 = r1.record
    result.read2 = r2.record
    return

  result.keep = false
  if not r1.keep and not r2.keep:
    result.pairReason = prrBothFailed
  elif not r1.keep:
    result.pairReason = prrForwardFailed
  else:
    result.pairReason = prrReverseFailed

template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc filterRange(src: ptr UncheckedArray[FQRecord], dst: ptr UncheckedArray[FilteredRead], n: int, opts: FastqFilterOptions) =
  for i in 0..<n:
    dst[i] = applyFilter(src[i], opts)

proc filterPairRange(src: ptr UncheckedArray[FQPair], dst: ptr UncheckedArray[FilteredPair], n: int, opts: PairedFastqFilterOptions) =
  for i in 0..<n:
    dst[i] = applyPairFilter(src[i], opts)

proc filterChunk(records: var seq[FQRecord], opts: FastqFilterOptions): seq[FilteredRead] =
  result = newSeq[FilteredRead](records.len)

  if records.len == 0:
    return

  if opts.threads == 1 or records.len < opts.parallelMinChunk:
    for i, rec in records:
      result[i] = applyFilter(rec, opts)
    return

  let workerCount = min(opts.threads, records.len)
  let blockSize = (records.len + workerCount - 1) div workerCount

  var m = createMaster()
  m.awaitAll:
    var idx = 0
    while idx < records.len:
      let n = min(blockSize, records.len - idx)
      m.spawn filterRange(atPtr(records, idx), atPtr(result, idx), n, opts)
      idx += n

proc filterPairChunk(records: var seq[FQPair], opts: PairedFastqFilterOptions): seq[FilteredPair] =
  result = newSeq[FilteredPair](records.len)

  if records.len == 0:
    return

  let threads = max(opts.forward.threads, opts.reverse.threads)

  if threads == 1 or records.len < max(opts.forward.parallelMinChunk, opts.reverse.parallelMinChunk):
    for i, rec in records:
      result[i] = applyPairFilter(rec, opts)
    return

  let workerCount = min(threads, records.len)
  let blockSize = (records.len + workerCount - 1) div workerCount

  var m = createMaster()
  m.awaitAll:
    var idx = 0
    while idx < records.len:
      let n = min(blockSize, records.len - idx)
      m.spawn filterPairRange(atPtr(records, idx), atPtr(result, idx), n, opts)
      idx += n

proc accumulateSingleStats(stats: var FastqFilterStats; filtered: FilteredRead) =
  inc stats.inputReads
  if filtered.keep:
    inc stats.outputReads
  else:
    inc stats.discardedReads
    inc stats.rejectedByReason[filtered.reason]

proc fastqFilter*(inputPath, outputPath: string, opts: FastqFilterOptions = defaultFastqFilterOptions()): FastqFilterStats =
  validateOptions(opts)

  var writer = openWriter(outputPath)
  defer:
    closeWriter(writer)

  var chunk = newSeqOfCap[FQRecord](opts.chunkSize)

  for rec in readFQ(inputPath):
    chunk.add(rec)
    if chunk.len >= opts.chunkSize:
      var filtered = filterChunk(chunk, opts)
      var passing = newSeqOfCap[FQRecord](filtered.len)
      for item in filtered:
        accumulateSingleStats(result, item)
        if item.keep:
          passing.add(item.record)
      writeRecords(writer, passing)
      chunk.setLen(0)

  if chunk.len > 0:
    var filtered = filterChunk(chunk, opts)
    var passing = newSeqOfCap[FQRecord](filtered.len)
    for item in filtered:
      accumulateSingleStats(result, item)
      if item.keep:
        passing.add(item.record)
    writeRecords(writer, passing)

proc fastqPairedFilter*(
  inputPath1, inputPath2: string,
  outputPath1, outputPath2: string,
  opts: PairedFastqFilterOptions = defaultPairedFastqFilterOptions()
): FastqPairedFilterStats =
  validateOptions(opts)

  var writer1 = openWriter(outputPath1)
  var writer2 = openWriter(outputPath2)
  defer:
    closeWriter(writer1)
    closeWriter(writer2)

  let chunkSize = max(opts.forward.chunkSize, opts.reverse.chunkSize)
  var chunk = newSeqOfCap[FQPair](chunkSize)

  for pair in readFQPair(inputPath1, inputPath2, checkNames = opts.checkNames):
    chunk.add(pair)
    if chunk.len >= chunkSize:
      var filtered = filterPairChunk(chunk, opts)
      var pass1 = newSeqOfCap[FQRecord](filtered.len)
      var pass2 = newSeqOfCap[FQRecord](filtered.len)

      for item in filtered:
        inc result.inputPairs

        inc result.forwardStats.inputReads
        if item.r1Reason == frrNone:
          inc result.forwardStats.outputReads
        else:
          inc result.forwardStats.discardedReads
          inc result.forwardStats.rejectedByReason[item.r1Reason]

        inc result.reverseStats.inputReads
        if item.r2Reason == frrNone:
          inc result.reverseStats.outputReads
        else:
          inc result.reverseStats.discardedReads
          inc result.reverseStats.rejectedByReason[item.r2Reason]

        if item.keep:
          inc result.outputPairs
          pass1.add(item.read1)
          pass2.add(item.read2)
        else:
          inc result.discardedPairs
          inc result.rejectedByPairReason[item.pairReason]

      writeRecords(writer1, pass1)
      writeRecords(writer2, pass2)
      chunk.setLen(0)

  if chunk.len > 0:
    var filtered = filterPairChunk(chunk, opts)
    var pass1 = newSeqOfCap[FQRecord](filtered.len)
    var pass2 = newSeqOfCap[FQRecord](filtered.len)

    for item in filtered:
      inc result.inputPairs

      inc result.forwardStats.inputReads
      if item.r1Reason == frrNone:
        inc result.forwardStats.outputReads
      else:
        inc result.forwardStats.discardedReads
        inc result.forwardStats.rejectedByReason[item.r1Reason]

      inc result.reverseStats.inputReads
      if item.r2Reason == frrNone:
        inc result.reverseStats.outputReads
      else:
        inc result.reverseStats.discardedReads
        inc result.reverseStats.rejectedByReason[item.r2Reason]

      if item.keep:
        inc result.outputPairs
        pass1.add(item.read1)
        pass2.add(item.read2)
      else:
        inc result.discardedPairs
        inc result.rejectedByPairReason[item.pairReason]

    writeRecords(writer1, pass1)
    writeRecords(writer2, pass2)
