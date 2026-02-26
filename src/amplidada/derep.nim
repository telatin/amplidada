import std/[algorithm, cpuinfo, strformat, tables]
import readfx
import malebolgia

type
  DerepOrder* = enum
    dsoAbundance,
    dsoFirstSeen

  DerepFastqOptions* = object
    chunkSize*: int
    threads*: int
    parallelMinChunk*: int
    keepReadMap*: bool
    order*: DerepOrder

  DerepUnique* = object
    sequence*: string
    abundance*: int64
    meanQual*: seq[float64]
    firstSeenReadIndex*: int64

  DerepFastqResult* = object
    totalReads*: int64
    uniqueCount*: int
    maxReadLength*: int
    uniques*: seq[DerepUnique]
    readToUnique*: seq[int]


type
  PartialDerep = object
    startIndex: int64
    readCount: int
    uniqSeqs: seq[string]
    counts: seq[int64]
    qualSums: seq[seq[float64]]
    firstSeen: seq[int64]
    readLocalUnique: seq[int]


template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc defaultDerepFastqOptions*(): DerepFastqOptions =
  DerepFastqOptions(
    chunkSize: 20000,
    threads: max(1, countProcessors()),
    parallelMinChunk: 512,
    keepReadMap: true,
    order: dsoAbundance
  )

proc validateOptions*(opts: DerepFastqOptions) =
  if opts.chunkSize <= 0:
    raise newException(ValueError, "chunkSize must be > 0")
  if opts.threads <= 0:
    raise newException(ValueError, "threads must be > 0")
  if opts.parallelMinChunk <= 0:
    raise newException(ValueError, "parallelMinChunk must be > 0")

proc describe*(opts: DerepFastqOptions): string =
  validateOptions(opts)
  &"chunkSize={opts.chunkSize} threads={opts.threads} parallelMinChunk={opts.parallelMinChunk} " &
    &"keepReadMap={opts.keepReadMap} order={opts.order}"

proc processDerepSlice(
  src: ptr UncheckedArray[FQRecord],
  n: int,
  globalStart: int64,
  keepReadMap: bool
): PartialDerep =
  result.startIndex = globalStart
  result.readCount = n
  if keepReadMap:
    result.readLocalUnique = newSeq[int](n)

  var seqToIdx = initTable[string, int]()

  for i in 0..<n:
    let rec = src[i]

    if rec.sequence.len != rec.quality.len:
      raise newException(ValueError, "Malformed FASTQ: sequence and quality lengths differ for read '" & rec.name & "'")

    var localIdx: int
    if seqToIdx.hasKey(rec.sequence):
      localIdx = seqToIdx[rec.sequence]
    else:
      localIdx = result.uniqSeqs.len
      seqToIdx[rec.sequence] = localIdx
      result.uniqSeqs.add(rec.sequence)
      result.counts.add(0)
      result.firstSeen.add(globalStart + i.int64)
      result.qualSums.add(newSeq[float64](rec.sequence.len))

    inc result.counts[localIdx]

    for j in 0..<rec.quality.len:
      let q = ord(rec.quality[j]) - 33
      if q < 0:
        raise newException(ValueError, "Malformed FASTQ: invalid quality character in read '" & rec.name & "'")
      result.qualSums[localIdx][j] += q.float64

    if keepReadMap:
      result.readLocalUnique[i] = localIdx

proc processDerepChunk(
  records: var seq[FQRecord],
  chunkStart: int64,
  opts: DerepFastqOptions
): seq[PartialDerep] =
  if records.len == 0:
    return @[]

  if opts.threads == 1 or records.len < opts.parallelMinChunk:
    return @[processDerepSlice(atPtr(records, 0), records.len, chunkStart, opts.keepReadMap)]

  let workerCount = min(opts.threads, records.len)
  let blockSize = (records.len + workerCount - 1) div workerCount

  result = newSeq[PartialDerep](workerCount)

  var m = createMaster()
  m.awaitAll:
    var worker = 0
    var idx = 0
    while idx < records.len:
      let n = min(blockSize, records.len - idx)
      let globalStart = chunkStart + idx.int64
      m.spawn processDerepSlice(atPtr(records, idx), n, globalStart, opts.keepReadMap) -> result[worker]
      idx += n
      inc worker

proc finalizeResult(
  seqs: seq[string],
  counts: seq[int64],
  firstSeen: seq[int64],
  qualSums: seq[seq[float64]],
  readToUniqueTmp: var seq[int],
  order: DerepOrder,
  totalReads: int64,
  maxReadLength: int,
  keepReadMap: bool
): DerepFastqResult =
  result.totalReads = totalReads
  result.maxReadLength = maxReadLength
  result.uniqueCount = seqs.len
  result.readToUnique = readToUniqueTmp

  var ordering = newSeq[int](seqs.len)
  for i in 0..<seqs.len:
    ordering[i] = i

  ordering.sort(proc(a, b: int): int =
    case order
    of dsoAbundance:
      if counts[a] != counts[b]:
        return cmp(counts[b], counts[a])
      if firstSeen[a] != firstSeen[b]:
        return cmp(firstSeen[a], firstSeen[b])
      cmp(seqs[a], seqs[b])
    of dsoFirstSeen:
      if firstSeen[a] != firstSeen[b]:
        return cmp(firstSeen[a], firstSeen[b])
      if counts[a] != counts[b]:
        return cmp(counts[b], counts[a])
      cmp(seqs[a], seqs[b])
  )

  var oldToNew = newSeq[int](seqs.len)
  result.uniques = newSeq[DerepUnique](seqs.len)

  for newIdx, oldIdx in ordering:
    oldToNew[oldIdx] = newIdx

    var meanQual = newSeq[float64](qualSums[oldIdx].len)
    for i in 0..<meanQual.len:
      meanQual[i] = qualSums[oldIdx][i] / counts[oldIdx].float64

    result.uniques[newIdx] = DerepUnique(
      sequence: seqs[oldIdx],
      abundance: counts[oldIdx],
      meanQual: meanQual,
      firstSeenReadIndex: firstSeen[oldIdx]
    )

  if keepReadMap:
    for i in 0..<result.readToUnique.len:
      result.readToUnique[i] = oldToNew[result.readToUnique[i]]
  else:
    result.readToUnique.setLen(0)

proc derepFastq*(inputPath: string, opts: DerepFastqOptions = defaultDerepFastqOptions()): DerepFastqResult =
  validateOptions(opts)

  var globalSeqToIdx = initTable[string, int]()
  var globalSeqs: seq[string] = @[]
  var globalCounts: seq[int64] = @[]
  var globalFirstSeen: seq[int64] = @[]
  var globalQualSums: seq[seq[float64]] = @[]

  var readToUniqueTmp: seq[int] = @[]

  var chunk = newSeqOfCap[FQRecord](opts.chunkSize)
  var totalReads = 0.int64
  var maxReadLength = 0

  proc mergeChunk(records: var seq[FQRecord]) =
    if records.len == 0:
      return

    let chunkStart = totalReads
    if opts.keepReadMap:
      readToUniqueTmp.setLen(readToUniqueTmp.len + records.len)

    let partials = processDerepChunk(records, chunkStart, opts)

    for part in partials:
      for i in 0..<part.uniqSeqs.len:
        let seq = part.uniqSeqs[i]
        var globalIdx: int

        if globalSeqToIdx.hasKey(seq):
          globalIdx = globalSeqToIdx[seq]
        else:
          globalIdx = globalSeqs.len
          globalSeqToIdx[seq] = globalIdx
          globalSeqs.add(seq)
          globalCounts.add(0)
          globalFirstSeen.add(part.firstSeen[i])
          globalQualSums.add(newSeq[float64](part.qualSums[i].len))

        globalCounts[globalIdx] += part.counts[i]
        if part.firstSeen[i] < globalFirstSeen[globalIdx]:
          globalFirstSeen[globalIdx] = part.firstSeen[i]

        for q in 0..<part.qualSums[i].len:
          globalQualSums[globalIdx][q] += part.qualSums[i][q]

      if opts.keepReadMap:
        for i in 0..<part.readCount:
          let seq = part.uniqSeqs[part.readLocalUnique[i]]
          let globalIdx = globalSeqToIdx[seq]
          let globalReadIndex = (part.startIndex - chunkStart).int + i
          readToUniqueTmp[chunkStart.int + globalReadIndex] = globalIdx

    totalReads += records.len.int64

  for rec in readFQ(inputPath):
    if rec.sequence.len > maxReadLength:
      maxReadLength = rec.sequence.len
    chunk.add(rec)
    if chunk.len >= opts.chunkSize:
      mergeChunk(chunk)
      chunk.setLen(0)

  if chunk.len > 0:
    mergeChunk(chunk)

  finalizeResult(
    seqs = globalSeqs,
    counts = globalCounts,
    firstSeen = globalFirstSeen,
    qualSums = globalQualSums,
    readToUniqueTmp = readToUniqueTmp,
    order = opts.order,
    totalReads = totalReads,
    maxReadLength = maxReadLength,
    keepReadMap = opts.keepReadMap
  )
