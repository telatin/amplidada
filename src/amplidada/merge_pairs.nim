import std/[algorithm, cpuinfo, strformat, tables]
import malebolgia
import amplidada/[dada, derep]

type
  MergeRejectReason* = enum
    mrrNoForwardAsv,
    mrrNoReverseAsv,
    mrrUnsynchronizedReads,
    mrrNoAlignment,
    mrrTooManyMismatches,
    mrrOverhangDisallowed

  MergePairsOptions* = object
    minOverlap*: int
    maxMismatch*: int
    trimOverhang*: bool
    threads*: int
    parallelMinPairs*: int

  MergedSequence* = object
    sequence*: string
    abundance*: int64

  MergePairsStats* = object
    forwardReadCount*: int64
    reverseReadCount*: int64
    synchronizedPairs*: int64
    unsynchronizedReads*: int64
    candidatePairs*: int64
    mergedPairs*: int64
    rejectedPairs*: int64
    uniqueAsvPairs*: int
    uniqueMergedSequences*: int
    rejectedByReason*: array[MergeRejectReason, int64]

  MergePairsResult* = object
    stats*: MergePairsStats
    merged*: seq[MergedSequence]

type
  PairCombo = object
    forwardAsv: int
    reverseAsv: int
    abundance: int64

  PartialMerge = object
    mergedCounts: Table[string, int64]
    mergedPairs: int64
    rejectedPairs: int64
    rejectedByReason: array[MergeRejectReason, int64]

proc defaultMergePairsOptions*(): MergePairsOptions =
  MergePairsOptions(
    minOverlap: 12,
    maxMismatch: 0,
    trimOverhang: false,
    threads: max(1, countProcessors()),
    parallelMinPairs: 128
  )

proc validateOptions*(opts: MergePairsOptions) =
  if opts.minOverlap <= 0:
    raise newException(ValueError, "minOverlap must be > 0")
  if opts.maxMismatch < 0:
    raise newException(ValueError, "maxMismatch must be >= 0")
  if opts.threads <= 0:
    raise newException(ValueError, "threads must be > 0")
  if opts.parallelMinPairs <= 0:
    raise newException(ValueError, "parallelMinPairs must be > 0")

proc describe*(opts: MergePairsOptions): string =
  validateOptions(opts)
  &"minOverlap={opts.minOverlap} maxMismatch={opts.maxMismatch} trimOverhang={opts.trimOverhang} " &
    &"threads={opts.threads} parallelMinPairs={opts.parallelMinPairs}"

proc reverseComplement*(sequence: string): string =
  result = newString(sequence.len)
  var j = 0
  for i in countdown(sequence.len - 1, 0):
    result[j] =
      case sequence[i]
      of 'A', 'a': 'T'
      of 'C', 'c': 'G'
      of 'G', 'g': 'C'
      of 'T', 't': 'A'
      of 'N', 'n': 'N'
      else: 'N'
    inc j

proc overlapLength(lenForward, lenReverse, offset: int): int =
  let startForward = max(0, offset)
  let startReverse = max(0, -offset)
  min(lenForward - startForward, lenReverse - startReverse)

proc tryMergePair(
  forwardSeq: string,
  reverseSeqRc: string,
  opts: MergePairsOptions,
  mergedOut: var string,
  reason: var MergeRejectReason
): bool =
  let lenF = forwardSeq.len
  let lenR = reverseSeqRc.len
  if lenF == 0 or lenR == 0:
    reason = mrrNoAlignment
    return false

  var found = false
  var bestOffset = 0
  var bestOverlap = -1
  var bestMismatch = high(int)
  var sawTooManyMismatch = false
  var sawOverhang = false
  var sawOverlap = false

  let minOffset = -lenR + opts.minOverlap
  let maxOffset = lenF - opts.minOverlap

  for offset in minOffset..maxOffset:
    let ovLen = overlapLength(lenF, lenR, offset)
    if ovLen < opts.minOverlap:
      continue
    sawOverlap = true

    let startF = max(0, offset)
    let startR = max(0, -offset)

    var mismatches = 0
    for i in 0..<ovLen:
      if forwardSeq[startF + i] != reverseSeqRc[startR + i]:
        inc mismatches
        if mismatches > opts.maxMismatch:
          break

    if mismatches > opts.maxMismatch:
      sawTooManyMismatch = true
      continue

    let overhang = (offset < 0) or (offset + lenR > lenF)
    if overhang and not opts.trimOverhang:
      sawOverhang = true
      continue

    if (not found) or (ovLen > bestOverlap) or
       (ovLen == bestOverlap and mismatches < bestMismatch) or
       (ovLen == bestOverlap and mismatches == bestMismatch and offset < bestOffset):
      found = true
      bestOffset = offset
      bestOverlap = ovLen
      bestMismatch = mismatches

  if not found:
    if sawOverhang:
      reason = mrrOverhangDisallowed
    elif sawTooManyMismatch:
      reason = mrrTooManyMismatches
    elif sawOverlap:
      reason = mrrNoAlignment
    else:
      reason = mrrNoAlignment
    return false

  let mergedStart = min(0, bestOffset)
  let mergedEnd = max(lenF, bestOffset + lenR)
  var merged = newStringOfCap(mergedEnd - mergedStart)
  for pos in mergedStart..<mergedEnd:
    let hasF = pos >= 0 and pos < lenF
    let rPos = pos - bestOffset
    let hasR = rPos >= 0 and rPos < lenR

    if hasF and hasR:
      merged.add(forwardSeq[pos])
    elif hasF:
      merged.add(forwardSeq[pos])
    else:
      merged.add(reverseSeqRc[rPos])

  mergedOut = merged
  reason = mrrNoAlignment
  true

proc addCount(target: var Table[string, int64], sequence: string, n: int64) =
  if sequence in target:
    target[sequence] = target[sequence] + n
  else:
    target[sequence] = n

proc processComboSlice(
  combos: ptr UncheckedArray[PairCombo],
  n: int,
  forwardAsvSeq: seq[string],
  reverseAsvSeqRc: seq[string],
  opts: MergePairsOptions
): PartialMerge =
  result.mergedCounts = initTable[string, int64]()

  for i in 0..<n:
    let combo = combos[i]
    var merged = ""
    var reason = mrrNoAlignment
    if tryMergePair(
      forwardSeq = forwardAsvSeq[combo.forwardAsv],
      reverseSeqRc = reverseAsvSeqRc[combo.reverseAsv],
      opts = opts,
      mergedOut = merged,
      reason = reason
    ):
      addCount(result.mergedCounts, merged, combo.abundance)
      result.mergedPairs += combo.abundance
    else:
      result.rejectedPairs += combo.abundance
      result.rejectedByReason[reason] += combo.abundance

template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc mergePairs*(
  derepForward: DerepFastqResult,
  derepReverse: DerepFastqResult,
  dadaForward: DadaResult,
  dadaReverse: DadaResult,
  opts: MergePairsOptions = defaultMergePairsOptions()
): MergePairsResult =
  validateOptions(opts)

  if derepForward.readToUnique.len == 0 or derepReverse.readToUnique.len == 0:
    raise newException(ValueError, "mergePairs requires derepFastq readToUnique maps for both directions")

  result.stats.forwardReadCount = derepForward.readToUnique.len.int64
  result.stats.reverseReadCount = derepReverse.readToUnique.len.int64
  result.stats.synchronizedPairs = min(
    derepForward.readToUnique.len,
    derepReverse.readToUnique.len
  ).int64
  result.stats.unsynchronizedReads = abs(derepForward.readToUnique.len - derepReverse.readToUnique.len).int64
  if result.stats.unsynchronizedReads > 0:
    result.stats.rejectedByReason[mrrUnsynchronizedReads] = result.stats.unsynchronizedReads

  var pairCounts = initTable[uint64, int64]()

  for i in 0..<result.stats.synchronizedPairs.int:
    let fUniq = derepForward.readToUnique[i]
    let rUniq = derepReverse.readToUnique[i]

    if fUniq < 0 or fUniq >= dadaForward.uniqueToAsv.len:
      result.stats.rejectedPairs += 1
      result.stats.rejectedByReason[mrrNoForwardAsv] += 1
      continue
    if rUniq < 0 or rUniq >= dadaReverse.uniqueToAsv.len:
      result.stats.rejectedPairs += 1
      result.stats.rejectedByReason[mrrNoReverseAsv] += 1
      continue

    let fAsv = dadaForward.uniqueToAsv[fUniq]
    let rAsv = dadaReverse.uniqueToAsv[rUniq]
    if fAsv < 0 or fAsv >= dadaForward.asvs.len:
      result.stats.rejectedPairs += 1
      result.stats.rejectedByReason[mrrNoForwardAsv] += 1
      continue
    if rAsv < 0 or rAsv >= dadaReverse.asvs.len:
      result.stats.rejectedPairs += 1
      result.stats.rejectedByReason[mrrNoReverseAsv] += 1
      continue

    let key = (uint64(fAsv.uint32) shl 32) or uint64(rAsv.uint32)
    if key in pairCounts:
      pairCounts[key] = pairCounts[key] + 1
    else:
      pairCounts[key] = 1
    result.stats.candidatePairs += 1

  var combos: seq[PairCombo] = @[]
  combos.setLen(pairCounts.len)
  var comboIdx = 0
  for key, abundance in pairCounts:
    let fAsv = int((key shr 32) and 0xffffffff'u64)
    let rAsv = int(key and 0xffffffff'u64)
    combos[comboIdx] = PairCombo(forwardAsv: fAsv, reverseAsv: rAsv, abundance: abundance)
    inc comboIdx
  result.stats.uniqueAsvPairs = combos.len

  var forwardAsvSeq = newSeq[string](dadaForward.asvs.len)
  for i, asv in dadaForward.asvs:
    forwardAsvSeq[i] = asv.sequence

  var reverseAsvSeqRc = newSeq[string](dadaReverse.asvs.len)
  for i, asv in dadaReverse.asvs:
    reverseAsvSeqRc[i] = reverseComplement(asv.sequence)

  var mergedCounts = initTable[string, int64]()

  if combos.len > 0:
    if opts.threads == 1 or combos.len < opts.parallelMinPairs:
      let partial = processComboSlice(atPtr(combos, 0), combos.len, forwardAsvSeq, reverseAsvSeqRc, opts)
      result.stats.mergedPairs += partial.mergedPairs
      result.stats.rejectedPairs += partial.rejectedPairs
      for reason in MergeRejectReason:
        result.stats.rejectedByReason[reason] += partial.rejectedByReason[reason]
      for seq, abundance in partial.mergedCounts:
        addCount(mergedCounts, seq, abundance)
    else:
      let workerCount = min(opts.threads, combos.len)
      let blockSize = (combos.len + workerCount - 1) div workerCount
      var partials = newSeq[PartialMerge](workerCount)

      var m = createMaster()
      m.awaitAll:
        var idx = 0
        var worker = 0
        while idx < combos.len:
          let n = min(blockSize, combos.len - idx)
          m.spawn processComboSlice(
            combos = atPtr(combos, idx),
            n = n,
            forwardAsvSeq = forwardAsvSeq,
            reverseAsvSeqRc = reverseAsvSeqRc,
            opts = opts
          ) -> partials[worker]
          idx += n
          inc worker

      for partial in partials:
        result.stats.mergedPairs += partial.mergedPairs
        result.stats.rejectedPairs += partial.rejectedPairs
        for reason in MergeRejectReason:
          result.stats.rejectedByReason[reason] += partial.rejectedByReason[reason]
        for seq, abundance in partial.mergedCounts:
          addCount(mergedCounts, seq, abundance)

  result.merged = newSeq[MergedSequence](mergedCounts.len)
  var idx = 0
  for seq, abundance in mergedCounts:
    result.merged[idx] = MergedSequence(sequence: seq, abundance: abundance)
    inc idx
  result.merged.sort(proc(a, b: MergedSequence): int =
    if a.abundance != b.abundance:
      return cmp(b.abundance, a.abundance)
    cmp(a.sequence, b.sequence)
  )

  result.stats.uniqueMergedSequences = result.merged.len

proc writeMergedTsv*(path: string, res: MergePairsResult) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open merged output file: " & path)
  defer: close(f)

  f.writeLine("merged_id\tsequence\tabundance")
  for i, row in res.merged:
    f.writeLine(&"MERGED{i + 1}\t{row.sequence}\t{row.abundance}")
