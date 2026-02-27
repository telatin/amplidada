import std/[algorithm, cpuinfo, sequtils, strformat, strutils, tables]
import malebolgia

type
  BimeraCall* = enum
    bcNone,
    bcExact,
    bcOneOffMismatch,
    bcOneOffIndel

  BimeraMethod* = enum
    bmPooled,
    bmConsensus,
    bmPerSample

  BimeraOptions* = object
    minFoldParentOverAbundance*: float64
    minParentAbundance*: int64
    allowOneOff*: bool
    minOneOffParentDistance*: int
    threads*: int
    parallelMinCandidates*: int

  BimeraConsensusOptions* = object
    minSampleFraction*: float64
    ignoreNNegatives*: int

  UniqueSequenceAbundance* = object
    sequence*: string
    abundance*: int64

  SampleSequenceAbundance* = object
    sampleId*: string
    sequence*: string
    abundance*: int64

  BimeraDenovoResult* = object
    sequences*: seq[string]
    abundances*: seq[int64]
    flags*: seq[bool]
    calls*: seq[BimeraCall]
    nBimeras*: int
    exactBimeraCount*: int
    oneOffMismatchBimeraCount*: int
    oneOffIndelBimeraCount*: int

  BimeraTableResult* = object
    sequences*: seq[string]
    presentSamples*: seq[int]
    flaggedSamples*: seq[int]
    flags*: seq[bool]
    calls*: seq[BimeraCall]
    nBimeras*: int
    exactBimeraCount*: int
    oneOffMismatchBimeraCount*: int
    oneOffIndelBimeraCount*: int

  RemoveBimeraResult* = object
    mode*: BimeraMethod
    inputUniqueCount*: int
    removedUniqueCount*: int
    keptUniqueCount*: int
    pooled*: BimeraDenovoResult
    consensus*: BimeraTableResult
    removed*: seq[SampleSequenceAbundance]
    kept*: seq[SampleSequenceAbundance]
    exactBimeraCount*: int
    oneOffMismatchBimeraCount*: int
    oneOffIndelBimeraCount*: int

type
  RankedUnique = object
    sequence: string
    abundance: int64

  RankedSampleRow = object
    sampleId: string
    sequence: string
    abundance: int64

  ScoreTop2 = object
    bestIdx: int
    bestScore: int
    secondIdx: int
    secondScore: int

proc defaultBimeraOptions*(): BimeraOptions =
  BimeraOptions(
    minFoldParentOverAbundance: 2.0,
    minParentAbundance: 8,
    allowOneOff: false,
    minOneOffParentDistance: 4,
    threads: max(1, countProcessors()),
    parallelMinCandidates: 128
  )

proc defaultBimeraConsensusOptions*(): BimeraConsensusOptions =
  BimeraConsensusOptions(
    minSampleFraction: 0.9,
    ignoreNNegatives: 1
  )

proc validateOptions*(opts: BimeraOptions) =
  if opts.minFoldParentOverAbundance <= 0:
    raise newException(ValueError, "minFoldParentOverAbundance must be > 0")
  if opts.minParentAbundance < 0:
    raise newException(ValueError, "minParentAbundance must be >= 0")
  if opts.minOneOffParentDistance < 0:
    raise newException(ValueError, "minOneOffParentDistance must be >= 0")
  if opts.threads <= 0:
    raise newException(ValueError, "threads must be > 0")
  if opts.parallelMinCandidates <= 0:
    raise newException(ValueError, "parallelMinCandidates must be > 0")

proc validateOptions*(opts: BimeraConsensusOptions) =
  if opts.minSampleFraction <= 0.0 or opts.minSampleFraction > 1.0:
    raise newException(ValueError, "minSampleFraction must be in (0, 1]")
  if opts.ignoreNNegatives < 0:
    raise newException(ValueError, "ignoreNNegatives must be >= 0")

proc describe*(opts: BimeraOptions): string =
  validateOptions(opts)
  &"minFoldParentOverAbundance={opts.minFoldParentOverAbundance} " &
    &"minParentAbundance={opts.minParentAbundance} allowOneOff={opts.allowOneOff} " &
    &"minOneOffParentDistance={opts.minOneOffParentDistance} threads={opts.threads} " &
    &"parallelMinCandidates={opts.parallelMinCandidates}"

proc describe*(opts: BimeraConsensusOptions): string =
  validateOptions(opts)
  &"minSampleFraction={opts.minSampleFraction} ignoreNNegatives={opts.ignoreNNegatives}"

proc toDnaUpper(sequence: string): string =
  result = newString(sequence.len)
  for i, c in sequence:
    result[i] =
      case c
      of 'a'..'z': chr(ord(c) - 32)
      else: c

proc normalizeUniques(input: openArray[UniqueSequenceAbundance]): seq[RankedUnique] =
  var bySeq = initTable[string, int64]()
  for row in input:
    if row.sequence.len == 0 or row.abundance <= 0:
      continue
    let seqUpper = toDnaUpper(row.sequence)
    bySeq[seqUpper] = bySeq.getOrDefault(seqUpper) + row.abundance

  result = newSeq[RankedUnique](bySeq.len)
  var i = 0
  for sequence, abundance in bySeq:
    result[i] = RankedUnique(sequence: sequence, abundance: abundance)
    inc i

  result.sort(proc(a, b: RankedUnique): int =
    if a.abundance != b.abundance:
      return cmp(b.abundance, a.abundance)
    cmp(a.sequence, b.sequence)
  )

proc normalizeSampleRows(input: openArray[SampleSequenceAbundance]): seq[RankedSampleRow] =
  var byKey = initTable[string, int64]()
  for row in input:
    if row.sampleId.len == 0 or row.sequence.len == 0 or row.abundance <= 0:
      continue
    let sampleId = row.sampleId.strip()
    if sampleId.len == 0:
      continue
    let seqUpper = toDnaUpper(row.sequence)
    let key = sampleId & "\t" & seqUpper
    byKey[key] = byKey.getOrDefault(key) + row.abundance

  result = newSeq[RankedSampleRow](byKey.len)
  var i = 0
  for key, abundance in byKey:
    let p = key.find('\t')
    result[i] = RankedSampleRow(
      sampleId: key[0 ..< p],
      sequence: key[p + 1 .. ^1],
      abundance: abundance
    )
    inc i

  result.sort(proc(a, b: RankedSampleRow): int =
    if a.sampleId != b.sampleId:
      return cmp(a.sampleId, b.sampleId)
    if a.abundance != b.abundance:
      return cmp(b.abundance, a.abundance)
    cmp(a.sequence, b.sequence)
  )

proc commonPrefixLen(a, b: string): int {.inline.} =
  let n = min(a.len, b.len)
  while result < n and a[result] == b[result]:
    inc result

proc commonSuffixLen(a, b: string): int {.inline.} =
  let n = min(a.len, b.len)
  while result < n and a[a.len - 1 - result] == b[b.len - 1 - result]:
    inc result

proc hammingDistance(a, b: string): int =
  if a.len != b.len:
    return max(a.len, b.len)
  for i in 0..<a.len:
    if a[i] != b[i]:
      inc result

proc editDistanceAtMost(a, b: string, maxDist: int): int =
  if maxDist < 0:
    return maxDist + 1

  let la = a.len
  let lb = b.len
  if abs(la - lb) > maxDist:
    return maxDist + 1

  if la == 0:
    return if lb <= maxDist: lb else: maxDist + 1
  if lb == 0:
    return if la <= maxDist: la else: maxDist + 1

  var prev = newSeq[int](lb + 1)
  var curr = newSeq[int](lb + 1)
  for j in 0..lb:
    prev[j] = j

  for i in 1..la:
    curr[0] = i
    var rowMin = curr[0]
    for j in 1..lb:
      let subCost = if a[i - 1] == b[j - 1]: 0 else: 1
      var best = prev[j] + 1
      if curr[j - 1] + 1 < best:
        best = curr[j - 1] + 1
      if prev[j - 1] + subCost < best:
        best = prev[j - 1] + subCost
      curr[j] = best
      if best < rowMin:
        rowMin = best

    if rowMin > maxDist:
      return maxDist + 1
    swap(prev, curr)

  if prev[lb] <= maxDist:
    prev[lb]
  else:
    maxDist + 1

proc anyDistinctPair(aBest, aSecond, bBest, bSecond: int): bool {.inline.} =
  if aBest >= 0 and bBest >= 0 and aBest != bBest:
    return true
  if aBest >= 0 and bSecond >= 0 and aBest != bSecond:
    return true
  if aSecond >= 0 and bBest >= 0 and aSecond != bBest:
    return true
  if aSecond >= 0 and bSecond >= 0 and aSecond != bSecond:
    return true
  false

proc updateTop2Idx(bestIdx: var int; secondIdx: var int; idx: int) {.inline.} =
  if bestIdx < 0:
    bestIdx = idx
  elif bestIdx != idx and secondIdx < 0:
    secondIdx = idx

proc initScoreTop2(): ScoreTop2 {.inline.} =
  ScoreTop2(bestIdx: -1, bestScore: high(int), secondIdx: -1, secondScore: high(int))

proc updateTop2Score(slot: var ScoreTop2; idx: int; score: int) {.inline.} =
  if slot.bestIdx < 0 or score < slot.bestScore:
    if slot.bestIdx != idx:
      slot.secondIdx = slot.bestIdx
      slot.secondScore = slot.bestScore
    slot.bestIdx = idx
    slot.bestScore = score
  elif slot.bestIdx != idx and (slot.secondIdx < 0 or score < slot.secondScore):
    slot.secondIdx = idx
    slot.secondScore = score

proc exactBimeraFromParents(query: string; parentSeqs: openArray[string]): bool =
  let seqLen = query.len
  if seqLen < 2 or parentSeqs.len < 2:
    return false

  var leftBest = newSeq[int](seqLen + 1)
  var leftSecond = newSeq[int](seqLen + 1)
  var rightBest = newSeq[int](seqLen + 1)
  var rightSecond = newSeq[int](seqLen + 1)
  for i in 0..seqLen:
    leftBest[i] = -1
    leftSecond[i] = -1
    rightBest[i] = -1
    rightSecond[i] = -1

  for parentIdx, parentSeq in parentSeqs:
    if parentSeq.len != seqLen:
      continue

    let leftLen = commonPrefixLen(query, parentSeq)
    let rightLen = commonSuffixLen(query, parentSeq)

    for k in 1..leftLen:
      updateTop2Idx(leftBest[k], leftSecond[k], parentIdx)
    for s in 1..rightLen:
      updateTop2Idx(rightBest[s], rightSecond[s], parentIdx)

  for k in 1..<(seqLen):
    let s = seqLen - k
    if anyDistinctPair(leftBest[k], leftSecond[k], rightBest[s], rightSecond[s]):
      return true
  false

proc oneOffMismatchBimeraFromParents(
  query: string;
  parentSeqs: openArray[string];
  minOneOffParentDistance: int
): bool =
  let seqLen = query.len
  if seqLen < 2 or parentSeqs.len < 2:
    return false

  var leftScores = newSeq[ScoreTop2](seqLen + 1)
  var rightScores = newSeq[ScoreTop2](seqLen + 1)
  for i in 0..seqLen:
    leftScores[i] = initScoreTop2()
    rightScores[i] = initScoreTop2()

  var eligibleEqualLen: seq[string] = @[]

  let maxParentDistanceForReject = max(0, minOneOffParentDistance - 1)
  for parentSeq in parentSeqs:
    if parentSeq.len != seqLen:
      continue

    var nearParent = false
    if minOneOffParentDistance > 0:
      nearParent = editDistanceAtMost(query, parentSeq, maxParentDistanceForReject) <= maxParentDistanceForReject
    if nearParent:
      continue

    eligibleEqualLen.add(parentSeq)

  if eligibleEqualLen.len < 2:
    return false

  for parentIdx, parentSeq in eligibleEqualLen:
    var prefixMismatches = 0
    for k in 1..<(seqLen):
      if query[k - 1] != parentSeq[k - 1]:
        inc prefixMismatches
      updateTop2Score(leftScores[k], parentIdx, prefixMismatches)

    var suffixMismatches = 0
    for s in 1..<(seqLen):
      if query[seqLen - s] != parentSeq[seqLen - s]:
        inc suffixMismatches
      updateTop2Score(rightScores[s], parentIdx, suffixMismatches)

  for k in 1..<(seqLen):
    let s = seqLen - k
    let left = leftScores[k]
    let right = rightScores[s]

    if left.bestIdx >= 0 and right.bestIdx >= 0 and left.bestIdx != right.bestIdx and
       left.bestScore + right.bestScore <= 1:
      return true
    if left.bestIdx >= 0 and right.secondIdx >= 0 and left.bestIdx != right.secondIdx and
       left.bestScore + right.secondScore <= 1:
      return true
    if left.secondIdx >= 0 and right.bestIdx >= 0 and left.secondIdx != right.bestIdx and
       left.secondScore + right.bestScore <= 1:
      return true
    if left.secondIdx >= 0 and right.secondIdx >= 0 and left.secondIdx != right.secondIdx and
       left.secondScore + right.secondScore <= 1:
      return true
  false

proc oneOffIndelBimeraFromParents(
  query: string;
  parentSeqs: openArray[string];
  minOneOffParentDistance: int
): bool =
  let seqLen = query.len
  if seqLen < 2 or parentSeqs.len < 2:
    return false

  var eligibleLenMinus1: seq[string] = @[]
  var eligibleLenPlus1: seq[string] = @[]

  let maxParentDistanceForReject = max(0, minOneOffParentDistance - 1)
  for parentSeq in parentSeqs:
    let lenDelta = parentSeq.len - seqLen
    if lenDelta < -1 or lenDelta > 1 or lenDelta == 0:
      continue

    var nearParent = false
    if minOneOffParentDistance > 0:
      nearParent = editDistanceAtMost(query, parentSeq, maxParentDistanceForReject) <= maxParentDistanceForReject
    if nearParent:
      continue

    if lenDelta == -1:
      eligibleLenMinus1.add(parentSeq)
    else:
      eligibleLenPlus1.add(parentSeq)

  proc hasIndelOneOff(candidates: seq[string]): bool =
    if candidates.len < 2:
      return false
    let targetLen = candidates[0].len
    for i in 0..<candidates.len:
      for j in 0..<candidates.len:
        if i == j:
          continue
        let left = candidates[i]
        let right = candidates[j]
        for k in 1..<targetLen:
          let recombinant = left[0 ..< k] & right[k .. ^1]
          if recombinant.len == seqLen:
            continue
          if editDistanceAtMost(query, recombinant, 1) <= 1:
            return true
    false

  if hasIndelOneOff(eligibleLenMinus1):
    return true
  if hasIndelOneOff(eligibleLenPlus1):
    return true
  false

proc classifyBimeraFromParents(
  sequence: string,
  parents: openArray[string],
  allowOneOff = false,
  minOneOffParentDistance = 4
): BimeraCall =
  let query = toDnaUpper(sequence)

  if exactBimeraFromParents(query, parents):
    return bcExact
  if not allowOneOff:
    return bcNone

  if oneOffMismatchBimeraFromParents(query, parents, minOneOffParentDistance):
    return bcOneOffMismatch
  if oneOffIndelBimeraFromParents(query, parents, minOneOffParentDistance):
    return bcOneOffIndel
  bcNone

proc isBimera*(
  sequence: string,
  parents: openArray[string],
  allowOneOff = false,
  minOneOffParentDistance = 4
): bool =
  classifyBimeraFromParents(sequence, parents, allowOneOff, minOneOffParentDistance) != bcNone

proc isBimeraAgainstCandidates(records: seq[RankedUnique], idx: int, opts: BimeraOptions): BimeraCall =
  let query = records[idx]
  let parentAbundanceCutoff = max(
    opts.minParentAbundance.float64,
    opts.minFoldParentOverAbundance * query.abundance.float64
  )

  var parentSeqs: seq[string] = @[]
  for j, candidate in records:
    if j == idx:
      continue
    if candidate.abundance.float64 <= parentAbundanceCutoff:
      continue
    if candidate.sequence.len != query.sequence.len and
       (not opts.allowOneOff or abs(candidate.sequence.len - query.sequence.len) > 1):
      continue
    parentSeqs.add(candidate.sequence)

  if parentSeqs.len < 2:
    return bcNone

  classifyBimeraFromParents(
    sequence = query.sequence,
    parents = parentSeqs,
    allowOneOff = opts.allowOneOff,
    minOneOffParentDistance = opts.minOneOffParentDistance
  )

template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc detectBimeraRange(
  records: seq[RankedUnique],
  outCalls: ptr UncheckedArray[BimeraCall],
  startIdx: int,
  n: int,
  opts: BimeraOptions
) =
  for i in 0..<n:
    let globalIdx = startIdx + i
    outCalls[i] = isBimeraAgainstCandidates(records, globalIdx, opts)

proc isBimeraDenovo*(
  uniques: openArray[UniqueSequenceAbundance],
  opts: BimeraOptions = defaultBimeraOptions()
): BimeraDenovoResult =
  validateOptions(opts)
  let records = normalizeUniques(uniques)

  result.sequences = newSeq[string](records.len)
  result.abundances = newSeq[int64](records.len)
  result.flags = newSeq[bool](records.len)
  result.calls = newSeq[BimeraCall](records.len)
  for i, row in records:
    result.sequences[i] = row.sequence
    result.abundances[i] = row.abundance

  if records.len == 0:
    return

  if opts.threads == 1 or records.len < opts.parallelMinCandidates:
    for i in 0..<records.len:
      result.calls[i] = isBimeraAgainstCandidates(records, i, opts)
  else:
    let workerCount = min(opts.threads, records.len)
    let blockSize = (records.len + workerCount - 1) div workerCount

    var m = createMaster()
    m.awaitAll:
      var idx = 0
      while idx < records.len:
        let n = min(blockSize, records.len - idx)
        m.spawn detectBimeraRange(records, atPtr(result.calls, idx), idx, n, opts)
        idx += n

  for i, call in result.calls:
    result.flags[i] = call != bcNone
    case call
    of bcNone:
      discard
    of bcExact:
      inc result.nBimeras
      inc result.exactBimeraCount
    of bcOneOffMismatch:
      inc result.nBimeras
      inc result.oneOffMismatchBimeraCount
    of bcOneOffIndel:
      inc result.nBimeras
      inc result.oneOffIndelBimeraCount

proc pooledUniquesFromRows(rows: openArray[RankedSampleRow]): seq[UniqueSequenceAbundance] =
  var bySeq = initTable[string, int64]()
  for row in rows:
    bySeq[row.sequence] = bySeq.getOrDefault(row.sequence) + row.abundance
  result = newSeq[UniqueSequenceAbundance](bySeq.len)
  var i = 0
  for sequence, abundance in bySeq:
    result[i] = UniqueSequenceAbundance(sequence: sequence, abundance: abundance)
    inc i

proc pooledUniquesFromSampleRows(rows: openArray[SampleSequenceAbundance]): seq[UniqueSequenceAbundance] =
  var bySeq = initTable[string, int64]()
  for row in rows:
    if row.sequence.len == 0 or row.abundance <= 0:
      continue
    let seqUpper = toDnaUpper(row.sequence)
    bySeq[seqUpper] = bySeq.getOrDefault(seqUpper) + row.abundance

  result = newSeq[UniqueSequenceAbundance](bySeq.len)
  var i = 0
  for sequence, abundance in bySeq:
    result[i] = UniqueSequenceAbundance(sequence: sequence, abundance: abundance)
    inc i

proc sequenceTableVote(nFlag, nSample: int, opts: BimeraConsensusOptions): bool =
  if nSample <= 0:
    return false
  if nFlag >= nSample:
    return true
  if nFlag <= 0:
    return false
  let threshold = (nSample - opts.ignoreNNegatives).float64 * opts.minSampleFraction
  nFlag.float64 >= threshold

proc callPriority(call: BimeraCall): int {.inline.} =
  case call
  of bcExact:
    3
  of bcOneOffMismatch:
    2
  of bcOneOffIndel:
    1
  of bcNone:
    0

proc strongerCall(a, b: BimeraCall): BimeraCall {.inline.} =
  if callPriority(a) >= callPriority(b):
    a
  else:
    b

proc dominantCall(exactCount, mismatchCount, indelCount: int): BimeraCall =
  if exactCount >= mismatchCount and exactCount >= indelCount and exactCount > 0:
    return bcExact
  if mismatchCount >= indelCount and mismatchCount > 0:
    return bcOneOffMismatch
  if indelCount > 0:
    return bcOneOffIndel
  bcNone

proc isBimeraDenovoTable*(
  rows: openArray[SampleSequenceAbundance],
  bimeraOpts: BimeraOptions = defaultBimeraOptions(),
  consensusOpts: BimeraConsensusOptions = defaultBimeraConsensusOptions()
): BimeraTableResult =
  validateOptions(bimeraOpts)
  validateOptions(consensusOpts)

  let normalized = normalizeSampleRows(rows)
  if normalized.len == 0:
    return

  var bySample = initTable[string, seq[UniqueSequenceAbundance]]()
  for row in normalized:
    bySample.mgetOrPut(row.sampleId, @[]).add(
      UniqueSequenceAbundance(sequence: row.sequence, abundance: row.abundance)
    )

  var presentSamples = initTable[string, int]()
  var flaggedSamples = initTable[string, int]()
  var exactSamples = initTable[string, int]()
  var mismatchSamples = initTable[string, int]()
  var indelSamples = initTable[string, int]()

  for _, sampleUniques in bySample:
    let sampleRes = isBimeraDenovo(sampleUniques, bimeraOpts)
    for i, sequence in sampleRes.sequences:
      inc presentSamples.mgetOrPut(sequence, 0)
      if sampleRes.flags[i]:
        inc flaggedSamples.mgetOrPut(sequence, 0)
        case sampleRes.calls[i]
        of bcExact:
          inc exactSamples.mgetOrPut(sequence, 0)
        of bcOneOffMismatch:
          inc mismatchSamples.mgetOrPut(sequence, 0)
        of bcOneOffIndel:
          inc indelSamples.mgetOrPut(sequence, 0)
        of bcNone:
          discard

  result.sequences = toSeq(presentSamples.keys)
  result.sequences.sort(system.cmp[string])
  result.presentSamples = newSeq[int](result.sequences.len)
  result.flaggedSamples = newSeq[int](result.sequences.len)
  result.flags = newSeq[bool](result.sequences.len)
  result.calls = newSeq[BimeraCall](result.sequences.len)

  for i, sequence in result.sequences:
    let nSample = presentSamples.getOrDefault(sequence)
    let nFlag = flaggedSamples.getOrDefault(sequence)
    result.presentSamples[i] = nSample
    result.flaggedSamples[i] = nFlag
    result.flags[i] = sequenceTableVote(nFlag, nSample, consensusOpts)
    if result.flags[i]:
      let call = dominantCall(
        exactSamples.getOrDefault(sequence),
        mismatchSamples.getOrDefault(sequence),
        indelSamples.getOrDefault(sequence)
      )
      result.calls[i] = call
      inc result.nBimeras
      case call
      of bcExact:
        inc result.exactBimeraCount
      of bcOneOffMismatch:
        inc result.oneOffMismatchBimeraCount
      of bcOneOffIndel:
        inc result.oneOffIndelBimeraCount
      of bcNone:
        discard
    else:
      result.calls[i] = bcNone

proc removeBimeraDenovo*(
  uniques: openArray[UniqueSequenceAbundance],
  opts: BimeraOptions = defaultBimeraOptions()
): seq[UniqueSequenceAbundance] =
  let denovo = isBimeraDenovo(uniques, opts)
  result = @[]
  for i, sequence in denovo.sequences:
    if not denovo.flags[i]:
      result.add(UniqueSequenceAbundance(sequence: sequence, abundance: denovo.abundances[i]))

proc removeBimeraDenovo*(
  rows: openArray[SampleSequenceAbundance],
  mode: BimeraMethod = bmConsensus,
  bimeraOpts: BimeraOptions = defaultBimeraOptions(),
  consensusOpts: BimeraConsensusOptions = defaultBimeraConsensusOptions()
): RemoveBimeraResult =
  validateOptions(bimeraOpts)
  validateOptions(consensusOpts)

  let normalized = normalizeSampleRows(rows)
  result.mode = mode
  result.inputUniqueCount = normalizeUniques(pooledUniquesFromRows(normalized)).len

  if normalized.len == 0:
    result.kept = @[]
    result.removed = @[]
    return

  case mode
  of bmPooled:
    let pooled = isBimeraDenovo(pooledUniquesFromRows(normalized), bimeraOpts)
    result.pooled = pooled
    result.exactBimeraCount = pooled.exactBimeraCount
    result.oneOffMismatchBimeraCount = pooled.oneOffMismatchBimeraCount
    result.oneOffIndelBimeraCount = pooled.oneOffIndelBimeraCount

    var flagged = initTable[string, bool]()
    for i, sequence in pooled.sequences:
      if pooled.flags[i]:
        flagged[sequence] = true

    for row in normalized:
      let outRow = SampleSequenceAbundance(sampleId: row.sampleId, sequence: row.sequence, abundance: row.abundance)
      if flagged.getOrDefault(row.sequence, false):
        result.removed.add(outRow)
      else:
        result.kept.add(outRow)

  of bmConsensus:
    let tableRes = isBimeraDenovoTable(rows, bimeraOpts, consensusOpts)
    result.consensus = tableRes
    result.exactBimeraCount = tableRes.exactBimeraCount
    result.oneOffMismatchBimeraCount = tableRes.oneOffMismatchBimeraCount
    result.oneOffIndelBimeraCount = tableRes.oneOffIndelBimeraCount

    var flagged = initTable[string, bool]()
    for i, sequence in tableRes.sequences:
      if tableRes.flags[i]:
        flagged[sequence] = true

    for row in normalized:
      let outRow = SampleSequenceAbundance(sampleId: row.sampleId, sequence: row.sequence, abundance: row.abundance)
      if flagged.getOrDefault(row.sequence, false):
        result.removed.add(outRow)
      else:
        result.kept.add(outRow)

  of bmPerSample:
    var bySample = initTable[string, seq[UniqueSequenceAbundance]]()
    for row in normalized:
      bySample.mgetOrPut(row.sampleId, @[]).add(
        UniqueSequenceAbundance(sequence: row.sequence, abundance: row.abundance)
      )

    var sampleFlagged = initTable[string, BimeraCall]()
    var uniqueFlaggedCall = initTable[string, BimeraCall]()
    for sampleId, sampleUniques in bySample:
      let sampleRes = isBimeraDenovo(sampleUniques, bimeraOpts)
      for i, sequence in sampleRes.sequences:
        let call = sampleRes.calls[i]
        if call != bcNone:
          sampleFlagged[sampleId & "\t" & sequence] = call
          uniqueFlaggedCall[sequence] = strongerCall(uniqueFlaggedCall.getOrDefault(sequence, bcNone), call)

    for _, call in uniqueFlaggedCall:
      case call
      of bcExact:
        inc result.exactBimeraCount
      of bcOneOffMismatch:
        inc result.oneOffMismatchBimeraCount
      of bcOneOffIndel:
        inc result.oneOffIndelBimeraCount
      of bcNone:
        discard

    for row in normalized:
      let outRow = SampleSequenceAbundance(sampleId: row.sampleId, sequence: row.sequence, abundance: row.abundance)
      if sampleFlagged.getOrDefault(row.sampleId & "\t" & row.sequence, bcNone) != bcNone:
        result.removed.add(outRow)
      else:
        result.kept.add(outRow)

  let keptUniq = normalizeUniques(pooledUniquesFromSampleRows(result.kept))
  result.keptUniqueCount = keptUniq.len
  result.removedUniqueCount = max(0, result.inputUniqueCount - result.keptUniqueCount)

proc writeSampleSequenceAbundanceTsv*(
  path: string,
  rows: openArray[SampleSequenceAbundance]
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open output file: " & path)
  defer:
    close(f)

  f.writeLine("SampleID\tSequence\tAbundance")
  for row in rows:
    f.writeLine(row.sampleId & "\t" & row.sequence & "\t" & $row.abundance)

proc writeUniqueSequenceAbundanceTsv*(
  path: string,
  rows: openArray[UniqueSequenceAbundance]
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open output file: " & path)
  defer:
    close(f)

  f.writeLine("Sequence\tAbundance")
  for row in rows:
    f.writeLine(row.sequence & "\t" & $row.abundance)
