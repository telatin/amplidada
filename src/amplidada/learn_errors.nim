import std/[algorithm, math, os, sets, strformat, strutils, tables]
import amplidada/derep

const
  DnaBases* = ['A', 'C', 'G', 'T']

type
  LearnErrorsOptions* = object
    qMin*: int
    qMax*: int
    pseudocount*: float64
    zeroTotalMismatchProb*: float64

  LearnErrorsSelfConsistOptions* = object
    enabled*: bool
    minIterations*: int
    maxIterations*: int
    tol*: float64
    maxCentersPerLength*: int

  LearnErrorsResult* = object
    qMin*: int
    qMax*: int
    counts*: array[16, seq[int64]]
    probs*: array[16, seq[float64]]
    totalTransitions*: int64
    skippedPositions*: int64
    centerCountByLength*: int

  LearnErrorsSelfConsistResult* = object
    matrix*: LearnErrorsResult
    iterationsRun*: int
    converged*: bool
    maxAbsProbDelta*: float64
    centerChangesTotal*: int
    centersByLength*: seq[tuple[seqLen: int, center: string]]

  LearnErrorsBatchResult* = object
    matrix*: LearnErrorsResult
    filesRequested*: int
    filesUsed*: int
    readsUsed*: int64
    basesUsed*: int64
    uniqueCountSum*: int
    selfConsistEnabled*: bool
    iterationsRun*: int
    converged*: bool
    maxAbsProbDelta*: float64
    centerChangesTotal*: int

type
  UniqueObs = object
    sequence: string
    abundance: int64
    meanQual: seq[float64]
    firstSeenGlobal: int64

  CenterCandidate = object
    sequence: string
    abundance: int64
    firstSeenGlobal: int64

proc defaultLearnErrorsOptions*(): LearnErrorsOptions =
  LearnErrorsOptions(
    qMin: 0,
    qMax: 40,
    pseudocount: 1.0,
    zeroTotalMismatchProb: 1e-7
  )

proc defaultLearnErrorsSelfConsistOptions*(): LearnErrorsSelfConsistOptions =
  LearnErrorsSelfConsistOptions(
    enabled: false,
    minIterations: 1,
    maxIterations: 8,
    tol: 1e-8,
    maxCentersPerLength: 500
  )

proc validateOptions*(opts: LearnErrorsOptions) =
  if opts.qMin < 0:
    raise newException(ValueError, "qMin must be >= 0")
  if opts.qMax < opts.qMin:
    raise newException(ValueError, "qMax must be >= qMin")
  if opts.pseudocount < 0:
    raise newException(ValueError, "pseudocount must be >= 0")
  if opts.zeroTotalMismatchProb < 0 or opts.zeroTotalMismatchProb >= 1.0 / 3.0:
    raise newException(ValueError, "zeroTotalMismatchProb must be in [0, 1/3)")

proc validateOptions*(opts: LearnErrorsSelfConsistOptions) =
  if opts.minIterations <= 0:
    raise newException(ValueError, "minIterations must be > 0")
  if opts.maxIterations <= 0:
    raise newException(ValueError, "maxIterations must be > 0")
  if opts.minIterations > opts.maxIterations:
    raise newException(ValueError, "minIterations must be <= maxIterations")
  if opts.tol < 0:
    raise newException(ValueError, "tol must be >= 0")
  if opts.maxCentersPerLength <= 0:
    raise newException(ValueError, "maxCentersPerLength must be > 0")

proc describe*(opts: LearnErrorsOptions): string =
  validateOptions(opts)
  &"qMin={opts.qMin} qMax={opts.qMax} pseudocount={opts.pseudocount} " &
    &"zeroTotalMismatchProb={opts.zeroTotalMismatchProb}"

proc describe*(opts: LearnErrorsSelfConsistOptions): string =
  validateOptions(opts)
  &"enabled={opts.enabled} minIterations={opts.minIterations} maxIterations={opts.maxIterations} " &
    &"tol={opts.tol} maxCentersPerLength={opts.maxCentersPerLength}"

proc baseIndex*(c: char): int =
  case c
  of 'A', 'a': 0
  of 'C', 'c': 1
  of 'G', 'g': 2
  of 'T', 't': 3
  else: -1

proc transitionIndex*(fromBase, toBase: char): int =
  let f = baseIndex(fromBase)
  let t = baseIndex(toBase)
  if f < 0 or t < 0:
    return -1
  f * 4 + t

proc qualityToBinIndex(q: float64; qMin, qMax: int): int =
  var qi = int(floor(q + 0.5))
  if qi < qMin:
    qi = qMin
  elif qi > qMax:
    qi = qMax
  qi - qMin

proc getCount*(res: LearnErrorsResult, fromBase, toBase: char, q: int): int64 =
  let idx = transitionIndex(fromBase, toBase)
  if idx < 0:
    raise newException(ValueError, "Invalid DNA base for count query")
  if q < res.qMin or q > res.qMax:
    raise newException(ValueError, "Quality outside matrix range")
  res.counts[idx][q - res.qMin]

proc getProb*(res: LearnErrorsResult, fromBase, toBase: char, q: int): float64 =
  let idx = transitionIndex(fromBase, toBase)
  if idx < 0:
    raise newException(ValueError, "Invalid DNA base for probability query")
  if q < res.qMin or q > res.qMax:
    raise newException(ValueError, "Quality outside matrix range")
  res.probs[idx][q - res.qMin]

proc derepBaseCount*(derepRes: DerepFastqResult): int64 =
  for u in derepRes.uniques:
    result += u.abundance * u.sequence.len.int64

proc writeLearnErrorsTsv*(path: string, res: LearnErrorsResult) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open error matrix output file: " & path)
  defer: close(f)

  f.writeLine("from\tto\tquality\tcount\tprob")

  for fromIdx in 0..3:
    for toIdx in 0..3:
      let idx = fromIdx * 4 + toIdx
      for q in res.qMin..res.qMax:
        let qIdx = q - res.qMin
        f.writeLine(
          &"{DnaBases[fromIdx]}\t{DnaBases[toIdx]}\t{q}\t{res.counts[idx][qIdx]}\t" &
          $res.probs[idx][qIdx]
        )

proc writeLearnErrorsCountsTsv*(path: string, res: LearnErrorsResult) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open transition counts output file: " & path)
  defer: close(f)

  f.writeLine("from\tto\tquality\tcount")
  for fromIdx in 0..3:
    for toIdx in 0..3:
      let idx = fromIdx * 4 + toIdx
      for q in res.qMin..res.qMax:
        let qIdx = q - res.qMin
        f.writeLine(&"{DnaBases[fromIdx]}\t{DnaBases[toIdx]}\t{q}\t{res.counts[idx][qIdx]}")

proc initResult(qMin, qMax: int): LearnErrorsResult

proc readLearnErrorsTsv*(path: string): LearnErrorsResult =
  if not fileExists(path):
    raise newException(IOError, "Error matrix TSV does not exist: " & path)

  var loaded = false
  var qMin = high(int)
  var qMax = low(int)

  type MatrixRec = object
    fromIdx: int
    toIdx: int
    q: int
    count: int64
    prob: float64

  var rows: seq[MatrixRec] = @[]

  var lineNo = 0
  for line in lines(path):
    inc lineNo
    let t = line.strip()
    if t.len == 0:
      continue
    if lineNo == 1 and t == "from\tto\tquality\tcount\tprob":
      continue

    let parts = t.split('\t')
    if parts.len < 5:
      raise newException(ValueError, &"Malformed matrix line {lineNo}: expected 5 columns")

    let fromC = if parts[0].len > 0: parts[0][0] else: '\0'
    let toC = if parts[1].len > 0: parts[1][0] else: '\0'
    let fromIdx = baseIndex(fromC)
    let toIdx = baseIndex(toC)
    if fromIdx < 0 or toIdx < 0:
      raise newException(ValueError, &"Malformed matrix line {lineNo}: invalid DNA base")

    var q: int
    var count: int64
    var prob: float64
    try:
      q = parseInt(parts[2].strip())
      count = parseBiggestInt(parts[3].strip()).int64
      prob = parseFloat(parts[4].strip())
    except ValueError:
      raise newException(ValueError, &"Malformed numeric value in matrix line {lineNo}")

    rows.add(MatrixRec(fromIdx: fromIdx, toIdx: toIdx, q: q, count: count, prob: prob))
    qMin = min(qMin, q)
    qMax = max(qMax, q)
    loaded = true

  if not loaded:
    raise newException(ValueError, "Error matrix TSV is empty: " & path)

  result = initResult(qMin, qMax)

  for rec in rows:
    let qIdx = rec.q - qMin
    let rowIdx = rec.fromIdx * 4 + rec.toIdx
    result.counts[rowIdx][qIdx] = rec.count
    result.probs[rowIdx][qIdx] = rec.prob
    result.totalTransitions += rec.count

  result.centerCountByLength = 0

proc cmpCandidate(a, b: CenterCandidate): int =
  if a.abundance != b.abundance:
    return cmp(b.abundance, a.abundance)
  if a.firstSeenGlobal != b.firstSeenGlobal:
    return cmp(a.firstSeenGlobal, b.firstSeenGlobal)
  cmp(a.sequence, b.sequence)

proc isBetterCandidate(a, b: CenterCandidate): bool =
  cmpCandidate(a, b) < 0

proc initResult(qMin, qMax: int): LearnErrorsResult =
  result.qMin = qMin
  result.qMax = qMax
  let qBins = result.qMax - result.qMin + 1
  for i in 0..15:
    result.counts[i] = newSeq[int64](qBins)
    result.probs[i] = newSeq[float64](qBins)

proc fillProbabilities(res: var LearnErrorsResult, opts: LearnErrorsOptions) =
  let qBins = res.qMax - res.qMin + 1

  for fromIdx in 0..3:
    for qIdx in 0..<qBins:
      var totalAtQ = 0'i64
      for toIdx in 0..3:
        totalAtQ += res.counts[fromIdx * 4 + toIdx][qIdx]

      if totalAtQ == 0:
        for toIdx in 0..3:
          if toIdx == fromIdx:
            res.probs[fromIdx * 4 + toIdx][qIdx] = 1.0 - 3.0 * opts.zeroTotalMismatchProb
          else:
            res.probs[fromIdx * 4 + toIdx][qIdx] = opts.zeroTotalMismatchProb
      else:
        let denom = totalAtQ.float64 + 4.0 * opts.pseudocount
        for toIdx in 0..3:
          res.probs[fromIdx * 4 + toIdx][qIdx] =
            (res.counts[fromIdx * 4 + toIdx][qIdx].float64 + opts.pseudocount) / denom

proc collectUniqueObs(derepResults: openArray[DerepFastqResult]): seq[UniqueObs] =
  var globalReadOffset = 0'i64
  for derepRes in derepResults:
    for u in derepRes.uniques:
      if u.sequence.len != u.meanQual.len:
        raise newException(ValueError, "Derep unique has meanQual length mismatch")
      result.add(UniqueObs(
        sequence: u.sequence,
        abundance: u.abundance,
        meanQual: u.meanQual,
        firstSeenGlobal: globalReadOffset + u.firstSeenReadIndex
      ))
    globalReadOffset += derepRes.totalReads

proc buildLengthGroups(uniqueObs: openArray[UniqueObs]): Table[int, seq[int]] =
  result = initTable[int, seq[int]]()
  for i in 0..<uniqueObs.len:
    let seqLen = uniqueObs[i].sequence.len
    var bucket = result.mgetOrPut(seqLen, @[])
    bucket.add(i)

proc buildCandidatePools(
  uniqueObs: openArray[UniqueObs],
  maxCentersPerLength: int
): Table[int, seq[CenterCandidate]] =
  result = initTable[int, seq[CenterCandidate]]()
  var agg = initTable[(int, string), CenterCandidate]()
  var lengths = initHashSet[int]()

  for u in uniqueObs:
    let seqLen = u.sequence.len
    let key = (seqLen, u.sequence)
    lengths.incl(seqLen)

    if key in agg:
      var a = agg[key]
      a.abundance += u.abundance
      if u.firstSeenGlobal < a.firstSeenGlobal:
        a.firstSeenGlobal = u.firstSeenGlobal
      agg[key] = a
    else:
      agg[key] = CenterCandidate(
        sequence: u.sequence,
        abundance: u.abundance,
        firstSeenGlobal: u.firstSeenGlobal
      )

  for seqLen in lengths:
    var candidates: seq[CenterCandidate] = @[]
    for key, cand in agg:
      if key[0] == seqLen:
        candidates.add(cand)
    candidates.sort(cmpCandidate)
    if candidates.len > maxCentersPerLength:
      candidates.setLen(maxCentersPerLength)
    result[seqLen] = candidates

proc topCentersByLength(candidatePools: Table[int, seq[CenterCandidate]]): Table[int, string] =
  result = initTable[int, string]()
  for seqLen, candidates in candidatePools:
    if candidates.len > 0:
      result[seqLen] = candidates[0].sequence

proc singleCentersToSeq(centers: Table[int, string]): seq[tuple[seqLen: int, center: string]] =
  for seqLen, c in centers:
    result.add((seqLen: seqLen, center: c))
  result.sort(proc(a, b: tuple[seqLen: int, center: string]): int =
    if a.seqLen != b.seqLen:
      return cmp(a.seqLen, b.seqLen)
    cmp(a.center, b.center)
  )

proc hammingDistance(a, b: string): int =
  let n = min(a.len, b.len)
  for i in 0..<n:
    if a[i] != b[i]: inc result
  result += abs(a.len - b.len)

proc nearestCenterIndex(s: string, centers: seq[string]): int =
  var bestDist = high(int)
  for i, c in centers:
    let d = hammingDistance(s, c)
    if d < bestDist:
      bestDist = d
      result = i
      if bestDist == 0: break

proc mostLikelyCenterIndex(
  u: UniqueObs,
  centers: seq[string],
  matrix: LearnErrorsResult
): int =
  const minProb = 1e-300
  var bestScore = -Inf
  for i, c in centers:
    var score = 0.0
    let n = min(u.sequence.len, c.len)
    for pos in 0..<n:
      let fromIdx = baseIndex(c[pos])
      let toIdx = baseIndex(u.sequence[pos])
      let p =
        if fromIdx < 0 or toIdx < 0: minProb
        else:
          let qIdx = qualityToBinIndex(u.meanQual[pos], matrix.qMin, matrix.qMax)
          max(minProb, matrix.probs[fromIdx * 4 + toIdx][qIdx])
      score += ln(p)
    if score > bestScore:
      bestScore = score
      result = i

proc initialMultiCentersByLength(
  uniqueObs: openArray[UniqueObs],
  maxCentersPerLength: int
): Table[int, seq[string]] =
  result = initTable[int, seq[string]]()
  let pools = buildCandidatePools(uniqueObs, maxCentersPerLength)
  for seqLen, candidates in pools:
    var seqs = newSeq[string](candidates.len)
    for i, cand in candidates:
      seqs[i] = cand.sequence
    result[seqLen] = seqs

proc multiCentersToSeq(
  centers: Table[int, seq[string]]
): seq[tuple[seqLen: int, center: string]] =
  for seqLen, seqList in centers:
    for c in seqList:
      result.add((seqLen: seqLen, center: c))
  result.sort(proc(a, b: tuple[seqLen: int, center: string]): int =
    if a.seqLen != b.seqLen:
      return cmp(a.seqLen, b.seqLen)
    cmp(a.center, b.center)
  )

proc accumulateFromCenters(
  uniqueObs: openArray[UniqueObs],
  centersByLength: Table[int, string],
  opts: LearnErrorsOptions
): LearnErrorsResult =
  result = initResult(opts.qMin, opts.qMax)
  result.centerCountByLength = centersByLength.len

  for u in uniqueObs:
    let seqLen = u.sequence.len
    if seqLen notin centersByLength:
      continue
    let centerSeq = centersByLength[seqLen]

    for pos in 0..<seqLen:
      let fromIdx = baseIndex(centerSeq[pos])
      let toIdx = baseIndex(u.sequence[pos])
      if fromIdx < 0 or toIdx < 0:
        result.skippedPositions += u.abundance
        continue

      let qIdx = qualityToBinIndex(u.meanQual[pos], result.qMin, result.qMax)
      result.counts[fromIdx * 4 + toIdx][qIdx] += u.abundance
      result.totalTransitions += u.abundance

  fillProbabilities(result, opts)

proc accumulateFromMultiCenters(
  uniqueObs: openArray[UniqueObs],
  centersByLength: Table[int, seq[string]],
  opts: LearnErrorsOptions
): LearnErrorsResult =
  result = initResult(opts.qMin, opts.qMax)
  var totalCenters = 0
  for _, cs in centersByLength:
    totalCenters += cs.len
  result.centerCountByLength = totalCenters

  for u in uniqueObs:
    let seqLen = u.sequence.len
    if seqLen notin centersByLength:
      continue
    let centers = centersByLength[seqLen]
    if centers.len == 0:
      continue
    let cIdx = nearestCenterIndex(u.sequence, centers)
    let centerSeq = centers[cIdx]

    for pos in 0..<seqLen:
      let fromIdx = baseIndex(centerSeq[pos])
      let toIdx = baseIndex(u.sequence[pos])
      if fromIdx < 0 or toIdx < 0:
        result.skippedPositions += u.abundance
        continue
      let qIdx = qualityToBinIndex(u.meanQual[pos], result.qMin, result.qMax)
      result.counts[fromIdx * 4 + toIdx][qIdx] += u.abundance
      result.totalTransitions += u.abundance

  fillProbabilities(result, opts)

proc accumulateFromMultiCentersLikelihood(
  uniqueObs: openArray[UniqueObs],
  centersByLength: Table[int, seq[string]],
  matrix: LearnErrorsResult,
  opts: LearnErrorsOptions
): LearnErrorsResult =
  result = initResult(opts.qMin, opts.qMax)
  var totalCenters = 0
  for _, cs in centersByLength:
    totalCenters += cs.len
  result.centerCountByLength = totalCenters

  for u in uniqueObs:
    let seqLen = u.sequence.len
    if seqLen notin centersByLength:
      continue
    let centers = centersByLength[seqLen]
    if centers.len == 0:
      continue
    let cIdx = mostLikelyCenterIndex(u, centers, matrix)
    let centerSeq = centers[cIdx]

    for pos in 0..<seqLen:
      let fromIdx = baseIndex(centerSeq[pos])
      let toIdx = baseIndex(u.sequence[pos])
      if fromIdx < 0 or toIdx < 0:
        result.skippedPositions += u.abundance
        continue
      let qIdx = qualityToBinIndex(u.meanQual[pos], result.qMin, result.qMax)
      result.counts[fromIdx * 4 + toIdx][qIdx] += u.abundance
      result.totalTransitions += u.abundance

  fillProbabilities(result, opts)

proc selectCentersByLikelihood(
  uniqueObs: openArray[UniqueObs],
  lengthGroups: Table[int, seq[int]],
  candidatePools: Table[int, seq[CenterCandidate]],
  matrix: LearnErrorsResult
): Table[int, string] =
  const minProb = 1e-300
  result = initTable[int, string]()

  for seqLen, obsIndexes in lengthGroups:
    if seqLen notin candidatePools or candidatePools[seqLen].len == 0:
      continue

    let candidates = candidatePools[seqLen]
    var bestCand = candidates[0]
    var bestScore = -Inf
    var haveBest = false

    for cand in candidates:
      var score = 0.0
      for obsIndex in obsIndexes:
        let u = uniqueObs[obsIndex]
        for pos in 0..<seqLen:
          let fromIdx = baseIndex(cand.sequence[pos])
          let toIdx = baseIndex(u.sequence[pos])
          let p =
            if fromIdx < 0 or toIdx < 0:
              minProb
            else:
              let qIdx = qualityToBinIndex(u.meanQual[pos], matrix.qMin, matrix.qMax)
              max(minProb, matrix.probs[fromIdx * 4 + toIdx][qIdx])
          score += u.abundance.float64 * ln(p)

      if (not haveBest) or (score > bestScore) or
         (score == bestScore and isBetterCandidate(cand, bestCand)):
        bestCand = cand
        bestScore = score
        haveBest = true

    result[seqLen] = bestCand.sequence

proc maxProbDelta(a, b: LearnErrorsResult): float64 =
  if a.qMin != b.qMin or a.qMax != b.qMax:
    raise newException(ValueError, "Cannot compare matrices with different quality ranges")
  for i in 0..15:
    for j in 0..<a.probs[i].len:
      result = max(result, abs(a.probs[i][j] - b.probs[i][j]))

proc learnErrorsFromDereps*(
  derepResults: openArray[DerepFastqResult],
  opts: LearnErrorsOptions = defaultLearnErrorsOptions()
): LearnErrorsResult =
  validateOptions(opts)
  let uniqueObs = collectUniqueObs(derepResults)
  let pools = buildCandidatePools(uniqueObs, 1)
  let centers = topCentersByLength(pools)
  accumulateFromCenters(uniqueObs, centers, opts)

proc learnErrorsFromDerep*(
  derepRes: DerepFastqResult,
  opts: LearnErrorsOptions = defaultLearnErrorsOptions()
): LearnErrorsResult =
  learnErrorsFromDereps([derepRes], opts)

proc learnErrorsSelfConsistentFromDereps*(
  derepResults: openArray[DerepFastqResult],
  learnOpts: LearnErrorsOptions = defaultLearnErrorsOptions(),
  scOpts: LearnErrorsSelfConsistOptions = defaultLearnErrorsSelfConsistOptions()
): LearnErrorsSelfConsistResult =
  validateOptions(learnOpts)
  validateOptions(scOpts)

  let uniqueObs = collectUniqueObs(derepResults)
  let candidatePools = buildCandidatePools(uniqueObs, scOpts.maxCentersPerLength)
  let lengthGroups = buildLengthGroups(uniqueObs)

  if not scOpts.enabled:
    let initialCenters = topCentersByLength(candidatePools)
    result.matrix = accumulateFromCenters(uniqueObs, initialCenters, learnOpts)
    result.iterationsRun = (if uniqueObs.len == 0: 0 else: 1)
    result.converged = true
    result.maxAbsProbDelta = 0.0
    result.centersByLength = singleCentersToSeq(initialCenters)
    return

  var previousMatrix = initResult(learnOpts.qMin, learnOpts.qMax)
  var havePrevious = false
  var currentCenters = initTable[int, string]()

  for iter in 1..scOpts.maxIterations:
    currentCenters =
      if iter == 1:
        topCentersByLength(candidatePools)
      else:
        selectCentersByLikelihood(uniqueObs, lengthGroups, candidatePools, previousMatrix)

    let currentMatrix = accumulateFromCenters(uniqueObs, currentCenters, learnOpts)

    var delta = 0.0
    if havePrevious:
      delta = maxProbDelta(previousMatrix, currentMatrix)

    result.matrix = currentMatrix
    result.iterationsRun = iter
    result.maxAbsProbDelta = delta

    if iter >= scOpts.minIterations and delta <= scOpts.tol:
      result.converged = true
      break

    if iter == scOpts.maxIterations:
      break

    previousMatrix = currentMatrix
    havePrevious = true

  result.centersByLength = singleCentersToSeq(currentCenters)

proc learnErrorsFromFastqPaths*(
  inputPaths: openArray[string],
  derepOpts: DerepFastqOptions = defaultDerepFastqOptions(),
  learnOpts: LearnErrorsOptions = defaultLearnErrorsOptions(),
  nbases: int64 = 100_000_000,
  selfConsistOpts: LearnErrorsSelfConsistOptions = defaultLearnErrorsSelfConsistOptions()
): LearnErrorsBatchResult =
  validateOptions(derepOpts)
  validateOptions(learnOpts)
  validateOptions(selfConsistOpts)

  if inputPaths.len == 0:
    raise newException(ValueError, "learnErrors requires at least one input FASTQ path")

  if nbases == 0:
    raise newException(ValueError, "nbases must be > 0 or negative to disable cutoff")

  result.filesRequested = inputPaths.len
  result.selfConsistEnabled = selfConsistOpts.enabled

  var dereps: seq[DerepFastqResult] = @[]
  for inputPath in inputPaths:
    let d = derepFastq(inputPath, derepOpts)
    dereps.add(d)
    inc result.filesUsed
    result.readsUsed += d.totalReads
    result.basesUsed += derepBaseCount(d)
    result.uniqueCountSum += d.uniqueCount

    if nbases > 0 and result.basesUsed >= nbases:
      break

  let scRes = learnErrorsSelfConsistentFromDereps(dereps, learnOpts, selfConsistOpts)
  result.matrix = scRes.matrix
  result.iterationsRun = scRes.iterationsRun
  result.converged = scRes.converged
  result.maxAbsProbDelta = scRes.maxAbsProbDelta
  result.centerChangesTotal = scRes.centerChangesTotal
