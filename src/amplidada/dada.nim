import std/[algorithm, cpuinfo, math, sets, strformat, tables]
import malebolgia
import amplidada/[derep, learn_errors]

type
  DadaOptions* = object
    omegaA*: float64
    omegaC*: float64
    minAbundance*: int64
    minHamming*: int
    minFold*: float64
    maxClusters*: int
    maxIterations*: int
    maxShuffle*: int
    useQuality*: bool
    threads*: int
    parallelMinGroups*: int

  DadaSelfConsistOptions* = object
    enabled*: bool
    minIterations*: int
    maxIterations*: int
    tol*: float64
    forceOmegaCZero*: bool
    finalDenoiseWithLatestErrors*: bool

  AsvRecord* = object
    sequence*: string
    abundance*: int64
    clusterSize*: int
    centerUniqueIndex*: int
    firstSeenReadIndex*: int64

  DadaResult* = object
    totalReads*: int64
    inputUniqueCount*: int
    asvCount*: int
    asvs*: seq[AsvRecord]
    uniqueToAsv*: seq[int]
    uniqueToCenterUnique*: seq[int]
    groupCount*: int
    iterationsRun*: int
    splitsAccepted*: int
    correctedUniques*: int
    uncorrectedUniques*: int

  DadaSelfConsistResult* = object
    errorMatrix*: LearnErrorsResult
    dada*: DadaResult
    iterationsRun*: int
    converged*: bool
    maxAbsProbDelta*: float64
    assignmentChanges*: int

type
  GroupDenoiseResult = object
    localToGlobal: seq[int]
    centers: seq[int]
    assignments: seq[int]
    clusterTotals: seq[int64]
    iterationsRun: int
    splitsAccepted: int

proc defaultDadaOptions*(): DadaOptions =
  DadaOptions(
    omegaA: 1e-40,
    omegaC: 1e-40,
    minAbundance: 2,
    minHamming: 1,
    minFold: 1.0,
    maxClusters: 0,
    maxIterations: 256,
    maxShuffle: 16,
    useQuality: true,
    threads: max(1, countProcessors()),
    parallelMinGroups: 2
  )

proc defaultDadaSelfConsistOptions*(): DadaSelfConsistOptions =
  DadaSelfConsistOptions(
    enabled: false,
    minIterations: 1,
    maxIterations: 8,
    tol: 1e-8,
    forceOmegaCZero: true,
    finalDenoiseWithLatestErrors: true
  )

proc validateOptions*(opts: DadaOptions) =
  if opts.omegaA <= 0 or opts.omegaA > 1:
    raise newException(ValueError, "omegaA must be in (0, 1]")
  if opts.omegaC < 0 or opts.omegaC > 1:
    raise newException(ValueError, "omegaC must be in [0, 1]")
  if opts.minAbundance <= 0:
    raise newException(ValueError, "minAbundance must be > 0")
  if opts.minHamming < 0:
    raise newException(ValueError, "minHamming must be >= 0")
  if opts.minFold < 0:
    raise newException(ValueError, "minFold must be >= 0")
  if opts.maxClusters < 0:
    raise newException(ValueError, "maxClusters must be >= 0")
  if opts.maxIterations <= 0:
    raise newException(ValueError, "maxIterations must be > 0")
  if opts.maxShuffle <= 0:
    raise newException(ValueError, "maxShuffle must be > 0")
  if opts.threads <= 0:
    raise newException(ValueError, "threads must be > 0")
  if opts.parallelMinGroups <= 0:
    raise newException(ValueError, "parallelMinGroups must be > 0")

proc validateOptions*(opts: DadaSelfConsistOptions) =
  if opts.minIterations <= 0:
    raise newException(ValueError, "minIterations must be > 0")
  if opts.maxIterations <= 0:
    raise newException(ValueError, "maxIterations must be > 0")
  if opts.minIterations > opts.maxIterations:
    raise newException(ValueError, "minIterations must be <= maxIterations")
  if opts.tol < 0:
    raise newException(ValueError, "tol must be >= 0")

proc describe*(opts: DadaOptions): string =
  validateOptions(opts)
  &"omegaA={opts.omegaA} omegaC={opts.omegaC} minAbundance={opts.minAbundance} minHamming={opts.minHamming} " &
    &"minFold={opts.minFold} maxClusters={opts.maxClusters} maxIterations={opts.maxIterations} " &
    &"maxShuffle={opts.maxShuffle} useQuality={opts.useQuality} threads={opts.threads} " &
    &"parallelMinGroups={opts.parallelMinGroups}"

proc describe*(opts: DadaSelfConsistOptions): string =
  validateOptions(opts)
  &"enabled={opts.enabled} minIterations={opts.minIterations} maxIterations={opts.maxIterations} " &
    &"tol={opts.tol} forceOmegaCZero={opts.forceOmegaCZero} " &
    &"finalDenoiseWithLatestErrors={opts.finalDenoiseWithLatestErrors}"

proc hammingDistance(a, b: string): int =
  if a.len != b.len:
    return max(a.len, b.len)
  for i in 0..<a.len:
    if a[i] != b[i]:
      inc result

proc qBinIndex(q: float64, qMin, qMax: int): int =
  var qi = int(floor(q + 0.5))
  if qi < qMin:
    qi = qMin
  elif qi > qMax:
    qi = qMax
  qi - qMin

proc poissonTailAtLeast(k: int64, lambda: float64): float64 =
  if k <= 0:
    return 1.0
  if lambda <= 0:
    return 0.0

  let kk = int(k)
  if kk <= 512 and lambda <= 700:
    if lambda < float64(kk):
      # Direct sum P(X=k) + P(X=k+1) + ... avoids catastrophic cancellation
      # in 1 - CDF when lambda is small relative to k.
      let logBase = -lambda + float64(kk) * ln(lambda) - lgamma(float64(kk) + 1.0)
      if logBase < -740.0:
        return 0.0
      var term = exp(logBase)
      var total = term
      for j in kk + 1..kk + 300:
        term *= lambda / float64(j)
        total += term
        if term < total * 1e-15:
          break
      return min(1.0, max(0.0, total))
    else:
      var term = exp(-lambda)  # P(X=0)
      var cdf = term
      for i in 1..<kk:
        term *= lambda / float64(i)
        cdf += term
        if term == 0.0:
          break
      return min(1.0, max(0.0, 1.0 - cdf))

  let z = (k.float64 - 0.5 - lambda) / sqrt(lambda)
  let tail = 0.5 * erfc(z / sqrt(2.0))
  min(1.0, max(0.0, tail))

proc abundancePValue*(reads: int64, expected: float64, prior = false): float64 =
  if reads <= 0:
    return 1.0
  if expected <= 0.0:
    return 0.0

  let tail = poissonTailAtLeast(reads, expected)
  if prior:
    return tail

  var denom: float64
  if expected < 1e-5:
    denom = expected - 0.5 * expected * expected
  else:
    denom = 1.0 - exp(-expected)

  if denom <= 0.0:
    return 1.0
  min(1.0, max(0.0, tail / denom))

proc cmpGlobalUnique(aGlobal, bGlobal: int, uniques: seq[DerepUnique]): int =
  let a = uniques[aGlobal]
  let b = uniques[bGlobal]
  if a.abundance != b.abundance:
    return cmp(b.abundance, a.abundance)
  if a.firstSeenReadIndex != b.firstSeenReadIndex:
    return cmp(a.firstSeenReadIndex, b.firstSeenReadIndex)
  cmp(a.sequence, b.sequence)

proc betterCenterLocal(aLocal, bLocal: int, localToGlobal: seq[int], uniques: seq[DerepUnique]): bool =
  cmpGlobalUnique(localToGlobal[aLocal], localToGlobal[bLocal], uniques) < 0

proc clusterTotals(assignments: seq[int], nClusters: int, localToGlobal: seq[int], uniques: seq[DerepUnique]): seq[int64] =
  result = newSeq[int64](nClusters)
  for localIdx in 0..<assignments.len:
    let c = assignments[localIdx]
    if c >= 0 and c < nClusters:
      result[c] += uniques[localToGlobal[localIdx]].abundance

proc lambdaCacheKey(centerLocal, obsLocal: int): uint64 =
  (uint64(centerLocal.uint32) shl 32) or uint64(obsLocal.uint32)

proc lambdaCenterToObs(
  centerLocal, obsLocal: int,
  localToGlobal: seq[int],
  uniques: seq[DerepUnique],
  err: LearnErrorsResult,
  useQuality: bool,
  cache: var Table[uint64, float64]
): float64 =
  let key = lambdaCacheKey(centerLocal, obsLocal)
  if key in cache:
    return cache[key]

  let centerSeq = uniques[localToGlobal[centerLocal]].sequence
  let obs = uniques[localToGlobal[obsLocal]]
  if centerSeq.len != obs.sequence.len:
    cache[key] = 0.0
    return 0.0

  let qFixed = max(err.qMin, min(err.qMax, 40))
  let qFixedIdx = qFixed - err.qMin

  var logLambda = 0.0
  for pos in 0..<centerSeq.len:
    let fromIdx = baseIndex(centerSeq[pos])
    let toIdx = baseIndex(obs.sequence[pos])
    if fromIdx < 0 or toIdx < 0:
      cache[key] = 0.0
      return 0.0

    let qIdx =
      if useQuality:
        qBinIndex(obs.meanQual[pos], err.qMin, err.qMax)
      else:
        qFixedIdx

    let p = err.probs[fromIdx * 4 + toIdx][qIdx]
    if p <= 0.0:
      cache[key] = 0.0
      return 0.0
    logLambda += ln(p)

  let lam =
    if logLambda < -745.0:
      0.0
    else:
      exp(logLambda)

  cache[key] = lam
  lam

proc chooseBestCluster(
  localIdx: int,
  centers: seq[int],
  totals: seq[int64],
  localToGlobal: seq[int],
  uniques: seq[DerepUnique],
  err: LearnErrorsResult,
  opts: DadaOptions,
  cache: var Table[uint64, float64]
): int =
  var bestCluster = 0
  var bestExpected = -1.0
  var haveBest = false

  for cIdx in 0..<centers.len:
    let centerLocal = centers[cIdx]
    let lam = lambdaCenterToObs(centerLocal, localIdx, localToGlobal, uniques, err, opts.useQuality, cache)
    let expected = lam * totals[cIdx].float64

    if (not haveBest) or expected > bestExpected or
       (expected == bestExpected and betterCenterLocal(centerLocal, centers[bestCluster], localToGlobal, uniques)):
      bestCluster = cIdx
      bestExpected = expected
      haveBest = true

  if haveBest and bestExpected > 0.0:
    return bestCluster

  # Deterministic fallback when every expectation is exactly zero.
  for cIdx in 0..<centers.len:
    if totals[cIdx] > totals[bestCluster] or
       (totals[cIdx] == totals[bestCluster] and
         betterCenterLocal(centers[cIdx], centers[bestCluster], localToGlobal, uniques)):
      bestCluster = cIdx
  bestCluster

proc denoiseLengthGroup(
  localToGlobal: seq[int],
  uniques: seq[DerepUnique],
  err: LearnErrorsResult,
  opts: DadaOptions
): GroupDenoiseResult =
  result.localToGlobal = localToGlobal
  if localToGlobal.len == 0:
    return

  let n = localToGlobal.len
  var centers = @[0]
  var centerToCluster = newSeq[int](n)
  for i in 0..<n:
    centerToCluster[i] = -1
  centerToCluster[0] = 0

  var assignments = newSeq[int](n)
  for i in 0..<n:
    assignments[i] = 0

  var cache = initTable[uint64, float64]()

  for outer in 1..opts.maxIterations:
    result.iterationsRun = outer

    for _ in 1..opts.maxShuffle:
      let totals = clusterTotals(assignments, centers.len, localToGlobal, uniques)
      var changed = false

      for localIdx in 0..<n:
        if centerToCluster[localIdx] >= 0:
          let forced = centerToCluster[localIdx]
          if assignments[localIdx] != forced:
            assignments[localIdx] = forced
            changed = true
          continue

        let best = chooseBestCluster(localIdx, centers, totals, localToGlobal, uniques, err, opts, cache)
        if assignments[localIdx] != best:
          assignments[localIdx] = best
          changed = true

      if not changed:
        break
    for cIdx, centerLocal in centers:
      assignments[centerLocal] = cIdx

    let totals = clusterTotals(assignments, centers.len, localToGlobal, uniques)
    var bestLocal = -1
    var bestP = 1.0

    for localIdx in 0..<n:
      if centerToCluster[localIdx] >= 0:
        continue

      let obs = uniques[localToGlobal[localIdx]]
      if obs.abundance < opts.minAbundance:
        continue

      let cIdx = assignments[localIdx]
      let centerLocal = centers[cIdx]
      let center = uniques[localToGlobal[centerLocal]]
      let hd = hammingDistance(center.sequence, obs.sequence)
      if hd < opts.minHamming:
        continue

      # Sum expected abundance over ALL current clusters (R DADA2 semantics):
      # lambda = sum_j (cluster_total_j × P(obs | center_j))
      var expected = 0.0
      for cj in 0..<centers.len:
        let lamJ = lambdaCenterToObs(centers[cj], localIdx, localToGlobal, uniques, err, opts.useQuality, cache)
        expected += lamJ * totals[cj].float64

      if expected > 0.0 and obs.abundance.float64 < opts.minFold * expected:
        continue

      let p = abundancePValue(obs.abundance, expected, prior = false)
      let pCorr = min(1.0, p * n.float64)

      if pCorr < bestP:
        bestP = pCorr
        bestLocal = localIdx
      elif pCorr == bestP and bestLocal >= 0:
        let aGlobal = localToGlobal[localIdx]
        let bGlobal = localToGlobal[bestLocal]
        if cmpGlobalUnique(aGlobal, bGlobal, uniques) < 0:
          bestLocal = localIdx

    if bestLocal >= 0 and bestP < opts.omegaA and
       (opts.maxClusters == 0 or centers.len < opts.maxClusters):
      let newCluster = centers.len
      centers.add(bestLocal)
      centerToCluster[bestLocal] = newCluster
      assignments[bestLocal] = newCluster
      inc result.splitsAccepted
      continue

    break

  result.centers = centers
  result.assignments = assignments
  result.clusterTotals = clusterTotals(assignments, centers.len, localToGlobal, uniques)

proc denoiseGroupRange(
  groupLists: ptr UncheckedArray[seq[int]],
  outResults: ptr UncheckedArray[GroupDenoiseResult],
  n: int,
  uniques: seq[DerepUnique],
  err: LearnErrorsResult,
  opts: DadaOptions
) =
  for i in 0..<n:
    outResults[i] = denoiseLengthGroup(groupLists[i], uniques, err, opts)

template atPtr[T](data: var seq[T]; i: int): ptr UncheckedArray[T] =
  cast[ptr UncheckedArray[T]](unsafeAddr data[i])

proc writeAsvTsv*(path: string, res: DadaResult) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open ASV output file: " & path)
  defer: close(f)

  f.writeLine("asv_id\tsequence\tabundance\tcluster_size\tcenter_unique_index\tfirst_seen_read_index")
  for i, asv in res.asvs:
    f.writeLine(&"ASV{i + 1}\t{asv.sequence}\t{asv.abundance}\t{asv.clusterSize}\t{asv.centerUniqueIndex}\t{asv.firstSeenReadIndex}")

proc initLearnErrorMatrix(learnOpts: LearnErrorsOptions): LearnErrorsResult =
  result.qMin = learnOpts.qMin
  result.qMax = learnOpts.qMax
  let qBins = result.qMax - result.qMin + 1
  for i in 0..15:
    result.counts[i] = newSeq[int64](qBins)
    result.probs[i] = newSeq[float64](qBins)

proc fillLearnErrorProbabilities(res: var LearnErrorsResult, learnOpts: LearnErrorsOptions) =
  let qBins = res.qMax - res.qMin + 1
  for fromIdx in 0..3:
    for qIdx in 0..<qBins:
      var totalAtQ = 0'i64
      for toIdx in 0..3:
        totalAtQ += res.counts[fromIdx * 4 + toIdx][qIdx]

      if totalAtQ == 0:
        for toIdx in 0..3:
          if toIdx == fromIdx:
            res.probs[fromIdx * 4 + toIdx][qIdx] = 1.0 - 3.0 * learnOpts.zeroTotalMismatchProb
          else:
            res.probs[fromIdx * 4 + toIdx][qIdx] = learnOpts.zeroTotalMismatchProb
      else:
        let denom = totalAtQ.float64 + 4.0 * learnOpts.pseudocount
        for toIdx in 0..3:
          res.probs[fromIdx * 4 + toIdx][qIdx] =
            (res.counts[fromIdx * 4 + toIdx][qIdx].float64 + learnOpts.pseudocount) / denom

proc estimateErrorMatrixFromDada*(
  derepRes: DerepFastqResult,
  dadaRes: DadaResult,
  learnOpts: LearnErrorsOptions = defaultLearnErrorsOptions()
): LearnErrorsResult =
  validateOptions(learnOpts)

  if dadaRes.uniqueToCenterUnique.len != derepRes.uniques.len:
    raise newException(ValueError, "Dada result unique mapping size mismatch")

  result = initLearnErrorMatrix(learnOpts)

  var centerLengths = initHashSet[int]()
  for uniqIdx, uniq in derepRes.uniques:
    var centerIdx = dadaRes.uniqueToCenterUnique[uniqIdx]
    if centerIdx < 0 or centerIdx >= derepRes.uniques.len:
      centerIdx = uniqIdx

    let center = derepRes.uniques[centerIdx]
    if center.sequence.len != uniq.sequence.len:
      result.skippedPositions += uniq.abundance * uniq.sequence.len.int64
      continue

    centerLengths.incl(center.sequence.len)

    for pos in 0..<uniq.sequence.len:
      let fromIdx = baseIndex(center.sequence[pos])
      let toIdx = baseIndex(uniq.sequence[pos])
      if fromIdx < 0 or toIdx < 0:
        result.skippedPositions += uniq.abundance
        continue

      let qIdx = qBinIndex(uniq.meanQual[pos], result.qMin, result.qMax)
      result.counts[fromIdx * 4 + toIdx][qIdx] += uniq.abundance
      result.totalTransitions += uniq.abundance

  result.centerCountByLength = centerLengths.len
  fillLearnErrorProbabilities(result, learnOpts)

proc maxProbDelta(a, b: LearnErrorsResult): float64 =
  if a.qMin != b.qMin or a.qMax != b.qMax:
    raise newException(ValueError, "Cannot compare error matrices with different quality ranges")
  for i in 0..15:
    for qIdx in 0..<a.probs[i].len:
      result = max(result, abs(a.probs[i][qIdx] - b.probs[i][qIdx]))

proc assignmentsChanged(prev, curr: seq[int]): bool =
  if prev.len == 0:
    return false
  if prev.len != curr.len:
    return true
  for i in 0..<prev.len:
    if prev[i] != curr[i]:
      return true
  false

proc dadaDenoise*(
  derepRes: DerepFastqResult,
  err: LearnErrorsResult,
  opts: DadaOptions = defaultDadaOptions()
): DadaResult =
  validateOptions(opts)
  if derepRes.uniques.len == 0:
    return DadaResult(
      totalReads: 0,
      inputUniqueCount: 0,
      asvCount: 0,
      asvs: @[],
      uniqueToAsv: @[],
      uniqueToCenterUnique: @[],
      groupCount: 0,
      iterationsRun: 0,
      splitsAccepted: 0,
      correctedUniques: 0,
      uncorrectedUniques: 0
    )

  var byLength = initTable[int, seq[int]]()
  for i, u in derepRes.uniques:
    if u.sequence.len != u.meanQual.len:
      raise newException(ValueError, "Derep unique has meanQual length mismatch")
    if byLength.hasKey(u.sequence.len):
      var bucket = byLength[u.sequence.len]
      bucket.add(i)
      byLength[u.sequence.len] = bucket
    else:
      byLength[u.sequence.len] = @[i]

  var lengths: seq[int] = @[]
  for seqLen, _ in byLength:
    lengths.add(seqLen)
  lengths.sort(system.cmp[int])

  var groups: seq[seq[int]] = @[]
  for seqLen in lengths:
    var idxs = byLength[seqLen]
    idxs.sort(proc(a, b: int): int = cmpGlobalUnique(a, b, derepRes.uniques))
    groups.add(idxs)

  var groupResults = newSeq[GroupDenoiseResult](groups.len)
  if opts.threads == 1 or groups.len < opts.parallelMinGroups:
    for i in 0..<groups.len:
      groupResults[i] = denoiseLengthGroup(groups[i], derepRes.uniques, err, opts)
  else:
    let workerCount = min(opts.threads, groups.len)
    let blockSize = (groups.len + workerCount - 1) div workerCount

    var m = createMaster()
    m.awaitAll:
      var idx = 0
      while idx < groups.len:
        let n = min(blockSize, groups.len - idx)
        m.spawn denoiseGroupRange(
          atPtr(groups, idx),
          atPtr(groupResults, idx),
          n,
          derepRes.uniques,
          err,
          opts
        )
        idx += n

  result.totalReads = derepRes.totalReads
  result.inputUniqueCount = derepRes.uniqueCount
  result.groupCount = groups.len
  result.uniqueToAsv = newSeq[int](derepRes.uniques.len)
  result.uniqueToCenterUnique = newSeq[int](derepRes.uniques.len)
  for i in 0..<result.uniqueToAsv.len:
    result.uniqueToAsv[i] = -1
    result.uniqueToCenterUnique[i] = -1

  type TmpAsv = object
    sequence: string
    abundance: int64
    clusterSize: int
    centerUniqueIndex: int
    firstSeenReadIndex: int64
    members: seq[int]

  var tmpAsvs: seq[TmpAsv] = @[]

  for gRes in groupResults:
    result.iterationsRun += gRes.iterationsRun
    result.splitsAccepted += gRes.splitsAccepted

    var correctionCache = initTable[uint64, float64]()
    var correctedMembersByCluster = newSeq[seq[int]](gRes.centers.len)
    var uncorrectedMembers: seq[int] = @[]

    for localIdx, cIdx in gRes.assignments:
      let centerLocal = gRes.centers[cIdx]
      if localIdx == centerLocal:
        correctedMembersByCluster[cIdx].add(localIdx)
        continue

      let globalIdx = gRes.localToGlobal[localIdx]
      let obs = derepRes.uniques[globalIdx]

      # Sum expected abundance over ALL clusters (R DADA2 semantics)
      var expected = 0.0
      for cj in 0..<gRes.centers.len:
        let lamJ = lambdaCenterToObs(
          centerLocal = gRes.centers[cj],
          obsLocal = localIdx,
          localToGlobal = gRes.localToGlobal,
          uniques = derepRes.uniques,
          err = err,
          useQuality = opts.useQuality,
          cache = correctionCache
        )
        expected += lamJ * gRes.clusterTotals[cj].float64
      let pCorr = abundancePValue(obs.abundance, expected, prior = true)

      if pCorr >= opts.omegaC:
        correctedMembersByCluster[cIdx].add(localIdx)
      else:
        uncorrectedMembers.add(localIdx)

    var clusterToTmpAsv = newSeq[int](gRes.centers.len)
    for cIdx, centerLocal in gRes.centers:
      let centerGlobal = gRes.localToGlobal[centerLocal]
      var abundance = 0'i64
      for localIdx in correctedMembersByCluster[cIdx]:
        let globalIdx = gRes.localToGlobal[localIdx]
        abundance += derepRes.uniques[globalIdx].abundance

      clusterToTmpAsv[cIdx] = tmpAsvs.len
      tmpAsvs.add(TmpAsv(
        sequence: derepRes.uniques[centerGlobal].sequence,
        abundance: abundance,
        clusterSize: correctedMembersByCluster[cIdx].len,
        centerUniqueIndex: centerGlobal,
        firstSeenReadIndex: derepRes.uniques[centerGlobal].firstSeenReadIndex,
        members: @[]
      ))

    for cIdx in 0..<gRes.centers.len:
      let tmpIdx = clusterToTmpAsv[cIdx]
      for localIdx in correctedMembersByCluster[cIdx]:
        let globalIdx = gRes.localToGlobal[localIdx]
        tmpAsvs[tmpIdx].members.add(globalIdx)
        result.uniqueToCenterUnique[globalIdx] = tmpAsvs[tmpIdx].centerUniqueIndex
        inc result.correctedUniques

    for localIdx in uncorrectedMembers:
      let globalIdx = gRes.localToGlobal[localIdx]
      tmpAsvs.add(TmpAsv(
        sequence: derepRes.uniques[globalIdx].sequence,
        abundance: derepRes.uniques[globalIdx].abundance,
        clusterSize: 1,
        centerUniqueIndex: globalIdx,
        firstSeenReadIndex: derepRes.uniques[globalIdx].firstSeenReadIndex,
        members: @[globalIdx]
      ))
      result.uniqueToCenterUnique[globalIdx] = globalIdx
      inc result.uncorrectedUniques

  var ordering = newSeq[int](tmpAsvs.len)
  for i in 0..<tmpAsvs.len:
    ordering[i] = i
  ordering.sort(proc(a, b: int): int =
    if tmpAsvs[a].abundance != tmpAsvs[b].abundance:
      return cmp(tmpAsvs[b].abundance, tmpAsvs[a].abundance)
    if tmpAsvs[a].firstSeenReadIndex != tmpAsvs[b].firstSeenReadIndex:
      return cmp(tmpAsvs[a].firstSeenReadIndex, tmpAsvs[b].firstSeenReadIndex)
    cmp(tmpAsvs[a].sequence, tmpAsvs[b].sequence)
  )

  var oldToNew = newSeq[int](tmpAsvs.len)
  result.asvs = newSeq[AsvRecord](tmpAsvs.len)
  for newIdx, oldIdx in ordering:
    oldToNew[oldIdx] = newIdx
    let t = tmpAsvs[oldIdx]
    result.asvs[newIdx] = AsvRecord(
      sequence: t.sequence,
      abundance: t.abundance,
      clusterSize: t.clusterSize,
      centerUniqueIndex: t.centerUniqueIndex,
      firstSeenReadIndex: t.firstSeenReadIndex
    )

  for oldIdx, t in tmpAsvs:
    let asvIdx = oldToNew[oldIdx]
    for globalIdx in t.members:
      result.uniqueToAsv[globalIdx] = asvIdx

  result.asvCount = result.asvs.len

proc dadaSelfConsistentFromInitialErrors*(
  derepRes: DerepFastqResult,
  initialErr: LearnErrorsResult,
  learnOpts: LearnErrorsOptions = defaultLearnErrorsOptions(),
  dadaOpts: DadaOptions = defaultDadaOptions(),
  scOpts: DadaSelfConsistOptions = defaultDadaSelfConsistOptions()
): DadaSelfConsistResult =
  validateOptions(learnOpts)
  validateOptions(dadaOpts)
  validateOptions(scOpts)

  if not scOpts.enabled:
    result.errorMatrix = initialErr
    result.dada = dadaDenoise(derepRes, result.errorMatrix, dadaOpts)
    result.iterationsRun = 1
    result.converged = true
    result.maxAbsProbDelta = 0.0
    result.assignmentChanges = 0
    return

  var iterDadaOpts = dadaOpts
  if scOpts.forceOmegaCZero:
    iterDadaOpts.omegaC = 0.0

  var currentErr = initialErr
  var prevAssignments: seq[int] = @[]
  var lastDada = DadaResult()
  var assignmentChangeCount = 0

  for iter in 1..scOpts.maxIterations:
    let currDada = dadaDenoise(derepRes, currentErr, iterDadaOpts)
    let newErr = estimateErrorMatrixFromDada(derepRes, currDada, learnOpts)
    let delta = maxProbDelta(currentErr, newErr)
    let changed = assignmentsChanged(prevAssignments, currDada.uniqueToCenterUnique)
    if changed:
      inc assignmentChangeCount

    result.iterationsRun = iter
    result.maxAbsProbDelta = delta
    lastDada = currDada

    currentErr = newErr
    prevAssignments = currDada.uniqueToCenterUnique

    if iter >= scOpts.minIterations and (not changed) and delta <= scOpts.tol:
      result.converged = true
      break

  result.assignmentChanges = assignmentChangeCount
  result.errorMatrix = currentErr

  if scOpts.finalDenoiseWithLatestErrors:
    result.dada = dadaDenoise(derepRes, result.errorMatrix, dadaOpts)
  else:
    result.dada = lastDada

proc dadaSelfConsistent*(
  derepRes: DerepFastqResult,
  learnOpts: LearnErrorsOptions = defaultLearnErrorsOptions(),
  dadaOpts: DadaOptions = defaultDadaOptions(),
  scOpts: DadaSelfConsistOptions = defaultDadaSelfConsistOptions()
): DadaSelfConsistResult =
  let initialErr = learnErrorsFromDerep(derepRes, learnOpts)
  dadaSelfConsistentFromInitialErrors(derepRes, initialErr, learnOpts, dadaOpts, scOpts)
