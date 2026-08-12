## loessErrfun_mod4: Error estimation function for binned quality scores.
##
## This implements the `loessErrfun_mod4` function from R DADA2, designed
## for NextSeq 2000 / NovaSeq binned quality scores.
##
## Key differences from standard loessErrfun:
## - Uses linear LOESS (degree=1) with span=0.95
## - Uses log10(tot) as weights instead of tot
## - Applies X40 floor: no fitted error probability can fall below its
##   substitution-specific Q40 estimate
## - MAX_ERROR_RATE = 0.25, MIN_ERROR_RATE = 1e-7

import std/[algorithm, math, sequtils, strformat]
import amplidada/[learn_errors, loess]

const
  Mod4MaxErrorRate* = 0.25
  Mod4MinErrorRate* = 1e-7
  Mod4Span* = 0.95
  Mod4Degree* = 1

proc loessErrfunMod4*(counts: array[16, seq[int64]], qMin, qMax: int): seq[seq[float64]] =
  ## Estimate error probabilities using loessErrfun_mod4.
  ##
  ## Parameters:
  ##   counts: 16 transition count arrays (one per from→to combination)
  ##   qMin, qMax: quality score range
  ##
  ## Returns: 16×(qMax-qMin+1) matrix of error probabilities.
  ## Row i corresponds to transition index i (from*4+to).
  ##
  ## This mirrors R DADA2's loessErrfun_mod4 function:
  ##   - Fits log10((errs + 1) / tot) ~ q using LOESS
  ##   - Weights = log10(tot)
  ##   - Span = 0.95, degree = 1
  ##   - Flat extrapolation beyond data range
  ##   - X40 floor applied to all quality bins
  ##   - MAX_ERROR_RATE and MIN_ERROR_RATE clamping

  let qBins = qMax - qMin + 1
  var qq = newSeq[float64](qBins)
  for i in 0..<qBins:
    qq[i] = float64(qMin + i)

  # est matrix: 12 rows (non-self transitions) × qBins
  # Transition order: A→C, A→G, A→T, C→A, C→G, C→T, G→A, G→C, G→T, T→A, T→C, T→G
  const transitionOrder = [
    (0, 1), (0, 2), (0, 3),  # A→C, A→G, A→T
    (1, 0), (1, 2), (1, 3),  # C→A, C→G, C→T
    (2, 0), (2, 1), (2, 3),  # G→A, G→C, G→T
    (3, 0), (3, 1), (3, 2)   # T→A, T→C, T→G
  ]

  var est = newSeq[seq[float64]](12)
  for i in 0..<12:
    est[i] = newSeq[float64](qBins)

  # For X40 calculation: store Q40 estimate for each transition
  var x40Estimates = newSeq[float64](12)

  for ti in 0..<12:
    let (fromIdx, toIdx) = transitionOrder[ti]
    let countIdx = fromIdx * 4 + toIdx

    # Get error counts and totals for this transition
    var errs = newSeq[int64](qBins)
    var tots = newSeq[int64](qBins)

    for qi in 0..<qBins:
      errs[qi] = counts[countIdx][qi]
      # Total = sum of all transitions from this base
      for toJ in 0..3:
        tots[qi] += counts[fromIdx * 4 + toJ][qi]

    # Compute regularized log mismatch rates: log10((errs + 1) / tot)
    var rlogp = newSeq[float64](qBins)
    var hasData = false
    for qi in 0..<qBins:
      if tots[qi] > 0:
        rlogp[qi] = log10((float64(errs[qi]) + 1.0) / float64(tots[qi]))
        hasData = true
      else:
        rlogp[qi] = NaN

    if not hasData:
      # No data for this transition - use fallback
      for qi in 0..<qBins:
        est[ti][qi] = Mod4MinErrorRate
      x40Estimates[ti] = Mod4MinErrorRate
      continue

    # Weights = log10(tot)
    var weights = newSeq[float64](qBins)
    for qi in 0..<qBins:
      if tots[qi] > 0:
        weights[qi] = log10(float64(tots[qi]))
      else:
        weights[qi] = 0.0

    # LOESS fit: rlogp ~ qq with weights
    let predictions = loessSmoothWithFallback(qq, rlogp, weights, qq, Mod4Span, Mod4Degree)

    # Flat extrapolation: extend edge values
    var pred = predictions
    var validStart = -1
    var validEnd = -1

    for qi in 0..<qBins:
      if pred[qi] == pred[qi]:  # Not NaN
        if validStart < 0:
          validStart = qi
        validEnd = qi

    if validStart >= 0:
      # Extend left edge value to the left
      for qi in 0..<validStart:
        pred[qi] = pred[validStart]
      # Extend right edge value to the right
      for qi in (validEnd + 1)..<qBins:
        pred[qi] = pred[validEnd]

    # Convert back from log10 scale and clamp
    var transitionProbs = newSeq[float64](qBins)
    for qi in 0..<qBins:
      if pred[qi] == pred[qi]:  # Not NaN
        transitionProbs[qi] = pow(10.0, pred[qi])
      else:
        transitionProbs[qi] = Mod4MinErrorRate

      # Clamp to [MIN_ERROR_RATE, MAX_ERROR_RATE]
      if transitionProbs[qi] > Mod4MaxErrorRate:
        transitionProbs[qi] = Mod4MaxErrorRate
      if transitionProbs[qi] < Mod4MinErrorRate:
        transitionProbs[qi] = Mod4MinErrorRate

    est[ti] = transitionProbs

    # Store Q40 estimate for X40 floor
    let q40Idx = qMax - qMin  # Index of Q40 in our array
    if q40Idx >= 0 and q40Idx < qBins:
      x40Estimates[ti] = transitionProbs[q40Idx]
    else:
      # Use last available quality bin
      x40Estimates[ti] = transitionProbs[^1]

  # Apply X40 floor: no fitted error probability can fall below its
  # substitution-specific Q40 estimate
  for ti in 0..<12:
    let x40 = x40Estimates[ti]
    for qi in 0..<qBins:
      if est[ti][qi] < x40:
        est[ti][qi] = x40

  # Build full 16×qBins probability matrix
  # Rows: A→A, A→C, A→G, A→T, C→A, C→C, C→G, C→T, ...
  result = newSeq[seq[float64]](16)
  for i in 0..<16:
    result[i] = newSeq[float64](qBins)

  var estRow = 0
  for fromIdx in 0..3:
    for toIdx in 0..3:
      let idx = fromIdx * 4 + toIdx
      if fromIdx == toIdx:
        # Diagonal: 1 - sum of other transitions from this base
        for qi in 0..<qBins:
          var sumOffDiag = 0.0
          for toJ in 0..3:
            if toJ != toIdx:
              sumOffDiag += result[fromIdx * 4 + toJ][qi]
          result[idx][qi] = 1.0 - sumOffDiag
          # Ensure diagonal is reasonable
          if result[idx][qi] < Mod4MinErrorRate:
            result[idx][qi] = 1.0 - 3.0 * Mod4MaxErrorRate
          if result[idx][qi] > 1.0:
            result[idx][qi] = 1.0 - 3.0 * Mod4MinErrorRate
      else:
        # Off-diagonal: use est matrix
        result[idx] = est[estRow]
        inc estRow

proc loessErrfunMod4FromResult*(res: LearnErrorsResult): seq[seq[float64]] =
  ## Convenience wrapper: apply loessErrfun_mod4 to a LearnErrorsResult.
  loessErrfunMod4(res.counts, res.qMin, res.qMax)

proc applyMod4ToResult*(res: var LearnErrorsResult) =
  ## Replace probabilities in LearnErrorsResult with loessErrfun_mod4 estimates.
  let newProbs = loessErrfunMod4(res.counts, res.qMin, res.qMax)
  for i in 0..15:
    res.probs[i] = newProbs[i]

proc describeMod4*(): string =
  &"loessErrfun_mod4: span={Mod4Span} degree={Mod4Degree} " &
  &"maxErrorRate={Mod4MaxErrorRate} minErrorRate={Mod4MinErrorRate}"

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

when isMainModule:
  import std/unittest

  suite "loessErrfun_mod4":
    test "basic estimation":
      var counts: array[16, seq[int64]]
      let qMin = 0
      let qMax = 40
      let qBins = qMax - qMin + 1

      # Initialize counts
      for i in 0..15:
        counts[i] = newSeq[int64](qBins)

      # Add some data: A→A transitions dominate at high Q
      for q in 20..40:
        counts[0][q] = 1000  # A→A
        counts[1][q] = 10    # A→C
        counts[2][q] = 5     # A→G
        counts[3][q] = 5     # A→T

      let probs = loessErrfunMod4(counts, qMin, qMax)

      # Check dimensions
      check probs.len == 16
      for p in probs:
        check p.len == qBins

      # A→A (diagonal) should be high at Q40
      let a2aQ40 = probs[0][40]
      check a2aQ40 > 0.9

      # A→C should be low at Q40
      let a2cQ40 = probs[1][40]
      check a2cQ40 < 0.1

    test "X40 floor":
      var counts: array[16, seq[int64]]
      let qMin = 0
      let qMax = 40
      let qBins = qMax - qMin + 1

      for i in 0..15:
        counts[i] = newSeq[int64](qBins)

      # Add data only at low quality
      for q in 0..10:
        counts[0][q] = 100
        counts[1][q] = 10

      let probs = loessErrfunMod4(counts, qMin, qMax)

      # Q40 estimate for A→C
      let x40 = probs[1][40]

      # All A→C probabilities should be >= X40
      for q in 0..40:
        check probs[1][q] >= x40 - 1e-10

    test "clamping":
      var counts: array[16, seq[int64]]
      let qMin = 0
      let qMax = 40
      let qBins = qMax - qMin + 1

      for i in 0..15:
        counts[i] = newSeq[int64](qBins)

      # Extreme case: all errors
      for q in 0..40:
        counts[0][q] = 0  # A→A
        counts[1][q] = 1000  # A→C
        counts[2][q] = 0
        counts[3][q] = 0

      let probs = loessErrfunMod4(counts, qMin, qMax)

      # A→C should be clamped to MAX_ERROR_RATE
      check probs[1][20] <= Mod4MaxErrorRate + 0.01
