## Simple LOESS (locally estimated scatterplot smoothing) implementation.
##
## This module provides a minimal LOESS smoother for 1-D data, supporting
## degree=1 (linear) and degree=2 (quadratic) local fits. It is designed
## to approximate R's `loess()` function for the specific use case of
## DADA2 error model estimation.
##
## The implementation uses:
## - Tricube weight function
## - Linear/quadratic local regression
## - Flat extrapolation beyond data range

import std/[algorithm, math, sequtils]

const
  DefaultSpan* = 0.75
  DefaultMinSpan* = 0.01

type
  LoessOptions* = object
    span*: float64
    degree*: int
    surface*: string  ## "direct" or "interpolate"

proc defaultLoessOptions*(): LoessOptions =
  LoessOptions(span: DefaultSpan, degree: 1, surface: "direct")

proc tricubeWeight(d: float64): float64 =
  ## Tricube weight function: (1 - |d|^3)^3 for |d| < 1, else 0
  let ad = abs(d)
  if ad >= 1.0:
    return 0.0
  let t = 1.0 - ad * ad * ad
  t * t * t

proc loessFit1D(
  x: openArray[float64],
  y: openArray[float64],
  weights: openArray[float64],
  xq: float64,
  span: float64,
  degree: int
): float64 =
  ## Fit a local polynomial of given degree at query point xq.
  ## Returns the predicted value, or NaN if fit fails.
  let n = x.len
  if n == 0:
    return NaN
  if n == 1:
    return y[0]

  # Compute distances and find neighborhood
  var distances = newSeq[float64](n)
  for i in 0..<n:
    distances[i] = abs(x[i] - xq)

  # Sort distances to find span neighborhood
  var sortedDist = distances.sorted
  let nNeighborhood = max(2, int(floor(span * float64(n))))
  let h = sortedDist[min(nNeighborhood - 1, n - 1)]

  if h == 0.0:
    # Exact match - return weighted average of exact matches
    var sumW = 0.0
    var sumWY = 0.0
    for i in 0..<n:
      if distances[i] == 0.0:
        let w = weights[i]
        sumW += w
        sumWY += w * y[i]
    if sumW > 0.0:
      return sumWY / sumW
    return y[0]

  # Compute tricube weights
  var localWeights = newSeq[float64](n)
  var nEffective = 0
  for i in 0..<n:
    let d = distances[i] / h
    localWeights[i] = tricubeWeight(d) * weights[i]
    if localWeights[i] > 1e-12:
      inc nEffective

  if nEffective == 0:
    return NaN

  # Fit local polynomial
  case degree
  of 0:
    # Local constant (weighted average)
    var sumW = 0.0
    var sumWY = 0.0
    for i in 0..<n:
      let w = localWeights[i]
      sumW += w
      sumWY += w * y[i]
    if sumW > 0.0:
      return sumWY / sumW
    return NaN
  of 1:
    # Local linear: y = a + b*(x - xq)
    var sumW = 0.0
    var sumWX = 0.0
    var sumWXX = 0.0
    var sumWY = 0.0
    var sumWXY = 0.0

    for i in 0..<n:
      let w = localWeights[i]
      if w < 1e-12: continue
      let dx = x[i] - xq
      sumW += w
      sumWX += w * dx
      sumWXX += w * dx * dx
      sumWY += w * y[i]
      sumWXY += w * dx * y[i]

    let denom = sumW * sumWXX - sumWX * sumWX
    if abs(denom) < 1e-12:
      # Fall back to weighted average
      if sumW > 0.0:
        return sumWY / sumW
      return NaN

    let a = (sumWXX * sumWY - sumWX * sumWXY) / denom
    return a
  of 2:
    # Local quadratic: y = a + b*(x - xq) + c*(x - xq)^2
    var sumW = 0.0
    var sumWX = 0.0
    var sumWXX = 0.0
    var sumWXXX = 0.0
    var sumWXXXX = 0.0
    var sumWY = 0.0
    var sumWXY = 0.0
    var sumWXXY = 0.0

    for i in 0..<n:
      let w = localWeights[i]
      if w < 1e-12: continue
      let dx = x[i] - xq
      let dx2 = dx * dx
      sumW += w
      sumWX += w * dx
      sumWXX += w * dx2
      sumWXXX += w * dx * dx2
      sumWXXXX += w * dx2 * dx2
      sumWY += w * y[i]
      sumWXY += w * dx * y[i]
      sumWXXY += w * dx2 * y[i]

    # Solve 3x3 system using Cramer's rule
    # [sumW    sumWX   sumWXX  ] [a]   [sumWY  ]
    # [sumWX   sumWXX  sumWXXX ] [b] = [sumWXY ]
    # [sumWXX  sumWXXX sumWXXXX] [c]   [sumWXXY]
    let a00 = sumW;   let a01 = sumWX;   let a02 = sumWXX
    let a10 = sumWX;  let a11 = sumWXX;  let a12 = sumWXXX
    let a20 = sumWXX; let a21 = sumWXXX; let a22 = sumWXXXX

    let det = a00 * (a11 * a22 - a12 * a21) -
              a01 * (a10 * a22 - a12 * a20) +
              a02 * (a10 * a21 - a11 * a20)

    if abs(det) < 1e-12:
      # Fall back to linear fit
      let denom2 = sumW * sumWXX - sumWX * sumWX
      if abs(denom2) < 1e-12:
        if sumW > 0.0:
          return sumWY / sumW
        return NaN
      return (sumWXX * sumWY - sumWX * sumWXY) / denom2

    let b0 = sumWY; let b1 = sumWXY; let b2 = sumWXXY

    let detA = b0 * (a11 * a22 - a12 * a21) -
               a01 * (b1 * a22 - a12 * b2) +
               a02 * (b1 * a21 - a11 * b2)

    return detA / det
  else:
    return NaN

proc loessSmooth*(
  x: openArray[float64],
  y: openArray[float64],
  weights: openArray[float64],
  xq: openArray[float64],
  span: float64 = DefaultSpan,
  degree: int = 1
): seq[float64] =
  ## Smooth y values at query points xq using LOESS.
  ##
  ## Parameters:
  ##   x: x-coordinates of data points
  ##   y: y-coordinates of data points
  ##   weights: weights for each data point
  ##   xq: x-coordinates where to predict
  ##   span: smoothing span (fraction of data in each neighborhood)
  ##   degree: degree of local polynomial (0, 1, or 2)
  ##
  ## Returns predicted values at xq. NaN indicates fit failure.

  result = newSeq[float64](xq.len)

  if x.len == 0 or xq.len == 0:
    return

  if weights.len == 0:
    # Default: uniform weights
    var uniformWeights = newSeq[float64](x.len)
    for i in 0..<x.len:
      uniformWeights[i] = 1.0
    for i in 0..<xq.len:
      result[i] = loessFit1D(x, y, uniformWeights, xq[i], span, degree)
  else:
    for i in 0..<xq.len:
      result[i] = loessFit1D(x, y, weights, xq[i], span, degree)

proc loessSmoothWithFallback*(
  x: openArray[float64],
  y: openArray[float64],
  weights: openArray[float64],
  xq: openArray[float64],
  span: float64 = DefaultSpan,
  degree: int = 1
): seq[float64] =
  ## LOESS smoothing with fallback to lower degree if fit fails.
  ##
  ## If degree=2 fitting fails, falls back to degree=1.
  ## If degree=1 fitting fails, falls back to degree=0 (local mean).

  result = loessSmooth(x, y, weights, xq, span, degree)

  # Check for NaN values and try lower degree
  if degree > 0:
    var hasNaN = false
    for v in result:
      if v != v:  # NaN check
        hasNaN = true
        break

    if hasNaN:
      let fallback = loessSmoothWithFallback(x, y, weights, xq, span, degree - 1)
      for i in 0..<xq.len:
        if result[i] != result[i]:  # NaN
          result[i] = fallback[i]

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

when isMainModule:
  import std/unittest

  suite "LOESS":
    test "constant function":
      let x = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let y = @[2.0, 2.0, 2.0, 2.0, 2.0]
      let w = @[1.0, 1.0, 1.0, 1.0, 1.0]
      let xq = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let pred = loessSmooth(x, y, w, xq, 1.0, 0)
      for p in pred:
        check abs(p - 2.0) < 0.01

    test "linear function":
      let x = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let y = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let w = @[1.0, 1.0, 1.0, 1.0, 1.0]
      let xq = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let pred = loessSmooth(x, y, w, xq, 1.0, 1)
      for i, p in pred:
        check abs(p - xq[i]) < 0.01

    test "extrapolation":
      let x = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let y = @[1.0, 2.0, 3.0, 4.0, 5.0]
      let w = @[1.0, 1.0, 1.0, 1.0, 1.0]
      let xq = @[0.0, 6.0]  # Outside range
      let pred = loessSmooth(x, y, w, xq, 1.0, 1)
      # Should extrapolate reasonably
      check pred[0] < 1.0
      check pred[1] > 5.0
