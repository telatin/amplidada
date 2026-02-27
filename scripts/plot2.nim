## plot.nim — DADA2-style error rate plot using ggplotnim
## Build:
##   arch -arm64 nim c -d:release \
##     --passL:"-L$(brew --prefix cairo)/lib" \
##     --passL:"-Wl,-rpath,$(brew --prefix cairo)/lib" \
##     -r scripts/plot.nim

import ggplotnim
import std/[math, random]

const
  Bases = ["A", "C", "G", "T"]
  MaxQ  = 40
  # Per off-diagonal substitution scaling factors (12 pairs, row-major excluding diagonal)
  BaseRates = [
    0.08, 0.12, 0.05,   # A2C, A2G, A2T
    0.06, 0.04, 0.09,   # C2A, C2G, C2T
    0.07, 0.03, 0.05,   # G2A, G2C, G2T
    0.06, 0.10, 0.04    # T2A, T2C, T2G
  ]

proc kernelSmooth(xs, ys: seq[float]; bw = 6.0): seq[float] =
  result = newSeq[float](xs.len)
  for i in 0 ..< xs.len:
    var wSum, vSum: float
    for j in 0 ..< ys.len:
      let w = exp(-0.5 * ((xs[i] - xs[j]) / bw) ^ 2)
      wSum += w
      vSum += w * log10(max(ys[j], 1e-10))
    result[i] = pow(10.0, vSum / wSum)

proc main() =
  randomize(42)

  var
    obsQ, obsErr, obsSize : seq[float]
    obsSub                : seq[string]
    smQ, smErr            : seq[float]
    smSub                 : seq[string]
    refQ, refErr          : seq[float]
    refSub                : seq[string]

  var rateIdx = 0

  for fb in Bases:
    for tb in Bases:
      let sub    = fb & "2" & tb
      let isDiag = (fb == tb)

      let scaleFactor =
        if isDiag: 1.0
        else:
          let s = BaseRates[rateIdx mod 12]
          inc rateIdx
          s

      var rawQ, rawErr: seq[float]

      for q in 1 .. MaxQ:
        let qf       = float(q)
        let phred    = pow(10.0, -qf / 10.0)
        let trueRate = phred / 3.0 * (if isDiag: 1.0 else: scaleFactor * 3.0 / 0.07)

        # Observation count: peaks ~Q18, log10-scaled so geom_point sizes are ~1.7–3.5
        let nObs   = max(50.0, 3000.0 * exp(-(qf - 18.0)^2 / 200.0) + rand(200.0))
        let ptSize = log10(nObs)

        # Observed rate: add log-normal noise; variance grows with quality (fewer reads)
        let sd      = 0.10 + 0.30 * sqrt(qf / float(MaxQ))
        let obsRate = clamp(trueRate * exp(gauss(0.0, sd)), 1e-6, 0.5)

        obsQ.add qf;       obsErr.add obsRate;  obsSize.add ptSize
        obsSub.add sub

        rawQ.add qf;       rawErr.add obsRate
        refQ.add qf;       refErr.add trueRate
        refSub.add sub

      # Smooth observed rates in log-space via Gaussian kernel
      let sm = kernelSmooth(rawQ, rawErr, bw = 6.0)
      for i in 0 ..< rawQ.len:
        smQ.add rawQ[i];   smErr.add sm[i]
        smSub.add sub

  # ── Assemble DataFrames ───────────────────────────────────────────
  let dfObs = toDf({
    "Quality":      obsQ,
    "ErrorFreq":    obsErr,
    "PtSize":       obsSize,
    "Substitution": obsSub
  })
  let dfSm = toDf({
    "Quality":      smQ,
    "ErrorFreq":    smErr,
    "Substitution": smSub
  })
  let dfRef = toDf({
    "Quality":      refQ,
    "ErrorFreq":    refErr,
    "Substitution": refSub
  })

  ggplot(dfObs, aes("Quality", "ErrorFreq")) +
    geom_point(aes(size = "PtSize"),
               color = some(parseHex("606060")),
               alpha = some(0.55)) +
    geom_line(data = dfSm,
              aes = aes("Quality", "ErrorFreq"),
              color = some(parseHex("111111")),
              size  = some(1.4)) +
    geom_line(data = dfRef,
              aes = aes("Quality", "ErrorFreq"),
              color = some(parseHex("DD2222")),
              size  = some(0.9)) +
    facet_wrap("Substitution") +        # back to plain column name
    scale_y_log10() +
    xlim(0.0, 41.0) +
    xlab("Consensus quality score", margin = 0.5) +
    ylab("Error frequency (log10)", margin = 1.5) +
    hideLegend() +
    theme_opaque() +
    theme_font_scale(1.4) +
    Theme(facetHeaderFont: some(font(20.0, bold = true))) +
    ggsave("dada2_error_rates.png", width = 3200, height = 3000)

when isMainModule:
  main()
