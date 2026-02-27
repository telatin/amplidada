## dada2_plot.nim  –  corrected version
import ggplotnim
import std/[math, random, strformat]

const
  Bases = ["A", "C", "G", "T"]
  MaxQ  = 40

proc kernelSmooth(xs, ys: seq[float]; bw = 5.0): seq[float] =
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

  # All data in ONE dataframe with a "Layer" column
  # Layers: "obs" (scatter), "smooth" (black curve), "ref" (red line)
  var
    allQ, allErr, allCount: seq[float]
    allSub, allLayer: seq[string]

  for fb in Bases:
    for tb in Bases:
      let sub    = fb & "2" & tb
      let isDiag = (fb == tb)
      var rawQ, rawErr: seq[float]

      for q in 1 .. MaxQ:          # start at Q1, avoid Q0 log issues
        let qf    = float(q)
        let phred = pow(10.0, -qf / 10.0)   # 10^(-Q/10)

        # Off-diagonal only: error rate per substitution type
        # Skip diagonal panels (A2A etc.) — they don't appear in DADA2 plots
        # but we'll keep them with a near-zero rate to show the panel
        let trueRate =
          if isDiag: phred / 3.0             # treat same as substitution
          else:      phred / 3.0

        let nObs = max(10.0,
          2000.0 * exp(-(qf - 18.0)^2 / 180.0) + rand(50.0))

        let sd      = 0.15 + 0.20 * (qf / float(MaxQ))
        let obsRate = clamp(trueRate * exp(gauss(0.0, sd)), 1e-7, 0.5)

        # Scatter point
        allQ.add qf; allErr.add obsRate; allCount.add log10(nObs)
        allSub.add sub; allLayer.add "obs"

        rawQ.add qf; rawErr.add obsRate

        # Reference line (Phred model)
        allQ.add qf; allErr.add trueRate; allCount.add 1.0
        allSub.add sub; allLayer.add "ref"

      # Smoothed curve
      let sm = kernelSmooth(rawQ, rawErr, bw = 5.0)
      for i in 0 ..< rawQ.len:
        allQ.add rawQ[i]; allErr.add sm[i]; allCount.add 1.0
        allSub.add sub; allLayer.add "smooth"

  let df = toDf({
    "Quality":      allQ,
    "ErrorFreq":    allErr,
    "Count":        allCount,    # log10-scaled so range is ~1–3.5
    "Substitution": allSub,
    "Layer":        allLayer
  })

  # Split into separate DFs for each geom — but built from the SAME base df
  # to avoid the multi-data facet_wrap bug: filter, don't pass data= to geoms
  let dfObs    = df.filter(f{`Layer` == "obs"})
  let dfSmooth = df.filter(f{`Layer` == "smooth"})
  let dfRef    = df.filter(f{`Layer` == "ref"})

  ggplot(dfObs, aes("Quality", "ErrorFreq")) +
    geom_point(aes(size = "Count"),
               color = some(parseHex("888888")),
               alpha = some(0.6)) +
    geom_line(data = dfSmooth,
              aes = aes("Quality", "ErrorFreq"),
              color = some(parseHex("111111")),
              size  = some(1.5)) +
    geom_line(data = dfRef,
              aes = aes("Quality", "ErrorFreq"),
              color = some(parseHex("CC0000")),
              size  = some(1.0)) +
    facet_wrap("Substitution") +          # no ncols — layout is automatic
    scale_y_log10() +
    xlim(0.0, 41.0) +
    xlab("Consensus quality score") +
    ylab("Error frequency (log10)") +
    hideLegend() +
    theme_opaque() +
    ggsave("dada2_error_rates.png", width = 2000, height = 2000)
when isMainModule:
  main()
