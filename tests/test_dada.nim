import std/[algorithm, os, osproc, strutils, times, unittest]
import amplidada

type DadaFastqRec = tuple[name, sequence, quality: string]

proc makeDadaTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc writeDadaFastq(path: string; records: openArray[DadaFastqRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

proc compileDadaCli(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let cmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let res = execCmdEx(cmd)
  check res.exitCode == 0

proc makeMatrix(matchProb, mismatchProb: float64): LearnErrorsResult =
  result.qMin = 40
  result.qMax = 40
  for i in 0..15:
    result.counts[i] = @[0'i64]
    result.probs[i] = @[0.0]
  for fromIdx in 0..3:
    for toIdx in 0..3:
      let idx = fromIdx * 4 + toIdx
      result.probs[idx][0] = if fromIdx == toIdx: matchProb else: mismatchProb

suite "dada denoising":
  test "singleton does not split into an ASV by abundance p-value":
    let dir = makeDadaTempDir("amplidada-dada-singleton")
    let input = dir / "reads.fastq"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<20:
      records.add(("a" & $i, "AAAA", "IIII"))
    records.add(("v1", "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let d = derepFastq(input)
    let err = makeMatrix(0.97, 0.01)

    var opts = defaultDadaOptions()
    opts.minAbundance = 1
    opts.omegaA = 0.1
    opts.threads = 1

    let res = dadaDenoise(d, err, opts)
    check d.uniqueCount == 2
    check res.asvCount == 1
    check res.asvs[0].abundance == 21
    check res.uniqueToAsv.len == d.uniqueCount
    check res.uniqueToAsv[0] == 0
    check res.uniqueToAsv[1] == 0

  test "abundant variant splits into a separate ASV":
    let dir = makeDadaTempDir("amplidada-dada-split")
    let input = dir / "reads.fastq"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<30:
      records.add(("a" & $i, "AAAA", "IIII"))
    for i in 0..<8:
      records.add(("b" & $i, "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let d = derepFastq(input)
    let err = makeMatrix(0.997, 0.001)

    var opts = defaultDadaOptions()
    opts.minAbundance = 2
    opts.omegaA = 1e-3
    opts.threads = 1

    let res = dadaDenoise(d, err, opts)
    check res.asvCount == 2

    var seqs: seq[string] = @[]
    for asv in res.asvs:
      seqs.add(asv.sequence)
    seqs.sort(system.cmp[string])
    check seqs == @["AAAA", "AAAT"]

  test "omegaC gating can keep a non-center unique as its own ASV":
    let dir = makeDadaTempDir("amplidada-dada-omegac")
    let input = dir / "reads.fastq"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<20:
      records.add(("a" & $i, "AAAA", "IIII"))
    records.add(("v1", "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let d = derepFastq(input)
    let err = makeMatrix(0.97, 0.01)

    var opts = defaultDadaOptions()
    opts.minAbundance = 1
    opts.omegaA = 0.1
    opts.omegaC = 0.2
    opts.threads = 1

    let res = dadaDenoise(d, err, opts)
    check res.asvCount == 2
    check res.correctedUniques == 1
    check res.uncorrectedUniques == 1
    check res.uniqueToCenterUnique[0] == 0
    check res.uniqueToCenterUnique[1] == 1

  test "different read lengths denoise independently":
    let dir = makeDadaTempDir("amplidada-dada-lengths")
    let input = dir / "reads.fastq"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<5:
      records.add(("a" & $i, "AAAA", "IIII"))
    for i in 0..<4:
      records.add(("b" & $i, "AAAAA", "IIIII"))
    writeDadaFastq(input, records)

    let d = derepFastq(input)
    let err = makeMatrix(0.99, 0.003333333333)

    var opts = defaultDadaOptions()
    opts.minAbundance = 1
    opts.threads = 2

    let res = dadaDenoise(d, err, opts)
    check res.asvCount == 2
    check res.groupCount == 2

  test "canonical self-consistency loop runs and returns final dada + errors":
    let dir = makeDadaTempDir("amplidada-dada-self-consist")
    let input = dir / "reads.fastq"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<24:
      records.add(("a" & $i, "AAAA", "IIII"))
    for i in 0..<6:
      records.add(("b" & $i, "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let d = derepFastq(input)

    var dOpts = defaultDadaOptions()
    dOpts.threads = 1
    dOpts.minAbundance = 1
    dOpts.omegaA = 0.05

    var scOpts = defaultDadaSelfConsistOptions()
    scOpts.enabled = true
    scOpts.minIterations = 2
    scOpts.maxIterations = 5
    scOpts.tol = 0.0
    scOpts.forceOmegaCZero = true

    let sc = dadaSelfConsistent(d, defaultLearnErrorsOptions(), dOpts, scOpts)
    check sc.iterationsRun >= 2
    check sc.iterationsRun <= 5
    check sc.errorMatrix.qMin == 0
    check sc.errorMatrix.qMax == 40
    check sc.dada.totalReads == d.totalReads
    check sc.dada.asvCount >= 1

  test "dada2 CLI runs with provided error matrix and writes outputs":
    let dir = makeDadaTempDir("amplidada-dada-cli")
    let input = dir / "reads.fastq"
    let matrixPath = dir / "err.tsv"
    let outPath = dir / "asv.tsv"
    let summaryPath = dir / "summary.txt"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<12:
      records.add(("a" & $i, "AAAA", "IIII"))
    for i in 0..<5:
      records.add(("b" & $i, "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let err = makeMatrix(0.999, 0.000333333333)
    writeLearnErrorsTsv(matrixPath, err)

    let binPath = dir / "dada2"
    compileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in=" & quoteShell(input) &
      " --out=" & quoteShell(outPath) &
      " --err-matrix=" & quoteShell(matrixPath) &
      " --summary=" & quoteShell(summaryPath) &
      " --omega-a=0.01 --min-abundance=2 --max-iter=64 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0
    check fileExists(outPath)
    check fileExists(summaryPath)

    let outLines = readFile(outPath).splitLines()
    check outLines.len >= 2
    check outLines[0] == "asv_id\tsequence\tabundance\tcluster_size\tcenter_unique_index\tfirst_seen_read_index"

    let summary = readFile(summaryPath)
    check summary.contains("err_loaded_from_file=true")
    check summary.contains("asv_count=2")

  test "dada2 CLI canonical self-consistency flags run in single-end mode":
    let dir = makeDadaTempDir("amplidada-dada-cli-self-consist")
    let input = dir / "reads.fastq"
    let outPath = dir / "asv.tsv"
    let summaryPath = dir / "summary.txt"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<16:
      records.add(("a" & $i, "AAAA", "IIII"))
    for i in 0..<4:
      records.add(("b" & $i, "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let binPath = dir / "dada2"
    compileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in " & quoteShell(input) &
      " --out " & quoteShell(outPath) &
      " --summary " & quoteShell(summaryPath) &
      " --dada-self-consist --dada-self-min-iter 2 --dada-self-max-iter 4 --dada-self-tol 0 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0

    let summary = readFile(summaryPath)
    check summary.contains("dada_self_consist_enabled=true")
    check summary.contains("dada_self_consist_ran=true")

  test "dada2 CLI omega-c can preserve non-center unique":
    let dir = makeDadaTempDir("amplidada-dada-cli-omegac")
    let input = dir / "reads.fastq"
    let matrixPath = dir / "err.tsv"
    let outPath = dir / "asv.tsv"
    let summaryPath = dir / "summary.txt"

    var records: seq[DadaFastqRec] = @[]
    for i in 0..<20:
      records.add(("a" & $i, "AAAA", "IIII"))
    records.add(("v1", "AAAT", "IIII"))
    writeDadaFastq(input, records)

    let err = makeMatrix(0.97, 0.01)
    writeLearnErrorsTsv(matrixPath, err)

    let binPath = dir / "dada2"
    compileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in " & quoteShell(input) &
      " --out " & quoteShell(outPath) &
      " --err-matrix " & quoteShell(matrixPath) &
      " --summary " & quoteShell(summaryPath) &
      " --omega-a 0.1 --omega-c 0.2 --min-abundance 1 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0

    let summary = readFile(summaryPath)
    check summary.contains("asv_count=2")
    check summary.contains("corrected_uniques=1")
    check summary.contains("uncorrected_uniques=1")
