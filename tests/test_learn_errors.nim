import std/[os, osproc, strutils, times, unittest]
import amplidada

type LearnRec = tuple[name, sequence, quality: string]

proc makeLearnTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc writeLearnFastq(path: string; records: openArray[LearnRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

proc almostEqLE(a, b: float64; eps = 1e-12): bool =
  abs(a - b) <= eps

proc compileLearnCli(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let cmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let res = execCmdEx(cmd)
  check res.exitCode == 0

suite "learnErrors phase 1":
  test "counts and probabilities are computed from derep centers":
    let dir = makeLearnTempDir("amplidada-learn-basic")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "AAAA", "IIII"),
      ("r3", "AACA", "IIII")
    ])

    var dOpts = defaultDerepFastqOptions()
    dOpts.threads = 1
    dOpts.keepReadMap = false

    let d = derepFastq(input, dOpts)

    var lOpts = defaultLearnErrorsOptions()
    lOpts.qMin = 0
    lOpts.qMax = 40
    lOpts.pseudocount = 1.0

    let e = learnErrorsFromDerep(d, lOpts)

    check e.centerCountByLength == 1
    check e.getCount('A', 'A', 40) == 11
    check e.getCount('A', 'C', 40) == 1

    check almostEqLE(e.getProb('A', 'A', 40), 12.0 / 16.0)
    check almostEqLE(e.getProb('A', 'C', 40), 2.0 / 16.0)
    check almostEqLE(e.getProb('A', 'G', 40), 1.0 / 16.0)
    check almostEqLE(e.getProb('A', 'T', 40), 1.0 / 16.0)

  test "quality bins are clamped to configured range":
    let dir = makeLearnTempDir("amplidada-learn-clamp")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "A", "!"),
      ("r2", "C", "h")
    ])

    var dOpts = defaultDerepFastqOptions()
    dOpts.threads = 1
    let d = derepFastq(input, dOpts)

    var lOpts = defaultLearnErrorsOptions()
    lOpts.qMin = 5
    lOpts.qMax = 40
    lOpts.pseudocount = 0.0

    let e = learnErrorsFromDerep(d, lOpts)
    check e.getCount('A', 'A', 5) == 1
    check e.getCount('A', 'C', 40) == 1

  test "zero-data quality bins use fallback identity probabilities":
    let dir = makeLearnTempDir("amplidada-learn-fallback")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "AAAA", "IIII")
    ])

    let d = derepFastq(input)

    var lOpts = defaultLearnErrorsOptions()
    lOpts.qMin = 0
    lOpts.qMax = 40
    lOpts.zeroTotalMismatchProb = 1e-6

    let e = learnErrorsFromDerep(d, lOpts)

    check almostEqLE(e.getProb('A', 'A', 10), 1.0 - 3.0e-6)
    check almostEqLE(e.getProb('A', 'C', 10), 1.0e-6)
    check almostEqLE(
      e.getProb('A', 'A', 10) + e.getProb('A', 'C', 10) + e.getProb('A', 'G', 10) + e.getProb('A', 'T', 10),
      1.0
    )

  test "learnErrors CLI accepts both --opt value and --opt=value forms":
    let dir = makeLearnTempDir("amplidada-learn-cli")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "AACA", "IIII")
    ])

    let binPath = dir / "learnErrors"
    compileLearnCli("src/cli/learnErrors.nim", binPath)

    let outA = dir / "err-a.tsv"
    let countsA = dir / "counts-a.tsv"
    let cmdA = quoteShell(binPath) &
      " --in " & quoteShell(input) &
      " --out " & quoteShell(outA) &
      " --counts-out " & quoteShell(countsA) &
      " --threads 1 --chunk-size 2"
    let resA = execCmdEx(cmdA)
    check resA.exitCode == 0
    check fileExists(outA)
    check fileExists(countsA)

    let outB = dir / "err-b.tsv"
    let cmdB = quoteShell(binPath) &
      " --in=" & quoteShell(input) &
      " --out=" & quoteShell(outB) &
      " --threads=1 --chunk-size=2"
    let resB = execCmdEx(cmdB)
    check resB.exitCode == 0
    check fileExists(outB)

    let lines = readFile(outB).splitLines()
    check lines.len > 1
    check lines[0] == "from\tto\tquality\tcount\tprob"

  test "learnErrorsFromFastqPaths combines files and respects nbases":
    let dir = makeLearnTempDir("amplidada-learn-multi")
    let inputA = dir / "a.fastq"
    let inputB = dir / "b.fastq"

    writeLearnFastq(inputA, @[
      ("a1", "AAAA", "IIII"),
      ("a2", "AAAA", "IIII")
    ])
    writeLearnFastq(inputB, @[
      ("b1", "CCCC", "IIII"),
      ("b2", "CCCC", "IIII")
    ])

    var dOpts = defaultDerepFastqOptions()
    dOpts.threads = 1

    var lOpts = defaultLearnErrorsOptions()
    lOpts.pseudocount = 0.0

    let batch = learnErrorsFromFastqPaths(@[inputA, inputB], dOpts, lOpts, nbases = 8)

    check batch.filesRequested == 2
    check batch.filesUsed == 1
    check batch.readsUsed == 2
    check batch.basesUsed == 8
    check batch.matrix.getCount('A', 'A', 40) == 8
    check batch.matrix.getCount('C', 'C', 40) == 0

  test "learnErrors CLI supports repeated --in and nbases cutoff":
    let dir = makeLearnTempDir("amplidada-learn-cli-multi")
    let inputA = dir / "a.fastq"
    let inputB = dir / "b.fastq"

    writeLearnFastq(inputA, @[
      ("a1", "AAAA", "IIII"),
      ("a2", "AAAA", "IIII")
    ])
    writeLearnFastq(inputB, @[
      ("b1", "CCCC", "IIII"),
      ("b2", "CCCC", "IIII")
    ])

    let binPath = dir / "learnErrors"
    compileLearnCli("src/cli/learnErrors.nim", binPath)

    let outTsv = dir / "err.tsv"
    let summary = dir / "summary.txt"
    let cmd = quoteShell(binPath) &
      " --in " & quoteShell(inputA) &
      " --in " & quoteShell(inputB) &
      " --nbases 8" &
      " --out " & quoteShell(outTsv) &
      " --summary " & quoteShell(summary) &
      " --threads 1 --chunk-size 8 --quiet"
    let res = execCmdEx(cmd)
    check res.exitCode == 0
    check fileExists(outTsv)
    check fileExists(summary)

    let summaryText = readFile(summary)
    check summaryText.contains("input_files_requested=2")
    check summaryText.contains("input_files_used=1")
    check summaryText.contains("input_reads=2")
    check summaryText.contains("input_bases=8")

  test "self-consistent learning honors min iterations and converges":
    let dir = makeLearnTempDir("amplidada-learn-self-consist")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "AAAA", "IIII"),
      ("r3", "AAAT", "IIII"),
      ("r4", "AAAT", "IIII")
    ])

    let d = derepFastq(input)
    var lOpts = defaultLearnErrorsOptions()
    lOpts.pseudocount = 1.0

    var scOpts = defaultLearnErrorsSelfConsistOptions()
    scOpts.enabled = true
    scOpts.minIterations = 2
    scOpts.maxIterations = 5
    scOpts.tol = 0.0
    scOpts.maxCentersPerLength = 8

    let scRes = learnErrorsSelfConsistentFromDereps([d], lOpts, scOpts)

    check scRes.iterationsRun >= 2
    check scRes.iterationsRun <= 5
    check scRes.converged
    check scRes.matrix.centerCountByLength == 1
    check scRes.centersByLength.len == 1
    check scRes.centersByLength[0].seqLen == 4
    check scRes.maxAbsProbDelta == 0.0

  test "self-consistent options validation rejects invalid settings":
    var scOpts = defaultLearnErrorsSelfConsistOptions()
    scOpts.enabled = true
    scOpts.minIterations = 3
    scOpts.maxIterations = 2

    expect ValueError:
      validateOptions(scOpts)

  test "learnErrors CLI supports self-consistency flags and summary fields":
    let dir = makeLearnTempDir("amplidada-learn-cli-self-consist")
    let input = dir / "reads.fastq"
    writeLearnFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "AAAA", "IIII"),
      ("r3", "AAAT", "IIII"),
      ("r4", "AAAT", "IIII")
    ])

    let binPath = dir / "learnErrors"
    compileLearnCli("src/cli/learnErrors.nim", binPath)

    let outTsv = dir / "err.tsv"
    let summary = dir / "summary.txt"
    let cmd = quoteShell(binPath) &
      " --in " & quoteShell(input) &
      " --out " & quoteShell(outTsv) &
      " --summary " & quoteShell(summary) &
      " --self-consist --min-iter 2 --max-iter 4 --tol 0 --max-centers-per-length 8 --quiet"
    let res = execCmdEx(cmd)
    check res.exitCode == 0
    check fileExists(outTsv)
    check fileExists(summary)

    let summaryText = readFile(summary)
    check summaryText.contains("self_consist_enabled=true")
    check summaryText.contains("self_consist_iterations=")
    check summaryText.contains("self_consist_converged=true")
