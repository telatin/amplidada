import std/[os, osproc, strutils, times, unittest]
import amplidada

type MpFastqRec = tuple[name, sequence, quality: string]

proc mpMakeTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc mpWriteFastq(path: string; records: openArray[MpFastqRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

proc mpMakeErrMatrix(matchProb, mismatchProb: float64): LearnErrorsResult =
  result.qMin = 40
  result.qMax = 40
  for i in 0..15:
    result.counts[i] = @[0'i64]
    result.probs[i] = @[0.0]
  for fromIdx in 0..3:
    for toIdx in 0..3:
      let idx = fromIdx * 4 + toIdx
      result.probs[idx][0] = if fromIdx == toIdx: matchProb else: mismatchProb

proc mpDadaOpts(): DadaOptions =
  result = defaultDadaOptions()
  result.threads = 1
  result.minAbundance = 1
  result.omegaA = 0.1

proc mpDerepOpts(): DerepFastqOptions =
  result = defaultDerepFastqOptions()
  result.threads = 1
  result.keepReadMap = true

proc mpCompileDadaCli(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let cmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let res = execCmdEx(cmd)
  check res.exitCode == 0

suite "mergePairs":
  test "reverse complement helper":
    check reverseComplement("ACGTTGCA") == "TGCAACGT"

  test "perfect overlap merges paired denoised reads":
    let dir = mpMakeTempDir("amplidada-merge-perfect")
    let fwd = dir / "R1.fastq"
    let rev = dir / "R2.fastq"

    mpWriteFastq(fwd, @[
      ("r1", "ACGTTGCA", "IIIIIIII"),
      ("r2", "ACGTTGCA", "IIIIIIII"),
      ("r3", "ACGTTGCA", "IIIIIIII")
    ])
    mpWriteFastq(rev, @[
      ("r1", "TGCAACGT", "IIIIIIII"),
      ("r2", "TGCAACGT", "IIIIIIII"),
      ("r3", "TGCAACGT", "IIIIIIII")
    ])

    let dF = derepFastq(fwd, mpDerepOpts())
    let dR = derepFastq(rev, mpDerepOpts())
    let err = mpMakeErrMatrix(0.999, 0.000333333333)

    let asvF = dadaDenoise(dF, err, mpDadaOpts())
    let asvR = dadaDenoise(dR, err, mpDadaOpts())

    var mOpts = defaultMergePairsOptions()
    mOpts.minOverlap = 4
    mOpts.maxMismatch = 0
    mOpts.trimOverhang = false
    mOpts.threads = 1
    let merged = mergePairs(dF, dR, asvF, asvR, mOpts)

    check merged.stats.mergedPairs == 3
    check merged.stats.rejectedPairs == 0
    check merged.merged.len == 1
    check merged.merged[0].sequence == "ACGTTGCA"
    check merged.merged[0].abundance == 3

  test "overhang handling respects trimOverhang flag":
    let dir = mpMakeTempDir("amplidada-merge-overhang")
    let fwd = dir / "R1.fastq"
    let rev = dir / "R2.fastq"

    mpWriteFastq(fwd, @[
      ("r1", "GTTG", "IIII"),
      ("r2", "GTTG", "IIII")
    ])
    mpWriteFastq(rev, @[
      ("r1", "TGCAAC", "IIIIII"),
      ("r2", "TGCAAC", "IIIIII")
    ])

    let dF = derepFastq(fwd, mpDerepOpts())
    let dR = derepFastq(rev, mpDerepOpts())
    let err = mpMakeErrMatrix(0.999, 0.000333333333)

    let asvF = dadaDenoise(dF, err, mpDadaOpts())
    let asvR = dadaDenoise(dR, err, mpDadaOpts())

    # trimOverhang=false (R DADA2 default): include tails → full merged "GTTGCA"
    var keepTails = defaultMergePairsOptions()
    keepTails.minOverlap = 4
    keepTails.maxMismatch = 0
    keepTails.trimOverhang = false
    keepTails.threads = 1

    let withTails = mergePairs(dF, dR, asvF, asvR, keepTails)
    check withTails.stats.mergedPairs == 2
    check withTails.merged.len == 1
    check withTails.merged[0].sequence == "GTTGCA"

    # trimOverhang=true: clip to forward span → merged = "GTTG"
    var trimOpts = keepTails
    trimOpts.trimOverhang = true
    let trimmed = mergePairs(dF, dR, asvF, asvR, trimOpts)
    check trimmed.stats.mergedPairs == 2
    check trimmed.merged.len == 1
    check trimmed.merged[0].sequence == "GTTG"

  test "dada2 CLI paired mode runs and writes merged output":
    let dir = mpMakeTempDir("amplidada-merge-cli")
    let fwd = dir / "R1.fastq"
    let rev = dir / "R2.fastq"
    let errPath = dir / "err.tsv"
    let outPath = dir / "merged.tsv"
    let summaryPath = dir / "summary.txt"
    let binPath = dir / "dada2"

    mpWriteFastq(fwd, @[
      ("r1", "ACGTTGCA", "IIIIIIII"),
      ("r2", "ACGTTGCA", "IIIIIIII"),
      ("r3", "ACGTTGCA", "IIIIIIII")
    ])
    mpWriteFastq(rev, @[
      ("r1", "TGCAACGT", "IIIIIIII"),
      ("r2", "TGCAACGT", "IIIIIIII"),
      ("r3", "TGCAACGT", "IIIIIIII")
    ])

    writeLearnErrorsTsv(errPath, mpMakeErrMatrix(0.999, 0.000333333333))
    mpCompileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in-forward " & quoteShell(fwd) &
      " --in-reverse " & quoteShell(rev) &
      " --out " & quoteShell(outPath) &
      " --summary " & quoteShell(summaryPath) &
      " --err-matrix " & quoteShell(errPath) &
      " --min-overlap 4 --max-mismatch 0 --min-abundance 1 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0
    check fileExists(outPath)
    check fileExists(summaryPath)

    let outLines = readFile(outPath).splitLines()
    check outLines.len >= 2
    check outLines[0] == "merged_id\tsequence\tabundance"
    check outLines[1].contains("\tACGTTGCA\t3")

    let summary = readFile(summaryPath)
    check summary.contains("mode=paired")
    check summary.contains("merged_count=1")
    check summary.contains("merged_pairs=3")

  test "dada2 CLI paired mode supports canonical self-consistency flags":
    let dir = mpMakeTempDir("amplidada-merge-cli-self-consist")
    let fwd = dir / "R1.fastq"
    let rev = dir / "R2.fastq"
    let outPath = dir / "merged.tsv"
    let summaryPath = dir / "summary.txt"
    let binPath = dir / "dada2"

    mpWriteFastq(fwd, @[
      ("r1", "ACGTTGCA", "IIIIIIII"),
      ("r2", "ACGTTGCA", "IIIIIIII"),
      ("r3", "ACGTTGCA", "IIIIIIII"),
      ("r4", "ACGTTGCA", "IIIIIIII")
    ])
    mpWriteFastq(rev, @[
      ("r1", "TGCAACGT", "IIIIIIII"),
      ("r2", "TGCAACGT", "IIIIIIII"),
      ("r3", "TGCAACGT", "IIIIIIII"),
      ("r4", "TGCAACGT", "IIIIIIII")
    ])

    mpCompileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in-forward " & quoteShell(fwd) &
      " --in-reverse " & quoteShell(rev) &
      " --out " & quoteShell(outPath) &
      " --summary " & quoteShell(summaryPath) &
      " --dada-self-consist --dada-self-min-iter 2 --dada-self-max-iter 4 --dada-self-tol 0 " &
      " --min-overlap 4 --max-mismatch 0 --min-abundance 1 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0
    check fileExists(outPath)
    check fileExists(summaryPath)

    let summary = readFile(summaryPath)
    check summary.contains("mode=paired")
    check summary.contains("dada_self_consist_enabled=true")
    check summary.contains("forward_dada_self_consist_ran=true")
    check summary.contains("reverse_dada_self_consist_ran=true")
    check summary.contains("forward_dada_self_consist_iterations=")
    check summary.contains("reverse_dada_self_consist_iterations=")

  test "dada2 CLI paired mode supports remove-bimera flags":
    let dir = mpMakeTempDir("amplidada-merge-cli-bimera")
    let fwd = dir / "R1.fastq"
    let rev = dir / "R2.fastq"
    let errPath = dir / "err.tsv"
    let outPath = dir / "merged.tsv"
    let summaryPath = dir / "summary.txt"
    let binPath = dir / "dada2"

    mpWriteFastq(fwd, @[
      ("r1", "ACGTTGCA", "IIIIIIII"),
      ("r2", "ACGTTGCA", "IIIIIIII"),
      ("r3", "ACGTTGCA", "IIIIIIII")
    ])
    mpWriteFastq(rev, @[
      ("r1", "TGCAACGT", "IIIIIIII"),
      ("r2", "TGCAACGT", "IIIIIIII"),
      ("r3", "TGCAACGT", "IIIIIIII")
    ])

    writeLearnErrorsTsv(errPath, mpMakeErrMatrix(0.999, 0.000333333333))
    mpCompileDadaCli("src/cli/dada2.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --in-forward " & quoteShell(fwd) &
      " --in-reverse " & quoteShell(rev) &
      " --out " & quoteShell(outPath) &
      " --summary " & quoteShell(summaryPath) &
      " --err-matrix " & quoteShell(errPath) &
      " --remove-bimera --bimera-method pooled --min-overlap 4 --max-mismatch 0 --min-abundance 1 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0
    check fileExists(outPath)
    check fileExists(summaryPath)

    let summary = readFile(summaryPath)
    check summary.contains("bimera_enabled=true")
    check summary.contains("bimera_mode=bmPooled")
    check summary.contains("merged_count_pre_bimera=")
    check summary.contains("merged_count=")
