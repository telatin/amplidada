import std/[os, osproc, times, unittest]
import amplidada

type BatchRec = tuple[name, sequence, quality: string]

proc makeBatchTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc writeBatchFastq(path: string; records: openArray[BatchRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

suite "filterAndTrim batch API":
  test "discover paired samples from standard naming":
    let dir = makeBatchTempDir("amplidada-discover")
    writeBatchFastq(dir / "S1_R1_001.fastq", @[("a", "ACGT", "IIII")])
    writeBatchFastq(dir / "S1_R2_001.fastq", @[("a", "TGCA", "IIII")])
    writeBatchFastq(dir / "S2_R1_001.fastq", @[("b", "ACGT", "IIII")])
    writeBatchFastq(dir / "S2_R2_001.fastq", @[("b", "TGCA", "IIII")])

    let samples = discoverPairedSamples(dir)
    check samples.len == 2
    check samples[0].sampleId == "S1"
    check samples[1].sampleId == "S2"
    check samples[0].readsForward.endsWith("S1_R1_001.fastq")
    check samples[0].readsReverse.endsWith("S1_R2_001.fastq")

  test "read sample sheet with relative paths":
    let dir = makeBatchTempDir("amplidada-sheet")
    let readsDir = dir / "reads"
    createDir(readsDir)

    writeBatchFastq(readsDir / "x_fwd.fastq", @[("r1", "ACGT", "IIII")])
    writeBatchFastq(readsDir / "x_rev.fastq", @[("r1", "TGCA", "IIII")])

    let sheet = dir / "samples.csv"
    var f: File
    doAssert open(f, sheet, fmWrite)
    f.writeLine("SampleID,reads_forward,reads_reverse")
    f.writeLine("sampleX,reads/x_fwd.fastq,reads/x_rev.fastq")
    close(f)

    let samples = readPairedSampleSheet(sheet)
    check samples.len == 1
    check samples[0].sampleId == "sampleX"
    check samples[0].readsForward.endsWith("reads/x_fwd.fastq")
    check samples[0].readsReverse.endsWith("reads/x_rev.fastq")

  test "build jobs creates expected output paths":
    let samples = @[
      PairedSampleInput(sampleId: "Mock A", readsForward: "/tmp/a_R1.fastq.gz", readsReverse: "/tmp/a_R2.fastq.gz")
    ]
    let jobs = buildPairedJobs(samples, "/tmp/out", "_F.fq.gz", "_R.fq.gz")
    check jobs.len == 1
    check jobs[0].outForward == "/tmp/out/Mock_A_F.fq.gz"
    check jobs[0].outReverse == "/tmp/out/Mock_A_R.fq.gz"

  test "batch filter runs in parallel and aggregates per-sample outputs":
    let dir = makeBatchTempDir("amplidada-batch")
    let inputDir = dir / "in"
    let outputDir = dir / "out"
    createDir(inputDir)
    createDir(outputDir)

    writeBatchFastq(inputDir / "Good_R1.fastq", @[("g1", "ACGT", "IIII")])
    writeBatchFastq(inputDir / "Good_R2.fastq", @[("g1", "TGCA", "IIII")])

    writeBatchFastq(inputDir / "Bad_R1.fastq", @[("b1", "ANNN", "IIII")])
    writeBatchFastq(inputDir / "Bad_R2.fastq", @[("b1", "TGCA", "IIII")])

    let samples = discoverPairedSamples(inputDir, r1Token = "_R1", r2Token = "_R2")
    let jobs = buildPairedJobs(samples, outputDir)

    var filterOpts = defaultPairedFastqFilterOptions()
    filterOpts.forward.minLen = 1
    filterOpts.reverse.minLen = 1
    filterOpts.forward.minQ = 0
    filterOpts.reverse.minQ = 0
    filterOpts.forward.maxEE = 100.0
    filterOpts.reverse.maxEE = 100.0
    filterOpts.forward.maxN = 0
    filterOpts.reverse.maxN = 0
    filterOpts.forward.threads = 1
    filterOpts.reverse.threads = 1

    var batchOpts = defaultBatchFilterOptions()
    batchOpts.sampleThreads = 2
    batchOpts.continueOnError = true

    let batch = filterAndTrimPaired(jobs, filterOpts, batchOpts)
    check batch.samples.len == 2
    check batch.failedSamples == 0

    var keptPairs = 0'i64
    for s in batch.samples:
      keptPairs += s.stats.outputPairs
    check keptPairs == 1
    check fileExists(outputDir / "Good_R1.filtered.fastq.gz")
    check fileExists(outputDir / "Good_R2.filtered.fastq.gz")

  test "filterAndTrim CLI accepts both --opt value and --opt=value forms":
    let dir = makeBatchTempDir("amplidada-cli-args")
    let inputDir = dir / "reads"
    let outDirA = dir / "out-a"
    let outDirB = dir / "out-b"
    createDir(inputDir)
    createDir(outDirA)
    createDir(outDirB)

    writeBatchFastq(inputDir / "S1_R1.fastq", @[("r1", "ACGTACGT", "IIIIIIII")])
    writeBatchFastq(inputDir / "S1_R2.fastq", @[("r1", "TGCATGCA", "IIIIIIII")])

    let binPath = dir / "filterAndTrim"
    compileBinary("src/cli/filterAndTrim.nim", binPath)

    let cmdA = quoteShell(binPath) &
      " --input-dir " & quoteShell(inputDir) &
      " --output-dir " & quoteShell(outDirA) &
      " --min-len 1 --min-q 0 --max-ee 100 --threads 1 --sample-threads 1"
    let resA = execCmdEx(cmdA)
    check resA.exitCode == 0
    check fileExists(outDirA / "S1_R1.filtered.fastq.gz")
    check fileExists(outDirA / "S1_R2.filtered.fastq.gz")

    let cmdB = quoteShell(binPath) &
      " --input-dir=" & quoteShell(inputDir) &
      " --output-dir=" & quoteShell(outDirB) &
      " --min-len=1 --min-q=0 --max-ee=100 --threads=1 --sample-threads=1"
    let resB = execCmdEx(cmdB)
    check resB.exitCode == 0
    check fileExists(outDirB / "S1_R1.filtered.fastq.gz")
    check fileExists(outDirB / "S1_R2.filtered.fastq.gz")

  test "filterAndTrim CLI --verbose prints per-sample progress":
    let dir = makeBatchTempDir("amplidada-cli-verbose")
    let inputDir = dir / "reads"
    let outDir = dir / "out"
    createDir(inputDir)
    createDir(outDir)

    writeBatchFastq(inputDir / "V1_R1.fastq", @[("r1", "ACGTACGT", "IIIIIIII")])
    writeBatchFastq(inputDir / "V1_R2.fastq", @[("r1", "TGCATGCA", "IIIIIIII")])

    let binPath = dir / "filterAndTrim"
    compileBinary("src/cli/filterAndTrim.nim", binPath)

    let cmd = quoteShell(binPath) &
      " --input-dir " & quoteShell(inputDir) &
      " --output-dir " & quoteShell(outDir) &
      " --min-len 1 --min-q 0 --max-ee 100 --threads 1 --sample-threads 1 --verbose"
    let res = execCmdEx(cmd)
    check res.exitCode == 0
    check "[filterAndTrim] Starting sample 1/1: V1" in res.output
    check "[filterAndTrim] Completed sample 1/1: V1" in res.output
