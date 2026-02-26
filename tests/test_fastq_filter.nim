import std/[os, strformat, times, unittest]
import readfx
import amplidada

type TestRec = tuple[name, sequence, quality: string]

proc makeTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc writeFastq(path: string; records: openArray[TestRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

proc readFastq(path: string): seq[FQRecord] =
  for rec in readFQ(path):
    result.add(rec)

proc readNames(path: string): seq[string] =
  for rec in readFQ(path):
    result.add(rec.name)

suite "fastqFilter rules and behavior":
  test "trim and truncLen change sequence as expected":
    let dir = makeTempDir("amplidada-trim")
    let input = dir / "in.fastq"
    let output = dir / "out.fastq"
    writeFastq(input, @[("r1", "AACCGGTT", "IIIIIIII")])

    var opts = defaultFastqFilterOptions()
    opts.trimLeft = 1
    opts.trimRight = 2
    opts.truncLen = 4
    opts.minLen = 1
    opts.minQ = 0
    opts.maxEE = 100.0
    opts.maxN = 0
    opts.threads = 1

    let stats = fastqFilter(input, output, opts)
    let outReads = readFastq(output)

    check stats.inputReads == 1
    check stats.outputReads == 1
    check outReads.len == 1
    check outReads[0].sequence == "ACCG"

  test "truncQ truncates at first low-quality base":
    let dir = makeTempDir("amplidada-truncq")
    let input = dir / "in.fastq"
    let output = dir / "out.fastq"
    writeFastq(input, @[("r1", "AACCGG", "II!III")])

    var opts = defaultFastqFilterOptions()
    opts.truncQ = 2
    opts.minLen = 1
    opts.minQ = 0
    opts.maxEE = 100.0
    opts.threads = 1

    let stats = fastqFilter(input, output, opts)
    let outReads = readFastq(output)

    check stats.outputReads == 1
    check outReads.len == 1
    check outReads[0].sequence == "AA"

  test "maxN rejects reads with too many ambiguous bases":
    let dir = makeTempDir("amplidada-maxn")
    let input = dir / "in.fastq"
    let output = dir / "out.fastq"
    writeFastq(input, @[("r1", "AANN", "IIII")])

    var opts = defaultFastqFilterOptions()
    opts.maxN = 1
    opts.minLen = 1
    opts.minQ = 0
    opts.maxEE = 100.0

    let stats = fastqFilter(input, output, opts)

    check stats.inputReads == 1
    check stats.outputReads == 0
    check stats.discardedReads == 1
    check stats.rejectedByReason[frrTooManyN] == 1

  test "minQ rejects when minimum quality is <= threshold":
    let dir = makeTempDir("amplidada-minq")
    let input = dir / "in.fastq"
    let output = dir / "out.fastq"
    writeFastq(input, @[("r1", "AAAA", "I#II")])

    var opts = defaultFastqFilterOptions()
    opts.minQ = 2
    opts.minLen = 1
    opts.maxEE = 100.0

    let stats = fastqFilter(input, output, opts)

    check stats.outputReads == 0
    check stats.rejectedByReason[frrLowMinQ] == 1

  test "maxEE rejects reads with high expected errors":
    let dir = makeTempDir("amplidada-maxee")
    let input = dir / "in.fastq"
    let output = dir / "out.fastq"
    writeFastq(input, @[("r1", "AAAA", "####")])

    var opts = defaultFastqFilterOptions()
    opts.minQ = 0
    opts.minLen = 1
    opts.maxEE = 2.0

    let stats = fastqFilter(input, output, opts)

    check stats.outputReads == 0
    check stats.rejectedByReason[frrTooManyExpectedErrors] == 1

  test "paired filter keeps only pairs where both reads pass":
    let dir = makeTempDir("amplidada-paired")
    let in1 = dir / "r1.fastq"
    let in2 = dir / "r2.fastq"
    let out1 = dir / "o1.fastq"
    let out2 = dir / "o2.fastq"

    writeFastq(in1, @[
      ("p1", "ACGT", "IIII"),
      ("p2", "ANGT", "IIII")
    ])
    writeFastq(in2, @[
      ("p1", "TGCA", "IIII"),
      ("p2", "TGCA", "IIII")
    ])

    var popts = defaultPairedFastqFilterOptions()
    popts.forward.minLen = 1
    popts.reverse.minLen = 1
    popts.forward.minQ = 0
    popts.reverse.minQ = 0
    popts.forward.maxEE = 100.0
    popts.reverse.maxEE = 100.0
    popts.forward.maxN = 0
    popts.reverse.maxN = 0
    popts.forward.threads = 2
    popts.reverse.threads = 2

    let stats = fastqPairedFilter(in1, in2, out1, out2, popts)

    check stats.inputPairs == 2
    check stats.outputPairs == 1
    check stats.discardedPairs == 1
    check stats.rejectedByPairReason[prrForwardFailed] == 1
    check stats.forwardStats.rejectedByReason[frrTooManyN] == 1

    check readNames(out1) == @["p1"]
    check readNames(out2) == @["p1"]

  test "multithreaded filtering is deterministic and matches single-thread output":
    let dir = makeTempDir("amplidada-threads")
    let input = dir / "in.fastq"
    let out1 = dir / "out-1.fastq"
    let out4 = dir / "out-4.fastq.gz"

    var records: seq[TestRec] = @[]
    for i in 0..<160:
      let name = &"r{i:03}"
      if i mod 5 == 0:
        records.add((name, "ACNTACGT", "IIIIIIII"))
      else:
        records.add((name, "ACGTACGT", "IIIIIIII"))
    writeFastq(input, records)

    var o1 = defaultFastqFilterOptions()
    o1.maxN = 0
    o1.minLen = 1
    o1.minQ = 0
    o1.maxEE = 100.0
    o1.chunkSize = 17
    o1.parallelMinChunk = 1
    o1.threads = 1

    var o4 = o1
    o4.threads = 4

    let s1 = fastqFilter(input, out1, o1)
    let s4 = fastqFilter(input, out4, o4)

    check s1.inputReads == s4.inputReads
    check s1.outputReads == s4.outputReads
    check s1.discardedReads == s4.discardedReads
    check s1.rejectedByReason == s4.rejectedByReason

    let n1 = readNames(out1)
    let n4 = readNames(out4)
    check n1 == n4
