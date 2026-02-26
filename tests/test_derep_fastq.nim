import std/[os, osproc, times, unittest]
import amplidada

type DerepRec = tuple[name, sequence, quality: string]

proc makeDerepTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc writeDerepFastq(path: string; records: openArray[DerepRec]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer: close(f)
  for rec in records:
    f.write("@", rec.name, "\n")
    f.write(rec.sequence, "\n")
    f.write("+\n")
    f.write(rec.quality, "\n")

proc almostEq(a, b: float64; eps = 1e-9): bool =
  abs(a - b) <= eps

proc compileCli(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let cmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let res = execCmdEx(cmd)
  check res.exitCode == 0

suite "derepFastq":
  test "counts, mean qualities and map are correct":
    let dir = makeDerepTempDir("amplidada-derep-basic")
    let input = dir / "reads.fastq"
    writeDerepFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "CCCC", "!!!!"),
      ("r3", "AAAA", "####"),
      ("r4", "AAAA", "IIII")
    ])

    var opts = defaultDerepFastqOptions()
    opts.chunkSize = 2
    opts.threads = 1
    opts.keepReadMap = true
    opts.order = dsoAbundance

    let res = derepFastq(input, opts)
    check res.totalReads == 4
    check res.uniqueCount == 2
    check res.uniques.len == 2
    check res.uniques[0].sequence == "AAAA"
    check res.uniques[0].abundance == 3
    check res.uniques[1].sequence == "CCCC"
    check res.uniques[1].abundance == 1
    check res.readToUnique == @[0, 1, 0, 0]

    for q in res.uniques[0].meanQual:
      check almostEq(q, (40.0 + 2.0 + 40.0) / 3.0)
    for q in res.uniques[1].meanQual:
      check almostEq(q, 0.0)

  test "ordering by first seen is deterministic":
    let dir = makeDerepTempDir("amplidada-derep-first")
    let input = dir / "reads.fastq"
    writeDerepFastq(input, @[
      ("r1", "CCCC", "IIII"),
      ("r2", "AAAA", "IIII"),
      ("r3", "AAAA", "IIII")
    ])

    var opts = defaultDerepFastqOptions()
    opts.order = dsoFirstSeen
    opts.keepReadMap = true

    let res = derepFastq(input, opts)
    check res.uniques[0].sequence == "CCCC"
    check res.uniques[1].sequence == "AAAA"
    check res.readToUnique == @[0, 1, 1]

  test "no map mode omits readToUnique":
    let dir = makeDerepTempDir("amplidada-derep-nomap")
    let input = dir / "reads.fastq"
    writeDerepFastq(input, @[("r1", "AAAA", "IIII")])

    var opts = defaultDerepFastqOptions()
    opts.keepReadMap = false

    let res = derepFastq(input, opts)
    check res.totalReads == 1
    check res.readToUnique.len == 0

  test "invalid quality character raises":
    let dir = makeDerepTempDir("amplidada-derep-malformed")
    let input = dir / "reads.fastq"
    writeDerepFastq(input, @[("r1", "AAAA", "III ")])

    expect(ValueError):
      discard derepFastq(input)

  test "multithreaded result matches single-thread result":
    let dir = makeDerepTempDir("amplidada-derep-mt")
    let input = dir / "reads.fastq"

    var records: seq[DerepRec] = @[]
    for i in 0..<240:
      if i mod 3 == 0:
        records.add(("r" & $i, "AAAACCCC", "IIIIIIII"))
      elif i mod 3 == 1:
        records.add(("r" & $i, "GGGGTTTT", "####IIII"))
      else:
        records.add(("r" & $i, "CCCCAAAA", "IIII####"))
    writeDerepFastq(input, records)

    var o1 = defaultDerepFastqOptions()
    o1.threads = 1
    o1.chunkSize = 19
    o1.parallelMinChunk = 1

    var o4 = o1
    o4.threads = 4

    let a = derepFastq(input, o1)
    let b = derepFastq(input, o4)

    check a.totalReads == b.totalReads
    check a.uniqueCount == b.uniqueCount
    check a.readToUnique == b.readToUnique
    check a.uniques.len == b.uniques.len

    for i in 0..<a.uniques.len:
      check a.uniques[i].sequence == b.uniques[i].sequence
      check a.uniques[i].abundance == b.uniques[i].abundance
      check a.uniques[i].meanQual.len == b.uniques[i].meanQual.len
      for q in 0..<a.uniques[i].meanQual.len:
        check almostEq(a.uniques[i].meanQual[q], b.uniques[i].meanQual[q])

  test "derepFastq CLI accepts both --opt value and --opt=value":
    let dir = makeDerepTempDir("amplidada-derep-cli")
    let input = dir / "reads.fastq"
    writeDerepFastq(input, @[
      ("r1", "AAAA", "IIII"),
      ("r2", "AAAA", "IIII"),
      ("r3", "CCCC", "IIII")
    ])

    let binPath = dir / "derepFastq"
    compileCli("src/cli/derepFastq.nim", binPath)

    let outA = dir / "uniques-a.tsv"
    let mapA = dir / "map-a.tsv"
    let cmdA = quoteShell(binPath) &
      " --in " & quoteShell(input) &
      " --out " & quoteShell(outA) &
      " --map-out " & quoteShell(mapA)
    let resA = execCmdEx(cmdA)
    check resA.exitCode == 0
    check fileExists(outA)
    check fileExists(mapA)

    let outB = dir / "uniques-b.tsv"
    let mapB = dir / "map-b.tsv"
    let cmdB = quoteShell(binPath) &
      " --in=" & quoteShell(input) &
      " --out=" & quoteShell(outB) &
      " --map-out=" & quoteShell(mapB)
    let resB = execCmdEx(cmdB)
    check resB.exitCode == 0
    check fileExists(outB)
    check fileExists(mapB)
