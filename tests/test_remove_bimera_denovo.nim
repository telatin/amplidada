import std/[os, osproc, sequtils, strutils, tables, times, unittest]
import amplidada

type
  TsvLongRow = tuple[sampleId, sequence: string, abundance: int64]

proc rbMakeTempDir(prefix: string): string =
  result = getTempDir() / (prefix & "-" & $cast[int64](epochTime() * 1_000_000))
  createDir(result)

proc rbWriteLongTsv(path: string; rows: openArray[TsvLongRow]) =
  var f: File
  doAssert open(f, path, fmWrite)
  defer:
    close(f)
  f.writeLine("SampleID\tSequence\tAbundance")
  for row in rows:
    f.writeLine(row.sampleId & "\t" & row.sequence & "\t" & $row.abundance)

proc rbCompileCli(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let cmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let res = execCmdEx(cmd)
  check res.exitCode == 0

suite "removeBimeraDenovo":
  let parentA = "AAAAACCCCC"
  let parentB = "GGGGGTTTTT"
  let exactChimera = "AAAAATTTTT"
  let oneOffChimera = "AAAACTTTTT"
  let oneOffDeletion = "AAAATTTTT"
  let oneOffInsertion = "AAAAACTTTTT"

  test "isBimera identifies exact and one-off bimeras":
    check isBimera(exactChimera, @[parentA, parentB])
    check not isBimera(oneOffChimera, @[parentA, parentB])
    check isBimera(oneOffChimera, @[parentA, parentB], allowOneOff = true, minOneOffParentDistance = 4)
    check not isBimera(oneOffDeletion, @[parentA, parentB], allowOneOff = false)
    check isBimera(oneOffDeletion, @[parentA, parentB], allowOneOff = true, minOneOffParentDistance = 4)
    check not isBimera(oneOffInsertion, @[parentA, parentB], allowOneOff = false)
    check isBimera(oneOffInsertion, @[parentA, parentB], allowOneOff = true, minOneOffParentDistance = 4)

  test "isBimeraDenovo pooled flags exact bimera":
    var opts = defaultBimeraOptions()
    opts.threads = 1
    opts.minFoldParentOverAbundance = 2.0
    opts.minParentAbundance = 8
    opts.allowOneOff = false

    let pooled = @[
      UniqueSequenceAbundance(sequence: parentA, abundance: 120),
      UniqueSequenceAbundance(sequence: parentB, abundance: 100),
      UniqueSequenceAbundance(sequence: exactChimera, abundance: 20),
      UniqueSequenceAbundance(sequence: "CCCCCAAAAA", abundance: 30)
    ]

    let res = isBimeraDenovo(pooled, opts)
    var bySeq = initTable[string, bool]()
    for i, sequence in res.sequences:
      bySeq[sequence] = res.flags[i]

    check bySeq.getOrDefault(exactChimera, false)
    check not bySeq.getOrDefault(parentA, false)
    check not bySeq.getOrDefault(parentB, false)
    check res.exactBimeraCount == 1
    check res.oneOffMismatchBimeraCount == 0
    check res.oneOffIndelBimeraCount == 0

  test "isBimeraDenovo pooled can flag indel one-off bimera":
    var opts = defaultBimeraOptions()
    opts.threads = 1
    opts.minFoldParentOverAbundance = 2.0
    opts.minParentAbundance = 8
    opts.allowOneOff = true
    opts.minOneOffParentDistance = 4

    let pooled = @[
      UniqueSequenceAbundance(sequence: parentA, abundance: 120),
      UniqueSequenceAbundance(sequence: parentB, abundance: 100),
      UniqueSequenceAbundance(sequence: oneOffDeletion, abundance: 20)
    ]

    let res = isBimeraDenovo(pooled, opts)
    var bySeq = initTable[string, bool]()
    for i, sequence in res.sequences:
      bySeq[sequence] = res.flags[i]

    check bySeq.getOrDefault(oneOffDeletion, false)
    check not bySeq.getOrDefault(parentA, false)
    check not bySeq.getOrDefault(parentB, false)
    check res.exactBimeraCount == 0
    check res.oneOffMismatchBimeraCount == 0
    check res.oneOffIndelBimeraCount == 1

  test "removeBimeraDenovo consensus and per-sample behave differently":
    let rows = @[
      SampleSequenceAbundance(sampleId: "S1", sequence: parentA, abundance: 50),
      SampleSequenceAbundance(sampleId: "S1", sequence: parentB, abundance: 40),
      SampleSequenceAbundance(sampleId: "S1", sequence: exactChimera, abundance: 10),
      SampleSequenceAbundance(sampleId: "S2", sequence: parentA, abundance: 4),
      SampleSequenceAbundance(sampleId: "S2", sequence: parentB, abundance: 3),
      SampleSequenceAbundance(sampleId: "S2", sequence: exactChimera, abundance: 2)
    ]

    var opts = defaultBimeraOptions()
    opts.threads = 1
    opts.minFoldParentOverAbundance = 1.5
    opts.minParentAbundance = 2

    let consensusRes = removeBimeraDenovo(rows, bmConsensus, opts, defaultBimeraConsensusOptions())
    check consensusRes.removedUniqueCount == 1
    check consensusRes.kept.allIt(it.sequence != exactChimera)
    check consensusRes.exactBimeraCount == 1
    check consensusRes.oneOffMismatchBimeraCount == 0
    check consensusRes.oneOffIndelBimeraCount == 0

    let perSampleRes = removeBimeraDenovo(rows, bmPerSample, opts, defaultBimeraConsensusOptions())
    check perSampleRes.kept.anyIt(it.sampleId == "S2" and it.sequence == exactChimera and it.abundance == 2)
    check perSampleRes.kept.allIt(not (it.sampleId == "S1" and it.sequence == exactChimera))
    check perSampleRes.exactBimeraCount == 1

  test "removeBimerDenovo CLI processes long table":
    let dir = rbMakeTempDir("amplidada-remove-bimera-cli")
    let inPath = dir / "asv_long.tsv"
    let outPath = dir / "asv_long.nobim.tsv"
    let summaryPath = dir / "summary.txt"
    let binPath = dir / "removeBimerDenovo"

    rbWriteLongTsv(inPath, @[
      ("S1", parentA, 50'i64),
      ("S1", parentB, 40'i64),
      ("S1", exactChimera, 10'i64),
      ("S2", parentA, 4'i64),
      ("S2", parentB, 3'i64),
      ("S2", exactChimera, 2'i64)
    ])

    rbCompileCli("src/cli/removeBimerDenovo.nim", binPath)
    let cmd = quoteShell(binPath) &
      " --input " & quoteShell(inPath) &
      " --output " & quoteShell(outPath) &
      " --summary " & quoteShell(summaryPath) &
      " --method consensus --threads 1 --quiet"
    let run = execCmdEx(cmd)
    check run.exitCode == 0
    check fileExists(outPath)
    check fileExists(summaryPath)

    let outText = readFile(outPath)
    check not outText.contains(exactChimera)
    let summaryText = readFile(summaryPath)
    check summaryText.contains("method=bmConsensus")
    check summaryText.contains("removed_unique_sequences=1")
    check summaryText.contains("bimera_exact_count=1")
    check summaryText.contains("bimera_one_off_mismatch_count=0")
    check summaryText.contains("bimera_one_off_indel_count=0")
