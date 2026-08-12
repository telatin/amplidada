import std/[unittest, strformat]
import amplidada/[sequence_table, dada, merge_pairs]

suite "SequenceTable":

  test "makeSequenceTable from ASVs - basic":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "TGCA", abundance: 50, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      ),
      SampleAsvs(
        name: "sample2",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 80, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "GGGG", abundance: 30, clusterSize: 1, centerUniqueIndex: 2, firstSeenReadIndex: 2)
        ]
      )
    ]

    let table = makeSequenceTable(samples)

    check table.samples == @["sample1", "sample2"]
    check table.sequences.len == 3  # ACGT, TGCA, GGGG

    # ACGT should be first (total abundance = 180)
    check table.sequences[0] == "ACGT"

    # Find indices
    var acgtIdx = -1
    var tgcaIdx = -1
    var ggggIdx = -1
    for i, seq in table.sequences:
      case seq
      of "ACGT": acgtIdx = i
      of "TGCA": tgcaIdx = i
      of "GGGG": ggggIdx = i

    check table.counts[0][acgtIdx] == 100  # sample1, ACGT
    check table.counts[0][tgcaIdx] == 50   # sample1, TGCA
    check table.counts[0][ggggIdx] == 0    # sample1, GGGG (not present)
    check table.counts[1][acgtIdx] == 80   # sample2, ACGT
    check table.counts[1][ggggIdx] == 30   # sample2, GGGG

  test "makeSequenceTable from ASVs - ordering by abundance":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 10, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "CCCC", abundance: 100, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      )
    ]

    let table = makeSequenceTable(samples, SequenceTableOptions(orderBy: obAbundance))
    check table.sequences[0] == "CCCC"  # Most abundant first
    check table.sequences[1] == "AAAA"

  test "makeSequenceTable from ASVs - ordering by NSamples":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 10, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "CCCC", abundance: 100, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      ),
      SampleAsvs(
        name: "sample2",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 5, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      )
    ]

    let table = makeSequenceTable(samples, SequenceTableOptions(orderBy: obNSamples))
    check table.sequences[0] == "AAAA"  # Present in 2 samples
    check table.sequences[1] == "CCCC"  # Present in 1 sample

  test "makeSequenceTable from ASVs - minAbundance filter":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "CCCC", abundance: 5, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      )
    ]

    let table = makeSequenceTable(samples, SequenceTableOptions(minAbundance: 10))
    check table.sequences.len == 1
    check table.sequences[0] == "AAAA"

  test "makeSequenceTable from ASVs - minSamples filter":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "CCCC", abundance: 50, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      ),
      SampleAsvs(
        name: "sample2",
        asvs: @[
          AsvRecord(sequence: "AAAA", abundance: 80, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      )
    ]

    let table = makeSequenceTable(samples, SequenceTableOptions(minSamples: 2))
    check table.sequences.len == 1
    check table.sequences[0] == "AAAA"

  test "makeSequenceTable from merged sequences":
    let samples = @[
      SampleMerged(
        name: "sample1",
        merged: @[
          MergedSequence(sequence: "ACGTACGT", abundance: 100),
          MergedSequence(sequence: "TGCA", abundance: 50)
        ]
      ),
      SampleMerged(
        name: "sample2",
        merged: @[
          MergedSequence(sequence: "ACGTACGT", abundance: 80),
          MergedSequence(sequence: "GGGG", abundance: 30)
        ]
      )
    ]

    let table = makeSequenceTable(samples)
    check table.samples == @["sample1", "sample2"]
    check table.sequences.len == 3

  test "writeTsv - basic output":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      )
    ]

    let table = makeSequenceTable(samples)
    let path = "tests/tmp_seqtab.tsv"

    writeTsv(path, table)

    # Read back and verify format
    var f: File
    check open(f, path, fmRead)
    var lines: seq[string] = @[]
    var line: string
    while f.readLine(line):
      lines.add(line)
    close(f)

    check lines.len == 2  # header + 1 data row
    check lines[0].startsWith("SampleID\t")
    check lines[1].startsWith("sample1\t")

    # Cleanup
    removeFile(path)

  test "totalReads":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0),
          AsvRecord(sequence: "TGCA", abundance: 50, clusterSize: 1, centerUniqueIndex: 1, firstSeenReadIndex: 1)
        ]
      ),
      SampleAsvs(
        name: "sample2",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 80, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      )
    ]

    let table = makeSequenceTable(samples)
    check table.totalReads() == 230  # 100 + 50 + 80

  test "samplesPresent":
    let samples = @[
      SampleAsvs(
        name: "sample1",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 100, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      ),
      SampleAsvs(
        name: "sample2",
        asvs: @[
          AsvRecord(sequence: "ACGT", abundance: 80, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      ),
      SampleAsvs(
        name: "sample3",
        asvs: @[]  # No ACGT in sample3
      )
    ]

    let table = makeSequenceTable(samples)
    check table.sequences.len == 1
    check table.samplesPresent(0) == 2  # ACGT present in sample1 and sample2

  test "empty input":
    let samples: seq[SampleAsvs] = @[]
    let table = makeSequenceTable(samples)
    check table.samples.len == 0
    check table.sequences.len == 0
    check table.counts.len == 0

  test "single sample single ASV":
    let samples = @[
      SampleAsvs(
        name: "only_sample",
        asvs: @[
          AsvRecord(sequence: "ACGTACGTACGT", abundance: 42, clusterSize: 1, centerUniqueIndex: 0, firstSeenReadIndex: 0)
        ]
      )
    ]

    let table = makeSequenceTable(samples)
    check table.samples == @["only_sample"]
    check table.sequences == @["ACGTACGTACGT"]
    check table.counts == @[@[42'i64]]
