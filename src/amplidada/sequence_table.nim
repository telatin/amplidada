import std/[algorithm, hashes, strformat, tables]
import amplidada/[dada, merge_pairs]

type
  OrderBy* = enum
    obAbundance,   ## Order columns by decreasing total abundance
    obNSamples,    ## Order columns by decreasing number of samples present in
    obNone         ## Keep first-seen order

  SequenceTableOptions* = object
    orderBy*: OrderBy
    minAbundance*: int64  ## Filter: exclude sequences with total abundance below this
    minSamples*: int      ## Filter: exclude sequences present in fewer samples

  SampleAsvs* = object
    name*: string
    asvs*: seq[AsvRecord]

  SampleMerged* = object
    name*: string
    merged*: seq[MergedSequence]

  SequenceTable* = object
    samples*: seq[string]
    sequences*: seq[string]
    counts*: seq[seq[int64]]  ## counts[sampleIdx][seqIdx]

proc defaultSequenceTableOptions*(): SequenceTableOptions =
  SequenceTableOptions(
    orderBy: obAbundance,
    minAbundance: 0,
    minSamples: 0
  )

proc validateOptions*(opts: SequenceTableOptions) =
  if opts.minAbundance < 0:
    raise newException(ValueError, "minAbundance must be >= 0")
  if opts.minSamples < 0:
    raise newException(ValueError, "minSamples must be >= 0")

proc describe*(opts: SequenceTableOptions): string =
  validateOptions(opts)
  let orderStr = case opts.orderBy
    of obAbundance: "abundance"
    of obNSamples: "nsamples"
    of obNone: "none"
  &"orderBy={orderStr} minAbundance={opts.minAbundance} minSamples={opts.minSamples}"

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

proc md5Hex*(s: string): string =
  ## Simple hash for sequence identifiers (not cryptographically secure,
  ## but sufficient for ASV IDs). Uses Nim's built-in hash function.
  let h = hash(s)
  result = &"{h:08x}"

proc makeTableFromAsvs(
  samples: seq[SampleAsvs],
  opts: SequenceTableOptions
): SequenceTable =
  ## Build a sequence table from denoised ASV records (single-end or pooled).
  ##
  ## Each sample contributes a set of (sequence, abundance) pairs. The result
  ## is a matrix with samples as rows and unique sequences as columns.

  # First pass: collect all unique sequences and their per-sample counts
  var seqIndex: Table[string, int] = initTable[string, int]()
  var sequences: seq[string] = @[]

  for sample in samples:
    for asv in sample.asvs:
      if not seqIndex.hasKey(asv.sequence):
        seqIndex[asv.sequence] = sequences.len
        sequences.add(asv.sequence)

  let nSamples = samples.len
  let nSeqs = sequences.len

  if nSeqs == 0:
    result = SequenceTable(
      samples: @[],
      sequences: @[],
      counts: @[]
    )
    return

  # Build count matrix
  var counts = newSeqOfCap[seq[int64]](nSamples)
  for i in 0..<nSamples:
    var row = newSeqOfCap[int64](nSeqs)
    for j in 0..<nSeqs:
      row.add(0'i64)
    counts.add(row)

  for i, sample in samples:
    for asv in sample.asvs:
      let j = seqIndex[asv.sequence]
      counts[i][j] += asv.abundance

  # Apply filters
  var keepSeqs = newSeqOfCap[bool](nSeqs)
  for j in 0..<nSeqs:
    keepSeqs.add(true)

  if opts.minAbundance > 0:
    for j in 0..<nSeqs:
      var total = 0'i64
      for i in 0..<nSamples:
        total += counts[i][j]
      if total < opts.minAbundance:
        keepSeqs[j] = false

  if opts.minSamples > 0:
    for j in 0..<nSeqs:
      var present = 0
      for i in 0..<nSamples:
        if counts[i][j] > 0:
          inc present
      if present < opts.minSamples:
        keepSeqs[j] = false

  # Determine column order
  var colOrder: seq[int]
  case opts.orderBy
  of obNone:
    colOrder = @[]
    for j in 0..<nSeqs:
      if keepSeqs[j]:
        colOrder.add(j)
  of obAbundance:
    var totals: seq[(int, int64)] = @[]
    for j in 0..<nSeqs:
      if not keepSeqs[j]: continue
      var total = 0'i64
      for i in 0..<nSamples:
        total += counts[i][j]
      totals.add((j, total))
    totals.sort(proc(a, b: (int, int64)): int = cmp(b[1], a[1]))
    colOrder = newSeqOfCap[int](totals.len)
    for t in totals:
      colOrder.add(t[0])
  of obNSamples:
    var presents: seq[(int, int)] = @[]
    for j in 0..<nSeqs:
      if not keepSeqs[j]: continue
      var present = 0
      for i in 0..<nSamples:
        if counts[i][j] > 0:
          inc present
      presents.add((j, present))
    presents.sort(proc(a, b: (int, int)): int = cmp(b[1], a[1]))
    colOrder = newSeqOfCap[int](presents.len)
    for p in presents:
      colOrder.add(p[0])

  # Build final table
  result.samples = newSeqOfCap[string](nSamples)
  for sample in samples:
    result.samples.add(sample.name)

  result.sequences = newSeqOfCap[string](colOrder.len)
  for j in colOrder:
    result.sequences.add(sequences[j])

  result.counts = newSeqOfCap[seq[int64]](nSamples)
  for i in 0..<nSamples:
    var row = newSeqOfCap[int64](colOrder.len)
    for j in colOrder:
      row.add(counts[i][j])
    result.counts.add(row)

proc makeTableFromMerged(
  samples: seq[SampleMerged],
  opts: SequenceTableOptions
): SequenceTable =
  ## Build a sequence table from merged paired-end ASVs.
  ##
  ## Same logic as makeTableFromAsvs but operates on MergedSequence records.

  var seqIndex: Table[string, int] = initTable[string, int]()
  var sequences: seq[string] = @[]

  for sample in samples:
    for merged in sample.merged:
      if not seqIndex.hasKey(merged.sequence):
        seqIndex[merged.sequence] = sequences.len
        sequences.add(merged.sequence)

  let nSamples = samples.len
  let nSeqs = sequences.len

  if nSeqs == 0:
    result = SequenceTable(
      samples: @[],
      sequences: @[],
      counts: @[]
    )
    return

  var counts = newSeqOfCap[seq[int64]](nSamples)
  for i in 0..<nSamples:
    var row = newSeqOfCap[int64](nSeqs)
    for j in 0..<nSeqs:
      row.add(0'i64)
    counts.add(row)

  for i, sample in samples:
    for merged in sample.merged:
      let j = seqIndex[merged.sequence]
      counts[i][j] += merged.abundance

  var keepSeqs = newSeqOfCap[bool](nSeqs)
  for j in 0..<nSeqs:
    keepSeqs.add(true)

  if opts.minAbundance > 0:
    for j in 0..<nSeqs:
      var total = 0'i64
      for i in 0..<nSamples:
        total += counts[i][j]
      if total < opts.minAbundance:
        keepSeqs[j] = false

  if opts.minSamples > 0:
    for j in 0..<nSeqs:
      var present = 0
      for i in 0..<nSamples:
        if counts[i][j] > 0:
          inc present
      if present < opts.minSamples:
        keepSeqs[j] = false

  var colOrder: seq[int]
  case opts.orderBy
  of obNone:
    colOrder = @[]
    for j in 0..<nSeqs:
      if keepSeqs[j]:
        colOrder.add(j)
  of obAbundance:
    var totals: seq[(int, int64)] = @[]
    for j in 0..<nSeqs:
      if not keepSeqs[j]: continue
      var total = 0'i64
      for i in 0..<nSamples:
        total += counts[i][j]
      totals.add((j, total))
    totals.sort(proc(a, b: (int, int64)): int = cmp(b[1], a[1]))
    colOrder = newSeqOfCap[int](totals.len)
    for t in totals:
      colOrder.add(t[0])
  of obNSamples:
    var presents: seq[(int, int)] = @[]
    for j in 0..<nSeqs:
      if not keepSeqs[j]: continue
      var present = 0
      for i in 0..<nSamples:
        if counts[i][j] > 0:
          inc present
      presents.add((j, present))
    presents.sort(proc(a, b: (int, int)): int = cmp(b[1], a[1]))
    colOrder = newSeqOfCap[int](presents.len)
    for p in presents:
      colOrder.add(p[0])

  result.samples = newSeqOfCap[string](nSamples)
  for sample in samples:
    result.samples.add(sample.name)

  result.sequences = newSeqOfCap[string](colOrder.len)
  for j in colOrder:
    result.sequences.add(sequences[j])

  result.counts = newSeqOfCap[seq[int64]](nSamples)
  for i in 0..<nSamples:
    var row = newSeqOfCap[int64](colOrder.len)
    for j in colOrder:
      row.add(counts[i][j])
    result.counts.add(row)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

proc makeSequenceTable*(
  samples: seq[SampleAsvs],
  opts: SequenceTableOptions = defaultSequenceTableOptions()
): SequenceTable =
  ## Build a sequence table from denoised ASV records.
  ##
  ## Parameters:
  ##   samples: One SampleAsvs per input sample, with name and ASV records
  ##   opts: Ordering and filtering options
  ##
  ## Returns a SequenceTable with samples as rows and sequences as columns.
  ## By default, columns are ordered by decreasing total abundance.
  validateOptions(opts)
  result = makeTableFromAsvs(samples, opts)

proc makeSequenceTable*(
  samples: seq[SampleMerged],
  opts: SequenceTableOptions = defaultSequenceTableOptions()
): SequenceTable =
  ## Build a sequence table from merged paired-end ASV records.
  ##
  ## Parameters:
  ##   samples: One SampleMerged per input sample, with name and merged records
  ##   opts: Ordering and filtering options
  ##
  ## Returns a SequenceTable with samples as rows and sequences as columns.
  ## By default, columns are ordered by decreasing total abundance.
  validateOptions(opts)
  result = makeTableFromMerged(samples, opts)

# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

proc writeTsv*(path: string, table: SequenceTable) =
  ## Write a sequence table to a TSV file.
  ##
  ## Format:
  ##   SampleID\tseq1\tseq2\t...
  ##   sample1\tcount11\tcount12\t...
  ##   sample2\tcount21\tcount22\t...
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open sequence table output file: " & path)
  defer:
    close(f)

  # Header row
  f.write("SampleID")
  for seq in table.sequences:
    f.write("\t" & seq)
  f.write("\n")

  # Data rows
  for i, sample in table.samples:
    f.write(sample)
    for count in table.counts[i]:
      f.write("\t" & $count)
    f.write("\n")

proc writeTsvWithIds*(path: string, table: SequenceTable) =
  ## Write a sequence table to a TSV file with ASV IDs (md5-based).
  ##
  ## Format:
  ##   SampleID\tASV_ID_1\tASV_ID_2\t...
  ##   sample1\tcount11\tcount12\t...
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open sequence table output file: " & path)
  defer:
    close(f)

  # Header row with ASV IDs
  f.write("SampleID")
  for seq in table.sequences:
    f.write("\t" & md5Hex(seq))
  f.write("\n")

  # Data rows
  for i, sample in table.samples:
    f.write(sample)
    for count in table.counts[i]:
      f.write("\t" & $count)
    f.write("\n")

proc writeFasta*(path: string, table: SequenceTable) =
  ## Write sequences from a sequence table to a FASTA file.
  ##
  ## Each sequence is written with an ASV ID header and total abundance.
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open FASTA output file: " & path)
  defer:
    close(f)

  for j, seq in table.sequences:
    var total = 0'i64
    for i in 0..<table.samples.len:
      total += table.counts[i][j]
    f.writeLine(">" & md5Hex(seq) & " totalcounts=" & $total)
    f.writeLine(seq)

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

proc totalReads*(table: SequenceTable): int64 =
  ## Total number of reads in the table.
  result = 0'i64
  for row in table.counts:
    for count in row:
      result += count

proc totalUniqueSequences*(table: SequenceTable): int =
  ## Number of unique sequences (columns).
  result = table.sequences.len

proc samplesPresent*(table: SequenceTable, seqIdx: int): int =
  ## Number of samples where a sequence is present (count > 0).
  result = 0
  for i in 0..<table.samples.len:
    if table.counts[i][seqIdx] > 0:
      inc result
