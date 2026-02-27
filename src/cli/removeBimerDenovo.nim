import std/[os, parseopt, sets, strformat, strutils]
import amplidada

type
  InputFormat = enum
    ifAuto,
    ifLong,
    ifUnique

  ParsedLongInput = object
    rows: seq[SampleSequenceAbundance]
    uniqueCount: int

  ParsedUniqueInput = object
    rows: seq[UniqueSequenceAbundance]

proc usage() =
  echo "Usage:"
  echo "  removeBimerDenovo --input <input.tsv> --output <output.tsv> [options]"
  echo ""
  echo "Input format:"
  echo "  --input-format <auto|long|unique>  Input TSV format (default: auto)"
  echo "    long:   SampleID, Sequence, Abundance"
  echo "    unique: Sequence, Abundance (or any TSV with these columns)"
  echo ""
  echo "Algorithm:"
  echo "  --method <pooled|consensus|per-sample>   Default: consensus"
  echo "  --min-fold-parent-over-abundance <float>"
  echo "  --min-parent-abundance <int>"
  echo "  --allow-one-off"
  echo "  --min-one-off-parent-distance <int>"
  echo "  --min-sample-fraction <float>            For consensus mode"
  echo "  --ignore-n-negatives <int>               For consensus mode"
  echo "  --threads <int>"
  echo "  --parallel-min-candidates <int>"
  echo ""
  echo "Reporting:"
  echo "  --summary <path>               Write summary file"
  echo "  --verbose                      Print progress"
  echo "  --quiet                        Suppress summary output"
  echo ""
  echo "General:"
  echo "  -h, --help"
  echo "  -v, --version"

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Error: " & msg)
  quit(2)

proc parseIntFlag(flag, value: string): int =
  try:
    parseInt(value)
  except ValueError:
    raise newException(ValueError, "Invalid integer for --" & flag & ": " & value)

proc parseFloatFlag(flag, value: string): float64 =
  try:
    parseFloat(value)
  except ValueError:
    raise newException(ValueError, "Invalid float for --" & flag & ": " & value)

proc normalizeArgs(args: seq[string]): seq[string] =
  let longNeedsValue = [
    "--input", "--output", "--input-format", "--method", "--summary",
    "--min-fold-parent-over-abundance", "--min-parent-abundance",
    "--min-one-off-parent-distance", "--min-sample-fraction", "--ignore-n-negatives",
    "--threads", "--parallel-min-candidates"
  ]

  var i = 0
  while i < args.len:
    let a = args[i]
    if a in longNeedsValue and i + 1 < args.len and (args[i + 1].len == 0 or args[i + 1][0] != '-'):
      result.add(a & "=" & args[i + 1])
      inc i
    else:
      result.add(a)
    inc i

proc splitTsv(line: string): seq[string] =
  line.split('\t')

proc toDnaUpperLocal(sequence: string): string =
  result = newString(sequence.len)
  for i, c in sequence:
    result[i] =
      case c
      of 'a'..'z': chr(ord(c) - 32)
      else: c

proc findColumn(header: seq[string], names: openArray[string]): int =
  for i, h in header:
    let col = h.strip().toLowerAscii()
    for n in names:
      if col == n:
        return i
  -1

proc detectInputFormat(path: string): InputFormat =
  var f: File
  if not open(f, path, fmRead):
    raise newException(IOError, "Unable to open input: " & path)
  defer:
    close(f)

  if endOfFile(f):
    raise newException(ValueError, "Input file is empty: " & path)
  let header = splitTsv(readLine(f))
  let idxSample = findColumn(header, ["sampleid", "sample", "sample_id"])
  let idxSeq = findColumn(header, ["sequence", "seq"])
  let idxAbd = findColumn(header, ["abundance", "count", "counts"])

  if idxSeq < 0 or idxAbd < 0:
    raise newException(ValueError, "Input must contain Sequence and Abundance columns")
  if idxSample >= 0:
    return ifLong
  ifUnique

proc parseLongInput(path: string): ParsedLongInput =
  var f: File
  if not open(f, path, fmRead):
    raise newException(IOError, "Unable to open input: " & path)
  defer:
    close(f)

  if endOfFile(f):
    raise newException(ValueError, "Input file is empty: " & path)
  let header = splitTsv(readLine(f))
  let idxSample = findColumn(header, ["sampleid", "sample", "sample_id"])
  let idxSeq = findColumn(header, ["sequence", "seq"])
  let idxAbd = findColumn(header, ["abundance", "count", "counts"])
  if idxSample < 0 or idxSeq < 0 or idxAbd < 0:
    raise newException(ValueError, "Long format requires SampleID, Sequence, Abundance columns")

  var seen = initHashSet[string]()
  while not endOfFile(f):
    let raw = readLine(f).strip()
    if raw.len == 0:
      continue
    let cols = splitTsv(raw)
    if max(idxSample, max(idxSeq, idxAbd)) >= cols.len:
      raise newException(ValueError, "Malformed TSV row: " & raw)
    let sampleId = cols[idxSample].strip()
    let sequence = cols[idxSeq].strip()
    let abundance = parseBiggestInt(cols[idxAbd].strip()).int64
    if sampleId.len == 0 or sequence.len == 0 or abundance <= 0:
      continue
    result.rows.add(SampleSequenceAbundance(sampleId: sampleId, sequence: sequence, abundance: abundance))
    seen.incl(toDnaUpperLocal(sequence))

  result.uniqueCount = seen.len

proc parseUniqueInput(path: string): ParsedUniqueInput =
  var f: File
  if not open(f, path, fmRead):
    raise newException(IOError, "Unable to open input: " & path)
  defer:
    close(f)

  if endOfFile(f):
    raise newException(ValueError, "Input file is empty: " & path)
  let header = splitTsv(readLine(f))
  let idxSeq = findColumn(header, ["sequence", "seq"])
  let idxAbd = findColumn(header, ["abundance", "count", "counts"])
  if idxSeq < 0 or idxAbd < 0:
    raise newException(ValueError, "Unique format requires Sequence and Abundance columns")

  while not endOfFile(f):
    let raw = readLine(f).strip()
    if raw.len == 0:
      continue
    let cols = splitTsv(raw)
    if max(idxSeq, idxAbd) >= cols.len:
      raise newException(ValueError, "Malformed TSV row: " & raw)
    let sequence = cols[idxSeq].strip()
    let abundance = parseBiggestInt(cols[idxAbd].strip()).int64
    if sequence.len == 0 or abundance <= 0:
      continue
    result.rows.add(UniqueSequenceAbundance(sequence: sequence, abundance: abundance))

proc writeSummary(
  path: string,
  mode: BimeraMethod,
  inputFormat: InputFormat,
  inputRows: int,
  inputUnique: int,
  outputRows: int,
  outputUnique: int,
  removedUnique: int,
  exactBimeraCount: int,
  oneOffMismatchBimeraCount: int,
  oneOffIndelBimeraCount: int,
  bimeraOpts: BimeraOptions,
  consensusOpts: BimeraConsensusOptions,
  uniqueRes: BimeraDenovoResult,
  tableRes: BimeraTableResult
) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open summary output: " & path)
  defer:
    close(f)

  f.writeLine("input_rows=" & $inputRows)
  f.writeLine("input_unique_sequences=" & $inputUnique)
  f.writeLine("output_rows=" & $outputRows)
  f.writeLine("output_unique_sequences=" & $outputUnique)
  f.writeLine("removed_unique_sequences=" & $removedUnique)
  f.writeLine("input_format=" & $inputFormat)
  f.writeLine("method=" & $mode)
  f.writeLine("bimera_options=" & bimeraOpts.describe())
  f.writeLine("consensus_options=" & consensusOpts.describe())
  f.writeLine("bimera_exact_count=" & $exactBimeraCount)
  f.writeLine("bimera_one_off_mismatch_count=" & $oneOffMismatchBimeraCount)
  f.writeLine("bimera_one_off_indel_count=" & $oneOffIndelBimeraCount)
  f.writeLine("bimera_total_count=" & $(exactBimeraCount + oneOffMismatchBimeraCount + oneOffIndelBimeraCount))
  f.writeLine("pooled_bimera_count=" & $uniqueRes.nBimeras)
  f.writeLine("consensus_bimera_count=" & $tableRes.nBimeras)

proc uniqueCountFromRows(rows: openArray[SampleSequenceAbundance]): int =
  var seen = initHashSet[string]()
  for row in rows:
    if row.sequence.len > 0:
      seen.incl(toDnaUpperLocal(row.sequence))
  seen.len

proc uniqueCountFromUniques(rows: openArray[UniqueSequenceAbundance]): int =
  var seen = initHashSet[string]()
  for row in rows:
    if row.sequence.len > 0:
      seen.incl(toDnaUpperLocal(row.sequence))
  seen.len

proc main() =
  var inputPath = ""
  var outputPath = ""
  var summaryPath = ""
  var quiet = false
  var verbose = false
  var format = ifAuto
  var mode = bmConsensus

  var opts = defaultBimeraOptions()
  var consensusOpts = defaultBimeraConsensusOptions()

  var minFoldSet = false
  var minParentSet = false

  var parser = initOptParser(normalizeArgs(commandLineParams()))
  try:
    for kind, key, value in parser.getopt():
      case kind
      of cmdLongOption, cmdShortOption:
        let flag = key.toLowerAscii()
        case flag
        of "help", "h":
          usage()
          return
        of "version", "v":
          echo AmpliDadaVersion
          return
        of "input":
          inputPath = value
        of "output":
          outputPath = value
        of "summary":
          summaryPath = value
        of "input-format":
          case value.toLowerAscii()
          of "auto":
            format = ifAuto
          of "long":
            format = ifLong
          of "unique":
            format = ifUnique
          else:
            fail("Invalid value for --input-format: " & value)
        of "method":
          case value.toLowerAscii()
          of "pooled":
            mode = bmPooled
          of "consensus":
            mode = bmConsensus
          of "per-sample", "per_sample", "persample":
            mode = bmPerSample
          else:
            fail("Invalid value for --method: " & value)
        of "min-fold-parent-over-abundance":
          opts.minFoldParentOverAbundance = parseFloatFlag(flag, value)
          minFoldSet = true
        of "min-parent-abundance":
          opts.minParentAbundance = parseIntFlag(flag, value).int64
          minParentSet = true
        of "allow-one-off":
          opts.allowOneOff = true
        of "min-one-off-parent-distance":
          opts.minOneOffParentDistance = parseIntFlag(flag, value)
        of "min-sample-fraction":
          consensusOpts.minSampleFraction = parseFloatFlag(flag, value)
        of "ignore-n-negatives":
          consensusOpts.ignoreNNegatives = parseIntFlag(flag, value)
        of "threads":
          opts.threads = parseIntFlag(flag, value)
        of "parallel-min-candidates":
          opts.parallelMinCandidates = parseIntFlag(flag, value)
        of "verbose":
          verbose = true
        of "quiet":
          quiet = true
        else:
          fail("Unknown option: --" & key)
      of cmdArgument:
        fail("Unexpected positional argument: " & key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  if inputPath.len == 0:
    fail("Missing required option: --input")
  if outputPath.len == 0:
    fail("Missing required option: --output")
  if quiet and verbose:
    fail("--quiet and --verbose are mutually exclusive")

  if not minFoldSet or not minParentSet:
    if mode == bmPooled:
      if not minFoldSet:
        opts.minFoldParentOverAbundance = 2.0
      if not minParentSet:
        opts.minParentAbundance = 8
    else:
      if not minFoldSet:
        opts.minFoldParentOverAbundance = 1.5
      if not minParentSet:
        opts.minParentAbundance = 2

  validateOptions(opts)
  validateOptions(consensusOpts)

  var resolvedFormat = format
  if resolvedFormat == ifAuto:
    resolvedFormat = detectInputFormat(inputPath)
    if verbose:
      echo &"[removeBimerDenovo] Auto-detected input format: {resolvedFormat}"

  try:
    case resolvedFormat
    of ifLong:
      let parsed = parseLongInput(inputPath)
      if parsed.rows.len == 0:
        writeSampleSequenceAbundanceTsv(outputPath, @[])
        if summaryPath.len > 0:
          writeSummary(
            summaryPath, mode, resolvedFormat, 0, 0, 0, 0, 0, 0, 0, 0,
            opts, consensusOpts, BimeraDenovoResult(), BimeraTableResult()
          )
        if not quiet:
          echo "No input rows after parsing."
        return

      if verbose:
        echo &"[removeBimerDenovo] Parsed rows={parsed.rows.len} unique={parsed.uniqueCount}"

      let removeRes = removeBimeraDenovo(parsed.rows, mode, opts, consensusOpts)
      writeSampleSequenceAbundanceTsv(outputPath, removeRes.kept)

      let outputUnique = uniqueCountFromRows(removeRes.kept)
      if summaryPath.len > 0:
        writeSummary(
          path = summaryPath,
          mode = mode,
          inputFormat = resolvedFormat,
          inputRows = parsed.rows.len,
          inputUnique = parsed.uniqueCount,
          outputRows = removeRes.kept.len,
          outputUnique = outputUnique,
          removedUnique = removeRes.removedUniqueCount,
          exactBimeraCount = removeRes.exactBimeraCount,
          oneOffMismatchBimeraCount = removeRes.oneOffMismatchBimeraCount,
          oneOffIndelBimeraCount = removeRes.oneOffIndelBimeraCount,
          bimeraOpts = opts,
          consensusOpts = consensusOpts,
          uniqueRes = removeRes.pooled,
          tableRes = removeRes.consensus
        )

      if not quiet:
        echo &"Method: {mode}"
        echo &"Input rows: {parsed.rows.len} unique sequences: {parsed.uniqueCount}"
        echo &"Output rows: {removeRes.kept.len} unique sequences: {outputUnique}"
        echo &"Removed unique sequences: {removeRes.removedUniqueCount}"
        echo &"Bimeras: exact={removeRes.exactBimeraCount} one_off_mismatch={removeRes.oneOffMismatchBimeraCount} one_off_indel={removeRes.oneOffIndelBimeraCount}"

    of ifUnique:
      if mode == bmPerSample:
        fail("--method=per-sample requires long format input (SampleID, Sequence, Abundance)")

      let parsed = parseUniqueInput(inputPath)
      if parsed.rows.len == 0:
        writeUniqueSequenceAbundanceTsv(outputPath, @[])
        if summaryPath.len > 0:
          writeSummary(
            summaryPath, mode, resolvedFormat, 0, 0, 0, 0, 0, 0, 0, 0,
            opts, consensusOpts, BimeraDenovoResult(), BimeraTableResult()
          )
        if not quiet:
          echo "No input rows after parsing."
        return

      if verbose:
        echo &"[removeBimerDenovo] Parsed pooled rows={parsed.rows.len}"

      let denovo = isBimeraDenovo(parsed.rows, opts)
      var kept: seq[UniqueSequenceAbundance] = @[]
      for i, sequence in denovo.sequences:
        if not denovo.flags[i]:
          kept.add(UniqueSequenceAbundance(sequence: sequence, abundance: denovo.abundances[i]))
      writeUniqueSequenceAbundanceTsv(outputPath, kept)

      let inUnique = uniqueCountFromUniques(parsed.rows)
      let outUnique = uniqueCountFromUniques(kept)
      if summaryPath.len > 0:
        writeSummary(
          path = summaryPath,
          mode = bmPooled,
          inputFormat = resolvedFormat,
          inputRows = parsed.rows.len,
          inputUnique = inUnique,
          outputRows = kept.len,
          outputUnique = outUnique,
          removedUnique = max(0, inUnique - outUnique),
          exactBimeraCount = denovo.exactBimeraCount,
          oneOffMismatchBimeraCount = denovo.oneOffMismatchBimeraCount,
          oneOffIndelBimeraCount = denovo.oneOffIndelBimeraCount,
          bimeraOpts = opts,
          consensusOpts = consensusOpts,
          uniqueRes = denovo,
          tableRes = BimeraTableResult()
        )

      if not quiet:
        echo "Method: pooled"
        echo &"Input rows: {parsed.rows.len} unique sequences: {inUnique}"
        echo &"Output rows: {kept.len} unique sequences: {outUnique}"
        echo &"Removed unique sequences: {max(0, inUnique - outUnique)}"
        echo &"Bimeras: exact={denovo.exactBimeraCount} one_off_mismatch={denovo.oneOffMismatchBimeraCount} one_off_indel={denovo.oneOffIndelBimeraCount}"

    of ifAuto:
      fail("Internal error: unresolved input format")
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
