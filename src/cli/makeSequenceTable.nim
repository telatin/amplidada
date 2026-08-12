import std/[algorithm, os, parseopt, strformat, strutils, tables]
import amplidada/[dada, merge_pairs, sequence_table, version]

proc usage() =
  echo "Usage:"
  echo "  makeSequenceTable --in <dir> --out <table.tsv> [options]"
  echo "  makeSequenceTable --in-file <path>... --out <table.tsv> [options]"
  echo ""
  echo "Build a sample × sequence count table from per-sample ASV outputs."
  echo ""
  echo "Inputs:"
  echo "  --in <dir>                       Directory with per-sample TSV files"
  echo "  --in-file <path>                 TSV file (can be repeated)"
  echo "  --input-format <format>          Format: dada2 (default)"
  echo ""
  echo "Options:"
  echo "  --order <abundance|nsamples|none>  Column ordering (default: abundance)"
  echo "  --min-abundance <int>            Filter: min total abundance (default: 0)"
  echo "  --min-samples <int>              Filter: min samples present (default: 0)"
  echo "  --with-ids                       Write ASV IDs (hashes) instead of sequences"
  echo "  --fasta <path>                   Also write sequences to FASTA file"
  echo "  --verbose                        Print progress"
  echo ""
  echo "General:"
  echo "  -h, --help"
  echo "  -v, --version"

proc fail(msg: string) {.noreturn.} =
  echo "Error: " & msg
  quit(1)

proc parseSampleName(path: string): string =
  ## Extract sample name from filename like "sample1.merged.tsv" or "sample1.asv.tsv"
  let p = splitFile(path)
  var name = p.name
  # Remove common suffixes
  if name.endsWith(".merged"):
    name = name[0..^8]
  elif name.endsWith(".asv"):
    name = name[0..^5]
  result = name

proc readAsvsFromTsv(path: string): tuple[name: string, asvs: seq[AsvRecord]] =
  ## Read ASV records from a dada2-style TSV file.
  ## Expected format: asv_id\tsequence\tabundance\t...
  result.name = parseSampleName(path)

  var f: File
  if not open(f, path, fmRead):
    raise newException(IOError, "Cannot open file: " & path)
  defer: close(f)

  var line: string
  var lineNum = 0
  while f.readLine(line):
    inc lineNum
    if lineNum == 1:
      continue  # Skip header
    let parts = line.split('\t')
    if parts.len >= 3:
      let seq = parts[1]
      let abundance = parseInt(parts[2]).int64
      result.asvs.add(AsvRecord(
        sequence: seq,
        abundance: abundance,
        clusterSize: 1,
        centerUniqueIndex: 0,
        firstSeenReadIndex: 0
      ))

proc readMergedFromTsv(path: string): tuple[name: string, merged: seq[MergedSequence]] =
  ## Read merged sequences from a dada2-style TSV file.
  ## Expected format: merged_id\tsequence\tabundance
  result.name = parseSampleName(path)

  var f: File
  if not open(f, path, fmRead):
    raise newException(IOError, "Cannot open file: " & path)
  defer: close(f)

  var line: string
  var lineNum = 0
  while f.readLine(line):
    inc lineNum
    if lineNum == 1:
      continue  # Skip header
    let parts = line.split('\t')
    if parts.len >= 3:
      let seq = parts[1]
      let abundance = parseInt(parts[2]).int64
      result.merged.add(MergedSequence(
        sequence: seq,
        abundance: abundance
      ))

proc main() =
  var inputDir = ""
  var inputFiles: seq[string] = @[]
  var outPath = ""
  var orderBy = obAbundance
  var minAbundance = 0'i64
  var minSamples = 0
  var withIds = false
  var fastaPath = ""
  var verbose = false

  var parser = initOptParser(commandLineParams())

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
        of "in":
          inputDir = value
        of "in-file":
          inputFiles.add(value)
        of "out":
          outPath = value
        of "order":
          case value.toLowerAscii()
          of "abundance": orderBy = obAbundance
          of "nsamples", "n-samples": orderBy = obNSamples
          of "none": orderBy = obNone
          else: fail("Invalid value for --order: " & value)
        of "min-abundance":
          minAbundance = parseInt(value).int64
        of "min-samples":
          minSamples = parseInt(value)
        of "with-ids":
          withIds = true
        of "fasta":
          fastaPath = value
        of "verbose":
          verbose = true
        else:
          fail("Unknown option: --" & key)
      of cmdArgument:
        fail("Unexpected positional argument: " & key)
      of cmdEnd:
        discard
  except ValueError as e:
    fail(e.msg)

  if outPath.len == 0:
    fail("Missing required option: --out")

  # Collect input files
  var tsvFiles: seq[string] = @[]

  if inputDir.len > 0:
    if not inputDir.dirExists():
      fail("Input directory does not exist: " & inputDir)
    for path in walkFiles(inputDir / "*.tsv"):
      tsvFiles.add(path)
    if tsvFiles.len == 0:
      fail("No TSV files found in: " & inputDir)
    tsvFiles.sort()
    if verbose:
      echo &"[makeSequenceTable] Found {tsvFiles.len} files in {inputDir}"

  for f in inputFiles:
    if not f.fileExists():
      fail("Input file does not exist: " & f)
    tsvFiles.add(f)

  if tsvFiles.len == 0:
    fail("No input files provided. Use --in <dir> or --in-file <path>")

  # Read all samples
  var samples: seq[SampleAsvs] = @[]

  for path in tsvFiles:
    try:
      let (name, asvs) = readAsvsFromTsv(path)
      if verbose:
        echo &"[makeSequenceTable] {name}: {asvs.len} ASVs from {path}"
      samples.add(SampleAsvs(name: name, asvs: asvs))
    except:
      # Try merged format
      try:
        let (name, merged) = readMergedFromTsv(path)
        # Convert merged to ASVs for uniform handling
        var asvs: seq[AsvRecord] = @[]
        for m in merged:
          asvs.add(AsvRecord(
            sequence: m.sequence,
            abundance: m.abundance,
            clusterSize: 1,
            centerUniqueIndex: 0,
            firstSeenReadIndex: 0
          ))
        if verbose:
          echo &"[makeSequenceTable] {name}: {merged.len} merged from {path}"
        samples.add(SampleAsvs(name: name, asvs: asvs))
      except:
        echo &"[WARN] Could not read {path}: {getCurrentException().msg}"

  if samples.len == 0:
    fail("No valid samples found in input files")

  let opts = SequenceTableOptions(
    orderBy: orderBy,
    minAbundance: minAbundance,
    minSamples: minSamples
  )

  let table = makeSequenceTable(samples, opts)

  if verbose:
    echo &"[makeSequenceTable] Building table: {table.samples.len} samples × {table.sequences.len} sequences"
    echo &"[makeSequenceTable] Total reads: {table.totalReads()}"

  if withIds:
    writeTsvWithIds(outPath, table)
  else:
    writeTsv(outPath, table)

  if fastaPath.len > 0:
    writeFasta(fastaPath, table)

  echo &"[makeSequenceTable] Done. Output: {outPath}"
  echo &"[makeSequenceTable] {table.samples.len} samples, {table.sequences.len} sequences, {table.totalReads()} total reads"

when isMainModule:
  main()
