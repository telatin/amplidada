import std/[os, parseopt, sets, strformat, strutils]
import amplidada

proc usage() =
  echo "Usage: derepFastq --in <reads.fastq[.gz]> --out <uniques.tsv> [options]"
  echo ""
  echo "Core options:"
  echo "  --in <path>                   Input FASTQ/FASTQ.GZ"
  echo "  --out <path>                  Output TSV with unique sequences"
  echo "  --map-out <path>              Optional read-to-unique mapping TSV"
  echo "  --order <abundance|first>     Ordering of unique sequences"
  echo "  --no-map                      Skip readToUnique map computation"
  echo ""
  echo "Performance:"
  echo "  --threads <int>               Worker threads"
  echo "  --chunk-size <int>            Reads per chunk"
  echo "  --parallel-min-chunk <int>    Minimum chunk size for multithreading"
  echo ""
  echo "Output/reporting:"
  echo "  --summary <path>              Optional run summary"
  echo "  --verbose                     Print progress details"
  echo "  --quiet                       Suppress summary output"
  echo ""
  echo "General:"
  echo "  -h, --help"
  echo "  -v, --version"

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Error: " & msg)
  quit(2)

proc parseIntFlag(flag, value: string): int =
  try:
    result = parseInt(value)
  except ValueError:
    raise newException(ValueError, "Invalid integer for --" & flag & ": " & value)

proc normalizeArgs(args: seq[string]): seq[string] =
  let longNeedsValue = toHashSet([
    "--in", "--out", "--map-out", "--order", "--threads",
    "--chunk-size", "--parallel-min-chunk", "--summary"
  ])

  var i = 0
  while i < args.len:
    let a = args[i]
    if a in longNeedsValue and i + 1 < args.len and (args[i + 1].len == 0 or args[i + 1][0] != '-'):
      result.add(a & "=" & args[i + 1])
      inc i
    else:
      result.add(a)
    inc i

proc meanQualToString(v: seq[float64]): string =
  for i, x in v:
    if i > 0:
      result.add(',')
    result.add(formatFloat(x, ffDecimal, 6))

proc writeUniques(path: string, res: DerepFastqResult) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open output file: " & path)
  defer: close(f)

  f.writeLine("index\tsequence\tabundance\tfirst_seen\tmean_qual")
  for i, u in res.uniques:
    f.writeLine(&"{i}\t{u.sequence}\t{u.abundance}\t{u.firstSeenReadIndex}\t{meanQualToString(u.meanQual)}")

proc writeMap(path: string, readToUnique: seq[int]) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open map output file: " & path)
  defer: close(f)

  f.writeLine("read_index\tunique_index")
  for i, u in readToUnique:
    f.writeLine(&"{i}\t{u}")

proc writeSummary(path: string, res: DerepFastqResult, opts: DerepFastqOptions) =
  var f: File
  if not open(f, path, fmWrite):
    raise newException(IOError, "Unable to open summary file: " & path)
  defer: close(f)

  f.writeLine("input_reads=" & $res.totalReads)
  f.writeLine("unique_sequences=" & $res.uniqueCount)
  f.writeLine("max_read_length=" & $res.maxReadLength)
  f.writeLine("map_size=" & $res.readToUnique.len)
  f.writeLine("options=" & opts.describe())

proc printSummary(res: DerepFastqResult, opts: DerepFastqOptions) =
  echo "Mode: derepFastq"
  echo &"Input reads: {res.totalReads}"
  echo &"Unique sequences: {res.uniqueCount}"
  echo &"Max read length: {res.maxReadLength}"
  echo &"Read->unique map size: {res.readToUnique.len}"
  echo &"Options: {opts.describe()}"

proc main() =
  var inputPath = ""
  var outputPath = ""
  var mapOutPath = ""
  var summaryPath = ""
  var quiet = false
  var verbose = false

  var opts = defaultDerepFastqOptions()

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
        of "in":
          inputPath = value
        of "out":
          outputPath = value
        of "map-out":
          mapOutPath = value
        of "summary":
          summaryPath = value
        of "threads":
          opts.threads = parseIntFlag(flag, value)
        of "chunk-size":
          opts.chunkSize = parseIntFlag(flag, value)
        of "parallel-min-chunk":
          opts.parallelMinChunk = parseIntFlag(flag, value)
        of "order":
          case value.toLowerAscii()
          of "abundance":
            opts.order = dsoAbundance
          of "first", "firstseen", "first-seen":
            opts.order = dsoFirstSeen
          else:
            fail("Invalid value for --order: " & value)
        of "no-map":
          opts.keepReadMap = false
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
    fail("Missing required option: --in")
  if outputPath.len == 0:
    fail("Missing required option: --out")
  if quiet and verbose:
    fail("--quiet and --verbose are mutually exclusive")
  if mapOutPath.len > 0 and not opts.keepReadMap:
    fail("--map-out cannot be used with --no-map")

  try:
    validateOptions(opts)

    if verbose:
      echo &"[derepFastq] Starting: {inputPath}"
      echo &"[derepFastq] Options: {opts.describe()}"

    let res = derepFastq(inputPath, opts)

    if verbose:
      echo &"[derepFastq] Writing unique table: {outputPath}"
    writeUniques(outputPath, res)

    if mapOutPath.len > 0:
      if verbose:
        echo &"[derepFastq] Writing read map: {mapOutPath}"
      writeMap(mapOutPath, res.readToUnique)

    if summaryPath.len > 0:
      writeSummary(summaryPath, res, opts)

    if not quiet:
      printSummary(res, opts)
  except CatchableError as e:
    fail(e.msg)

when isMainModule:
  main()
