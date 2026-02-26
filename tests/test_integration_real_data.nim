import std/[os, osproc, strformat, strutils, tables, times]
import amplidada

const
  EnvRunSlow = "AMPLIDADA_RUN_SLOW_INTEGRATION"
  EnvMaxSamples = "AMPLIDADA_INTEGRATION_MAX_SAMPLES"

proc envBool(name: string; default = false): bool =
  let raw = getEnv(name, if default: "1" else: "0").strip().toLowerAscii()
  raw in ["1", "true", "yes", "on"]

proc envInt(name: string; default: int): int =
  let raw = getEnv(name, "").strip()
  if raw.len == 0:
    return default
  try:
    result = parseInt(raw)
  except ValueError:
    result = default

proc fail(msg: string) {.noreturn.} =
  stderr.writeLine("Integration test failure: " & msg)
  quit(2)

proc parseSummary(path: string): Table[string, string] =
  for line in lines(path):
    let t = line.strip()
    if t.len == 0:
      continue
    let p = t.find('=')
    if p <= 0:
      continue
    result[t[0 ..< p]] = t[p + 1 .. ^1]

proc requireInt(summary: Table[string, string], key: string): int64 =
  if key notin summary:
    fail("Missing key in summary: " & key)
  try:
    parseBiggestInt(summary[key]).int64
  except ValueError:
    fail("Invalid integer in summary for key " & key & ": " & summary[key])

proc main() =
  if not envBool(EnvRunSlow):
    echo &"Skipping optional real-data integration test. Set {EnvRunSlow}=1 to run."
    quit(0)

  let dataDir = "data"
  if not dirExists(dataDir):
    fail("Data directory not found: " & dataDir)

  let samples = discoverPairedSamples(
    inputDir = dataDir,
    r1Token = "_R1_001",
    r2Token = "_R2_001",
    recursive = false
  )

  if samples.len == 0:
    fail("No paired samples discovered in " & dataDir)

  let maxSamples = max(1, min(samples.len, envInt(EnvMaxSamples, samples.len)))
  let selected = samples[0 ..< maxSamples]

  let outDir = getTempDir() / ("amplidada-integration-" & $epochTime().int64)
  createDir(outDir)

  let binPath = outDir / "dada2"
  let nimCachePath = outDir / ".nimcache"
  let compileCmd =
    "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(binPath) & " " & quoteShell("src/cli/dada2.nim")
  let compileRes = execCmdEx(compileCmd)
  if compileRes.exitCode != 0:
    fail("Failed to compile dada2 CLI:\n" & compileRes.output)

  echo &"Running optional integration on {selected.len} paired sample(s) from {dataDir}"

  var totalCandidate = 0'i64
  var totalMerged = 0'i64

  for i, sample in selected:
    let mergedOut = outDir / (sample.sampleId & ".merged.tsv")
    let summaryOut = outDir / (sample.sampleId & ".summary.txt")
    let runCmd =
      quoteShell(binPath) &
      " --in-forward " & quoteShell(sample.readsForward) &
      " --in-reverse " & quoteShell(sample.readsReverse) &
      " --out " & quoteShell(mergedOut) &
      " --summary " & quoteShell(summaryOut) &
      " --threads 2 --merge-threads 2 --min-overlap 12 --max-mismatch 2 --trim-overhang " &
      " --max-iter 64 --max-shuffle 12 --quiet"

    echo &"  [{i + 1}/{selected.len}] {sample.sampleId}"
    let runRes = execCmdEx(runCmd)
    if runRes.exitCode != 0:
      fail(&"dada2 failed on sample {sample.sampleId}:\n{runRes.output}")
    if not fileExists(summaryOut):
      fail("Missing summary output for sample " & sample.sampleId)
    if not fileExists(mergedOut):
      fail("Missing merged output for sample " & sample.sampleId)

    let summary = parseSummary(summaryOut)
    if summary.getOrDefault("mode") != "paired":
      fail("Unexpected mode in summary for sample " & sample.sampleId)

    let candidate = requireInt(summary, "merge_candidate_pairs")
    let merged = requireInt(summary, "merged_pairs")
    let rejected = requireInt(summary, "merge_rejected_pairs")

    if candidate <= 0:
      fail("No merge candidates for sample " & sample.sampleId)
    if merged < 0 or rejected < 0:
      fail("Invalid merged/rejected counts for sample " & sample.sampleId)
    if merged + rejected != candidate:
      fail("Merged + rejected does not equal candidate for sample " & sample.sampleId)

    totalCandidate += candidate
    totalMerged += merged

  if totalCandidate <= 0:
    fail("No candidates across selected samples")
  if totalMerged <= 0:
    fail("No merged pairs across selected samples")

  echo &"Integration OK. samples={selected.len} candidate_pairs={totalCandidate} merged_pairs={totalMerged}"

when isMainModule:
  main()
