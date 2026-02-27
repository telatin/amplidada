import std/[os, osproc, strutils, unittest]
import amplidada

proc compileBinary(srcPath, outPath: string) =
  let nimCachePath = outPath & ".nimcache"
  let compileCmd = "nim c --threads:on --hints:off --verbosity:0 --nimcache:" & quoteShell(nimCachePath) &
    " --out:" & quoteShell(outPath) & " " & quoteShell(srcPath)
  let result = execCmdEx(compileCmd)
  check result.exitCode == 0
  if result.exitCode != 0:
    checkpoint "Compile output:\n" & result.output

proc runHelp(binPath: string) =
  let result = execCmdEx(quoteShell(binPath) & " --help")
  check result.exitCode == 0
  check "Usage:" in result.output

suite "amplidada smoke tests":
  test "library import and version":
    check AmpliDadaVersion.len > 0
    check libVersion() == AmpliDadaVersion

    let defaultFilter = defaultFastqFilterOptions()
    check defaultFilter.truncLen == 0
    check defaultFilter.maxEE == 2.0

  test "CLI startup --help":
    let binDir = getTempDir() / "amplidada-smoke-bin"
    createDir(binDir)

    let fastqFilterBin = binDir / "fastqFilter"
    let learnErrorsBin = binDir / "learnErrors"
    let derepFastqBin = binDir / "derepFastq"
    let filterAndTrimBin = binDir / "filterAndTrim"
    let dada2Bin = binDir / "dada2"
    let removeBimerDenovoBin = binDir / "removeBimerDenovo"

    compileBinary("src/cli/fastqFilter.nim", fastqFilterBin)
    compileBinary("src/cli/learnErrors.nim", learnErrorsBin)
    compileBinary("src/cli/derepFastq.nim", derepFastqBin)
    compileBinary("src/cli/filterAndTrim.nim", filterAndTrimBin)
    compileBinary("src/cli/dada2.nim", dada2Bin)
    compileBinary("src/cli/removeBimerDenovo.nim", removeBimerDenovoBin)

    runHelp(fastqFilterBin)
    runHelp(learnErrorsBin)
    runHelp(derepFastqBin)
    runHelp(filterAndTrimBin)
    runHelp(dada2Bin)
    runHelp(removeBimerDenovoBin)
