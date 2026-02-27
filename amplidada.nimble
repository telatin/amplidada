# Package

version       = "0.1.0"
author        = "Andrea Telatin"
description   = "Nim library for DADA2 denoising."
license       = "GPL-3.0"
srcDir        = "src"
binDir        = "bin"
bin           = @[
  "cli/dada2=dada2",
  "cli/derepFastq=derepFastq",
  "cli/fastqFilter=fastqFilter",
  "cli/filterAndTrim=filterAndTrim",
  "cli/learnErrors=learnErrors",
  "cli/removeBimerDenovo=removeBimerDenovo"
]

# Dependencies

requires "nim >= 2.0.0"
requires "readfx >= 0.3.1"
requires "argparse >= 4.0.0"
requires "malebolgia >= 1.3.0"

task docs, "Generate HTML documentation":
  exec "nim doc --project --index:on --nimcache:docs/.nimcache --outdir:docs src/amplidada.nim"

task test, "Run smoke tests":
  exec "nim c -r --threads:on --nimcache:.nimcache --hints:off --verbosity:0 tests/test_all.nim"

task integration, "Run optional slow real-data integration test (set AMPLIDADA_RUN_SLOW_INTEGRATION=1)":
  exec "nim c -r --threads:on --nimcache:.nimcache-integration --hints:off --verbosity:0 tests/test_integration_real_data.nim"

task tools, "Build CLI binaries into ./bin":
  exec "mkdir -p bin && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/dada2.nim && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/derepFastq.nim && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/fastqFilter.nim && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/filterAndTrim.nim && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/learnErrors.nim && nim c --threads:on -d:release --nimcache:.nimcache-tools --outdir:bin src/cli/removeBimerDenovo.nim"
