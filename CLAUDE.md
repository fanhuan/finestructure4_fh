# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fineSTRUCTURE v4 (`fs`) is a C++/C bioinformatics tool for Bayesian population structure inference from dense haplotype data. It integrates three components: ChromoPainter (haplotype painting), ChromoCombine (chromosome combination), and fineSTRUCTURE (MCMC clustering). A secondary binary `mixPainter` supports the badMIXTURE analysis.

## Build System

Uses GNU Autotools (`autoconf`/`automake`). Standard build:

```sh
./configure
make
make install   # optional
```

To compile `mixPainter` after `make`:
```sh
make mixPainter
```

If `configure` is missing or stale, regenerate it (requires `automake`, `autoreconf`):
```sh
./autogen.sh
```

On Mac, `clang` lacks OpenMP support. Use a GCC compiler explicitly:
```sh
./configure CXX=g++-11   # substitute your version
# or use the convenience script:
./configure.mac
```

There are no external library dependencies as of v4.1.0 — zlib is bundled in `cp/`.

## Code Architecture

The codebase is organized into four logical layers:

### 1. `fs` wrapper (root directory)
`fs.cpp` is the entry point. It dispatches to one of three modes — project mode (`FsProject`), classic finestructure mode, or classic chromopainter mode. The wrapper layer handles the **project file** (`.cp` file), which stores all parameters and pipeline state across stages.

- `fsproject.*` — `FsProject` class: reads/writes the `.cp` project file, orchestrates stage execution, manages HPC command generation
- `fscmds.*` — `FsPar`/`FsCmd` classes: parameter and command definitions
- `fssettings.*` — `FsSettings`/`FsSettingsValue`: key-value settings storage
- `fsparam.*` — parameter parsing and defaults
- `fsdonor.*` — donor population management
- `fsutils.*` — shared utilities
- `fsconstants.h` — compile-time constants

### 2. `finestructure/` — Bayesian clustering core (C++)
Contains the MCMC inference engine and related statistics. Key files:
- `fines.h/cpp` — top-level entry points for the finestructure algorithm
- `data.h/cpp` — data structures for the chunk-count coancestry matrix
- `state.h/cpp` — MCMC state representation (population assignments, tree)
- `node.h/cpp` — tree node structure
- `prior.h/cpp` — Dirichlet process prior
- `inf1.h/cpp` — core inference (split/merge moves)
- `infmcmc.h/cpp` — MCMC driver
- `infadmixture.*`, `infconcordance.*` — admixture and concordance analysis
- `infextract*.h/cpp` (1–5, donors) — output extraction at different summarization levels
- `rng.*` — random number generator
- `fsxml.*` — XML output format

### 3. `cp/` — ChromoPainter (C)
Modified ChromoPainter v1 (by Garrett Hellenthal). Performs haplotype copying model inference (Li & Stephens model) to produce the coancestry chunk-count matrix consumed by finestructure. The bundled zlib files (`gz*.c`, `deflate.c`, etc.) handle compressed I/O.

### 4. `chromocombine/` — ChromoCombine (C++)
Combines ChromoPainter output across multiple chromosomes into a single coancestry matrix.

## Pipeline Stages

The `fs` tool runs a 4-stage pipeline stored in the `.cp` project file:

| Stage | Description |
|-------|-------------|
| s1    | ChromoPainter EM (parameter estimation) |
| s2    | ChromoCombine (merge chromosomes) |
| s3    | fineSTRUCTURE MCMC |
| s4    | Tree building (maximization) |

Run all stages automatically with `-go`. For large datasets use `-hpc <N>` to generate parallelizable shell commands.

## Input Data Format

Inputs are in CHROMOPAINTERv2 format:
- `.phase` — phased haplotype data
- `.recombfile` — recombination map
- `.ids` — individual IDs and population labels

Conversion scripts in `scripts/`:
- `vcf2cp.pl` — VCF to ChromoPainter format
- `impute2chromopainter.pl` — IMPUTE2 format
- `beagle2chromopainter.pl` — Beagle format
- `makeuniformrecfile.pl` — create uniform recombination map

## Running Examples

```sh
cd examples/example1
fs example_cp.cp -n -phasefiles example_cp.phase -recombfiles example_cp.recombfile -idfile example_cp.ids -go
```

Example 2 demonstrates multi-chromosome analysis with European data. Example 4 contains `FinestructureLibrary.R` / `FinestructureExample.R` for downstream R analysis.

## R Analysis

`FinestructureLibrary.R` (root and `examples/example4/`) provides functions for reading and plotting fineSTRUCTURE XML output. Not part of the compiled binary.
