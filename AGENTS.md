# AGENTS.md

> **DEPRECATED / ARCHIVE CANDIDATE.** This project is no longer actively
> developed. The functionality is better served by a Python rewrite (pandas /
> NumPy / SciPy / matplotlib); see "Python migration note" below and the README.
> Prefer retiring this repository over extending it.

## Project Overview

An [ImageJ](https://imagej.net) plugin (Maven project) for automated analysis of
time-series microscopy data (e.g. actin/adaptor dynamics). It aggregates,
normalises, resamples, and analyses CSV data tracks and detection-map TIFFs
produced by upstream acquisition/segmentation plugins.

- **Language:** Java 11
- **Build tool:** Maven (parent `pom-scijava:35.0.0`, SciJava plugin framework)
- **Group/artifact:** `net.calm:adaptdataprocessing` (package root `net.calm.adaptdataprocessing`)
- **License:** header comments say GPLv3, but `pom.xml` declares Simplified BSD
  (see "Gotchas")

## Commands

No Makefile, no tests, no lint config exist.

- **Build:** `mvn verify`
- **Build (with custom settings for GitHub Packages auth):**
  `mvn --batch-mode --update-snapshots -Dinternal.repo.password="$PAT" --settings mvn_settings.xml verify`
- **CI:** GitHub Actions (`.github/workflows/maven.yml`) runs the above on
  every push using JDK 11 (AdoptOpenJDK) and the `PAT` secret.

The single dependency `com.github.djpbarry:IAClassLibrary` is resolved from
jitpack.io / Maven Central (declared in `pom.xml` repositories). When building
locally you may need `mvn_settings.xml` to authenticate against
`maven.pkg.github.com/djpbarry/*` via `$PAT`.

## Architecture

Four classes, all in `src/main/java/net/calm/adaptdataprocessing/DataProcessing/`:

| Class | Role |
|-------|------|
| `DataAnalytics` | Plain data holder (bean) for per-track statistics: peak index, zero crossing, +/-, pos/neg means, min/max. Getters/setters only. |
| `DataFileAverager` | Main processing engine. Reads a directory of `.csv` tracks (via `FileReader` from `IAClassLibrary`), truncates, normalises, computes mean data, optionally plots (`ij.gui.Plot`), and writes aggregated outputs (`Aggregated_Data`, `Collated_Data`, `mean_data.csv`, `file_list.txt`). |
| `DataResampler` | `PlugIn` that resamples track data to a new sample rate (uses `DSPProcessor.upScale`). |
| `DetectionMapAnalyser` | `PlugIn` that analyses TIFF detection maps and produces correlation maps + per-column SD signal files. |

**Control/data flow:**
1. A `PlugIn` (`DataResampler`, `DetectionMapAnalyser`) or the `DataFileAverager`
   entry point is invoked by ImageJ with a directory prompt (`Utilities.getFolder`).
2. `DataFileAverager.run(dir)` is the workhorse: lists `.csv` files →
   `FileReader.getParamList`/`readData` → truncate → normalise → `calcMeanData`
   → optional plots → `outputFileList`.
3. Statistical math (min/max/mean/stddev by percentile) is delegated to
   `net.calm.iaclasslibrary.IAClasses.DataStatistics`, and file I/O to
   `net.calm.iaclasslibrary.IO.FileReader` / `UtilClasses.*`.

The `IAClassLibrary` dependency supplies most non-trivial logic
(`DataStatistics`, `FileReader`, `GenUtils`, `Utilities`, `DSPProcessor`).
It is an external repo published under `com.github.djpbarry`.

## Conventions & Style

- **Naming:** classes use `DataProcessing` package; the package/artifact name is
  `DataProcessing` (capitalised, no underscore) despite the repo/project being
  "AdaptDataProcessing". Filenames match class names exactly.
- **Java style:** older-style codebase — package-private access on helper methods
  (e.g. `getParamIndices`, `getMapSDSig`, `generateFile`), fields often `final`
  for constants, `@author barry05` in Javadoc headers.
- **GPLv3 header** is copy-pasted at the top of every source file.
- Heavy use of commented-out code blocks (entire old implementations left in
  place). Do not assume commented code is dead to be deleted; it documents
  prior behaviour.
- No `static` `main` methods in active code — entry is via ImageJ's `PlugIn`
  interface (`run(String arg)`) or direct construction of `DataFileAverager`.

## Gotchas

- **License mismatch:** every source file carries a GPLv3 notice, and the repo
  ships a GPLv3 `LICENSE` file, but `pom.xml` declares `bsd_2` and organization
  "Francis Crick Institute". Don't "fix" one without resolving this discrepancy
  explicitly.
- **Hardcoded file/folder names:** `DataFileAverager` relies on literal names
  `file_list.txt`, `Aggregated_Data`, `Collated_Data`, `mean_data.csv`;
  `DetectionMapAnalyser` expects `SignalMap.tif`, `Maps/`, `PlotDataFiles/`
  on disk. These are implicit output/input contracts, not config.
- **Constant flags:** `DataFileAverager` has `NORM`, `TRUNC`, `COLLATE` all
  hardcoded `true`, with the corresponding code paths partially commented out
  (e.g. `aggregateData`, `colateData`, selection dialogs). Changing these flags
  may reference methods that are commented out.
- **Directory delimiter** is platform-dependent via `GenUtils.getDelimiter()`,
  not hardcoded `/` or `\`.
- **Dependency resolution** requires jitpack/central repositories (in `pom.xml`)
  and, for GitHub Packages snapshots, `mvn_settings.xml` + a `PAT`. The recent
  commits are all "pom dependencies" bumps.
- No unit tests exist; verification is manual/visual inside ImageJ.

## Testing

There is no test framework or test source. To validate changes, build with
`mvn verify` and, for behavioural changes, run the plugin inside an ImageJ/Fiji
installation against real `.csv` track directories or map folders.

## Python migration note

Most of this codebase is redundant and should not be extended. The *meaningful*
logic is a small slice of `DataFileAverager` (read CSV tracks → truncate →
normalise → per-timepoint mean/σ → write `mean_data.csv` / `file_list.txt`) and
the `DetectionMapAnalyser` correlation-map step. The rest is:

- **Dead code** — `DataResampler.run()` is fully commented out; `aggregateData`,
  `colateData`, `showSelectionDialog`, `performAnalysis`, `dataAnalysis`, and
  `getFrameRate` are present but uncalled/deactivated.
- **Reimplemented stdlib** — hand-rolled statistics, NaN handling, and CSV
  parsing (delegated to the external `IAClassLibrary`, itself homegrown).
- **Build/UX ceremony** — the Maven/SciJava/ImageJ wrapper exists only to run
  ~three analysis classes inside Java ImageJ.

A Python rewrite (pandas `read_csv`/`groupby().agg(['mean','std'])`, NumPy
clipping/zero-crossing, `scipy.signal` resampling, matplotlib for plots) would
collapse this to roughly ~200 lines and is the recommended path forward — keep
an ImageJ entry point via `pyimagej`/`scyjava` only if the in-ImageJ UX is
actually required.
