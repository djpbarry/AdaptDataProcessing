# Adapt Data Processing

> **⚠️ Deprecated — archive candidate.** This project is no longer actively
> developed. The analysis it performs is better served by a Python rewrite
> (pandas / NumPy / SciPy / matplotlib). See [Recommendation](#recommendation)
> below.

An [ImageJ](https://imagej.net) plugin (Java/Maven) for automated analysis of
time-series microscopy data, originally built at the Francis Crick Institute
(2014) to measure actin/adaptor dynamics.

## What it does

The plugin reads tracking output produced by upstream acquisition/segmentation
plugins and derives summary statistics.

- **`DataFileAverager`** — the main engine. Reads a directory of `.csv` tracks,
  truncates each track at the first NaN/negative-velocity reversal, normalises
  selected parameters to `[0, 1]`, then computes the **per-timepoint mean and
  standard deviation** across files. Writes `mean_data.csv` and `file_list.txt`
  (and optionally displays per-parameter plots via ImageJ's `ij.gui.Plot`).
- **`DataResampler`** — intended to resample track data to a new sample rate;
  its core logic (`run()`) is currently commented out, so it is effectively
  non-functional.
- **`DetectionMapAnalyser`** — analyses TIFF detection maps (`SignalMap.tif`
  plus `Maps/` and `PlotDataFiles/` folders), computing per-column signal
  standard-deviation and correlation maps.

Statistical helpers (mean/stddev/extrema, zero-crossings) and file I/O are
supplied by the external `com.github.djpbarry:IAClassLibrary` dependency.

## Building

```bash
mvn verify
```

The single dependency is resolved from jitpack.io / Maven Central. CI
(`.github/workflows/maven.yml`) additionally authenticates against
`maven.pkg.github.com/djpbarry/*` using `mvn_settings.xml` and a `PAT` secret.

## Recommendation

Most of this codebase is redundant: dead/commented-out code, hand-rolled
statistics and CSV parsing, and an ImageJ/Maven wrapper around a small amount of
actual analysis. The meaningful work can be expressed in a few hundred lines of
Python:

| Current Java | Python equivalent |
|---|---|
| custom `FileReader` + nested `ArrayList`s | `pandas.read_csv` → DataFrame |
| truncate / normalise | `df.clip()`, `(x-min)/(max-min)` |
| per-timepoint mean/σ across files | `df.groupby(...).agg(['mean','std'])` |
| resampling | `scipy.signal.resample_poly` |
| zero-crossing / velocity stats | a few lines of NumPy |
| correlation maps | `np.correlate` / `scipy.ndimage` |
| plotting | `matplotlib` |

A self-contained version of the core `DataFileAverager` pipeline (read tracks →
truncate → normalise → per-timepoint mean/σ) is roughly:

```python
import sys
from pathlib import Path

import pandas as pd


def average_tracks(directory, norm_columns=(), velocity_col="Velocity"):
    """Core of `DataFileAverager` in Python.

    Reads every `.csv` track in `directory` (one header row), truncates each
    track at the first NaN or negative-velocity reversal, normalises
    `norm_columns` to [0, 1] per file, then writes per-timepoint mean and
    standard deviation to `mean_data.csv` plus a `file_list.txt`.
    """
    files = sorted(Path(directory).glob("*.csv"))
    frames, names = [], []

    for path in files:
        df = pd.read_csv(path)
        vel = df[velocity_col].to_numpy()

        # Truncate at the first NaN, or the first time velocity returns to
        # >= 0 after having gone negative.
        nan_mask = df.isna().any(axis=1).to_numpy()
        pos = neg = False
        cut = len(df)
        for i, v in enumerate(vel):
            if nan_mask[i] or (neg and v >= 0.0):
                cut = i
                break
            if pos and v < 0.0:
                neg = True
            elif v >= 0.0:
                pos = True

        df = df.iloc[:cut].copy()
        for col in norm_columns:
            lo = df[col].where(df[col] >= 0.0).min()
            hi = df[col].max()
            df[col] = (df[col] - lo) / (hi - lo)

        frames.append(df.reset_index(drop=True))
        names.append(path.name)

    # Per-timepoint mean/std across files (tracks may differ in length).
    data = pd.concat(frames, keys=names, names=["file", "t"])
    grouped = data.groupby("t")

    out = grouped.mean(numeric_only=True).copy()
    out.insert(0, "N", grouped.size())
    out.to_csv(Path(directory) / "mean_data.csv")

    Path(directory, "file_list.txt").write_text("\n".join(names) + "\n")
    return out


if __name__ == "__main__":
    average_tracks(sys.argv[1], norm_columns=sys.argv[2:])
```

For new work, prefer a Python implementation over extending this plugin. If the
in-ImageJ user experience is genuinely required, retain a thin ImageJ entry
point and call Python via `pyimagej`/`scyjava` rather than porting the analysis
logic back into Java.
