# easy-matrix-method

A Python framework for analyzing **Imaging Atmospheric Cherenkov Telescope (IACT)** data using the **Matrix Background Method**.
It provides an automated pipeline for generating matched OFF runs, modeling and subtracting background, detecting sources, measuring spectra, and producing diagnostic plots.

---

## Overview

This repository implements the **Matrix Method** for background estimation in VERITAS (and other IACT) data analyses.
The method improves upon traditional reflected-region or ring-background techniques, especially for **extended or complex sources**, by constructing a matrix that relates ON and OFF regions in camera coordinates.

The package automates the key steps in a high-level IACT analysis pipeline:

1. Generate matched OFF run lists for each target field.
2. Build background models using exposure and acceptance corrections.
3. Perform source detection and significance calculations.
4. Extract energy spectra and light curves.
5. Produce publication-quality maps and diagnostic plots.

---

## Features

* **Matched OFF run selection** based on zenith angle, azimuth, and observation conditions.
* **Background modeling** using the matrix method for extended source analysis.
* **Excess and significance map generation** with adaptive smoothing.
* **Spectral reconstruction** (including effective area correction and energy binning).
* **Diagnostic visualization** for quality control and inspection.
* **Support for both point-like and extended source analyses.**

---

## Installation

Clone this repository and set up the environment:

```bash
git clone https://github.com/ruoyushang/easy-matrix-method.git
cd easy-matrix-method
```

It is recommended to use a conda environment with ROOT and the required Python dependencies:

```bash
conda create -n ematrix python=3.10 numpy scipy matplotlib astropy uproot
conda activate ematrix
```

Make sure your **ROOT** installation (with Python bindings) is available in the environment.

---

## Dependencies

The main dependencies are:

* **ROOT** (for reading and handling VERITAS DST / event files)
* **NumPy**, **SciPy** — numerical operations and optimization
* **Astropy** — coordinate transformations and FITS file handling
* **Matplotlib** — plotting and diagnostics
* **uproot** — for ROOT I/O in pure Python

---

## Workflow

A typical analysis using this package involves:

1. **Prepare run lists**
   Create a text file (e.g., `RunList_Crab.txt`) listing run numbers and corresponding source names.

2. **Generate matched OFF runs**
   The pipeline matches OFF runs with similar observing conditions for each ON run.

3. **Run background modeling and subtraction**
   The main script `EasyMatrixMethod.py` constructs the background model and computes excess counts.

4. **Detect sources and compute significance**
   Apply the Li & Ma (1983) significance calculation per sky bin.

5. **Extract spectra**
   Use the reconstructed events and effective area to produce differential flux measurements.

6. **Generate diagnostic plots**
   Visualize run-by-run statistics, acceptance curves, and residual maps.

---

## Usage Example

Example command:

```bash
python EasyMatrixMethod.py --runlist RunList_Crab.txt --output ./results/Crab/
```

Typical arguments:

| Option          | Description                                |
| --------------- | ------------------------------------------ |
| `--runlist`     | Text file listing ON runs and source names |
| `--output`      | Output directory for results               |
| `--energy-bins` | Optional custom energy binning             |
| `--smooth`      | Enable adaptive smoothing for sky maps     |
| `--verbose`     | Print detailed diagnostic information      |

---

## Outputs

After running the pipeline, the following outputs are produced:

| Output Type             | Description                                       |
| ----------------------- | ------------------------------------------------- |
| `sky_map.fits`          | Sky map of excess counts                          |
| `bkg_map.fits`          | Background model map                              |
| `significance_map.fits` | Significance map using Li & Ma formula            |
| `spectrum.txt`          | Measured energy spectrum                          |
| `diagnostic_plots/`     | QA plots for acceptance, alpha distribution, etc. |
| `logs/`                 | Run logs and matching summary                     |

---

## Acknowledgements

This code is developed for **VERITAS IACT analysis** and may be adapted for other Cherenkov telescope arrays.
It builds on methods used in standard VERITAS background modeling workflows.

Special thanks to the VERITAS Collaboration supporting this work.

---

*Maintained by Ruo-Yu Shang, Barnard College, Columbia University.*

