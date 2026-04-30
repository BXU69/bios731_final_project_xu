# BIOS 731 Final Project: Latent-State Matrix-Logarithmic Covariance Regression

---

## Overview

This project develops a **Latent-State Matrix-Logarithmic Covariance Regression (LS-MLCR)** model for mapping resting-state functional connectivity to task-state functional connectivity in neuroimaging data. The model applies a matrix-logarithm transformation to symmetric positive-definite (SPD) covariance matrices, fits an edge-wise linear regression in the resulting Euclidean space, and discovers latent subgroups of subjects with distinct rest-to-task mapping profiles via a finite Gaussian mixture estimated by the EM algorithm. Inference is carried out via nonparametric bootstrap parallelized as a SLURM job array on the Emory RSPH high-performance computing cluster.

A simulation study across 18 factorial scenarios (n ∈ {50, 100, 200}, K ∈ {2, 3, 4}, p ∈ {6, 15}, 500 replicates each) evaluates parameter recovery, BIC model selection accuracy, and clustering accuracy (Adjusted Rand Index).

---

## Repository Structure

```
bios731_final_project/
│
├── R/                          # Core R functions
│   ├── em_functions.R          # EM algorithm, E/M-step, model selection (BIC),
│   │                           #   label alignment, and matrix utilities (mat_log, vech, etc.)
│   ├── simulate_data.R         # Data generating function sim_lcr_data() for simulation study
│   ├── run_bootstrap.R         # Nonparametric bootstrap: each SLURM task runs a batch
│   │                           #   of replicates and saves results as .RDA files
│   └── combine_results.R       # Aggregates bootstrap .RDA outputs and computes
│                               #   standard errors and confidence intervals
│
├── simulations/                # Simulation study scripts
│   ├── run_sim.R               # Main simulation script; takes SLURM array task ID as argument,
│   │                           #   runs 500 Monte Carlo replicates for the corresponding scenario
│   └── run_sim.sh              # SLURM batch script (18-task array, one task per scenario)
│
├── results/
│   └── 20260414/               # Simulation results (one .RDA file per scenario, 1-18)
│       ├── 1.RDA
│       ├── 2.RDA
│       └── ...
│
├── figures_and_tables.Rmd      # Reproduces all figures and tables in the report
├── figures_and_tables.html     # Rendered output of figures_and_tables.Rmd
│
├── bios_731_final_project_report.pdf   # Final project report (PDF)
└── bios731_final_project.Rproj         # RStudio project file
```

---

## Reproducing the Results

### Prerequisites

Install the required R packages:

```r
install.packages(c("tidyverse", "here", "mclust", "tictoc", "kableExtra", "scales"))
```

### Step 1: Run the simulation study

The simulation study is designed to run on a SLURM cluster. From the project root directory:

```bash
sbatch simulations/run_sim.sh
```

This submits an 18-task array job (`#SBATCH --array=1-18`). Each task runs 500 Monte Carlo replicates for one scenario and saves its output to `results/20260414/<task_id>.RDA`. Pre-computed results are already included in the repository under `results/20260414/`.

To run a single scenario locally (e.g., scenario 1) for testing:

```bash
Rscript simulations/run_sim.R 1
```

### Step 2: Reproduce figures and tables

Open `figures_and_tables.Rmd` in RStudio and knit the document, or run:

```r
rmarkdown::render("figures_and_tables.Rmd")
```

This will load the pre-computed simulation results from `results/20260414/`, reproduce all figures and tables in the report, and render them to `figures_and_tables.html`.

### Step 3: Bootstrap inference (future work)

Bootstrap inference on real fMRI data is designed to run via SLURM job arrays using `R/run_bootstrap.R`. After all bootstrap tasks complete, aggregate results with:

```r
source("R/combine_results.R")
```

This step requires a real fMRI dataset (`results/full_fit.RDA`) not included in this repository.
