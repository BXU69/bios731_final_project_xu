# BIOS 731 Final Project: Latent-State Matrix-Logarithmic Covariance Regression

---

## Overview

This project develops a **Latent-State Matrix-Logarithmic Covariance Regression (LS-MLCR)** model for mapping resting-state functional connectivity to task-state functional connectivity in neuroimaging data. The model applies a matrix-logarithm transformation to symmetric positive-definite (SPD) covariance matrices, fits an edge-wise linear regression in the resulting Euclidean space, and discovers latent subgroups of subjects with distinct rest-to-task mapping profiles via a finite Gaussian mixture estimated by the EM algorithm. Inference is carried out via nonparametric bootstrap parallelized as a SLURM job array on the Emory RSPH high-performance computing cluster.

A simulation study across 18 factorial scenarios (n ∈ {50, 100, 200}, K ∈ {2, 3, 4}, p ∈ {6, 15}, 500 replicates each) evaluates parameter recovery, BIC model selection accuracy, and clustering accuracy (Adjusted Rand Index).

---

## Methods

### Data Transformation

For each subject $i = 1, \ldots, n$, we observe two $p \times p$ symmetric positive-definite (SPD) covariance matrices: $C_i^{(\text{rest})}$ and $C_i^{(\text{task})}$. Since SPD matrices lie on a curved Riemannian manifold, we apply the **matrix logarithm** to map them into an unconstrained Euclidean space:

$$L_i^{(\cdot)} = \log(C_i^{(\cdot)}) = U \, \text{diag}(\log \lambda_1, \ldots, \log \lambda_p) \, U^\top$$

where $C = U \Lambda U^\top$ is the eigendecomposition. We then half-vectorize the upper triangle (including diagonal) via $\text{vech}(\cdot)$, yielding $q = p(p+1)/2$-dimensional vectors $r_i$ (resting-state, predictor) and $s_i$ (task-state, response).

### Model

We assume $K$ latent states with mixing proportions $\pi_k > 0$, $\sum_k \pi_k = 1$, and latent assignment $z_i \in \{1, \ldots, K\}$. Conditional on $z_i = k$, the task connectivity vector follows an **edge-wise Gaussian regression**:

$$s_i \mid z_i = k \;\sim\; \mathcal{N}_q \left(\beta_0^{(k)} + \text{diag}(\beta_1^{(k)}) \, r_i, \;\text{diag}(\sigma^{2(k)})\right)$$

where $\beta_0^{(k)}, \beta_1^{(k)} \in \mathbb{R}^q$ are state-specific intercept and slope vectors, and $\sigma^{2(k)} \in \mathbb{R}^q$ is a vector of positive residual variances. The diagonal covariance structure means each edge is modeled by an independent univariate regression. The marginal distribution of $s_i$ is a $K$-component Gaussian mixture:

$$p(s_i \mid r_i, \theta) = \sum_{k=1}^K \pi_k \, \mathcal{N}_q\left(s_i;\; \beta_0^{(k)} + \text{diag}(\beta_1^{(k)}) r_i,\; \text{diag}(\sigma^{2(k)})\right)$$

where $\theta = \{ \pi_k, \beta_0^{(k)}, \beta_1^{(k)}, \sigma^{2(k)} \}_{k=1}^K$ collects all model parameters.

### EM Algorithm

The model is estimated via the EM algorithm, iterating between:

**E-step.** Compute posterior responsibilities for each subject $i$ and state $k$:

$$\gamma_{ik} = \frac{\pi_k \, \phi_q(s_i;\, \mu_i^{(k)},\, \text{diag}(\sigma^{2(k)}))}{\sum_{\ell=1}^K \pi_\ell \, \phi_q(s_i;\, \mu_i^{(\ell)},\, \text{diag}(\sigma^{2(\ell)}))}$$

where $\mu_i^{(k)} = \beta_0^{(k)} + \text{diag}(\beta_1^{(k)}) r_i$ and $\phi_q$ denotes the $q$-dimensional Gaussian density. Computed via the log-sum-exp trick for numerical stability.

**M-step.** With $N_k = \sum_i \gamma_{ik}$, $\bar{r}_{kj} = N_k^{-1} \sum_i \gamma_{ik} r_{ij}$, $\bar{s}_{kj} = N_k^{-1} \sum_i \gamma_{ik} s_{ij}$, update in closed form for each state $k$ and edge $j$:

$$\hat{\pi}_k = \frac{N_k}{n}, \qquad \hat{\beta}_{1j}^{(k)} = \frac{\sum_i \gamma_{ik}(r_{ij} - \bar{r}_{kj})(s_{ij} - \bar{s}_{kj})}{\sum_i \gamma_{ik}(r_{ij} - \bar{r}_{kj})^2}$$

$$\hat{\beta}_{0j}^{(k)} = \bar{s}_{kj} - \hat{\beta}_{1j}^{(k)} \bar{r}_{kj}, \qquad \hat{\sigma}_j^{2(k)} = \max \left(\frac{\sum_i \gamma_{ik}(s_{ij} - \hat{\beta}_{0j}^{(k)} - \hat{\beta}_{1j}^{(k)} r_{ij})^2}{N_k},\; 10^{-6}\right)$$

The algorithm is run with $n_{\text{init}} = 20$ random restarts (K-means initialization), returning the solution with the highest final log-likelihood. The number of states $K$ is selected via BIC: $\text{BIC}(K) = -2\ell(\hat{\theta}_K) + d_K \log n$, where $d_K = K - 1 + 2qK$.

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
├── figures_and_tables.pdf      # Rendered output of figures_and_tables.Rmd
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

This will load the pre-computed simulation results from `results/20260414/`, reproduce all figures and tables in the report, and render them to `figures_and_tables.pdf`.

### Step 3: Bootstrap inference (future work)

Bootstrap inference on real fMRI data is designed to run via SLURM job arrays using `R/run_bootstrap.R`. After all bootstrap tasks complete, aggregate results with:

```r
source("R/combine_results.R")
```

This step requires a real fMRI dataset (`results/full_fit.RDA`) not included in this repository.
