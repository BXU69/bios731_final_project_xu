####################################################################
# Baijia Xu
# April 2026
#
# Simulation study for Latent-State Matrix-Log Covariance Regression
# Evaluates parameter recovery, model selection (BIC), and
# clustering accuracy (adjusted Rand index) across scenarios
####################################################################


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(tictoc))


wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################

source(here::here("R", "em_functions.R"))
source(here::here("R", "simulate_data.R"))

###############################################################
## set simulation design elements
###############################################################

nsim = 500

n = c(50, 100, 200)
K = c(2, 3, 4)
p = c(6, 15)

params = expand.grid(n = n, K = K, p = p)

## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 1
  params = params[scenario,]
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
  params = params[scenario,]
}

cat(sprintf("Scenario %d: n = %d, K = %d, p = %d\n", scenario, params$n, params$K, params$p))

###############################################################
## start simulation code
###############################################################

p = params$p
seed = floor(runif(nsim, 1, 10000))
results = as.list(rep(NA, nsim))

for(i in 1:nsim){
  set.seed(seed[i])

  ####################
  # simulate data
  simdata = sim_lcr_data(n = params$n, p = p, K = params$K, seed = seed[i])

  ####################
  # fit EM at true K and K +/- 1 for model selection evaluation
  K_true = params$K
  K_candidates = unique(pmax(1, (K_true - 1):(K_true + 1)))

  tic()
  sel = select_K(simdata$s_mat, simdata$r_mat,
                 K_max = max(K_candidates),
                 n_init = 10, tol = 1e-4, max_iter = 300)
  time_sel = toc(quiet = TRUE)

  # restrict table to candidate K values
  sel_table = sel$table[sel$table$K %in% K_candidates, ]

  # BIC-selected K
  K_bic = sel_table$K[which.min(sel_table$BIC)]

  ####################
  # extract fit at true K
  fit_true_K = sel$fits[[K_true]]

  # align estimated parameters to true parameters
  aligned = align_labels(fit_true_K$params, simdata$true_params, K_true)

  # parameter estimation error: Frobenius norm of beta1 difference
  beta1_err = sqrt(sum((aligned$beta1 - simdata$true_params$beta1)^2))

  # clustering accuracy via adjusted Rand index
  # assign each subject to highest-responsibility state
  gamma = em_estep(simdata$s_mat, simdata$r_mat, fit_true_K$params)
  z_hat = apply(gamma, 1, which.max)

  # need to align z_hat labels before computing ARI
  # ARI is permutation-invariant, so no alignment needed

  ari = adjustedRandIndex(simdata$true_z, z_hat)

  ####################
  # store results
  estimates = tibble(
    sim = i,
    n = params$n,
    p = p,
    K_true = K_true,
    K_bic = K_bic,
    bic_correct = (K_bic == K_true),
    beta1_frob_err = beta1_err,
    ari = ari,
    converged = fit_true_K$converged,
    n_iter = fit_true_K$n_iter,
    time = time_sel$toc - time_sel$tic,
    seed = seed[i]
  )

  results[[i]] = estimates

  if(i %% 10 == 0) cat(sprintf("  completed sim %d/%d\n", i, nsim))
}


## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("results"), Date), showWarnings = FALSE)

filename = paste0(here::here("results", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

cat(sprintf("Saved results to %s\n", filename))

###############################################################
## end sim
###############################################################
