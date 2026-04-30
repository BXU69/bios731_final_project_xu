####################################################################
# Baijia Xu
# April 2026
#
# Nonparametric bootstrap for Latent-State Matrix-Log Covariance
# Regression. Each SLURM array task runs batch_size replicates.
####################################################################


suppressPackageStartupMessages(library(here))


wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## source functions
###############################################################

source(here::here("R", "em_functions.R"))

###############################################################
## set bootstrap parameters
###############################################################

batch_size = 10

if(doLocal) {
  task_id = 1
}else{
  task_id = as.numeric(commandArgs(trailingOnly=TRUE))
}

###############################################################
## load data and full-data fit
###############################################################

load(here::here("results", "processed_data.RDA"))
# expects: s_mat, r_mat

load(here::here("results", "full_fit.RDA"))
# expects: full_fit (output of em_lcr)

n = nrow(s_mat)
K = full_fit$K

###############################################################
## run bootstrap replicates
###############################################################

boot_results = as.list(rep(NA, batch_size))

for(b in 1:batch_size) {
  rep_id = (task_id - 1) * batch_size + b

  tryCatch({
    set.seed(rep_id * 2749 + 31)

    # resample subjects with replacement
    idx = sample(1:n, n, replace = TRUE)
    s_boot = s_mat[idx, , drop = FALSE]
    r_boot = r_mat[idx, , drop = FALSE]

    # fit EM with warm start from full-data estimate
    fit_boot = em_lcr(s_boot, r_boot, K = K,
                      n_init = 1,
                      init_params = full_fit$params,
                      tol = 1e-6,
                      max_iter = 500)

    boot_results[[b]] = list(
      rep_id = rep_id,
      params = fit_boot$params,
      converged = fit_boot$converged,
      n_iter = fit_boot$n_iter,
      loglik_final = fit_boot$loglik[fit_boot$n_iter]
    )
  }, error = function(e) {
    boot_results[[b]] <<- list(
      rep_id = rep_id,
      error = conditionMessage(e)
    )
  })

  if(b %% 5 == 0) cat(sprintf("  task %d: completed replicate %d/%d\n", task_id, b, batch_size))
}


###############################################################
## save results
###############################################################

dir.create(here::here("results", "boot"), showWarnings = FALSE)
filename = paste0(here::here("results", "boot"), "/boot_task", task_id, ".RDA")
save(boot_results, file = filename)
cat(sprintf("Saved %s\n", filename))

###############################################################
## end bootstrap task
###############################################################
