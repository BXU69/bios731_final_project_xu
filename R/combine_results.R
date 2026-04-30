####################################################################
# Baijia Xu
# April 2026
#
# Combine bootstrap results and compute inference summaries
# Run after all bootstrap SLURM jobs have completed
####################################################################


suppressPackageStartupMessages(library(here))

source(here::here("R", "em_functions.R"))

###############################################################
## load full-data fit (reference for label alignment)
###############################################################

load(here::here("results", "full_fit.RDA"))
# expects: full_fit

K = full_fit$K
q = ncol(full_fit$params$beta0)

###############################################################
## load all bootstrap result files
###############################################################

boot_dir = here::here("results", "boot")
boot_files = list.files(boot_dir, pattern = "^boot_task.*\\.RDA$", full.names = TRUE)

if (length(boot_files) == 0) stop("No bootstrap result files found in results/boot/")

cat("Loading", length(boot_files), "bootstrap files...\n")

all_params = list()
n_errors = 0

for (f in boot_files) {
  load(f)  # loads boot_results

  for (b in seq_along(boot_results)) {
    res = boot_results[[b]]

    if (!is.null(res$error)) {
      n_errors = n_errors + 1
      next
    }

    if (!res$converged) next

    # align labels to full-data fit
    aligned = align_labels(res$params, full_fit$params, K)
    all_params = c(all_params, list(aligned))
  }
}

cat(sprintf("Collected %d successful bootstrap replicates (%d errors)\n",
            length(all_params), n_errors))

###############################################################
## extract bootstrap distributions
###############################################################

B = length(all_params)

# stack parameter arrays: B x K x q for beta0, beta1, sigma2; B x K for pi
pi_boot = matrix(NA, B, K)
beta0_boot = array(NA, dim = c(B, K, q))
beta1_boot = array(NA, dim = c(B, K, q))
sigma2_boot = array(NA, dim = c(B, K, q))

for (b in 1:B) {
  pi_boot[b, ] = all_params[[b]]$pi
  beta0_boot[b, , ] = all_params[[b]]$beta0
  beta1_boot[b, , ] = all_params[[b]]$beta1
  sigma2_boot[b, , ] = all_params[[b]]$sigma2
}

###############################################################
## compute bootstrap SE and 95% percentile CI
###############################################################

# mixing proportions
pi_est = full_fit$params$pi
pi_se = apply(pi_boot, 2, sd)
pi_ci_lo = apply(pi_boot, 2, quantile, probs = 0.025)
pi_ci_hi = apply(pi_boot, 2, quantile, probs = 0.975)

pi_summary = data.frame(
  state = 1:K,
  estimate = pi_est,
  se = pi_se,
  ci_lower = pi_ci_lo,
  ci_upper = pi_ci_hi
)

# regression parameters (beta0, beta1, sigma2) per state and edge
make_summary = function(est_mat, boot_arr, param_name) {
  rows = list()
  for (k in 1:K) {
    for (j in 1:q) {
      boot_vals = boot_arr[, k, j]
      rows = c(rows, list(data.frame(
        param = param_name,
        state = k,
        edge = j,
        estimate = est_mat[k, j],
        se = sd(boot_vals),
        ci_lower = quantile(boot_vals, 0.025),
        ci_upper = quantile(boot_vals, 0.975)
      )))
    }
  }
  do.call(rbind, rows)
}

beta0_summary = make_summary(full_fit$params$beta0, beta0_boot, "beta0")
beta1_summary = make_summary(full_fit$params$beta1, beta1_boot, "beta1")
sigma2_summary = make_summary(full_fit$params$sigma2, sigma2_boot, "sigma2")

param_summary = rbind(beta0_summary, beta1_summary, sigma2_summary)
rownames(param_summary) = NULL

###############################################################
## save and print
###############################################################

save(pi_summary, param_summary, pi_boot, beta0_boot, beta1_boot, sigma2_boot,
     file = here::here("results", "bootstrap_inference.RDA"))

cat("\n=== Mixing Proportions ===\n")
print(pi_summary, digits = 4)

cat("\n=== Parameter Summary (first 20 rows) ===\n")
print(head(param_summary, 20), digits = 4)

cat(sprintf("\nFull summary saved to %s\n",
            here::here("results", "bootstrap_inference.RDA")))

###############################################################
## end
###############################################################
