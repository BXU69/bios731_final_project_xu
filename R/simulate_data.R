#####################################################################################
# Baijia Xu
# April 2026
#
# Data simulation for Latent-State Matrix-Log Covariance Regression
# Generates synthetic data in the tangent space (vech of log-covariance)
#####################################################################################


#' Simulate data for the Latent-State Matrix-Log Covariance Regression
#'
#' @param n number of subjects
#' @param p number of brain regions (default 6)
#' @param K number of latent states
#' @param pi K-vector of mixing proportions (default: equal)
#' @param seed random seed
#' @return list with s_mat, r_mat, true_z, true_params, p, q, K
sim_lcr_data = function(n, p = 6, K = 2, pi = NULL, seed = 42) {
  set.seed(seed)

  q = p * (p + 1) / 2

  # default mixing proportions: equal
  if (is.null(pi)) {
    pi = rep(1 / K, K)
  }

  # generate state-specific parameters
  beta0 = matrix(NA, K, q)
  beta1 = matrix(NA, K, q)
  sigma2 = matrix(NA, K, q)

  for (k in 1:K) {
    # intercepts: different baselines per state, spread apart for identifiability
    beta0[k, ] = rnorm(q, mean = (k - 1) * 1.5, sd = 0.5)

    # slopes: in [0.3, 0.8], distinct across states
    # state 1: lower slopes, state 2: higher slopes, etc.
    slope_center = 0.3 + (k - 1) * (0.5 / max(K - 1, 1))
    beta1[k, ] = runif(q, min = slope_center - 0.1, max = slope_center + 0.1)

    # residual variances: around 0.1, slight variation across edges
    sigma2[k, ] = runif(q, min = 0.05, max = 0.15)
  }

  true_params = list(
    pi = pi,
    beta0 = beta0,
    beta1 = beta1,
    sigma2 = sigma2
  )

  # assign subjects to states via multinomial
  true_z = sample(1:K, size = n, replace = TRUE, prob = pi)

  # generate predictor (resting-state, in tangent space)
  r_mat = matrix(rnorm(n * q), n, q)

  # generate response (task-state, in tangent space)
  s_mat = matrix(NA, n, q)

  for (i in 1:n) {
    k = true_z[i]
    mu_i = beta0[k, ] + beta1[k, ] * r_mat[i, ]
    s_mat[i, ] = rnorm(q, mean = mu_i, sd = sqrt(sigma2[k, ]))
  }

  list(
    s_mat = s_mat,
    r_mat = r_mat,
    true_z = true_z,
    true_params = true_params,
    p = p,
    q = q,
    K = K
  )
}
