#####################################################################################
# Baijia Xu
# April 2026
#
# EM algorithm for Latent-State Matrix-Logarithmic Covariance Regression
# Functions for matrix utilities, EM steps, model selection, and label alignment
#####################################################################################


#####################################################################################
# Matrix Utilities
#####################################################################################

#' Matrix logarithm via eigendecomposition
#' Adds eps * I if any eigenvalue is non-positive (to ensure log is defined)
mat_log = function(C, eps = 1e-8) {
  eig = eigen(C, symmetric = TRUE)
  vals = eig$values
  vecs = eig$vectors

  # shift if any eigenvalue is non-positive

  if (any(vals <= 0)) {
    shift = abs(min(vals)) + eps
    vals = vals + shift
  }

  vecs %*% diag(log(vals)) %*% t(vecs)
}


#' Matrix exponential via eigendecomposition
mat_exp = function(A) {
  eig = eigen(A, symmetric = TRUE)
  eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
}


#' Upper-triangular vectorization (including diagonal)
#' Extracts upper triangle column-by-column (R default)
vech = function(M) {
  M[upper.tri(M, diag = TRUE)]
}


#' Inverse vech: reconstruct symmetric matrix from vech vector
vech_inv = function(v, p) {
  M = matrix(0, p, p)
  M[upper.tri(M, diag = TRUE)] = v
  M = M + t(M)
  diag(M) = diag(M) / 2  # diagonal was doubled by adding transpose
  M
}


#####################################################################################
# EM Components
#####################################################################################

#' Log of multivariate normal density under diagonal covariance
#' Returns n-vector of log densities
#' s_mat: n x q matrix of responses
#' mu_mat: n x q matrix of means
#' sigma2_vec: q-vector of variances (diagonal entries)
log_dnorm_diag = function(s_mat, mu_mat, sigma2_vec) {
  n = nrow(s_mat)
  q = ncol(s_mat)

  # log density = -q/2 log(2pi) - 1/2 sum log(sigma2_j) - 1/2 sum (s_j - mu_j)^2/sigma2_j
  log_const = -q / 2 * log(2 * pi) - 0.5 * sum(log(sigma2_vec))

  # residuals scaled by variance: (s - mu)^2 / sigma2, summed over q dimensions
  # sweep sigma2_vec across columns
  resid_sq = (s_mat - mu_mat)^2
  resid_scaled = sweep(resid_sq, 2, sigma2_vec, FUN = "/")
  log_const - 0.5 * rowSums(resid_scaled)
}


#' E-step: compute n x K responsibility matrix
#' Uses log-sum-exp trick for numerical stability
em_estep = function(s_mat, r_mat, params) {
  n = nrow(s_mat)
  K = length(params$pi)

  # log_prob[i, k] = log(pi_k) + log N(s_i | mu_ik, Sigma_k)
  log_prob = matrix(NA, n, K)

  for (k in 1:K) {
    # mean matrix: beta0_k + beta1_k * r_i (element-wise)
    mu_mat = sweep(sweep(r_mat, 2, params$beta1[k, ], FUN = "*"), 2, params$beta0[k, ], FUN = "+")
    log_prob[, k] = log(params$pi[k]) + log_dnorm_diag(s_mat, mu_mat, params$sigma2[k, ])
  }

  # log-sum-exp for numerical stability
  log_max = apply(log_prob, 1, max)
  log_sum = log_max + log(rowSums(exp(log_prob - log_max)))

  gamma = exp(log_prob - log_sum)

  # safety: clamp to avoid exact 0 or 1
  gamma = pmax(gamma, 1e-300)
  gamma = gamma / rowSums(gamma)

  gamma
}


#' M-step: update parameters from responsibilities
#' gamma: n x K responsibility matrix
em_mstep = function(s_mat, r_mat, gamma) {
  n = nrow(s_mat)
  q = ncol(s_mat)
  K = ncol(gamma)

  # effective cluster sizes
  N_k = colSums(gamma)
  # floor to prevent degenerate clusters
  N_k = pmax(N_k, 2)

  pi_new = N_k / n

  beta0_new = matrix(NA, K, q)
  beta1_new = matrix(NA, K, q)
  sigma2_new = matrix(NA, K, q)

  for (k in 1:K) {
    w = gamma[, k]

    for (j in 1:q) {
      # weighted least squares: s_ij ~ beta0_kj + beta1_kj * r_ij
      # with weights w_i = gamma_ik
      sw = sum(w)
      wr = sum(w * r_mat[, j])
      ws = sum(w * s_mat[, j])
      wrr = sum(w * r_mat[, j]^2)
      wrs = sum(w * r_mat[, j] * s_mat[, j])

      # WLS closed form for simple linear regression
      denom = sw * wrr - wr^2
      if (abs(denom) < 1e-12) {
        # degenerate case: fall back to intercept-only
        beta0_new[k, j] = ws / sw
        beta1_new[k, j] = 0
      } else {
        beta1_new[k, j] = (sw * wrs - wr * ws) / denom
        beta0_new[k, j] = (ws - beta1_new[k, j] * wr) / sw
      }

      # weighted residual variance
      resid = s_mat[, j] - beta0_new[k, j] - beta1_new[k, j] * r_mat[, j]
      sigma2_new[k, j] = sum(w * resid^2) / sw
    }
  }

  # variance floor
  sigma2_new = pmax(sigma2_new, 1e-6)

  list(
    pi = pi_new,
    beta0 = beta0_new,
    beta1 = beta1_new,
    sigma2 = sigma2_new
  )
}


#' Observed data log-likelihood (log-sum-exp)
log_lik_obs = function(s_mat, r_mat, params) {
  n = nrow(s_mat)
  K = length(params$pi)

  log_prob = matrix(NA, n, K)

  for (k in 1:K) {
    mu_mat = sweep(sweep(r_mat, 2, params$beta1[k, ], FUN = "*"), 2, params$beta0[k, ], FUN = "+")
    log_prob[, k] = log(params$pi[k]) + log_dnorm_diag(s_mat, mu_mat, params$sigma2[k, ])
  }

  log_max = apply(log_prob, 1, max)
  sum(log_max + log(rowSums(exp(log_prob - log_max))))
}


#' Label alignment: permute state labels to best match a reference fit
#' Minimizes sum_k ||beta1_new[perm[k],] - beta1_ref[k,]||^2
align_labels = function(params_new, params_ref, K) {
  if (K == 1) return(params_new)

  if (K <= 4) {
    # enumerate all K! permutations
    perms = combinat_perms(K)
    best_cost = Inf
    best_perm = 1:K

    for (p in perms) {
      cost = sum((params_new$beta1[p, , drop = FALSE] - params_ref$beta1)^2)
      if (cost < best_cost) {
        best_cost = cost
        best_perm = p
      }
    }
  } else {
    # greedy matching for K > 4
    # compute cost matrix
    cost_mat = matrix(0, K, K)
    for (i in 1:K) {
      for (j in 1:K) {
        cost_mat[i, j] = sum((params_new$beta1[i, ] - params_ref$beta1[j, ])^2)
      }
    }
    best_perm = rep(0, K)
    used = logical(K)
    for (step in 1:K) {
      # find unassigned (i, j) pair with minimum cost
      min_val = Inf
      min_i = min_j = 0
      for (i in 1:K) {
        if (best_perm[i] != 0) next
        for (j in 1:K) {
          if (used[j]) next
          if (cost_mat[i, j] < min_val) {
            min_val = cost_mat[i, j]
            min_i = i
            min_j = j
          }
        }
      }
      best_perm[min_j] = min_i
      used[min_j] = TRUE
    }
    # best_perm[j] = i means new's cluster i maps to ref's cluster j
    # we need the inverse: for ref cluster k, which new cluster to use
    inv_perm = rep(0, K)
    for (j in 1:K) inv_perm[j] = best_perm[j]
    best_perm = inv_perm
  }

  # apply permutation
  params_aligned = list(
    pi = params_new$pi[best_perm],
    beta0 = params_new$beta0[best_perm, , drop = FALSE],
    beta1 = params_new$beta1[best_perm, , drop = FALSE],
    sigma2 = params_new$sigma2[best_perm, , drop = FALSE]
  )
  params_aligned
}


#' Generate all permutations of 1:K (for K <= 4)
combinat_perms = function(K) {
  if (K == 1) return(list(1))
  perms = list()
  sub = combinat_perms(K - 1)
  for (s in sub) {
    for (pos in 1:K) {
      new_perm = append(s, K, after = pos - 1)
      perms = c(perms, list(new_perm))
    }
  }
  perms
}


#####################################################################################
# Main EM Function
#####################################################################################

#' Fit Latent-State Matrix-Log Covariance Regression via EM
#'
#' @param s_mat n x q matrix of vectorized log task covariances
#' @param r_mat n x q matrix of vectorized log rest covariances
#' @param K number of latent states
#' @param n_init number of random restarts
#' @param tol convergence tolerance on log-likelihood
#' @param max_iter maximum EM iterations
#' @param verbose print progress
#' @return list with params, loglik trajectory, convergence info
em_lcr = function(s_mat, r_mat, K, n_init = 20, tol = 1e-4, max_iter = 500,
                  min_iter = 20, init_method = "kmeans",
                  verbose = FALSE, init_params = NULL) {
  n = nrow(s_mat)
  q = ncol(s_mat)

  best_ll = -Inf
  best_result = NULL

  n_runs = if (!is.null(init_params)) 1 else n_init

  for (init in 1:n_runs) {

    # initialize parameters
    if (!is.null(init_params) && init == 1) {
      params = init_params
    } else {
      params = initialize_params(s_mat, r_mat, K, seed = init * 137 + 7,
                                 method = init_method)
    }

    # EM iterations
    ll_vec = rep(NA, max_iter)
    converged = FALSE

    for (iter in 1:max_iter) {
      # E-step
      gamma = em_estep(s_mat, r_mat, params)

      # M-step
      params = em_mstep(s_mat, r_mat, gamma)

      # observed log-likelihood
      ll_vec[iter] = log_lik_obs(s_mat, r_mat, params)

      # convergence check: absolute change in log-likelihood, after min_iter
      if (iter >= min_iter) {
        abs_change = abs(ll_vec[iter] - ll_vec[iter - 1])
        if (abs_change < tol) {
          converged = TRUE
          break
        }
      }

      if (verbose && iter %% 50 == 0) {
        cat(sprintf("  init %d, iter %d, ll = %.4f\n", init, iter, ll_vec[iter]))
      }
    }

    n_iter = iter
    ll_vec = ll_vec[1:n_iter]

    if (ll_vec[n_iter] > best_ll) {
      best_ll = ll_vec[n_iter]
      best_result = list(
        params = params,
        loglik = ll_vec,
        converged = converged,
        n_iter = n_iter,
        K = K,
        gamma = gamma
      )
    }
  }

  best_result
}


#' Initialize EM parameters
#' method = "kmeans": K-means on s_mat then OLS within clusters (smart, fast convergence)
#' method = "random": random responsibilities drawn from Dirichlet (slower, shows full path)
initialize_params = function(s_mat, r_mat, K, seed = 42, method = "kmeans") {
  set.seed(seed)
  n = nrow(s_mat)
  q = ncol(s_mat)

  if (method == "random") {
    # draw random responsibilities from a flat Dirichlet, then do one M-step
    gamma_init = matrix(runif(n * K), n, K)
    gamma_init = gamma_init / rowSums(gamma_init)
    return(em_mstep(s_mat, r_mat, gamma_init))
  }

  # K-means on s_mat to get initial cluster assignments
  km = tryCatch(
    kmeans(s_mat, centers = K, nstart = 5, iter.max = 50),
    error = function(e) NULL
  )

  if (is.null(km)) {
    z_init = sample(1:K, n, replace = TRUE)
  } else {
    z_init = km$cluster
  }

  # ensure all clusters have at least 2 members
  for (k in 1:K) {
    if (sum(z_init == k) < 2) {
      largest = which.max(table(z_init))
      idx = which(z_init == largest)
      z_init[idx[1:2]] = k
    }
  }

  pi = rep(NA, K)
  beta0 = matrix(NA, K, q)
  beta1 = matrix(NA, K, q)
  sigma2 = matrix(NA, K, q)

  for (k in 1:K) {
    idx = which(z_init == k)
    pi[k] = length(idx) / n

    for (j in 1:q) {
      s_j = s_mat[idx, j]
      r_j = r_mat[idx, j]

      if (length(idx) >= 3 && var(r_j) > 1e-10) {
        fit = lm.fit(cbind(1, r_j), s_j)
        beta0[k, j] = fit$coefficients[1]
        beta1[k, j] = fit$coefficients[2]
        sigma2[k, j] = max(mean(fit$residuals^2), 1e-6)
      } else {
        beta0[k, j] = mean(s_j)
        beta1[k, j] = 0.5
        sigma2[k, j] = var(s_j) + 0.1
      }
    }

    # add small perturbation to break symmetry across restarts
    beta0[k, ] = beta0[k, ] + rnorm(q, 0, 0.1)
    beta1[k, ] = beta1[k, ] + rnorm(q, 0, 0.05)
  }

  list(pi = pi, beta0 = beta0, beta1 = beta1, sigma2 = sigma2)
}


#####################################################################################
# Model Selection
#####################################################################################

#' Fit for K = 1:K_max and return BIC/AIC table
#' @param s_mat n x q response matrix
#' @param r_mat n x q predictor matrix
#' @param K_max maximum number of states to try
#' @param ... additional arguments passed to em_lcr
select_K = function(s_mat, r_mat, K_max = 4, ...) {
  n = nrow(s_mat)
  q = ncol(s_mat)

  results = data.frame(
    K = 1:K_max,
    loglik = NA,
    n_params = NA,
    AIC = NA,
    BIC = NA,
    converged = NA,
    n_iter = NA
  )

  fits = list()

  for (K in 1:K_max) {
    fit = em_lcr(s_mat, r_mat, K, ...)  # inherits tol/min_iter from ...

    # number of free parameters:
    # (K-1) mixing proportions + K*q intercepts + K*q slopes + K*q variances
    n_params = (K - 1) + 3 * K * q

    ll = fit$loglik[fit$n_iter]

    results$loglik[K] = ll
    results$n_params[K] = n_params
    results$AIC[K] = -2 * ll + 2 * n_params
    results$BIC[K] = -2 * ll + log(n) * n_params
    results$converged[K] = fit$converged
    results$n_iter[K] = fit$n_iter

    fits[[K]] = fit
  }

  list(table = results, fits = fits)
}
