#' Pure-R Reference IES Update — Chen & Oliver (2013) eq. 12
#'
#' A textbook, pure-R implementation of one Iterative Ensemble Smoother
#' (IES) parameter upgrade as published in Chen & Oliver (2013) eq. 12.
#' This function is independent of [ensemble_solution()] (the C++
#' kernel) and exists for two purposes:
#'
#' 1. **Independent numerical validation.** PESTO's C++ kernel can be
#'    cross-checked against this reference at machine precision (the
#'    package's regression test
#'    `tests/testthat/test-ensemble-solution-sign.R` and the
#'    paired-seed protocol in `inst/scripts/i2_paired_seed_check.R`
#'    both rely on this).
#' 2. **Self-contained `pestpp-ies` comparison.** Vignettes can compare
#'    PESTO native IES to this textbook reference without requiring the
#'    upstream `pestpp-ies` binary, making the comparison story
#'    reproducible on CRAN check farms and other binary-free
#'    environments.
#'
#' The implementation is deliberately pedagogical — readability over
#' speed. For production work use [ensemble_solution()] (the C++
#' kernel; \eqn{\sim 35\times} faster on production-scale ensembles).
#'
#' @section Sign convention:
#' This reference uses the textbook sign `obs_resid = obs - sim` with a
#' positive leading sign in the upgrade. PESTO's C++ kernel uses the
#' equivalent `obs_resid = sim - obs` with a leading negative sign
#' (see `?ensemble_solution`). Both yield identical upgrades; the
#' difference is purely cosmetic. The two are cross-validated to
#' machine precision in
#' `inst/scripts/i2_paired_seed_check.R`.
#'
#' @param par_ensemble Numeric matrix `(n_par x n_real)`. Current
#'   parameter ensemble (each column is one realisation).
#' @param obs_ensemble Numeric matrix `(n_obs x n_real)`. Simulated
#'   observations from the current ensemble.
#' @param obs_target Numeric vector `(n_obs)`. Target (measured)
#'   observation values.
#' @param weights Numeric vector `(n_obs)`. Observation weights
#'   (typically `1 / sd(obs_noise)`).
#' @param lambda Numeric scalar. Marquardt damping parameter. `1.0` is
#'   the textbook GLM update; larger values dampen more aggressively.
#' @return A numeric matrix `(n_par x n_real)`: the parameter upgrade
#'   (add to `par_ensemble` to get the next iterate).
#'
#' @references
#' Chen, Y. & Oliver, D. S. (2013). Levenberg-Marquardt forms of the
#' iterative ensemble smoother for efficient history matching and
#' uncertainty quantification. *Computational Geosciences*, 17(4),
#' 689--703. \doi{10.1007/s10596-013-9351-5}
#'
#' Evensen, G. (2018). Analysis of iterative ensemble smoothers for
#' solving inverse problems. *Computational Geosciences*, 22(3),
#' 885--908. \doi{10.1007/s10596-018-9731-y}
#'
#' @seealso [ensemble_solution()] for the production C++ kernel.
#'
#' @examples
#' # Simple linear inverse problem: identify the true theta given noisy obs.
#' set.seed(20260425L)
#' n_par <- 4L; n_obs <- 12L; n_real <- 30L
#' G <- matrix(rnorm(n_obs * n_par) / sqrt(n_par), n_obs, n_par)
#' theta_true <- rnorm(n_par)
#' y_obs <- as.numeric(G %*% theta_true) + rnorm(n_obs, sd = 0.05)
#' weights <- rep(1 / 0.05, n_obs)
#'
#' par_ens <- matrix(rnorm(n_par * n_real, sd = 1), n_par, n_real)
#' obs_ens <- G %*% par_ens
#'
#' # One textbook IES upgrade.
#' upg <- pesto_reference_ies(
#'   par_ensemble = par_ens,
#'   obs_ensemble = obs_ens,
#'   obs_target   = y_obs,
#'   weights      = weights,
#'   lambda       = 1.0
#' )
#' dim(upg)
#'
#' # Apply the upgrade and check phi reduces.
#' par_next <- par_ens + upg
#' obs_next <- G %*% par_next
#' Y <- matrix(rep(y_obs, n_real), n_obs, n_real)
#' phi0 <- mean(compute_phi(Y - obs_ens,  weights))
#' phi1 <- mean(compute_phi(Y - obs_next, weights))
#' phi1 < phi0
#'
#' @export
pesto_reference_ies <- function(par_ensemble,
                                obs_ensemble,
                                obs_target,
                                weights,
                                lambda = 1.0) {
  .check_reference_ies_inputs(
    par_ensemble, obs_ensemble, obs_target, weights, lambda
  )

  # Ensemble geometry ----------------------------------------------------
  n_real <- ncol(par_ensemble)
  scale  <- 1 / sqrt(n_real - 1)
  w_diag <- diag(weights)

  par_mean <- rowMeans(par_ensemble)
  obs_mean <- rowMeans(obs_ensemble)
  par_diff <- par_ensemble - par_mean
  obs_diff <- obs_ensemble - obs_mean

  # Residuals (textbook obs - sim sign convention) -----------------------
  y_resid <- matrix(rep(obs_target, n_real), nrow = length(obs_target)) -
             obs_ensemble

  # SVD of the scaled, weighted observation difference matrix ------------
  svd_res  <- svd(scale * w_diag %*% obs_diff)
  sigma    <- svd_res$d
  inv_term <- 1 / (sigma^2 + (lambda + 1))

  upgrade_factor <-
    svd_res$v %*% diag(sigma * inv_term, nrow = length(sigma)) %*%
    t(svd_res$u) %*% (w_diag %*% y_resid)

  upgrade <- (scale * par_diff) %*% upgrade_factor
  upgrade
}


# Internal helpers ------------------------------------------------------

#' Validate inputs to `pesto_reference_ies()`
#'
#' Checks matrix shapes, vector-length agreement, and `lambda`
#' admissibility. Called from the function preamble so the body reads
#' as the IES update equation.
#'
#' @noRd
#' @keywords internal
.check_reference_ies_inputs <- function(par_ensemble, obs_ensemble,
                                        obs_target, weights, lambda) {
  .assert_matrix(par_ensemble, "par_ensemble")
  .assert_matrix(obs_ensemble, "obs_ensemble")
  .assert_same_ncol(par_ensemble, obs_ensemble,
    "par_ensemble", "obs_ensemble"
  )
  if (length(obs_target) != nrow(obs_ensemble)) {
    stop(
      sprintf(
        "`obs_target` length (%d) must equal `nrow(obs_ensemble)` (%d).",
        length(obs_target), nrow(obs_ensemble)
      ),
      call. = FALSE
    )
  }
  if (length(weights) != length(obs_target)) {
    stop(
      sprintf(
        "`weights` (length %d) must match `obs_target` (length %d).",
        length(weights), length(obs_target)
      ),
      call. = FALSE
    )
  }
  .assert_nonneg_scalar(lambda, "lambda")
  invisible(TRUE)
}
