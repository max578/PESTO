# Parameter-recovery (across seeds) and algorithm-performance (data-size
# scaling) tests for the IES driver.

#' @srrstats {G5.6} **Parameter recovery.** On data simulated from a known
#'   parameter vector, the IES posterior mean recovers the generating
#'   parameters (corroborated by the linear-Gaussian, crop-growth and SEIR
#'   recoveries in test-ies-callback.R / test-ode-templates.R).
#' @srrstats {G5.6a} Recovery is asserted within a tolerance, never as exact
#'   equality.
#' @srrstats {G5.6b} Recovery is confirmed across multiple random seeds (both
#'   the simulated data and the prior ensemble are stochastic): the regular
#'   suite checks three seeds; the extended suite checks eight.
#' @srrstats {G5.7} **Algorithm performance.** Recovery error decreases as the
#'   number of observations grows -- the expected behaviour of a consistent
#'   estimator -- and the filter driver's posterior spread tightens as windows
#'   accrue (test-ies-filter.R).
#' @noRd
NULL

.linear_recovery_rmse <- function(seed, npar = 3L, nobs = 9L, nreal = 100L,
                                  sigma = 0.05, noptmax = 8L) {
  set.seed(seed)
  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0, 0.3, -1.1)[seq_len(npar)]
  y   <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = sigma)
  obs <- stats::setNames(y, paste0("o", seq_len(nobs)))
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  set.seed(seed + 500L)
  fit <- pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                            obs_sd = sigma, noptmax = noptmax, verbose = FALSE)
  sqrt(mean((colMeans(as.matrix(fit$par_ensemble[, -1L])) - theta_true)^2))
}

test_that("parameter recovery succeeds across several seeds", {
  errs <- vapply(1:3, .linear_recovery_rmse, numeric(1L))
  expect_true(
    all(errs < 0.05),
    info = sprintf("recovery RMSE per seed: %s",
                   paste(signif(errs, 3L), collapse = ", "))
  )
})

test_that("parameter recovery holds across many seeds (extended)", {
  skip_if_not_extended()
  errs <- vapply(1:8, .linear_recovery_rmse, numeric(1L))
  expect_true(
    all(errs < 0.05),
    info = sprintf("recovery RMSE per seed: %s",
                   paste(signif(errs, 3L), collapse = ", "))
  )
})

test_that("recovery error decreases as the number of observations grows", {
  scaling_rmse <- function(nobs, seed) {
    .linear_recovery_rmse(seed, npar = 3L, nobs = nobs, nreal = 150L,
                          sigma = 0.1, noptmax = 10L)
  }
  # Average over seeds so the trend is robust, not a single-draw artefact.
  small <- mean(vapply(1:4, function(s) scaling_rmse(4L,  s), numeric(1L)))
  large <- mean(vapply(1:4, function(s) scaling_rmse(40L, s), numeric(1L)))
  expect_lt(large, small)           # more data -> better recovery
  expect_lt(large, 0.5 * small)     # and materially so
})
