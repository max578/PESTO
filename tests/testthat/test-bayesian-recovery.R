# Posterior- and prior-recovery standards for the Bayesian & Monte Carlo
# category (srr BS7), plus posterior-estimate validation (BS4.2).

#' @srrstats {BS4.2} Posterior estimates are validated: on problems with a known
#'   answer the posterior mean recovers the truth / the closed-form posterior
#'   (here and in test-correctness-analytic.R), and the spread-ESS / coverage
#'   diagnostics validate the dispersion.
#' @srrstats {BS7.0} Recovery of the parametric estimates of a prior: the
#'   empirical moments of a prior ensemble drawn from N(mu, Sigma) recover
#'   mu and Sigma within tolerance.
#' @srrstats {BS7.1} Recovery of a prior in the absence of additional data: with
#'   non-informative observations (obs_sd = Inf -> zero weight) the posterior
#'   ensemble equals the prior.
#' @srrstats {BS7.2} Recovery of the expected posterior given a prior and data:
#'   on a linear-Gaussian problem the posterior mean matches the closed-form
#'   Gaussian-conjugate posterior (the primary test is in
#'   test-correctness-analytic.R; a compact check is repeated here).
#' @srrstats {BS7.4} Fitted (simulated-observation) values are returned on the
#'   same scale as the input observations.
#' @srrstats {BS7.4a} The scale assumption is tested explicitly: with
#'   observations on a non-zero, non-unit scale the posterior-predictive
#'   ensemble is recovered on that same scale (not standardised).
#' @noRd
NULL

test_that("empirical prior moments recover the generating moments (BS7.0)", {
  set.seed(20L)
  mu    <- c(0.3, 0.1, 1.5)
  sds   <- c(0.5, 0.2, 0.8)
  nreal <- 4000L
  prior <- sapply(seq_along(mu), function(j) stats::rnorm(nreal, mu[j], sds[j]))
  expect_equal(colMeans(prior), mu, tolerance = 0.05)
  expect_equal(apply(prior, 2L, stats::sd), sds, tolerance = 0.05)
})

test_that("non-informative data leaves the posterior equal to the prior (BS7.1)", {
  set.seed(21L)
  npar <- 3L; nobs <- 6L; nreal <- 150L
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  obs <- stats::setNames(stats::rnorm(nobs), paste0("o", seq_len(nobs)))
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                            obs_sd = Inf, noptmax = 4L, verbose = FALSE)
  # Zero observation weight => no update => posterior is exactly the prior.
  expect_equal(as.matrix(fit$par_ensemble[, -1L]), prior, ignore_attr = TRUE)
})

test_that("posterior mean matches the closed-form conjugate posterior (BS7.2)", {
  set.seed(22L)
  npar <- 3L; nobs <- 24L; nreal <- 400L; sigma <- 0.03
  G   <- matrix(stats::rnorm(nobs * npar), nobs, npar) / sqrt(npar)
  tt  <- c(0.8, -1.2, 1.5)
  y   <- as.numeric(G %*% tt) + stats::rnorm(nobs, sd = sigma)
  obs <- stats::setNames(y, paste0("o", seq_len(nobs)))
  prior <- matrix(stats::rnorm(nreal * npar, sd = 4), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  set.seed(7L)
  fit <- pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                            obs_sd = sigma, noptmax = 30L, verbose = FALSE)
  post_mean <- colMeans(as.matrix(fit$par_ensemble[, -1L]))
  mu0    <- colMeans(prior)
  P0_inv <- diag(1 / apply(prior, 2L, stats::var))
  prec   <- P0_inv + crossprod(G) / sigma^2
  mu_pos <- drop(solve(prec, P0_inv %*% mu0 + crossprod(G, y) / sigma^2))
  expect_equal(unname(post_mean), unname(mu_pos), tolerance = 1e-2)
})

test_that("fitted values are returned on the observation scale (BS7.4/7.4a)", {
  set.seed(23L)
  npar <- 3L; nobs <- 8L; nreal <- 120L
  G  <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  tt <- c(1.0, -0.5, 2.0)
  # Observations on a deliberately non-zero, non-unit scale.
  offset <- 100
  y   <- as.numeric(G %*% tt) + offset + stats::rnorm(nobs, sd = 0.1)
  obs <- stats::setNames(y, paste0("o", seq_len(nobs)))
  fn  <- function(theta) theta %*% t(G) + offset
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_callback(fn, prior, obs, obs_sd = 0.1,
                            noptmax = 8L, verbose = FALSE)
  fitted_mean <- colMeans(as.matrix(fit$obs_ensemble[, -1L]))
  # On the same (offset ~100) scale as the observations, and a close fit.
  expect_gt(min(fitted_mean), offset - 10)
  expect_lt(max(abs(fitted_mean - y)), 1)
})
