# Correctness tests against independent oracles: a closed-form analytic
# solution and a fixed-version external implementation (pestpp-ies 5.2.16).
# The reference is never the implementation under test (Independent Oracle
# Principle) -- the analytic posterior is derived from first principles and the
# pestpp-ies golden was produced by separate software.

#' @srrstats {G5.0} Tests use inverse problems with analytically known
#'   properties as the field-standard analogue of canonical reference data sets:
#'   a linear-Gaussian problem with a known generating `theta_true` and a
#'   closed-form least-squares / Gaussian-conjugate posterior (here), plus ODEs
#'   with closed-form solutions (test-ode-templates.R).
#' @srrstats {G5.4} **Correctness.** On a linear-Gaussian inverse problem the
#'   IES posterior mean converges to the closed-form solution, matching both the
#'   ordinary-least-squares estimate and the Gaussian-conjugate posterior mean
#'   to within 1e-3.
#' @srrstats {G5.4a} The reference is an independent analytic oracle (OLS /
#'   conjugate posterior), not the implementation under test. Correctness is
#'   corroborated by cross-implementation agreement elsewhere: the localised vs
#'   standard GLM kernel agree to 1e-10 (test-inflation-localisation.R) and the
#'   RK4 vs deSolve ODE solvers agree (test-ode-templates.R).
#' @srrstats {G5.4b} **Comparison with a previous implementation.** The IES
#'   posterior mean matches that of fixed-version pestpp-ies 5.2.16 on the frozen
#'   `linear_p20_n50` benchmark problem, using a stored golden fixture
#'   (`fixtures/pestpp_ies_tier1_golden.rds`; provenance in its `$provenance`).
#'   The full reproducible multi-tool benchmark vs PEST 18.25 + pestpp-ies 5.2.16
#'   underlies G1.5 / G1.6.
#' @srrstats {G5.5} Every correctness test fixes the RNG via set.seed().
#' @noRd
NULL

test_that("IES posterior mean matches the closed-form least-squares solution", {
  set.seed(2026L)
  npar <- 4L; nobs <- 30L; nreal <- 500L; sigma <- 0.02
  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar) / sqrt(npar)
  theta_true <- c(1.0, -0.5, 2.0, 0.25)
  y   <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = sigma)
  obs <- stats::setNames(y, paste0("o", seq_len(nobs)))

  # Broad prior so the data dominates: the posterior mean approaches the
  # generalised-least-squares solution independently of the prior-covariance
  # convention, giving a convention-robust analytic oracle.
  prior <- matrix(stats::rnorm(nreal * npar, sd = 5), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))

  set.seed(7L)
  fit <- pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                            obs_sd = sigma, noptmax = 30L, verbose = FALSE)
  post_mean <- colMeans(as.matrix(fit$par_ensemble[, -1L]))

  # (i) Ordinary least squares: (G'G)^{-1} G'y (R = sigma^2 I cancels).
  ols <- drop(solve(crossprod(G), crossprod(G, y)))
  expect_equal(unname(post_mean), unname(ols), tolerance = 1e-3)

  # (ii) Gaussian-conjugate posterior mean from the empirical prior moments and
  # the observation covariance R = sigma^2 I.
  mu0     <- colMeans(prior)
  P0_inv  <- diag(1 / apply(prior, 2L, stats::var))
  GtRinvG <- crossprod(G) / sigma^2
  mu_post <- drop(solve(P0_inv + GtRinvG,
                        P0_inv %*% mu0 + crossprod(G, y) / sigma^2))
  expect_equal(unname(post_mean), unname(mu_post), tolerance = 1e-3)
})

test_that("IES posterior mean matches fixed-version pestpp-ies 5.2.16 (golden)", {
  g <- readRDS(test_path("fixtures", "pestpp_ies_tier1_golden.rds"))
  obs <- stats::setNames(g$obs, sprintf("o%02d", seq_along(g$obs)))

  set.seed(g$rng_seed)
  fit <- pesto_ies_callback(function(theta) theta %*% t(g$X), g$prior, obs,
                            obs_sd = g$obs_sd, noptmax = g$noptmax,
                            lambda = 1.0, verbose = FALSE)
  post_mean <- colMeans(as.matrix(fit$par_ensemble[, -1L]))

  # Agreement with the independently-produced pestpp-ies posterior mean
  # (observed RMSE ~0.005, maxabs ~0.011 on parameters of scale ~1.3).
  expect_lt(sqrt(mean((post_mean - g$pestpp_posterior_mean)^2)), 0.05)
  expect_lt(max(abs(post_mean - g$pestpp_posterior_mean)), 0.05)

  # And PESTO sits within the benchmarked parity band relative to the truth.
  expect_lt(sqrt(mean((post_mean - g$truth)^2)), 0.05)
})
