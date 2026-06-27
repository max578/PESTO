# Noise-susceptibility / stochastic-stability tests for the IES driver.

#' @srrstats {G5.9} Stochastic behaviour is tested for stability:
#' @srrstats {G5.9a} adding trivial noise (at the scale of `.Machine$double.eps`)
#'   to the observations does not meaningfully change the posterior mean;
#' @srrstats {G5.9b} running under different random seeds / prior-ensemble draws
#'   does not meaningfully change the recovered posterior mean.
#' @noRd
NULL

test_that("G5.9a trivial observation noise does not change the posterior", {
  set.seed(3L)
  npar <- 3L; nobs <- 8L; nreal <- 60L
  G   <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  y   <- as.numeric(G %*% c(1, -0.5, 2)) + stats::rnorm(nobs, sd = 0.05)
  obs <- stats::setNames(y, paste0("o", seq_len(nobs)))
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))

  run <- function(o) {
    set.seed(99L)
    colMeans(as.matrix(
      pesto_ies_callback(function(theta) theta %*% t(G), prior, o,
                         obs_sd = 0.05, noptmax = 6L,
                         verbose = FALSE)$par_ensemble[, -1L]
    ))
  }
  m0 <- run(obs)
  m1 <- run(obs + .Machine$double.eps * abs(obs))
  expect_lt(max(abs(m0 - m1)), 1e-6)
})

test_that("G5.9b different prior-ensemble seeds give a stable posterior mean", {
  set.seed(1L)
  npar <- 3L; nobs <- 8L
  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1, -0.5, 2)
  obs <- stats::setNames(
    as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05),
    paste0("o", seq_len(nobs))
  )

  mean_for_seed <- function(s) {
    set.seed(s)
    prior <- matrix(stats::rnorm(200L * npar), 200L, npar,
                    dimnames = list(NULL, paste0("p", seq_len(npar))))
    colMeans(as.matrix(
      pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                         obs_sd = 0.05, noptmax = 8L,
                         verbose = FALSE)$par_ensemble[, -1L]
    ))
  }
  means  <- vapply(1:5, mean_for_seed, numeric(npar))
  spread <- max(apply(means, 1L, function(z) diff(range(z))))
  expect_lt(spread, 0.02)            # observed ~5e-4; truth scale ~1.3
})
