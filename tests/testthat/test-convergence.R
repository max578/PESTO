# Convergence-checker standards (srr BS1.4, BS4.3, BS4.4, BS4.6, BS4.7) for the
# opt-in phi_tol early-stopping rule of pesto_ies_callback().

#' @srrstats {BS1.4} The convergence checker is opt-in: `phi_tol = NULL`
#'   (default) runs the full `noptmax` with no checker, while a non-NULL
#'   `phi_tol` enables the phi-reduction stopping rule. Both paths are tested.
#' @srrstats {BS4.3} A convergence checker is implemented -- the phi-reduction
#'   stopping rule of White (2018), referenced in `?pesto_ies_callback`:
#'   iteration stops once the relative reduction in mean phi falls below `phi_tol`.
#' @srrstats {BS4.4} Computation can be stopped on convergence (`converged =
#'   TRUE`, `n_iterations < noptmax`), but not by default (`phi_tol = NULL`).
#' @srrstats {BS4.6} A checker-stopped run is identical to a fixed-iteration run
#'   of the same length without the checker -- the checker only halts the loop,
#'   it does not alter the update.
#' @srrstats {BS4.7} The threshold parameter behaves monotonically: a smaller
#'   (stricter) `phi_tol` yields at least as many iterations as a larger one.
#' @noRd
NULL

.conv_linear <- function(seed = 1L, npar = 3L, nobs = 9L, nreal = 100L) {
  set.seed(seed)
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  obs <- stats::setNames(
    as.numeric(G %*% c(1, -0.5, 2)[seq_len(npar)]) +
      stats::rnorm(nobs, sd = 0.05),
    paste0("o", seq_len(nobs))
  )
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  list(fn = function(theta) theta %*% t(G), obs = obs, prior = prior)
}

.conv_crop <- function(seed = 101L, nreal = 120L) {
  set.seed(seed)
  times <- seq(0, 120, by = 15)
  fm <- crop_growth_forward_model(times = times)
  truth_v <- c(r = 0.06, b_max = 1400, b0 = 20)
  biomass <- as.numeric(pesto_evaluate(
    fm, matrix(truth_v, nrow = 1L, dimnames = list(NULL, names(truth_v)))
  ))
  obs <- stats::setNames(biomass + stats::rnorm(length(biomass), sd = 20),
                         paste0("t", seq_along(biomass)))
  prior <- cbind(r = stats::runif(nreal, 0.02, 0.12),
                 b_max = stats::runif(nreal, 900, 2000),
                 b0 = stats::runif(nreal, 5, 60))
  list(fm = fm, obs = obs, prior = prior)
}

test_that("phi_tol = NULL runs the full noptmax with no checker (BS1.4)", {
  p   <- .conv_linear()
  fit <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                            noptmax = 6L, verbose = FALSE)
  expect_false(fit$converged)
  expect_equal(fit$n_iterations, 6L)
  fit2 <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                             noptmax = 6L, phi_tol = NULL, verbose = FALSE)
  expect_identical(as.matrix(fit$par_ensemble[, -1L]),
                   as.matrix(fit2$par_ensemble[, -1L]))
})

test_that("the checker stops the run early on convergence (BS4.3/4.4)", {
  p   <- .conv_crop()
  fit <- pesto_ies_callback(p$fm, p$prior, p$obs, obs_sd = 20,
                            noptmax = 15L, phi_tol = 0.2, verbose = FALSE)
  expect_true(fit$converged)
  expect_lt(fit$n_iterations, 15L)
  expect_equal(length(fit$iterations), fit$n_iterations)   # records trimmed
})

test_that("a checker-stopped run equals the equivalent fixed-length run (BS4.6)", {
  p <- .conv_linear()
  stopped <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                                noptmax = 12L, phi_tol = 0.5, verbose = FALSE)
  fixed <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                              noptmax = stopped$n_iterations, verbose = FALSE)
  expect_identical(as.matrix(stopped$par_ensemble[, -1L]),
                   as.matrix(fixed$par_ensemble[, -1L]))
})

test_that("a stricter tolerance yields more iterations (BS4.7)", {
  p <- .conv_crop()
  loose  <- pesto_ies_callback(p$fm, p$prior, p$obs, obs_sd = 20,
                               noptmax = 15L, phi_tol = 0.5, verbose = FALSE)
  strict <- pesto_ies_callback(p$fm, p$prior, p$obs, obs_sd = 20,
                               noptmax = 15L, phi_tol = 0.05, verbose = FALSE)
  expect_gt(strict$n_iterations, loose$n_iterations)
})

test_that("phi_tol is validated", {
  p <- .conv_linear()
  expect_error(
    pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05, phi_tol = -1,
                       verbose = FALSE),
    "phi_tol"
  )
  expect_error(
    pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05, phi_tol = 0,
                       verbose = FALSE),
    "phi_tol"
  )
})
