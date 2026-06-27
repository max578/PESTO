# Edge-condition and error/return-validity tests for the IES driver and the
# forward-model contract.

#' @srrstats {G5.2} Error and warning behaviour is demonstrated through tests:
#'   malformed inputs raise errors with informative messages, exercised here and
#'   across the per-function validation tests (test-ies-callback.R,
#'   test-forward-model.R, test-multifidelity.R, test-ode-templates.R,
#'   test-manifest.R).
#' @srrstats {G5.2a} Condition messages are distinguishable to their source: a
#'   source scan of every stop()/warning()/message() call found no two full
#'   messages identical -- shared prefixes (e.g. "[apsim_callback] real <i>") are
#'   completed by stage-specific causes (edit / run / extractor failed).
#' @srrstats {G5.2b} Each error path is triggered explicitly and its message
#'   matched (the expect_error() calls here and in the validation tests).
#' @srrstats {G5.3} Return objects contain no missing or undefined values where
#'   none are expected: posterior finiteness is asserted here, and the
#'   ensemble-update kernels are checked for NA/Inf in test-ensemble-solution.R.
#' @srrstats {G5.8} **Edge conditions** produce clear errors or defined
#'   behaviour, covering:
#' @srrstats {G5.8a} zero-length data (empty observations; a zero-realisation
#'   prior ensemble);
#' @srrstats {G5.8b} data of unsupported types (a non-function forward model is
#'   rejected by the typed contract; a non-numeric ensemble raises a coercion
#'   warning rather than passing silently);
#' @srrstats {G5.8c} all-`NA` and all-identical fields (an all-`NA` observation
#'   vector is caught; an all-identical, zero-variance prior column is handled
#'   without producing NA in the posterior);
#' @srrstats {G5.8d} data outside the usual scope (an under-determined problem
#'   with more parameters than realisations runs and returns finite values via
#'   SVD truncation).
#' @noRd
NULL

.edge_problem <- function(seed = 1L, npar = 3L, nobs = 8L, nreal = 60L) {
  set.seed(seed)
  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0)[seq_len(npar)]
  obs <- stats::setNames(
    as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05),
    paste0("o", seq_len(nobs))
  )
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  list(fn = function(theta) theta %*% t(G), obs = obs, prior = prior)
}

test_that("G5.8a zero-length data raises clear errors", {
  p <- .edge_problem()
  expect_error(
    pesto_ies_callback(p$fn, p$prior, numeric(0), obs_sd = 0.05,
                       noptmax = 2L, verbose = FALSE),
    "n_obs|obs"
  )
  expect_error(
    pesto_ies_callback(p$fn, matrix(numeric(0), 0L, 3L), p$obs, obs_sd = 0.05,
                       noptmax = 2L, verbose = FALSE),
    "at least 2 realisations"
  )
})

test_that("G5.8b unsupported input types do not pass silently", {
  p <- .edge_problem()
  # The typed forward-model contract rejects a non-function with a clear error.
  expect_error(pesto_forward_model(fn = 42), "function")
  # A non-numeric ensemble coerces to NA with a warning (not silent success);
  # the run cannot then yield a valid fit, so the warning is the observable.
  char_prior <- matrix(letters[1:9], 3L, 3L,
                       dimnames = list(NULL, c("p1", "p2", "p3")))
  expect_warning(
    tryCatch(
      pesto_ies_callback(p$fn, char_prior, p$obs, obs_sd = 0.05,
                         noptmax = 2L, verbose = FALSE),
      error = function(e) invisible(NULL)
    ),
    "NA|coercion"
  )
})

test_that("G5.8c all-NA and all-identical fields behave correctly", {
  p <- .edge_problem(nreal = 40L)
  # All-NA observations are caught, not silently fit.
  expect_error(
    pesto_ies_callback(
      p$fn, p$prior,
      stats::setNames(rep(NA_real_, 8L), paste0("o", seq_len(8L))),
      obs_sd = 0.05, noptmax = 2L, verbose = FALSE
    ),
    "successful realisations|NA"
  )
  # An all-identical (zero-variance) prior column is handled gracefully and the
  # posterior contains no NA (G5.3).
  ident <- p$prior
  ident[, 2L] <- 7.0
  fit <- pesto_ies_callback(p$fn, ident, p$obs, obs_sd = 0.05,
                            noptmax = 3L, verbose = FALSE)
  expect_false(anyNA(as.matrix(fit$par_ensemble[, -1L])))
})

test_that("G5.8d under-determined problems (npar > nreal) run and stay finite", {
  set.seed(9L)
  npar <- 10L; nobs <- 8L; nreal <- 3L
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  obs <- stats::setNames(
    as.numeric(G %*% stats::rnorm(npar)) + stats::rnorm(nobs, sd = 0.05),
    paste0("o", seq_len(nobs))
  )
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_callback(function(theta) theta %*% t(G), prior, obs,
                            obs_sd = 0.05, noptmax = 3L, verbose = FALSE)
  est <- as.matrix(fit$par_ensemble[, -1L])
  expect_equal(dim(est), c(nreal, npar))
  expect_true(all(is.finite(est)))           # G5.3: no NA/Inf in the return
})
