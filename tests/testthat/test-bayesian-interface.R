# Interface / return-object standards for the Bayesian & Monte Carlo category:
# starting values, warm-starting, verbosity control, returned convergence
# statistics, and the print / plot methods.

#' @srrstats {BS2.7} Starting values are explicitly controlled: the
#'   `prior_ensemble` passed to the driver *is* the set of starting points (one
#'   row per realisation).
#' @srrstats {BS2.8} The result of a previous run can be used as the starting
#'   point of a subsequent run -- the posterior ensemble is itself a valid prior
#'   ensemble (warm start); exercised here and intrinsic to the sequential
#'   filter.
#' @srrstats {BS2.13} Messages / progress indicators are suppressible via
#'   `verbose = FALSE` while warnings and errors are retained; tested here.
#' @srrstats {BS5.1} The return object carries metadata on the inputs: parameter
#'   and observation names / dimensions (par_ensemble, obs_ensemble) and the
#'   assimilation inputs obs_target / obs_sd / weights.
#' @srrstats {BS5.3} Convergence statistics are returned: the per-realisation
#'   phi trace by iteration plus the per-iteration dispersion diagnostics.
#' @srrstats {BS5.5} Absence of convergence is diagnosable from the returned phi
#'   trace and the spread-ESS ratio (the ensemble-collapse diagnostic).
#' @noRd
NULL

.bi_problem <- function(seed = 1L, npar = 3L, nobs = 6L, nreal = 80L) {
  set.seed(seed)
  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0)[seq_len(npar)]
  obs <- stats::setNames(
    as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05),
    paste0("o", seq_len(nobs))
  )
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  list(fn = function(theta) theta %*% t(G), obs = obs, prior = prior,
       theta_true = theta_true)
}

test_that("the prior ensemble is the controlled set of starting values (BS2.7)", {
  p <- .bi_problem()
  fit <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                            noptmax = 3L, verbose = FALSE)
  # Same number of realisations in == out; parameter names preserved.
  expect_equal(nrow(fit$par_ensemble), nrow(p$prior))
  expect_equal(setdiff(names(fit$par_ensemble), "real_name"), colnames(p$prior))
})

test_that("a previous run's posterior warm-starts the next run (BS2.8)", {
  p <- .bi_problem()
  set.seed(11L)
  f1 <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                           noptmax = 2L, verbose = FALSE)
  warm <- as.matrix(f1$par_ensemble[, -1L])          # posterior as new prior
  set.seed(12L)
  f2 <- pesto_ies_callback(p$fn, warm, p$obs, obs_sd = 0.05,
                           noptmax = 2L, verbose = FALSE)
  e1 <- sqrt(mean((colMeans(as.matrix(f1$par_ensemble[, -1L])) - p$theta_true)^2))
  e2 <- sqrt(mean((colMeans(as.matrix(f2$par_ensemble[, -1L])) - p$theta_true)^2))
  expect_s3_class(f2, "pesto_ies_result")
  expect_lt(e2, e1 + 0.005)                           # warm start does not materially worsen
})

test_that("verbose toggles progress output without hiding errors (BS2.13)", {
  p <- .bi_problem()
  expect_message(
    pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05, noptmax = 2L,
                       verbose = TRUE),
    "pesto_ies_callback"
  )
  expect_no_message(
    pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05, noptmax = 2L,
                       verbose = FALSE)
  )
  # An error still surfaces even with verbose = FALSE.
  expect_error(
    pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = -1, verbose = FALSE),
    "obs_sd"
  )
})

test_that("the return object exposes input metadata and convergence stats (BS5.1/5.3/5.5)", {
  p <- .bi_problem()
  fit <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                            noptmax = 4L, verbose = FALSE)
  # BS5.1 input metadata
  expect_named(fit$obs_target, names(p$obs))
  expect_named(fit$weights, names(p$obs))
  expect_true(all(c("real_name", colnames(p$prior)) %in% names(fit$par_ensemble)))
  # BS5.3 convergence statistics returned
  expect_true(all(c("iteration", "phi") %in% names(fit$phi)))
  mphi <- tapply(fit$phi$phi, fit$phi$iteration, mean)
  expect_lt(mphi[length(mphi)], mphi[1L])             # phi descends
  # BS5.5 collapse diagnostic available per step
  ess <- vapply(fit$iterations, function(it) it$spread_ess_ratio, numeric(1L))
  expect_true(all(is.finite(ess)))
})

test_that("print and plot methods work on a fitted result (BS6.0/6.1)", {
  p   <- .bi_problem()
  fit <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                            noptmax = 3L, verbose = FALSE)
  # S3 dispatch via print()/plot() is unreliable under devtools::load_all() for
  # this S7-importing package (it registers cleanly on install / R CMD check);
  # call the methods directly, matching the package test convention.
  pr  <- getS3method("print", "pesto_ies_result")
  out <- utils::capture.output(res <- pr(fit))
  expect_identical(res, fit)                          # returns x invisibly
  expect_true(any(grepl("pesto_ies_result", out)))
  expect_true(any(grepl("mean phi", out)))
  pl <- getS3method("plot", "pesto_ies_result")
  expect_s3_class(pl(fit), "ggplot")
})
