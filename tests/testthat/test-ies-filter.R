# Tests for the sequential (filter-mode) IES driver: recovery, the
# in-season tightening property, window validation, the iterated filter,
# multi-fidelity scheduling, failure tolerance, and manifest support.

filter_problem <- function(seed = 1L, npar = 3L, nobs = 9L, nreal = 80L) {
  set.seed(seed)
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0)[seq_len(npar)]
  y <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05)
  names(y) <- paste0("o", seq_len(nobs))
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  list(G = G, fn = function(theta) theta %*% t(G), obs = y,
       prior = prior, theta_true = theta_true)
}

test_that("filter recovers the truth across windows", {
  p <- filter_problem()
  fit <- pesto_ies_filter(p$fn, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:3, 4:6, 7:9), verbose = FALSE)
  est <- unname(colMeans(as.matrix(fit$par_ensemble[, -1L])))
  expect_equal(est, p$theta_true, tolerance = 0.1)
  expect_s3_class(fit, "pesto_ies_filter_result")
  expect_s3_class(fit, "pesto_ies_result")
  expect_length(fit$windows, 3L)
})

test_that("posterior sd tightens as windows accrue", {
  p <- filter_problem()
  fit <- pesto_ies_filter(p$fn, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:3, 4:6, 7:9), verbose = FALSE)
  sd_trace <- vapply(fit$windows, function(w) mean(w$par_sd), numeric(1))
  # Non-increasing: each window assimilates more, never widening.
  expect_true(all(diff(sd_trace) <= 1e-8))
  # And it genuinely shrinks from first to last.
  expect_lt(sd_trace[length(sd_trace)], sd_trace[1L])
})

test_that("windows are validated (range, disjointness, type)", {
  p <- filter_problem(nobs = 6L, npar = 2L)
  args <- list(forward_model = p$fn, prior_ensemble = p$prior,
               obs = p$obs, obs_sd = 0.05, verbose = FALSE)
  expect_error(do.call(pesto_ies_filter, c(args, list(windows = 1:3))),
               "list")
  expect_error(
    do.call(pesto_ies_filter, c(args, list(windows = list(1:3, 4:9)))),
    "in \\[1, 6\\]"
  )
  expect_error(
    do.call(pesto_ies_filter, c(args, list(windows = list(1:3, 3:5)))),
    "disjoint"
  )
  expect_error(
    do.call(pesto_ies_filter, c(args, list(windows = list(c(1L, 1L, 2L))))),
    "duplicate"
  )
})

test_that("iterated filter (window_noptmax > 1) runs and recovers", {
  p <- filter_problem()
  fit <- pesto_ies_filter(p$fn, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:5, 6:9), window_noptmax = 3L,
                          verbose = FALSE)
  est <- unname(colMeans(as.matrix(fit$par_ensemble[, -1L])))
  expect_equal(est, p$theta_true, tolerance = 0.1)
})

test_that("multi-fidelity filter records provenance and honours schedule", {
  p <- filter_problem()
  cheap <- function(theta) theta %*% t(p$G) + 0.3
  mf <- pesto_multifidelity_model(
    levels = list(
      pesto_forward_model(fn = cheap, n_obs = 9L, fidelity = 0L),
      pesto_forward_model(fn = p$fn,  n_obs = 9L, fidelity = 1L)
    ),
    costs = c(1, 20)
  )
  fit <- pesto_ies_filter(mf, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:3, 4:6, 7:9),
                          fidelity_schedule = c(0L, 1L, 1L), verbose = FALSE)
  expect_identical(fit$fidelity$type, "multifidelity")
  expect_identical(fit$fidelity$schedule, c(0L, 1L, 1L))
  expect_identical(fit$fidelity$final_level, 1L)
})

test_that("on_failure='na' tolerates failed realisations", {
  p <- filter_problem()
  flaky <- function(theta) {
    r <- theta %*% t(p$G)
    r[1L, ] <- NA            # one persistent failure
    r
  }
  fit <- pesto_ies_filter(flaky, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:4, 5:9), on_failure = "na",
                          verbose = FALSE)
  expect_gt(fit$failure_rate, 0)
  expect_true(all(vapply(fit$windows, function(w) w$n_failures >= 1L,
                         logical(1))))
})

test_that("a filter result becomes an ies_filter manifest that verifies", {
  skip_if_not_installed("yaml")
  p <- filter_problem()
  fit <- pesto_ies_filter(p$fn, p$prior, p$obs, obs_sd = 0.05,
                          windows = list(1:3, 4:6, 7:9), verbose = FALSE)
  m <- as_manifest(fit, run_id = "filt_rt")
  expect_identical(m@method, "ies_filter")
  dir <- tempfile("filt_manifest_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  yaml_path <- file.path(dir, "m.yaml")
  write_manifest(m, yaml_path)
  back <- read_manifest(yaml_path)
  expect_identical(back@method, "ies_filter")
  expect_true(isTRUE(verify_manifest(back)$ok))
})
