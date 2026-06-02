# Tests for the pesto_forward_model contract object and its evaluation
# engine: construction / validation, coercion, shape + failure
# accounting, the bulk vs per-row paths, and parallel determinism.

linear_model <- function(seed = 1L, npar = 3L, nobs = 4L) {
  set.seed(seed)
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  list(G = G, fn = function(theta) theta %*% t(G), npar = npar, nobs = nobs)
}

test_that("pesto_forward_model coerces numeric integer fields", {
  fm <- pesto_forward_model(fn = identity, n_obs = 6, fidelity = 2)
  expect_true(S7::S7_inherits(fm, pesto_forward_model))
  expect_identical(fm@n_obs, 6L)
  expect_identical(fm@fidelity, 2L)
})

test_that("validator rejects malformed policy arguments", {
  expect_error(pesto_forward_model(fn = identity, on_failure = "skip"),
               "on_failure")
  expect_error(pesto_forward_model(fn = identity, parallel = "gpu"),
               "parallel")
  expect_error(pesto_forward_model(fn = identity, max_fail_frac = 2),
               "max_fail_frac")
  expect_error(pesto_forward_model(fn = identity, n_obs = 0),
               "n_obs")
  expect_error(pesto_forward_model(fn = identity, fidelity = -1L),
               "fidelity")
})

test_that("as_forward_model wraps functions and passes objects through", {
  fm1 <- as_forward_model(function(theta) theta, on_failure = "stop")
  expect_true(S7::S7_inherits(fm1, pesto_forward_model))
  expect_identical(fm1@on_failure, "stop")
  expect_identical(as_forward_model(fm1), fm1)
})

test_that("pesto_evaluate returns shape + failure attributes", {
  lm <- linear_model()
  fm <- pesto_forward_model(fn = lm$fn, n_obs = lm$nobs)
  theta <- matrix(stats::rnorm(5L * lm$npar), 5L, lm$npar)
  out <- pesto_evaluate(fm, theta)
  expect_equal(dim(out), c(5L, lm$nobs))
  expect_identical(attr(out, "n_failures"), 0L)
  expect_identical(attr(out, "fail_idx"), integer(0))
})

test_that("n_obs is inferred when left NA", {
  lm <- linear_model()
  fm <- pesto_forward_model(fn = lm$fn)            # n_obs = NA
  theta <- matrix(stats::rnorm(4L * lm$npar), 4L, lm$npar)
  expect_equal(ncol(pesto_evaluate(fm, theta)), lm$nobs)
})

test_that("param_names enforce + reorder theta columns", {
  fn <- function(theta) theta[, "a", drop = FALSE] - theta[, "b", drop = FALSE]
  fm <- pesto_forward_model(fn = fn, n_obs = 1L, param_names = c("a", "b"))
  # Supplied in the wrong order: engine reorders to (a, b).
  theta <- matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE,
                  dimnames = list(NULL, c("b", "a")))
  out <- pesto_evaluate(fm, theta)
  expect_equal(as.numeric(out), c(2 - 1, 4 - 3))
  bad <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(NULL, "a"))
  expect_error(pesto_evaluate(fm, bad), "missing parameters")
})

test_that("bulk and per-row paths agree numerically", {
  lm <- linear_model(seed = 7L)
  theta <- matrix(stats::rnorm(6L * lm$npar), 6L, lm$npar)
  bulk <- pesto_forward_model(fn = lm$fn, n_obs = lm$nobs)
  # map_fn forces the per-row path even under the serial strategy.
  perrow <- pesto_forward_model(fn = lm$fn, n_obs = lm$nobs, map_fn = lapply)
  expect_equal(pesto_evaluate(bulk, theta)[],
               pesto_evaluate(perrow, theta)[])
})

test_that("on_failure='stop' aborts on NA rows and on errors", {
  lm <- linear_model()
  theta <- matrix(stats::rnorm(4L * lm$npar), 4L, lm$npar)
  na_fn <- function(theta) {
    r <- theta %*% t(lm$G); r[1L, ] <- NA; r
  }
  fm_stop <- pesto_forward_model(fn = na_fn, n_obs = lm$nobs,
                                 on_failure = "stop")
  expect_error(pesto_evaluate(fm_stop, theta), "on_failure")

  err_fn <- function(theta) stop("boom")
  fm_err <- pesto_forward_model(fn = err_fn, n_obs = lm$nobs,
                                on_failure = "stop")
  expect_error(pesto_evaluate(fm_err, theta), "failed|boom")
})

test_that("max_fail_frac aborts past the threshold", {
  lm <- linear_model()
  theta <- matrix(stats::rnorm(4L * lm$npar), 4L, lm$npar)
  na_fn <- function(theta) {
    r <- theta %*% t(lm$G); r[1L, ] <- NA; r   # 1/4 = 25% failures
  }
  ok  <- pesto_forward_model(fn = na_fn, n_obs = lm$nobs, max_fail_frac = 0.5)
  bad <- pesto_forward_model(fn = na_fn, n_obs = lm$nobs, max_fail_frac = 0.1)
  expect_identical(attr(pesto_evaluate(ok, theta), "n_failures"), 1L)
  expect_error(pesto_evaluate(bad, theta), "max_fail_frac")
})

test_that("multicore evaluation matches serial and is reproducible", {
  skip_on_os("windows")              # mclapply is serial on Windows
  skip_on_cran()
  if (parallel::detectCores() < 2L) skip("needs >= 2 cores")
  lm <- linear_model(seed = 3L)
  theta <- matrix(stats::rnorm(8L * lm$npar), 8L, lm$npar)
  serial <- pesto_forward_model(fn = lm$fn, n_obs = lm$nobs)
  multi  <- pesto_forward_model(fn = lm$fn, n_obs = lm$nobs,
                                parallel = "multicore", n_cores = 2L)
  old <- RNGkind("L'Ecuyer-CMRG"); on.exit(RNGkind(old[[1L]]), add = TRUE)
  set.seed(42L); a <- pesto_evaluate(multi, theta)
  set.seed(42L); b <- pesto_evaluate(multi, theta)
  expect_equal(a[], b[])                       # stream-reproducible
  expect_equal(a[], pesto_evaluate(serial, theta)[])  # same answer
})
