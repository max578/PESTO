# Tests for pesto_ies_callback() — the in-process R-side IES driver
# introduced for D4 (apsimx native R callback path).

test_that("linear-Gaussian recovery: posterior mean moves toward truth", {
  set.seed(42)
  npar  <- 4L
  nobs  <- 8L
  nreal <- 120L
  sigma <- 0.05

  G          <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0, 0.25)
  y          <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = sigma)
  names(y)   <- paste0("o", seq_len(nobs))

  forward <- function(theta) theta %*% t(G)

  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  prior_mean <- colMeans(prior)
  prior_err  <- sqrt(mean((prior_mean - theta_true)^2))

  fit <- pesto_ies_callback(
    forward_model  = forward,
    prior_ensemble = prior,
    obs            = y,
    obs_sd         = sigma,
    noptmax        = 6L,
    lambda         = 1.0,
    verbose        = FALSE
  )

  expect_s3_class(fit, "pesto_ies_callback_result")
  expect_true(data.table::is.data.table(fit$par_ensemble))
  expect_true(data.table::is.data.table(fit$obs_ensemble))
  expect_equal(nrow(fit$par_ensemble), nreal)
  expect_equal(ncol(fit$par_ensemble), npar + 1L)        # +real_name
  expect_equal(setdiff(names(fit$par_ensemble), "real_name"),
               paste0("p", seq_len(npar)))

  post_mean <- colMeans(as.matrix(fit$par_ensemble[, -1L]))
  post_err  <- sqrt(mean((post_mean - theta_true)^2))

  # Posterior strictly closer to truth than the prior was
  expect_lt(post_err, prior_err * 0.5)

  # Phi decreases across iterations (mean phi)
  iter_phi <- vapply(fit$iterations, `[[`, numeric(1L), "mean_phi")
  expect_lt(iter_phi[length(iter_phi)], iter_phi[1L])

  # Bookkeeping
  expect_equal(fit$failure_rate, 0)
  expect_equal(fit$n_forward_evals, nreal * (6L + 1L))   # noptmax + final
})

test_that("on_failure='na' tolerates per-realisation NA returns", {
  set.seed(7)
  npar <- 2L; nobs <- 3L; nreal <- 30L
  G <- matrix(c(1, 0, 0, 1, 1, -1), nobs, npar)
  theta_true <- c(0.5, -0.25)
  y <- as.numeric(G %*% theta_true)
  names(y) <- paste0("o", seq_len(nobs))

  # Forward model that fails on a fixed minority of realisations
  forward <- function(theta) {
    out <- theta %*% t(G)
    bad <- seq_len(nrow(theta)) %% 7L == 0L   # ~14% failure
    out[bad, ] <- NA_real_
    out
  }

  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, c("a", "b")))

  fit <- pesto_ies_callback(
    forward_model  = forward,
    prior_ensemble = prior,
    obs            = y,
    obs_sd         = 0.1,
    noptmax        = 3L,
    on_failure     = "na",
    verbose        = FALSE
  )

  expect_gt(fit$failure_rate, 0)
  expect_lt(fit$failure_rate, 0.5)
  expect_equal(nrow(fit$par_ensemble), nreal)
})

test_that("on_failure='stop' aborts when forward_model errors", {
  set.seed(7)
  npar <- 2L; nobs <- 2L; nreal <- 10L
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, c("a", "b")))

  bad_forward <- function(theta) stop("simulated failure")

  expect_error(
    pesto_ies_callback(
      forward_model  = bad_forward,
      prior_ensemble = prior,
      obs            = c(o1 = 0.0, o2 = 0.0),
      obs_sd         = 0.1,
      noptmax        = 2L,
      on_failure     = "stop",
      verbose        = FALSE
    ),
    regexp = "forward_model failed"
  )
})

test_that("input validation rejects malformed arguments", {
  prior <- matrix(stats::rnorm(20), 10L, 2L,
                  dimnames = list(NULL, c("a", "b")))
  obs <- c(o1 = 1.0, o2 = 1.0)
  f <- function(theta) matrix(0, nrow(theta), 2L)

  expect_error(
    pesto_ies_callback("not_a_function", prior, obs, 0.1, verbose = FALSE),
    "forward_model"
  )
  expect_error(
    pesto_ies_callback(f, prior, obs, obs_sd = -1, verbose = FALSE),
    "obs_sd"
  )
  expect_error(
    pesto_ies_callback(f, prior, obs, 0.1, noptmax = 0L, verbose = FALSE),
    "noptmax"
  )
  expect_error(
    pesto_ies_callback(f, prior[1L, , drop = FALSE], obs, 0.1,
                       verbose = FALSE),
    "at least 2 realisations"
  )
  expect_error(
    pesto_ies_callback(f, prior, obs, 0.1,
                       parcov = c(-1, 1), verbose = FALSE),
    "parcov"
  )
})

test_that("data.table prior ensemble is accepted", {
  set.seed(11)
  prior_dt <- data.table::data.table(
    real_name = paste0("r", 1:20),
    a = stats::rnorm(20),
    b = stats::rnorm(20)
  )
  f <- function(theta) {
    matrix(c(theta[, "a"] + theta[, "b"],
             theta[, "a"] - theta[, "b"]),
           nrow = nrow(theta), ncol = 2L)
  }
  fit <- pesto_ies_callback(
    forward_model  = f,
    prior_ensemble = prior_dt,
    obs            = c(o1 = 0.3, o2 = -0.1),
    obs_sd         = 0.05,
    noptmax        = 3L,
    verbose        = FALSE
  )
  expect_equal(setdiff(names(fit$par_ensemble), "real_name"), c("a", "b"))
})
