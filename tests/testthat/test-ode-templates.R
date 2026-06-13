# Tests for the ODE / compartmental forward-model templates: the generic
# ode_forward_model() builder, the crop-growth and SEIR specialisations,
# the RK4 / deSolve solver agreement, input validation, the NA-failure
# contract, and -- the headline -- forward-simulate-then-invert recovery
# of known parameters through pesto_ies_callback().

# Generic builder -------------------------------------------------------

test_that("ode_forward_model integrates a scalar linear ODE accurately", {
  # dy/dt = -k y has the closed form y(t) = y0 exp(-k t).
  fm <- ode_forward_model(
    derivs      = function(t, y, theta) -theta[["k"]] * y,
    y0          = c(y = 1),
    times       = seq(0, 4, by = 1),
    param_names = "k"
  )
  expect_true(S7::S7_inherits(fm, pesto_forward_model))
  theta <- matrix(0.5, nrow = 1L, dimnames = list(NULL, "k"))
  out <- as.numeric(pesto_evaluate(fm, theta))
  expect_equal(out, exp(-0.5 * (1:4)), tolerance = 1e-4)
})

test_that("a theta-dependent initial state is honoured", {
  # y0 itself is a calibrated parameter; observe the state at t = 1.
  fm <- ode_forward_model(
    derivs      = function(t, y, theta) -0.3 * y,
    y0          = function(theta) c(y = theta[["y0"]]),
    times       = c(0, 1),
    param_names = "y0"
  )
  out <- as.numeric(pesto_evaluate(fm, matrix(c(2, 5), ncol = 1L,
                                             dimnames = list(NULL, "y0"))))
  expect_equal(out, c(2, 5) * exp(-0.3), tolerance = 1e-4)
})

test_that("the rk4 and desolve solvers agree on a smooth system", {
  skip_if_not_installed("deSolve")
  times <- seq(0, 120, by = 15)
  truth <- matrix(c(0.06, 1400, 20), nrow = 1L,
                  dimnames = list(NULL, c("r", "b_max", "b0")))
  rk4 <- crop_growth_forward_model(times = times, solver = "rk4")
  ds  <- crop_growth_forward_model(times = times, solver = "desolve")
  expect_equal(as.numeric(pesto_evaluate(rk4, truth)),
               as.numeric(pesto_evaluate(ds, truth)),
               tolerance = 1e-2)
})

# Validation ------------------------------------------------------------

test_that("ode_forward_model rejects malformed inputs", {
  ok_derivs <- function(t, y, theta) -y
  expect_error(ode_forward_model(ok_derivs, c(y = 1), times = 0),
               "strictly increasing")
  expect_error(ode_forward_model(ok_derivs, c(y = 1), times = c(0, 0, 1)),
               "strictly increasing")
  expect_error(ode_forward_model("not a function", c(y = 1), times = c(0, 1)),
               "derivs")
  expect_error(ode_forward_model(ok_derivs, "bad", times = c(0, 1)),
               "y0")
  expect_error(
    ode_forward_model(ok_derivs, c(y = 1), times = c(0, 1), n_steps = 0L),
    "n_steps"
  )
})

test_that("a failing realisation becomes an NA row, not an error", {
  # derivs errors only for the second realisation (k < 0 sentinel).
  fm <- ode_forward_model(
    derivs = function(t, y, theta) {
      if (theta[["k"]] < 0) stop("bad parameter")
      -theta[["k"]] * y
    },
    y0          = c(y = 1),
    times       = c(0, 1, 2),
    param_names = "k",
    on_failure  = "na"
  )
  theta <- matrix(c(0.5, -1), ncol = 1L, dimnames = list(NULL, "k"))
  out <- pesto_evaluate(fm, theta)
  expect_equal(dim(out), c(2L, 2L))
  expect_identical(attr(out, "n_failures"), 1L)
  expect_true(all(is.na(out[2L, ])))
  expect_true(all(is.finite(out[1L, ])))
})

# Crop-growth template --------------------------------------------------

test_that("crop_growth_forward_model produces a sigmoid that recovers", {
  set.seed(101L)
  times <- seq(0, 120, by = 15)
  fm    <- crop_growth_forward_model(times = times)
  expect_identical(fm@n_obs, length(times) - 1L)

  truth_v <- c(r = 0.06, b_max = 1400, b0 = 20)
  truth   <- matrix(truth_v, nrow = 1L, dimnames = list(NULL, names(truth_v)))
  biomass <- as.numeric(pesto_evaluate(fm, truth))
  # A logistic curve: monotone increasing, approaching but below b_max.
  expect_true(all(diff(biomass) > 0))
  expect_true(max(biomass) < truth_v[["b_max"]])

  # Forward-simulate-then-invert: the posterior mean is far closer to the
  # truth than the prior mean (relative parameter error contracts).
  nreal <- 120L
  obs        <- biomass + stats::rnorm(length(biomass), sd = 20)
  names(obs) <- paste0("t", seq_along(obs))
  prior <- cbind(
    r     = stats::runif(nreal, 0.02, 0.12),
    b_max = stats::runif(nreal, 900, 2000),
    b0    = stats::runif(nreal, 5, 60)
  )
  fit  <- pesto_ies_callback(fm, prior, obs, obs_sd = 20, noptmax = 8L,
                             verbose = FALSE)
  post <- colMeans(as.matrix(fit$par_ensemble[, -1L]))[names(truth_v)]

  rel <- function(p) sqrt(mean(((p - truth_v) / truth_v)^2))
  expect_lt(rel(post), 0.10)
  expect_lt(rel(post), rel(colMeans(prior)))
})

# SEIR template ---------------------------------------------------------

test_that("seir_forward_model produces an epidemic curve that recovers", {
  set.seed(202L)
  times <- seq(0, 60, by = 5)
  fm    <- seir_forward_model(times = times, n_pop = 1000, i0 = 1)
  expect_identical(fm@n_obs, length(times) - 1L)

  truth_v <- c(beta = 0.6, sigma = 0.2, gamma = 0.1)
  truth   <- matrix(truth_v, nrow = 1L, dimnames = list(NULL, names(truth_v)))
  prev    <- as.numeric(pesto_evaluate(fm, truth))
  # An epidemic curve: rises to a single peak then declines.
  expect_gt(which.max(prev), 1L)
  expect_lt(which.max(prev), length(prev))

  # Recover the parameters and, in particular, the reproduction number.
  nreal <- 200L
  obs        <- prev + stats::rnorm(length(prev), sd = 3)
  names(obs) <- paste0("d", seq_along(obs))
  prior <- cbind(
    beta  = stats::runif(nreal, 0.3, 0.9),
    sigma = stats::runif(nreal, 0.1, 0.4),
    gamma = stats::runif(nreal, 0.05, 0.2)
  )
  fit  <- pesto_ies_callback(fm, prior, obs, obs_sd = 3, noptmax = 10L,
                             verbose = FALSE)
  post <- colMeans(as.matrix(fit$par_ensemble[, -1L]))[names(truth_v)]

  r0_truth <- truth_v[["beta"]] / truth_v[["gamma"]]
  r0_post  <- post[["beta"]] / post[["gamma"]]
  expect_equal(r0_post, r0_truth, tolerance = 0.15)

  # The posterior-mean curve fits the truth to within the noise level.
  sim_post <- as.numeric(pesto_evaluate(
    fm, matrix(post, nrow = 1L, dimnames = list(NULL, names(post)))
  ))
  expect_lt(sqrt(mean((sim_post - prev)^2)), 5)
})

test_that("seir_forward_model rejects an infeasible seeding", {
  expect_error(seir_forward_model(times = c(0, 1), n_pop = 10, i0 = 10),
               "smaller than")
  expect_error(seir_forward_model(times = c(0, 1), n_pop = -1),
               "n_pop")
})

# Composition with the multi-fidelity stack -----------------------------

test_that("ODE templates compose into a multi-fidelity stack", {
  times <- seq(0, 60, by = 10)
  cheap <- crop_growth_forward_model(times = times, n_steps = 2L,
                                     fidelity = 0L)
  fine  <- crop_growth_forward_model(times = times, n_steps = 20L,
                                     fidelity = 1L)
  mf <- pesto_multifidelity_model(levels = list(cheap, fine),
                                  costs = c(1, 10))
  truth <- matrix(c(0.06, 1400, 20), nrow = 1L,
                  dimnames = list(NULL, c("r", "b_max", "b0")))
  lo <- as.numeric(pesto_evaluate(mf, truth, level = 0L))
  hi <- as.numeric(pesto_evaluate(mf, truth, level = 1L))
  expect_length(lo, length(times) - 1L)
  expect_equal(lo, hi, tolerance = 0.05 * max(hi))
})
