# Tests for the multi-fidelity bridge: container construction +
# validation, level evaluation, the affine control-variate primitive,
# and the fidelity-schedule resolver.

two_level_mf <- function() {
  set.seed(11L)
  G <- matrix(stats::rnorm(12L), 4L, 3L)
  truth <- function(theta) theta %*% t(G)
  biased <- function(theta) theta %*% t(G) + 0.5
  pesto_multifidelity_model(
    levels = list(
      pesto_forward_model(fn = biased, n_obs = 4L, fidelity = 0L),
      pesto_forward_model(fn = truth,  n_obs = 4L, fidelity = 1L)
    ),
    costs = c(1, 30)
  )
}

test_that("constructor coerces functions and defaults costs", {
  mf <- pesto_multifidelity_model(
    levels = list(function(theta) theta, function(theta) theta * 2)
  )
  expect_true(S7::S7_inherits(mf, pesto_multifidelity_model))
  expect_true(S7::S7_inherits(mf@levels[[1L]], pesto_forward_model))
  expect_equal(mf@costs, c(1, 2))
})

test_that("validator rejects malformed level stacks", {
  expect_error(pesto_multifidelity_model(levels = list()), "at least one")
  expect_error(
    pesto_multifidelity_model(
      levels = list(pesto_forward_model(fn = identity)),
      costs = c(1, 2)
    ),
    "one entry per"
  )
  expect_error(
    pesto_multifidelity_model(
      levels = list(pesto_forward_model(fn = identity)),
      costs = -1
    ),
    "positive"
  )
})

test_that("pesto_evaluate selects the requested fidelity level", {
  mf <- two_level_mf()
  theta <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2L, byrow = TRUE)
  cheap <- pesto_evaluate(mf, theta, level = 0L)
  exp_t <- pesto_evaluate(mf, theta, level = 1L)
  expect_equal(cheap - exp_t, matrix(0.5, nrow = 2L, ncol = 4L),
               ignore_attr = TRUE)
  # Default level is the highest fidelity.
  expect_equal(pesto_evaluate(mf, theta)[], exp_t[], ignore_attr = TRUE)
  expect_error(pesto_evaluate(mf, theta, level = 2L), "level")
  expect_error(pesto_evaluate(mf, theta, level = -1L), "level")
})

test_that("mf_control_variate recovers a noiseless affine map", {
  set.seed(5L)
  low_all <- matrix(stats::rnorm(60L), ncol = 3L)
  sub <- 1:8
  low_sub  <- low_all[sub, , drop = FALSE]
  a_true <- c(0.2, -1.0, 3.0); b_true <- c(1.5, 0.8, -0.4)
  high_sub <- sweep(sweep(low_sub, 2L, b_true, "*"), 2L, a_true, "+")
  cv <- mf_control_variate(low_all, high_sub, low_sub)
  expect_equal(attr(cv, "slope"), b_true, tolerance = 1e-8)
  expect_equal(attr(cv, "intercept"), a_true, tolerance = 1e-8)
  # Corrected full-ensemble prediction equals the true affine lift.
  expect_equal(cv[], sweep(sweep(low_all, 2L, b_true, "*"),
                           2L, a_true, "+"), ignore_attr = TRUE)
})

test_that("mf_control_variate degrades gracefully on zero-variance low", {
  low_all <- matrix(c(rep(2, 5L), stats::rnorm(5L)), ncol = 2L)
  sub <- 1:3
  low_sub  <- low_all[sub, , drop = FALSE]   # column 1 is constant
  high_sub <- matrix(c(7, 8, 9, 1, 2, 3), ncol = 2L)
  cv <- mf_control_variate(low_all, high_sub, low_sub)
  expect_identical(attr(cv, "slope")[1L], 0)
  expect_equal(attr(cv, "intercept")[1L], mean(high_sub[, 1L]))
  expect_true(is.na(attr(cv, "subset_cor")[1L]))
})

test_that("mf_control_variate enforces conformable shapes", {
  expect_error(
    mf_control_variate(matrix(1, 2L, 2L), matrix(1, 1L, 3L),
                       matrix(1, 1L, 3L)),
    "same number of columns"
  )
  expect_error(
    mf_control_variate(matrix(1, 2L, 2L), matrix(1, 2L, 2L),
                       matrix(1, 1L, 2L)),
    "same number of rows"
  )
})

test_that("fidelity schedule resolves, pads, and validates", {
  expect_equal(.resolve_fidelity_schedule(NULL, 2L, 4L), rep(1L, 4L))
  expect_equal(.resolve_fidelity_schedule(c(0L, 1L), 2L, 4L),
               c(0L, 1L, 1L, 1L))
  expect_error(.resolve_fidelity_schedule(c(0L, 2L), 2L, 4L), "in \\[0, 1\\]")
})

test_that("driver honours the fidelity schedule and recovers truth", {
  mf <- two_level_mf()
  set.seed(21L)
  prior <- matrix(stats::rnorm(40L * 3L), 40L, 3L,
                  dimnames = list(NULL, paste0("p", 1:3)))
  # Synthesise observations from the high-fidelity level at a known truth.
  y <- as.numeric(pesto_evaluate(mf, matrix(c(1, -0.5, 2), nrow = 1L),
                                 level = 1L))
  fit <- pesto_ies_callback(
    mf, prior, stats::setNames(y, paste0("o", 1:4)), obs_sd = 0.1,
    noptmax = 5L, fidelity_schedule = c(0L, 0L, 1L, 1L, 1L), verbose = FALSE
  )
  est <- colMeans(as.matrix(fit$par_ensemble[, -1]))
  expect_equal(unname(est), c(1, -0.5, 2), tolerance = 0.25)
})
