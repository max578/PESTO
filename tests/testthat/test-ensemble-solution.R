test_that("ensemble_solution GLM form produces correct dimensions", {
  set.seed(42)
  npar <- 10
  nobs <- 20
  nreal <- 30

  par_diff <- matrix(rnorm(npar * nreal), nrow = npar, ncol = nreal)
  obs_diff <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  obs_resid <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  par_resid <- matrix(rnorm(npar * nreal), nrow = npar, ncol = nreal)
  weights <- abs(rnorm(nobs)) + 0.1
  parcov_inv <- abs(rnorm(npar)) + 0.1
  Am <- matrix(rnorm(npar * (nreal - 1)), nrow = npar, ncol = nreal - 1)

  result <- ensemble_solution(
    par_diff = par_diff,
    obs_diff = obs_diff,
    obs_resid = obs_resid,
    par_resid = par_resid,
    weights = weights,
    parcov_inv = parcov_inv,
    Am = Am,
    cur_lam = 1.0,
    eigthresh = 1e-6,
    use_approx = TRUE,
    use_prior_scaling = FALSE,
    iter = 1
  )

  # Result should be nreal x npar

  expect_equal(nrow(result), nreal)
  expect_equal(ncol(result), npar)
  expect_false(any(is.na(result)))
  expect_false(any(is.infinite(result)))
})

test_that("ensemble_solution MDA form produces correct dimensions", {
  set.seed(123)
  npar <- 5
  nobs <- 10
  nreal <- 20

  par_diff <- matrix(rnorm(npar * nreal), nrow = npar, ncol = nreal)
  obs_diff <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  obs_resid <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  obs_err <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)

  result <- ensemble_solution_mda(
    par_diff = par_diff,
    obs_diff = obs_diff,
    obs_resid = obs_resid,
    obs_err = obs_err,
    cur_lam = 1.0,
    eigthresh = 1e-6
  )

  expect_equal(nrow(result), nreal)
  expect_equal(ncol(result), npar)
  expect_false(any(is.na(result)))
})

test_that("ensemble_solution with upgrade_2 (non-approx) works", {
  set.seed(99)
  npar <- 8
  nobs <- 15
  nreal <- 25

  par_diff <- matrix(rnorm(npar * nreal), nrow = npar, ncol = nreal)
  obs_diff <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  obs_resid <- matrix(rnorm(nobs * nreal), nrow = nobs, ncol = nreal)
  par_resid <- matrix(rnorm(npar * nreal), nrow = npar, ncol = nreal)
  weights <- abs(rnorm(nobs)) + 0.1
  parcov_inv <- abs(rnorm(npar)) + 0.1
  Am <- matrix(rnorm(npar * (nreal - 1)), nrow = npar, ncol = nreal - 1)

  result <- ensemble_solution(
    par_diff = par_diff,
    obs_diff = obs_diff,
    obs_resid = obs_resid,
    par_resid = par_resid,
    weights = weights,
    parcov_inv = parcov_inv,
    Am = Am,
    cur_lam = 0.5,
    eigthresh = 1e-6,
    use_approx = FALSE,
    use_prior_scaling = TRUE,
    iter = 3,
    reg_factor = 0.5
  )

  expect_equal(nrow(result), nreal)
  expect_equal(ncol(result), npar)
})

test_that("compute_phi produces correct values", {
  # Simple test: identity weights, known residuals
  residuals <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  weights <- c(1, 1, 1)

  phi <- compute_phi(residuals, weights)

  expect_equal(length(phi), 2)
  expect_equal(phi[1], sum(c(1, 2, 3)^2))
  expect_equal(phi[2], sum(c(4, 5, 6)^2))
})

test_that("compute_phi respects weights", {
  residuals <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)
  weights <- c(2, 3)

  phi <- compute_phi(residuals, weights)

  # phi = sum(w^2 * r^2) = 4*1 + 9*1 = 13

expect_equal(phi[1], 13)
  expect_equal(phi[2], 13)
})
