# Tests for the finite-ensemble pathology countermeasures: spread-ESS
# diagnostic, covariance inflation (RTPS / adaptive / multiplicative), and
# covariance localisation (correlation-based + Gaspari-Cohn), plus the
# explicit-gain localised kernel and the additive (NULL-default) wiring.

# -- ensemble_spread_ess ------------------------------------------------------

test_that("spread-ESS separates healthy from collapsed ensembles", {
  set.seed(1L)
  good <- matrix(rnorm(6L * 80L), 6L, 80L)
  d_good <- ensemble_spread_ess(good)
  expect_equal(d_good$r_max, 6)
  expect_gt(d_good$ess_ratio, 0.6)
  expect_lte(d_good$ess_ratio, 1)

  v <- rnorm(6L)
  bad <- outer(v, rnorm(80L)) + matrix(rnorm(6L * 80L, sd = 1e-3), 6L, 80L)
  d_bad <- ensemble_spread_ess(bad)
  expect_lt(d_bad$ess_ratio, 0.3)
  expect_lt(d_bad$ess, d_good$ess)
})

test_that("spread-ESS handles degenerate and too-small inputs", {
  zero <- matrix(0, 4L, 5L)
  expect_equal(ensemble_spread_ess(zero)$ess, 1)
  expect_error(ensemble_spread_ess(matrix(0, 4L, 1L)),
               "at least 2 columns")
})

# -- gaspari_cohn -------------------------------------------------------------

test_that("Gaspari-Cohn taper has the right shape and support", {
  d <- matrix(c(0, 0.5, 1, 1.5, 2, 2.5), nrow = 1L)
  g <- as.numeric(gaspari_cohn(d, radius = 1))
  expect_equal(g[1], 1)                 # G(0) = 1
  expect_equal(g[5], 0)                 # G(2c) = 0
  expect_equal(g[6], 0)                 # beyond support
  expect_true(all(diff(g) <= 1e-9))     # monotone non-increasing
  expect_true(all(g >= 0 & g <= 1))
})

test_that("Gaspari-Cohn validates its arguments", {
  expect_error(gaspari_cohn(matrix(1, 1L, 1L), radius = 0), "positive")
  expect_error(gaspari_cohn(matrix(-1, 1L, 1L), radius = 1), "non-negative")
})

# -- correlation_localisation -------------------------------------------------

test_that("correlation localisation keeps planted links, kills spurious", {
  set.seed(2L)
  npar <- 8L; nobs <- 5L; nreal <- 50L
  pd <- matrix(rnorm(npar * nreal), npar, nreal)
  od <- matrix(rnorm(nobs * nreal), nobs, nreal)
  od[1L, ] <- od[1L, ] + 2 * pd[1L, ]     # genuine par1 <-> obs1 link

  loc <- correlation_localisation(pd, od, taper = "hard")
  expect_equal(loc$rho[1L, 1L], 1)        # planted link retained
  expect_true(loc$frac_active < 0.5)      # most spurious entries damped
  expect_true(all(loc$rho %in% c(0, 1)))  # hard taper is an indicator
})

test_that("soft taper is a graded ramp; manual threshold is honoured", {
  set.seed(3L)
  pd <- matrix(rnorm(6L * 40L), 6L, 40L)
  od <- matrix(rnorm(4L * 40L), 4L, 40L)
  soft <- correlation_localisation(pd, od, threshold = 0.2, taper = "soft")
  expect_true(all(soft$rho >= 0 & soft$rho <= 1))
  expect_equal(soft$threshold, 0.2)
  expect_true(any(soft$rho > 0 & soft$rho < 1))  # genuinely graded
})

test_that("automatic floor is reproducible under set.seed", {
  pd <- matrix(stats::rnorm(6L * 40L), 6L, 40L)
  od <- matrix(stats::rnorm(4L * 40L), 4L, 40L)
  set.seed(99L); a <- correlation_localisation(pd, od, n_shuffle = 3L)
  set.seed(99L); b <- correlation_localisation(pd, od, n_shuffle = 3L)
  expect_identical(a$threshold, b$threshold)
  expect_identical(a$rho, b$rho)
})

test_that("correlation localisation validates shape", {
  expect_error(
    correlation_localisation(matrix(1, 3L, 5L), matrix(1, 2L, 4L)),
    "same number of"
  )
  expect_error(
    correlation_localisation(matrix(1, 3L, 2L), matrix(1, 2L, 2L)),
    "at least 3 realisations"
  )
})

# -- ensemble_solution_localised ----------------------------------------------

test_that("localised kernel with rho = 1 matches the standard GLM update", {
  set.seed(4L)
  npar <- 5L; nreal <- 30L; nobs <- 8L
  par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
  obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
  obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
  par_resid <- matrix(rnorm(npar * nreal, sd = 0.1), npar, nreal)
  weights   <- rep(1.3, nobs)
  parcov_inv <- rep(1, npar)
  Am <- matrix(0, 0, 0)

  u_std <- ensemble_solution(par_diff, obs_diff, obs_resid, par_resid,
                             weights, parcov_inv, Am, cur_lam = 2,
                             eigthresh = 1e-6, use_approx = TRUE)
  u_loc <- ensemble_solution_localised(par_diff, obs_diff, obs_resid,
                                       weights, matrix(1, npar, nobs),
                                       cur_lam = 2, eigthresh = 1e-6)
  expect_equal(u_loc, u_std, tolerance = 1e-10)
})

test_that("localised kernel with rho = 0 yields no update", {
  set.seed(5L)
  npar <- 4L; nreal <- 20L; nobs <- 6L
  u <- ensemble_solution_localised(
    matrix(rnorm(npar * nreal), npar, nreal),
    matrix(rnorm(nobs * nreal), nobs, nreal),
    matrix(rnorm(nobs * nreal), nobs, nreal),
    rep(1, nobs), matrix(0, npar, nobs), cur_lam = 1
  )
  expect_equal(dim(u), c(nreal, npar))
  expect_true(all(abs(u) < 1e-12))
  expect_error(
    ensemble_solution_localised(
      matrix(0, npar, nreal), matrix(0, nobs, nreal),
      matrix(0, nobs, nreal), rep(1, nobs), matrix(1, npar, nobs + 1L),
      cur_lam = 1
    ),
    "npar x nobs"
  )
})

# -- pesto_inflation / pesto_localisation specs -------------------------------

test_that("pesto_inflation validates and prints", {
  expect_s3_class(pesto_inflation(), "pesto_inflation")
  expect_equal(pesto_inflation()$method, "none")
  expect_error(pesto_inflation("rtps", alpha = 1.5), "alpha")
  expect_error(pesto_inflation("multiplicative", factor = 0.5), "factor")
  expect_error(pesto_inflation("adaptive", retention_floor = 0), "retention")
  # Resolve the method directly: dispatch via print() is unreliable under
  # devtools::load_all() for this S7-importing package, but registers cleanly
  # on install (exercised by R CMD check).
  pr <- getS3method("print", "pesto_inflation")
  expect_match(paste(utils::capture.output(pr(pesto_inflation("rtps"))),
                     collapse = "\n"), "method: rtps")
})

test_that("pesto_localisation validates and prints", {
  expect_s3_class(pesto_localisation(), "pesto_localisation")
  expect_error(pesto_localisation("distance"), "radius")
  expect_error(
    pesto_localisation("distance", radius = 1),
    "either `distances` or both"
  )
  loc <- pesto_localisation("distance", radius = 2,
                            distances = matrix(1, 3L, 2L))
  expect_equal(loc$radius, 2)
  pr <- getS3method("print", "pesto_localisation")
  expect_match(paste(utils::capture.output(pr(pesto_localisation("correlation"))),
                     collapse = "\n"), "automatic floor")
})

# -- .apply_inflation behaviour -----------------------------------------------

test_that("multiplicative inflation scales anomalies exactly, keeps the mean", {
  set.seed(6L)
  bg <- matrix(rnorm(40L * 3L), 40L, 3L)
  post <- bg * 0.5                                   # collapsed posterior
  out <- .apply_inflation(post, bg, pesto_inflation("multiplicative",
                                                    factor = 2))
  expect_equal(colMeans(out$par), colMeans(post))    # mean preserved
  expect_equal(.col_sd(out$par), .col_sd(post) * 2)  # spread doubled
  expect_equal(out$factor, 2)
})

test_that("RTPS re-expands collapsed spread toward the background", {
  set.seed(7L)
  bg <- matrix(rnorm(60L * 4L), 60L, 4L)
  post <- sweep(sweep(bg, 2L, colMeans(bg), "-") * 0.3, 2L, colMeans(bg), "+")
  ret_before <- mean(.col_sd(post) / .col_sd(bg))
  out <- .apply_inflation(post, bg, pesto_inflation("rtps", alpha = 0.8))
  ret_after <- mean(.col_sd(out$par) / .col_sd(bg))
  expect_gt(ret_after, ret_before)
  expect_equal(colMeans(out$par), colMeans(post))    # mean untouched
})

# -- Integration: additive default + inflation effect -------------------------

test_that("NULL inflation/localisation leaves the callback driver unchanged", {
  set.seed(8L)
  npar <- 3L; nobs <- 6L; nreal <- 60L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  y <- as.numeric(G %*% c(1, -0.5, 2)) + rnorm(nobs, sd = 0.05)
  f <- function(theta) theta %*% t(G)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  args <- list(forward_model = f, prior_ensemble = prior,
               obs = stats::setNames(y, paste0("o", seq_len(nobs))),
               obs_sd = 0.05, noptmax = 5L, verbose = FALSE)
  a <- do.call(pesto_ies_callback, args)
  b <- do.call(pesto_ies_callback,
               c(args, list(inflation = NULL, localisation = NULL)))
  expect_equal(as.matrix(a$par_ensemble[, -1]),
               as.matrix(b$par_ensemble[, -1]))
  # Diagnostics are recorded even with both countermeasures off.
  expect_true(is.finite(a$iterations[[1L]]$spread_ess_ratio))
  expect_equal(a$iterations[[1L]]$inflation_method, "none")
})

test_that("inflation increases retained posterior spread (under-dispersion)", {
  set.seed(42L)
  npar <- 6L; nobs <- 10L; nreal <- 24L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  y <- as.numeric(G %*% rnorm(npar)) + rnorm(nobs, sd = 0.3)
  f <- function(theta) theta %*% t(G)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  spread <- function(inf) {
    fit <- pesto_ies_callback(f, prior,
                              stats::setNames(y, paste0("o", seq_len(nobs))),
                              obs_sd = 0.3, noptmax = 12L,
                              inflation = inf, verbose = FALSE)
    mean(apply(as.matrix(fit$par_ensemble[, -1]), 2L, stats::sd))
  }
  s_none <- spread(NULL)
  s_rtps <- spread(pesto_inflation("rtps", alpha = 0.6))
  expect_gt(s_rtps, s_none * 1.3)        # materially more retained spread
})

test_that("localisation runs end to end and records diagnostics", {
  set.seed(9L)
  npar <- 5L; nobs <- 8L; nreal <- 40L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  y <- as.numeric(G %*% rnorm(npar)) + rnorm(nobs, sd = 0.1)
  f <- function(theta) theta %*% t(G)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_callback(f, prior,
                            stats::setNames(y, paste0("o", seq_len(nobs))),
                            obs_sd = 0.1, noptmax = 4L,
                            localisation = pesto_localisation("correlation"),
                            verbose = FALSE)
  d1 <- fit$iterations[[1L]]
  expect_equal(d1$localisation, "correlation")
  expect_true(is.finite(d1$loc_frac_active))
  expect_true(d1$loc_frac_active >= 0 && d1$loc_frac_active <= 1)
})

test_that("filter driver accepts inflation + localisation", {
  set.seed(10L)
  npar <- 3L; nobs <- 9L; nreal <- 50L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  y <- as.numeric(G %*% c(1, -0.5, 2)) + rnorm(nobs, sd = 0.05)
  f <- function(theta) theta %*% t(G)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_filter(
    f, prior, stats::setNames(y, paste0("o", seq_len(nobs))),
    obs_sd = 0.05, windows = list(1:3, 4:6, 7:9),
    inflation = pesto_inflation("rtps"),
    localisation = pesto_localisation("correlation"),
    verbose = FALSE
  )
  expect_s3_class(fit, "pesto_ies_filter_result")
  expect_equal(fit$windows[[1L]]$inflation_method, "rtps")
  expect_true(is.finite(fit$windows[[1L]]$spread_ess_ratio))
})

# -- .euclidean_distances -----------------------------------------------------

test_that("euclidean distance helper matches stats::dist", {
  pc <- matrix(c(0, 0, 1, 1), 2L, 2L, byrow = TRUE)
  oc <- matrix(c(0, 0, 0, 3), 2L, 2L, byrow = TRUE)
  d <- .euclidean_distances(pc, oc)
  expect_equal(d[1L, 1L], 0)
  expect_equal(d[1L, 2L], 3)
  expect_equal(d[2L, 1L], sqrt(2))
})
