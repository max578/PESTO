# Posterior-calibration gate for the native R iterative ensemble smoother.
#
# The oracle in this file is the closed-form Gaussian-conjugate posterior of a
# linear model, which is independent of every line of PESTO: for
#
#   d | theta ~ N(G theta, C_D),   theta ~ N(mu_0, C_0),
#
# the posterior is Gaussian with
#
#   C_post = (C_0^-1 + G' C_D^-1 G)^-1,
#   mu_post = C_post (C_0^-1 mu_0 + G' C_D^-1 d),
#
# (Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd edn, ch. 14;
# Rasmussen, C.E. & Williams, C.K.I. (2006). *Gaussian Processes for Machine
# Learning*, eq. 2.8). No expected value below is read off a PESTO run.
#
# What is graded here that is graded nowhere else in the suite: the posterior
# **second moment** and the **realised coverage** of the ensemble's central
# intervals. `test-bayesian-recovery.R` grades the posterior mean, and the
# spread-ESS diagnostic is a participation ratio -- invariant to a global
# rescaling of the anomalies -- so neither can see a uniformly collapsed
# spread. The unperturbed smoother
# (`obs_perturbation = "none"`) returns a spread 5--20 per cent of the
# analytic posterior on these problems; that is the failure this file gates.
#
#' @srrstats {BS4.2} Posterior estimates are validated against a closed-form
#'   oracle in both moments: the mean in test-bayesian-recovery.R /
#'   test-correctness-analytic.R, and the covariance plus the realised
#'   coverage of nominal 50 / 80 / 95 per cent intervals here.
#' @srrstats {BS7.2} Recovery of the expected posterior given a prior and
#'   data, graded on the full covariance matrix rather than the mean alone.
#' @noRd
NULL


# Build one linear-Gaussian inverse problem together with its analytic
# posterior. Everything returned is closed-form; nothing is run through PESTO.
.lin_gauss_problem <- function(seed, npar = 3L, nobs = 12L, sigma = 0.2,
                               prior_sd = 1.5) {
  set.seed(seed)
  G     <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta <- stats::rnorm(npar)
  y     <- as.numeric(G %*% theta) + stats::rnorm(nobs, sd = sigma)
  mu0   <- rep(0, npar)
  s0    <- rep(prior_sd, npar)

  c0_inv <- diag(1 / s0^2, nrow = npar)
  prec   <- c0_inv + crossprod(G) / sigma^2
  c_post <- solve(prec)
  mu_pos <- drop(c_post %*% (c0_inv %*% mu0 + crossprod(G, y) / sigma^2))

  list(
    G = G, y = y, sigma = sigma, npar = npar, nobs = nobs,
    mu0 = mu0, s0 = s0,
    obs = stats::setNames(y, paste0("o", seq_len(nobs))),
    fn  = function(t) t %*% t(G),
    c_post = c_post, mu_post = mu_pos, sd_post = sqrt(diag(c_post))
  )
}

.draw_prior <- function(p, nreal) {
  m <- vapply(seq_len(p$npar),
              function(j) stats::rnorm(nreal, p$mu0[j], p$s0[j]),
              numeric(nreal))
  colnames(m) <- paste0("p", seq_len(p$npar))
  m
}


test_that("the ensemble covariance recovers the analytic posterior covariance", {
  # Tolerance justification. Two error sources only: Monte-Carlo error in the
  # sample covariance of `nreal` draws, whose relative standard error on each
  # entry is O(sqrt(2 / nreal)) = 3.2 per cent at nreal = 2000; and the
  # ensemble-approximation error of the ES-MDA update itself, of the same
  # order. A 15 per cent relative Frobenius tolerance is therefore ~4 sample
  # standard errors -- loose enough not to flake, and roughly six times
  # tighter than the ~0.98 relative Frobenius error the unperturbed update
  # returns on this problem.
  p     <- .lin_gauss_problem(seed = 4011L)
  nreal <- 2000L

  set.seed(99L)
  prior <- .draw_prior(p, nreal)
  fit <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma,
                            noptmax = 4L, parcov = p$s0^2, verbose = FALSE)
  post <- as.matrix(fit$par_ensemble[, -1L])

  rel_frob <- norm(stats::cov(post) - p$c_post, "F") / norm(p$c_post, "F")
  expect_lt(rel_frob, 0.15)

  # Marginal spread, parameter by parameter: within 15 per cent of analytic.
  sd_ratio <- apply(post, 2L, stats::sd) / p$sd_post
  expect_equal(unname(sd_ratio), rep(1, p$npar), tolerance = 0.15)

  # The mean must not regress while the spread is fixed.
  expect_equal(unname(colMeans(post)), unname(p$mu_post),
               tolerance = 0.1 * mean(p$sd_post))
})


test_that("the unperturbed smoother is the one that collapses the spread", {
  # The same problem, same seeds, `obs_perturbation = "none"`: the mean is
  # still right and the spread is an order of magnitude too small. This is
  # the behaviour the default no longer has, pinned so it cannot come back
  # unnoticed.
  p     <- .lin_gauss_problem(seed = 4011L)
  nreal <- 2000L

  set.seed(99L)
  prior <- .draw_prior(p, nreal)
  fit <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma,
                            noptmax = 4L, parcov = p$s0^2,
                            obs_perturbation = "none", verbose = FALSE)
  post <- as.matrix(fit$par_ensemble[, -1L])

  sd_ratio <- apply(post, 2L, stats::sd) / p$sd_post
  expect_true(all(sd_ratio < 0.3))
  expect_equal(unname(colMeans(post)), unname(p$mu_post),
               tolerance = 0.2 * mean(p$sd_post))
})


test_that("nominal 50 / 80 / 95 per cent intervals realise their coverage", {
  # Coverage protocol. For each of `n_rep` independent problems (a fresh G, a
  # fresh data realisation, a fresh prior ensemble) the analytic posterior
  # gives a *known* truth to cover: a draw theta* ~ N(mu_post, C_post), which
  # by construction falls inside the exact central-`lev` interval with
  # probability `lev`. Counting how often the *ensemble's* interval contains
  # theta* over n_rep * npar independent trials estimates the realised
  # coverage of the smoother's intervals.
  #
  # Tolerance justification. With n_rep = 30 and npar = 3 there are 90
  # trials, so the binomial standard error is at most 0.053 (at lev = 0.5)
  # and 0.023 at lev = 0.95. The band below is +/- 0.12, i.e. between 2.3 and
  # 5 standard errors -- wide enough to be seed-robust, and far tighter than
  # the gap to the unperturbed smoother, whose realised coverage on these
  # problems is a few per cent at every nominal level.
  skip_on_cran()
  levels_nom <- c(0.50, 0.80, 0.95)
  n_rep      <- 30L
  nreal      <- 400L
  hits       <- matrix(FALSE, nrow = n_rep * 3L, ncol = length(levels_nom))

  row <- 0L
  for (r in seq_len(n_rep)) {
    p <- .lin_gauss_problem(seed = 7000L + r)
    set.seed(8000L + r)
    prior <- .draw_prior(p, nreal)
    fit <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma,
                              noptmax = 4L, parcov = p$s0^2, verbose = FALSE)
    post <- as.matrix(fit$par_ensemble[, -1L])

    # A truth drawn from the analytic posterior (closed form, not from PESTO).
    ch    <- chol(p$c_post)
    truth <- p$mu_post + drop(crossprod(ch, stats::rnorm(p$npar)))

    for (li in seq_along(levels_nom)) {
      lev <- levels_nom[li]
      qs  <- apply(post, 2L, stats::quantile,
                   probs = c((1 - lev) / 2, 1 - (1 - lev) / 2))
      hits[(row + 1L):(row + p$npar), li] <-
        qs[1L, ] <= truth & truth <= qs[2L, ]
    }
    row <- row + p$npar
  }

  realised <- colMeans(hits)
  for (li in seq_along(levels_nom)) {
    expect_gt(realised[li], levels_nom[li] - 0.12)
    expect_lt(realised[li], levels_nom[li] + 0.12)
  }
})


test_that("the ES-MDA inflation schedule is enforced, not assumed", {
  p <- .lin_gauss_problem(seed = 4012L, nobs = 8L)
  set.seed(3L)
  prior <- .draw_prior(p, 60L)

  # sum(1 / alpha) must be exactly 1: a uniform alpha = 2 over 4 steps
  # assimilates the likelihood twice and is refused.
  expect_error(
    pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 4L,
                       mda_alpha = rep(2, 4L), verbose = FALSE),
    "sum\\(1 / alpha\\)"
  )
  # `lambda` is read as alpha - 1 and validated the same way.
  expect_error(
    pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 4L,
                       lambda = 1.0, verbose = FALSE),
    "sum\\(1 / alpha\\)"
  )
  # A valid schedule is accepted, and the derived lambda schedule recorded.
  alpha <- c(2, 4, 8, 8)           # 1/2 + 1/4 + 1/8 + 1/8 = 1
  fit <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 4L,
                            mda_alpha = alpha, verbose = FALSE)
  expect_equal(fit$obs_perturbation$alpha, alpha)
  expect_equal(vapply(fit$iterations, function(it) it$lambda, numeric(1L)),
               alpha - 1)
  expect_identical(fit$obs_perturbation$scheme, "mda")
})


test_that("mda_alpha is refused under obs_perturbation = 'none'", {
  p <- .lin_gauss_problem(seed = 4013L, nobs = 8L)
  set.seed(4L)
  prior <- .draw_prior(p, 40L)
  expect_error(
    pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 2L,
                       obs_perturbation = "none", mda_alpha = c(2, 2),
                       verbose = FALSE),
    "only meaningful"
  )
})


test_that("the noise ensemble and seed are recorded on the result", {
  p <- .lin_gauss_problem(seed = 4014L, nobs = 8L)
  set.seed(5L)
  prior <- .draw_prior(p, 50L)
  fit <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 3L,
                            seed = 12345L, verbose = FALSE)
  op <- fit$obs_perturbation
  expect_identical(op$seed, 12345L)
  expect_equal(dim(op$obs_noise), c(p$nobs, 50L))
  expect_identical(rownames(op$obs_noise), names(p$obs))
  # The final step drew N(0, alpha_3 C_D): a chi-square check on the realised
  # scale, generous enough for 50 * 8 = 400 draws.
  realised <- stats::sd(as.numeric(op$obs_noise))
  expect_equal(realised, p$sigma * sqrt(op$alpha[3L]), tolerance = 0.15)

  # A seeded run reproduces exactly; the seed reaches the manifest.
  fit2 <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 3L,
                             seed = 12345L, verbose = FALSE)
  expect_equal(as.matrix(fit2$par_ensemble[, -1L]),
               as.matrix(fit$par_ensemble[, -1L]))
  expect_identical(as_manifest(fit)@seed, 12345L)
})


test_that("phi is reported against the unperturbed observations", {
  # Perturbation enters the update, not the reported objective function: a
  # run and its `obs_perturbation = "none"` twin must report the same phi at
  # iteration 1, where both start from the identical prior ensemble.
  p <- .lin_gauss_problem(seed = 4015L, nobs = 10L)
  set.seed(6L)
  prior <- .draw_prior(p, 80L)
  a <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 2L,
                          verbose = FALSE)
  b <- pesto_ies_callback(p$fn, prior, p$obs, obs_sd = p$sigma, noptmax = 2L,
                          obs_perturbation = "none", verbose = FALSE)
  expect_equal(a$phi[iteration == 1L]$phi, b$phi[iteration == 1L]$phi)
})


test_that("the filter also returns a calibrated posterior", {
  # Disjoint windows assimilate each observation exactly once, so the
  # sequential filter must land on the same analytic posterior as the batch
  # smoother. Same oracle, same tolerance rationale, fewer realisations.
  skip_on_cran()
  p <- .lin_gauss_problem(seed = 4016L, nobs = 12L)
  set.seed(11L)
  prior <- .draw_prior(p, 2000L)
  fit <- pesto_ies_filter(p$fn, prior, p$obs, obs_sd = p$sigma,
                          windows = list(1:4, 5:8, 9:12),
                          window_noptmax = 2L, parcov = p$s0^2,
                          verbose = FALSE)
  post <- as.matrix(fit$par_ensemble[, -1L])
  sd_ratio <- apply(post, 2L, stats::sd) / p$sd_post
  expect_equal(unname(sd_ratio), rep(1, p$npar), tolerance = 0.20)
  expect_identical(fit$obs_perturbation$scheme, "mda")
  expect_equal(fit$windows[[1L]]$mda_alpha, c(2, 2))
})


test_that("the filter refuses a lambda it would have to reinterpret", {
  p <- .lin_gauss_problem(seed = 4017L, nobs = 6L)
  set.seed(12L)
  prior <- .draw_prior(p, 40L)
  expect_error(
    pesto_ies_filter(p$fn, prior, p$obs, obs_sd = p$sigma,
                     windows = list(1:3, 4:6), lambda = 1.0, verbose = FALSE),
    "mda_alpha"
  )
})
