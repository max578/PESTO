# Typed abstention (P-08): the reliability gates that used to `stop()` or
# silently warn now return a classed `pesto_abstention`, visible to the
# orchestra's `is_orchestra_decline()` by the `_abstention` class-name
# convention (ORCHESTRA_dev/integration/refusal_contract.R).

test_that("pesto_abstention() constructs a classed, convention-conforming decline", {
  a <- pesto_abstention("degenerate_ensemble",
                        detail = "1 of 40 realisations succeeded",
                        diagnostics = list(iteration = 2L, n_ok = 1L))
  expect_true(is_pesto_abstention(a))
  expect_true(grepl("_abstention$", class(a)[1]))
  expect_true(a$abstained)
  expect_equal(a$reason, "degenerate_ensemble")
  expect_false(is_pesto_abstention(42))
  expect_output(print(a), "pesto_abstention")
})

test_that("pesto_ies_callback() abstains, typed, on a degenerate ensemble", {
  set.seed(1)
  npar <- 2L; nobs <- 3L; nreal <- 20L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, c("a", "b")))
  # Every realisation but one fails from the first evaluation onward --
  # fewer than 2 successful realisations at iteration 1.
  bad_forward <- function(theta) {
    out <- theta %*% t(G)
    out[-1L, ] <- NA_real_
    out
  }
  fit <- pesto_ies_callback(
    forward_model  = bad_forward,
    prior_ensemble = prior,
    obs            = stats::setNames(rep(0, nobs), paste0("o", seq_len(nobs))),
    obs_sd         = 0.1,
    noptmax        = 2L,
    on_failure     = "na",
    verbose        = FALSE
  )
  expect_true(is_pesto_abstention(fit))
  expect_equal(fit$reason, "degenerate_ensemble")
  expect_equal(fit$diagnostics$iteration, 1L)
})

test_that("pesto_ies_callback() still supports a highly-parameterised ensemble (npar > nreal)", {
  # No structural over-determination guard is added on this axis: see the
  # NOTE in R/pesto_run.R and PESTO_refusals.md. This mirrors
  # test-edge-conditions.R::G5.8d and asserts the capability is untouched.
  set.seed(2)
  npar <- 10L; nreal <- 3L; nobs <- 8L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  fit <- pesto_ies_callback(
    forward_model  = function(t) t %*% t(G),
    prior_ensemble = prior,
    obs            = stats::setNames(rep(0, nobs), paste0("o", seq_len(nobs))),
    obs_sd         = 0.1,
    noptmax        = 2L,
    verbose        = FALSE
  )
  expect_false(is_pesto_abstention(fit))
  expect_equal(dim(as.matrix(fit$par_ensemble[, -1L])), c(nreal, npar))
})

test_that("pesto_ies_filter() abstains, typed, on a degenerate ensemble in a window", {
  set.seed(3)
  npar <- 2L; nobs <- 4L; nreal <- 20L
  G <- matrix(rnorm(nobs * npar), nobs, npar)
  prior <- matrix(rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, c("a", "b")))
  bad_forward <- function(theta) {
    out <- theta %*% t(G)
    out[-1L, ] <- NA_real_
    out
  }
  fit <- pesto_ies_filter(
    forward_model  = bad_forward,
    prior_ensemble = prior,
    obs            = stats::setNames(rep(0, nobs), paste0("o", seq_len(nobs))),
    obs_sd         = 0.1,
    windows        = list(1:2, 3:4),
    on_failure     = "na",
    verbose        = FALSE
  )
  expect_true(is_pesto_abstention(fit))
  expect_equal(fit$reason, "degenerate_ensemble")
})

test_that("pesto_surrogate_ies() abstains, typed, off its favourable training regime", {
  set.seed(4)
  np <- 30L; no <- 12L; nr <- 30L   # ratio nr/np = 1, far below the floor of 5
  pe <- matrix(rnorm(nr * np), nr, np)
  oe <- matrix(rnorm(nr * no), nr, no)
  yt <- rnorm(no)

  fit <- suppressWarnings(pesto_surrogate_ies(
    par_ensemble = pe, obs_ensemble = oe, obs_target = yt,
    weights = rep(1, no), parcov_inv = rep(1, np)
  ))
  expect_true(is_pesto_abstention(fit))
  expect_equal(fit$reason, "surrogate_off_design")
})

test_that("check_surrogate_regime() keeps its own warn-and-return-logical contract", {
  # P-08 gates pesto_surrogate_ies(), not check_surrogate_regime() itself --
  # its existing contract (test-export-surface.R) must be unaffected.
  expect_true(check_surrogate_regime(n_params = 4, n_train = 100))
  expect_warning(check_surrogate_regime(n_params = 30, n_train = 30))
})

test_that("as_manifest() on a pesto_abstention emits a manifest with summary$abstained", {
  a <- pesto_abstention("over_determined", detail = "nreal <= npar",
                        diagnostics = list(nreal = 5L, npar = 5L))
  m <- as_manifest(a, method = "ies_callback")
  expect_s7_class(m, pesto_ensemble_manifest)
  expect_true(m@summary$abstained)
  expect_equal(m@summary$reason, "over_determined")
  expect_equal(nrow(m@params), 0L)
  expect_equal(nrow(m@outputs), 0L)

  # Round-trips through write/read.
  f <- tempfile(fileext = ".yaml")
  write_manifest(m, f)
  m2 <- read_manifest(f)
  expect_true(m2@summary$abstained)
  expect_equal(m2@summary$reason, "over_determined")
  expect_true(verify_manifest(m2)$ok)
})

test_that("as_manifest() on a normal callback result leaves summary NULL", {
  set.seed(5)
  G <- matrix(1, 2L, 2L)
  fit <- pesto_ies_callback(
    function(t) t %*% t(G),
    matrix(rnorm(20L), 10L, 2L, dimnames = list(NULL, c("a", "b"))),
    c(o1 = 0.1, o2 = -0.2), obs_sd = 0.1, noptmax = 2L, verbose = FALSE
  )
  expect_false(is_pesto_abstention(fit))
  m <- as_manifest(fit, seed = 1L)
  expect_null(m@summary)
})

test_that("is_orchestra_decline()-style regex recognises pesto_abstention", {
  # Mirrors ORCHESTRA_dev/integration/refusal_contract.R's
  # is_orchestra_decline() predicate without depending on the workspace file
  # (a coordination artefact, not installable package code).
  is_orchestra_decline_like <- function(x) {
    cls <- class(x)
    any(grepl("_(refusal|abstention)$", cls))
  }
  expect_true(is_orchestra_decline_like(pesto_abstention("degenerate_ensemble")))
})
