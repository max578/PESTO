# Contract tests: the same forward-model contract must produce the same
# IES posterior however it reaches the driver -- bare function, typed
# pesto_forward_model, or a single-level multi-fidelity container. This
# is the mode-invariance property that lets the native-callback and
# .pst-file adapter modes share one contract.

contract_problem <- function() {
  set.seed(101L)
  npar <- 3L; nobs <- 6L; nreal <- 60L
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- c(1.0, -0.5, 2.0)
  y <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05)
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  list(
    fn    = function(theta) theta %*% t(G),
    obs   = stats::setNames(y, paste0("o", seq_len(nobs))),
    prior = prior, theta_true = theta_true
  )
}

test_that("bare function and typed forward model give identical posteriors", {
  p <- contract_problem()
  args <- list(prior_ensemble = p$prior, obs = p$obs, obs_sd = 0.05,
               noptmax = 4L, verbose = FALSE)

  bare  <- do.call(pesto_ies_callback, c(list(forward_model = p$fn), args))
  typed <- do.call(pesto_ies_callback,
                   c(list(forward_model = pesto_forward_model(
                            fn = p$fn, n_obs = length(p$obs))), args))

  expect_equal(as.matrix(bare$par_ensemble[, -1]),
               as.matrix(typed$par_ensemble[, -1]))
  expect_equal(bare$obs_ensemble, typed$obs_ensemble)
})

test_that("single-level multi-fidelity equals the plain callback path", {
  p <- contract_problem()
  args <- list(prior_ensemble = p$prior, obs = p$obs, obs_sd = 0.05,
               noptmax = 4L, verbose = FALSE)

  bare <- do.call(pesto_ies_callback, c(list(forward_model = p$fn), args))
  mf1  <- pesto_multifidelity_model(
    levels = list(pesto_forward_model(fn = p$fn, n_obs = length(p$obs)))
  )
  via_mf <- do.call(pesto_ies_callback,
                    c(list(forward_model = mf1), args))

  expect_equal(as.matrix(bare$par_ensemble[, -1]),
               as.matrix(via_mf$par_ensemble[, -1]))
})

test_that("the contract recovers the truth within tolerance", {
  p <- contract_problem()
  fit <- pesto_ies_callback(p$fn, p$prior, p$obs, obs_sd = 0.05,
                            noptmax = 6L, verbose = FALSE)
  est <- unname(colMeans(as.matrix(fit$par_ensemble[, -1])))
  expect_equal(est, p$theta_true, tolerance = 0.15)
})

test_that("apsim_callback emits a closure honouring the contract", {
  skip_if_not_installed("apsimx")
  # Stub APSIM so no real binary is needed; assert the closure satisfies
  # the forward-model contract (nreal x nobs, finite, wrappable).
  tmpf <- tempfile(fileext = ".apsimx"); writeLines("stub", tmpf)
  on.exit(unlink(tmpf), add = TRUE)
  workdir <- tempfile("apsim_cb_"); on.exit(unlink(workdir, recursive = TRUE),
                                            add = TRUE)
  store <- new.env(parent = emptyenv()); store$v <- list()
  fm_closure <- apsim_callback(
    template          = tmpf,
    param_map         = list(a = "A", b = "B"),
    output_extractor  = function(sim) {
      utils::tail(vapply(store$v, `[[`, numeric(1L), "value"), 2L)
    },
    workdir           = workdir,
    param_writer      = function(file, src.dir, node, value) {
      store$v <- c(store$v, list(list(value = value))); invisible(TRUE)
    },
    simulation_runner = function(file, src.dir) list()
  )
  # Wrap the apsimx closure in the typed contract and evaluate.
  fm <- pesto_forward_model(fn = fm_closure, n_obs = 2L,
                            param_names = c("a", "b"))
  theta <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2L, byrow = TRUE,
                  dimnames = list(NULL, c("a", "b")))
  out <- pesto_evaluate(fm, theta)
  expect_equal(dim(out), c(2L, 2L))
  expect_true(all(is.finite(out)))
})
