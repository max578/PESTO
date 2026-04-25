# Sign-convention regression test for ensemble_solution().
#
# The C++ kernel implements:
#   Delta_theta = -Delta_theta' V Sigma (Sigma^2 + (lambda+1)I)^{-1} U^T . obs_resid
# Combined with the gradient of Phi = ||W(g(theta) - d_obs)||^2, the leading
# minus requires obs_resid = sim - obs for the update to act as a descent
# step. Passing obs - sim inverts the gradient and phi diverges geometrically.
#
# This test fixes the convention on a well-conditioned linear inverse problem
# where convergence is unambiguous. See INVESTIGATIONS.md item I1
# (closed 2026-04-25) for the full diagnosis.

test_that("ensemble_solution drives phi monotonically down with obs_resid = sim - obs", {
  set.seed(20260425L)
  npar  <- 6L
  nobs  <- 30L
  nreal <- 40L
  n_iter <- 4L

  # Linear forward operator with mild conditioning.
  G <- matrix(rnorm(nobs * npar) / sqrt(npar), nobs, npar)
  theta_true <- rnorm(npar)
  obs_noise_sd <- 0.05
  y_obs <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = obs_noise_sd)

  weights    <- rep(1 / obs_noise_sd, nobs)
  parcov_inv <- rep(1.0, npar)

  run_forward <- function(P) G %*% P

  # Common starting ensemble (so both runs see identical priors).
  set.seed(7L)
  par0 <- matrix(rnorm(npar * nreal, sd = 1.0), npar, nreal)
  obs0 <- run_forward(par0)
  Y    <- matrix(rep(y_obs, nreal), nobs, nreal)
  phi0 <- mean(compute_phi(Y - obs0, weights))

  ies_run <- function(sign_convention) {
    par_ens <- par0
    obs_ens <- obs0
    phi_trace <- numeric(n_iter + 1L)
    phi_trace[1L] <- phi0
    for (it in seq_len(n_iter)) {
      pm <- rowMeans(par_ens)
      om <- rowMeans(obs_ens)
      pd <- par_ens - pm
      od <- obs_ens - om
      obs_resid <- switch(sign_convention,
                          sim_minus_obs = obs_ens - Y,
                          obs_minus_sim = Y - obs_ens)
      Am <- matrix(rnorm(npar * (nreal - 1L)), npar, nreal - 1L)
      upg <- ensemble_solution(
        par_diff   = pd,
        obs_diff   = od,
        obs_resid  = obs_resid,
        par_resid  = pd,
        weights    = weights,
        parcov_inv = parcov_inv,
        Am         = Am,
        cur_lam    = 1.0,
        iter       = as.integer(it)
      )
      par_ens <- par_ens + t(upg)
      obs_ens <- run_forward(par_ens)
      phi_trace[it + 1L] <- mean(compute_phi(Y - obs_ens, weights))
    }
    phi_trace
  }

  # Documented convention: sim - obs.
  phi_correct <- ies_run("sim_minus_obs")

  # Strict monotone decrease across all iterations.
  expect_true(
    all(diff(phi_correct) < 0),
    info = sprintf("phi trace not monotone: %s",
                   paste(signif(phi_correct, 4), collapse = " -> "))
  )

  # Substantial reduction (well-conditioned linear problem converges fast).
  expect_lt(phi_correct[length(phi_correct)] / phi_correct[1L], 1e-2)

  # Inverted convention must diverge — protects against future docstring drift.
  phi_wrong <- ies_run("obs_minus_sim")
  expect_gt(phi_wrong[length(phi_wrong)], phi_wrong[1L])
})

test_that("compute_phi is sign-invariant in obs_resid (sanity)", {
  r <- matrix(rnorm(50), 10, 5)
  w <- abs(rnorm(10)) + 0.1
  expect_equal(compute_phi(r, w), compute_phi(-r, w))
})
