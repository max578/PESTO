# obs_perturbation.R -- Per-realisation observation perturbation for the
# native R iterative ensemble smoother.
#
# An ensemble smoother that assimilates the same, *unperturbed* observation
# vector into every realisation does not return a posterior sample. Every
# member feels an identical data pull, so the between-member spread is
# whatever prior variance the update failed to remove -- a cloud around the
# maximum-likelihood fit, not a draw from p(theta | d). On a linear-Gaussian
# problem the effect is large and does not shrink with ensemble size.
#
# The standard remedy is to give each realisation its own perturbed data
# vector,
#
#   d_j = d + e_j,   e_j ~ N(0, alpha C_D),
#
# and, because an iterative smoother assimilates the *same* data at every
# step, to inflate the measurement covariance by alpha_k at step k subject to
#
#   sum_k 1 / alpha_k = 1
#
# so that the total information assimilated is exactly one likelihood
# (Emerick & Reynolds 2013; Evensen 2018). Together these two make the
# update the ensemble-smoother-with-multiple-data-assimilation (ES-MDA)
# analysis step
#
#   m_j <- m_j + C_MD (C_DD + alpha_k C_D)^-1 (d + sqrt(alpha_k) e_j - g(m_j)).
#
# PESTO's GLM kernel already computes exactly that form: `ensemble_solution()`
# builds the inverse term (s^2 + (cur_lam + 1) I)^-1 on weight-scaled
# anomalies, and with C_D = W^-2 the identity
#
#   alpha_k = cur_lam_k + 1
#
# makes the Marquardt damping and the MDA inflation the same number. The
# drivers therefore derive the lambda schedule from the alpha schedule; no
# kernel change is required.
#
# References
#   Emerick, A.A. & Reynolds, A.C. (2013). Ensemble smoother with multiple
#     data assimilation. Computers & Geosciences, 55, 3--15.
#   Evensen, G. (2018). Analysis of iterative ensemble smoothers for solving
#     inverse problems. Computational Geosciences, 22(3), 885--908.
#   Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
#     iterative ensemble smoother. Computational Geosciences, 17(4), 689--703.


# Draw one nobs x nreal observation-noise ensemble, e_j ~ N(0, scale * C_D),
# with C_D = diag(obs_sd^2). Observations whose standard deviation is not
# finite carry zero weight in the update, so they are given zero perturbation
# rather than an infinite one (an infinite draw would propagate NaN through
# the weight-scaled residual instead of being annihilated by the zero weight).
.draw_obs_noise <- function(nobs, nreal, obs_sd, scale = 1) {
  sd_vec <- as.numeric(obs_sd) * sqrt(scale)
  keep   <- is.finite(sd_vec) & sd_vec > 0
  e      <- matrix(0.0, nrow = nobs, ncol = nreal)
  if (any(keep)) {
    nk <- sum(keep)
    e[keep, ] <- matrix(stats::rnorm(nk * nreal), nrow = nk, ncol = nreal) *
      sd_vec[keep]
  }
  e
}


# Resolve the ES-MDA inflation schedule for `n_steps` assimilations of the
# same observation block, and return it together with the matching Marquardt
# lambda schedule (alpha = lambda + 1). Three sources, in precedence order:
# an explicit `mda_alpha`; a user-supplied `lambda` (read as alpha - 1); the
# uniform default alpha_k = n_steps of Emerick & Reynolds (2013). The
# sum(1 / alpha) == 1 condition is the definition of the scheme, so a
# schedule that violates it is an error, not a warning.
.mda_alpha_schedule <- function(mda_alpha, lambda, lambda_supplied, n_steps,
                                lambda_arg = "lambda") {
  if (!is.null(mda_alpha)) {
    alpha  <- .recycle_to(as.numeric(mda_alpha), n_steps)
    source <- "`mda_alpha`"
  } else if (isTRUE(lambda_supplied)) {
    alpha  <- .recycle_to(as.numeric(lambda), n_steps) + 1.0
    source <- sprintf("`%s` (read as alpha - 1)", lambda_arg)
  } else {
    alpha  <- rep(as.numeric(n_steps), n_steps)
    source <- "the default uniform ES-MDA schedule"
  }
  if (anyNA(alpha) || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop(
      sprintf("The ES-MDA inflation schedule from %s must be finite and > 0.",
              source),
      call. = FALSE
    )
  }
  total <- sum(1.0 / alpha)
  if (abs(total - 1.0) > 1e-8) {
    stop(
      sprintf(
        paste0(
          "The ES-MDA inflation schedule from %s gives sum(1 / alpha) = %.6g, ",
          "but the scheme requires exactly 1 (Emerick & Reynolds 2013). ",
          "Supply `mda_alpha` summing correctly in reciprocal (the uniform ",
          "choice is rep(%d, %d)), or set `obs_perturbation = \"none\"` for ",
          "the unperturbed, maximum-likelihood-seeking update."
        ),
        source, total, n_steps, n_steps
      ),
      call. = FALSE
    )
  }
  list(alpha = alpha, lambda = alpha - 1.0, source = source)
}


# Validate + normalise the `obs_perturbation` argument.
.check_obs_perturbation <- function(obs_perturbation) {
  match.arg(obs_perturbation, c("mda", "none"))
}
