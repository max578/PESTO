# Finite-ensemble pathology countermeasures for the iterative ensemble
# smoother: covariance inflation (against under-dispersion) and covariance
# localisation (against spurious finite-sample correlations). The user-facing
# specifications [pesto_inflation()] and [pesto_localisation()] are
# lightweight, validated control objects; the loop-side machinery lives in the
# leading-dot helpers below and is shared by [pesto_ies_callback()] and
# [pesto_ies_filter()] through [.ies_apply_update()].


# -- Inflation ----------------------------------------------------------------

#' Covariance Inflation Specification for IES
#'
#' Builds a control object that tells [pesto_ies_callback()] and
#' [pesto_ies_filter()] how to counteract ensemble under-dispersion --- the
#' progressive collapse of posterior spread that an iterative ensemble smoother
#' suffers with a finite ensemble. Inflation re-expands the post-update
#' parameter spread each iteration; the default `method = "none"` leaves the
#' update byte-identical to the un-inflated smoother.
#'
#' Four methods are offered. `"rtps"` is relaxation-to-prior-spread (Whitaker &
#' Hamill 2012): each parameter's posterior anomalies are rescaled by
#' \eqn{\alpha(\sigma^{b} - \sigma^{a})/\sigma^{a} + 1}, where \eqn{\sigma^{b}}
#' and \eqn{\sigma^{a}} are the background (pre-update) and analysis
#' (post-update) standard deviations. Being per-parameter, it re-inflates the
#' directions that collapsed hardest, so it is the spectrally-aware workhorse.
#' `"adaptive"` is a global, magnitude-targeting scheme: it measures the mean
#' spread-retention ratio \eqn{q = \mathrm{mean}_j(\sigma^{a}_j/\sigma^{b}_j)}
#' and, when `q` falls below `retention_floor`, applies a single multiplicative
#' factor \eqn{\min(\texttt{max\_factor}, \texttt{retention\_floor}/q)} to
#' restore the lost variance magnitude. `"multiplicative"` applies a fixed
#' `factor` every iteration. `"none"` disables inflation.
#'
#' The companion *diagnostic* is the spectral spread-ESS
#' ([ensemble_spread_ess()]), recorded each iteration regardless of method: it
#' reports the effective number of variance-carrying directions and is what
#' detects directional collapse. Because that participation ratio is invariant
#' to a global rescaling, a global (`"multiplicative"` / `"adaptive"`) inflation
#' restores variance *magnitude* but not the spectral *shape*; `"rtps"` is the
#' method that reshapes the spectrum. The two compose well.
#'
#' @param method Character. One of `"none"` (default), `"rtps"`, `"adaptive"`,
#'   `"multiplicative"`.
#' @param alpha Numeric in \[0, 1\]. RTPS relaxation coefficient (default 0.5);
#'   used only when `method = "rtps"`.
#' @param factor Numeric \eqn{\ge} 1. Fixed inflation factor for
#'   `method = "multiplicative"` (default 1, i.e. no inflation).
#' @param retention_floor Numeric in (0, 1\]. Target floor on the mean
#'   spread-retention ratio for `method = "adaptive"` (default 0.5).
#' @param max_factor Numeric \eqn{\ge} 1. Upper bound on any single-iteration
#'   inflation factor for `"rtps"` and `"adaptive"` (default 5).
#' @return An object of class `"pesto_inflation"`.
#' @references
#' Whitaker, J.S. & Hamill, T.M. (2012). Evaluating methods to account for
#' system errors in ensemble data assimilation. *Monthly Weather Review*,
#' 140(9), 3078--3089.
#' @seealso [pesto_localisation()], [ensemble_spread_ess()],
#'   [pesto_ies_callback()].
#' @examples
#' inf <- pesto_inflation("rtps", alpha = 0.5)
#' inf
#' @export
pesto_inflation <- function(method = c("none", "rtps", "adaptive",
                                       "multiplicative"),
                            alpha = 0.5,
                            factor = 1.0,
                            retention_floor = 0.5,
                            max_factor = 5.0) {
  method <- match.arg(method)
  .check_scalar_in(alpha, "alpha", 0, 1)
  .check_scalar_ge(factor, "factor", 1)
  .check_scalar_in(retention_floor, "retention_floor", 0, 1,
                   include_lower = FALSE)
  .check_scalar_ge(max_factor, "max_factor", 1)

  structure(
    list(
      method          = method,
      alpha           = as.numeric(alpha),
      factor          = as.numeric(factor),
      retention_floor = as.numeric(retention_floor),
      max_factor      = as.numeric(max_factor)
    ),
    class = "pesto_inflation"
  )
}

#' @export
print.pesto_inflation <- function(x, ...) {
  cat("<pesto_inflation>\n")
  cat("  method:", x$method, "\n")
  if (x$method == "rtps") {
    cat("  alpha: ", x$alpha, " (max_factor ", x$max_factor, ")\n", sep = "")
  } else if (x$method == "adaptive") {
    cat("  retention_floor: ", x$retention_floor,
        " (max_factor ", x$max_factor, ")\n", sep = "")
  } else if (x$method == "multiplicative") {
    cat("  factor:", x$factor, "\n")
  }
  invisible(x)
}

# Apply the inflation spec to a post-update ensemble. `par_post` and
# `par_background` are nreal x npar matrices (post- and pre-update). Returns the
# inflated ensemble plus the realised factor and spread-retention ratio.
.apply_inflation <- function(par_post, par_background, inflation) {
  npar     <- ncol(par_post)
  mean_a   <- colMeans(par_post)
  anom     <- sweep(par_post, 2L, mean_a, "-")

  if (inflation$method == "multiplicative") {
    fac <- inflation$factor
    par_new <- sweep(anom * fac, 2L, mean_a, "+")
    return(list(par = par_new, factor = fac, retention = NA_real_))
  }

  sd_b <- .col_sd(par_background)
  sd_a <- .col_sd(par_post)

  if (inflation$method == "rtps") {
    fac <- rep(1.0, npar)
    safe <- sd_a > 0
    fac[safe] <- inflation$alpha * (sd_b[safe] - sd_a[safe]) / sd_a[safe] + 1.0
    fac <- pmin(pmax(fac, 0.0), inflation$max_factor)
    par_new <- sweep(sweep(anom, 2L, fac, "*"), 2L, mean_a, "+")
    safe_ret <- sd_b > 0 & sd_a > 0
    q <- if (any(safe_ret)) mean(sd_a[safe_ret] / sd_b[safe_ret]) else NA_real_
    return(list(par = par_new, factor = fac, retention = q))
  }

  # method == "adaptive": global magnitude-targeting inflation.
  safe <- sd_b > 0 & sd_a > 0
  q <- if (any(safe)) mean(sd_a[safe] / sd_b[safe]) else 1.0
  fac <- if (is.finite(q) && q > 0 && q < inflation$retention_floor) {
    min(inflation$max_factor, inflation$retention_floor / q)
  } else {
    1.0
  }
  par_new <- sweep(anom * fac, 2L, mean_a, "+")
  list(par = par_new, factor = fac, retention = q)
}


# -- Localisation -------------------------------------------------------------

#' Covariance Localisation Specification for IES
#'
#' Builds a control object that tells [pesto_ies_callback()] and
#' [pesto_ies_filter()] how to taper the ensemble Kalman gain, suppressing the
#' spurious long-range parameter-observation correlations that a finite
#' ensemble manufactures. Localisation is applied as a Schur (elementwise)
#' product on the explicit gain inside [ensemble_solution_localised()]; the
#' default `method = "none"` leaves the standard SVD update untouched.
#'
#' `"correlation"` is the iterative-ensemble-smoother-native automatic
#' localisation of Luo & Bhakta (2020): it needs no parameter or observation
#' coordinates, estimating a noise floor from the ensemble itself and damping
#' sample correlations that fall below it (see [correlation_localisation()]).
#' This is the recommended default for parameter-estimation problems whose
#' parameters carry no spatial metric. `"distance"` is classical
#' distance-based localisation: a Gaspari-Cohn taper ([gaspari_cohn()]) of a
#' parameter-to-observation distance matrix, for problems where such a metric
#' exists --- supply either `distances` directly or `par_coords` + `obs_coords`
#' (Euclidean distances are then computed), together with `radius`.
#'
#' @param method Character. One of `"none"` (default), `"correlation"`,
#'   `"distance"`.
#' @param taper Character. `"hard"` (default) or `"soft"`; passed to
#'   [correlation_localisation()] for `method = "correlation"`.
#' @param threshold Numeric. Correlation noise floor; negative (default -1)
#'   triggers automatic per-iteration estimation. `method = "correlation"`.
#' @param n_shuffle Integer \eqn{\ge} 1. Permutation replicates for the
#'   automatic floor (default 1). `method = "correlation"`.
#' @param quantile Numeric in (0, 1). Quantile of the spurious-correlation
#'   distribution used as the floor (default 0.95). `method = "correlation"`.
#' @param distances Matrix (npar x nobs) or `NULL`. Precomputed
#'   parameter-to-observation distances for `method = "distance"`.
#' @param par_coords,obs_coords Matrices (npar x d, nobs x d) or `NULL`.
#'   Parameter / observation coordinates; Euclidean distances are derived when
#'   `distances` is `NULL`. `method = "distance"`.
#' @param radius Numeric (> 0) or `NULL`. Gaspari-Cohn localisation radius;
#'   required for `method = "distance"`.
#' @return An object of class `"pesto_localisation"`.
#' @references
#' Luo, X. & Bhakta, T. (2020). Automatic and adaptive localization for
#' ensemble-based history matching. *Journal of Petroleum Science and
#' Engineering*, 184, 106559.
#' @seealso [pesto_inflation()], [correlation_localisation()], [gaspari_cohn()].
#' @examples
#' loc <- pesto_localisation("correlation", taper = "soft")
#' loc
#' @export
pesto_localisation <- function(method = c("none", "correlation", "distance"),
                               taper = c("hard", "soft"),
                               threshold = -1.0,
                               n_shuffle = 1L,
                               quantile = 0.95,
                               distances = NULL,
                               par_coords = NULL,
                               obs_coords = NULL,
                               radius = NULL) {
  method <- match.arg(method)
  taper  <- match.arg(taper)
  .check_scalar_in(quantile, "quantile", 0, 1, include_lower = FALSE,
                   include_upper = FALSE)
  n_shuffle <- as.integer(n_shuffle)
  if (length(n_shuffle) != 1L || is.na(n_shuffle) || n_shuffle < 1L) {
    stop("`n_shuffle` must be a positive integer scalar.", call. = FALSE)
  }

  if (method == "distance") {
    if (is.null(radius) || !is.numeric(radius) || length(radius) != 1L ||
        radius <= 0) {
      stop("`method = \"distance\"` requires a positive scalar `radius`.",
           call. = FALSE)
    }
    has_dist   <- !is.null(distances)
    has_coords <- !is.null(par_coords) && !is.null(obs_coords)
    if (!has_dist && !has_coords) {
      stop(
        paste0("`method = \"distance\"` requires either `distances` or both ",
               "`par_coords` and `obs_coords`."),
        call. = FALSE
      )
    }
    if (has_dist) {
      distances <- as.matrix(distances)
      if (any(distances < 0)) {
        stop("`distances` must be non-negative.", call. = FALSE)
      }
    }
  }

  structure(
    list(
      method     = method,
      taper      = taper,
      threshold  = as.numeric(threshold),
      n_shuffle  = n_shuffle,
      quantile   = as.numeric(quantile),
      distances  = distances,
      par_coords = par_coords,
      obs_coords = obs_coords,
      radius     = if (is.null(radius)) NULL else as.numeric(radius)
    ),
    class = "pesto_localisation"
  )
}

#' @export
print.pesto_localisation <- function(x, ...) {
  cat("<pesto_localisation>\n")
  cat("  method:", x$method, "\n")
  if (x$method == "correlation") {
    cat("  taper: ", x$taper,
        if (x$threshold < 0) "  (automatic floor)" else
          paste0("  (threshold ", x$threshold, ")"),
        "\n", sep = "")
  } else if (x$method == "distance") {
    cat("  radius:", x$radius, "\n")
  }
  invisible(x)
}

# Build the npar x nobs localisation taper for one update step from the
# current anomalies. Returns a list(rho, threshold, frac_active) or NULL when
# localisation is inactive. `par_diff` / `obs_diff` are npar x nreal and
# nobs x nreal (the anomalies already formed by the caller).
.localisation_rho <- function(localisation, par_diff, obs_diff) {
  if (is.null(localisation) || localisation$method == "none") {
    return(NULL)
  }
  npar <- nrow(par_diff)
  nobs <- nrow(obs_diff)

  if (localisation$method == "correlation") {
    rb <- correlation_localisation(
      par_diff  = par_diff,
      obs_diff  = obs_diff,
      threshold = localisation$threshold,
      taper     = localisation$taper,
      n_shuffle = localisation$n_shuffle,
      quantile  = localisation$quantile
    )
    return(list(rho = rb$rho, threshold = rb$threshold,
                frac_active = rb$frac_active))
  }

  # method == "distance"
  dmat <- localisation$distances
  if (is.null(dmat)) {
    dmat <- .euclidean_distances(localisation$par_coords,
                                 localisation$obs_coords)
  }
  if (nrow(dmat) != npar || ncol(dmat) != nobs) {
    stop(
      sprintf(paste0("Localisation distance matrix is %d x %d but the update ",
                     "needs %d x %d (npar x nobs)."),
              nrow(dmat), ncol(dmat), npar, nobs),
      call. = FALSE
    )
  }
  rho <- gaspari_cohn(dmat, localisation$radius)
  list(rho = rho, threshold = NA_real_,
       frac_active = mean(rho > 0))
}

# Pairwise Euclidean distances between rows of par_coords (npar x d) and rows
# of obs_coords (nobs x d), returning a npar x nobs matrix.
.euclidean_distances <- function(par_coords, obs_coords) {
  par_coords <- as.matrix(par_coords)
  obs_coords <- as.matrix(obs_coords)
  if (ncol(par_coords) != ncol(obs_coords)) {
    stop("`par_coords` and `obs_coords` must have the same number of columns.",
         call. = FALSE)
  }
  np <- nrow(par_coords)
  no <- nrow(obs_coords)
  d  <- matrix(0.0, nrow = np, ncol = no)
  for (j in seq_len(no)) {
    diff <- sweep(par_coords, 2L, obs_coords[j, ], "-")
    d[, j] <- sqrt(rowSums(diff * diff))
  }
  d
}


# -- Shared update core -------------------------------------------------------

# One IES GLM update step with optional localisation and inflation, shared by
# the callback and filter drivers. `par_ok` / `obs_ok` are the
# complete-case nreal_ok x npar and nreal_ok x nobs matrices for this step;
# `obs_vec` is the length-nobs target; `prior_mean` the (window- or run-level)
# prior parameter mean. Returns the updated ensemble block, the per-realisation
# phi, and a diagnostics list (spread-ESS + inflation + localisation).
.ies_apply_update <- function(par_ok, obs_ok, obs_vec, weights, parcov_inv,
                              prior_mean, lambda, eigthresh, use_approx, iter,
                              inflation = NULL, localisation = NULL) {
  nreal_ok <- nrow(par_ok)
  npar     <- ncol(par_ok)
  nobs     <- ncol(obs_ok)

  par_mean  <- colMeans(par_ok)
  obs_mean  <- colMeans(obs_ok)
  par_diff  <- t(sweep(par_ok, 2L, par_mean, "-"))        # npar x nreal_ok
  obs_diff  <- t(sweep(obs_ok, 2L, obs_mean, "-"))        # nobs x nreal_ok
  obs_resid <- matrix(obs_vec, nrow = nobs, ncol = nreal_ok) - t(obs_ok)
  par_resid <- t(sweep(par_ok, 2L, prior_mean, "-"))

  phi_vec <- compute_phi(obs_resid, weights)

  # Localisation: build the taper (per-iteration for the correlation method).
  loc_block <- .localisation_rho(localisation, par_diff, obs_diff)
  loc_label <- if (is.null(loc_block)) "none" else localisation$method
  loc_thresh <- if (is.null(loc_block)) NA_real_ else loc_block$threshold
  loc_frac   <- if (is.null(loc_block)) NA_real_ else loc_block$frac_active

  if (!is.null(loc_block)) {
    if (!use_approx) {
      warning(
        paste0("Localisation uses the approximate (upgrade_1) update; the ",
               "`use_approx = FALSE` null-space correction is ignored while ",
               "localisation is active."),
        call. = FALSE
      )
    }
    upgrade <- ensemble_solution_localised(
      par_diff  = par_diff,
      obs_diff  = obs_diff,
      obs_resid = obs_resid,
      weights   = weights,
      rho       = loc_block$rho,
      cur_lam   = lambda,
      eigthresh = eigthresh
    )
  } else {
    Am_k <- matrix(0.0, nrow = npar, ncol = max(nreal_ok - 1L, 1L))
    upgrade <- ensemble_solution(
      par_diff          = par_diff,
      obs_diff          = obs_diff,
      obs_resid         = obs_resid,
      par_resid         = par_resid,
      weights           = weights,
      parcov_inv        = parcov_inv,
      Am                = Am_k,
      cur_lam           = lambda,
      eigthresh         = eigthresh,
      use_approx        = use_approx,
      use_prior_scaling = FALSE,
      iter              = as.integer(iter),
      reg_factor        = -1.0
    )
  }

  par_post <- par_ok - upgrade

  inf_label  <- "none"
  inf_factor <- 1.0
  retention  <- NA_real_
  if (!is.null(inflation) && inflation$method != "none") {
    ai <- .apply_inflation(par_post, par_ok, inflation)
    par_post   <- ai$par
    inf_label  <- inflation$method
    inf_factor <- if (length(ai$factor) > 1L) mean(ai$factor) else ai$factor
    retention  <- ai$retention
  }

  # Spread-ESS collapse diagnostic on the post-update (post-inflation)
  # ensemble. Computed even when both countermeasures are off, so the manifest
  # always carries the dispersion trace.
  post_mean     <- colMeans(par_post)
  par_diff_post <- t(sweep(par_post, 2L, post_mean, "-"))
  se <- ensemble_spread_ess(par_diff_post)

  list(
    par_ok_new = par_post,
    phi        = phi_vec,
    diag       = list(
      spread_ess       = se$ess,
      spread_ess_ratio = se$ess_ratio,
      inflation_method = inf_label,
      inflation_factor = inf_factor,
      retention        = retention,
      localisation     = loc_label,
      loc_threshold    = loc_thresh,
      loc_frac_active  = loc_frac
    )
  )
}


# -- Small validators ---------------------------------------------------------

# Column-wise standard deviation of a matrix (nreal x npar) -> length-npar.
.col_sd <- function(m) {
  apply(m, 2L, stats::sd)
}

.check_scalar_in <- function(x, name, lower, upper,
                             include_lower = TRUE, include_upper = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a numeric scalar.", name), call. = FALSE)
  }
  lo_ok <- if (include_lower) x >= lower else x > lower
  hi_ok <- if (include_upper) x <= upper else x < upper
  if (!lo_ok || !hi_ok) {
    stop(
      sprintf("`%s` must lie in %s%g, %g%s.", name,
              if (include_lower) "[" else "(", lower, upper,
              if (include_upper) "]" else ")"),
      call. = FALSE
    )
  }
  invisible(x)
}

.check_scalar_ge <- function(x, name, lower) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower) {
    stop(sprintf("`%s` must be a numeric scalar >= %g.", name, lower),
         call. = FALSE)
  }
  invisible(x)
}

# Driver-side guards: a countermeasure argument must be NULL or the matching
# specification object (built by its constructor, so already self-validated).
.check_inflation <- function(inflation) {
  if (!is.null(inflation) && !inherits(inflation, "pesto_inflation")) {
    stop("`inflation` must be NULL or a `pesto_inflation()` object.",
         call. = FALSE)
  }
  invisible(inflation)
}

.check_localisation <- function(localisation) {
  if (!is.null(localisation) &&
      !inherits(localisation, "pesto_localisation")) {
    stop("`localisation` must be NULL or a `pesto_localisation()` object.",
         call. = FALSE)
  }
  invisible(localisation)
}
