# pesto_filter.R -- Sequential (filter-mode) Iterative Ensemble Smoother.
#
# Where pesto_ies_callback() is a *smoother* (assimilate all observations
# in one batch), pesto_ies_filter() is a *filter*: it assimilates
# time-ordered observation windows one after another, the posterior of
# each window becoming the prior of the next. Calibration parameters are
# static; what evolves is the information assimilated. The use case is
# in-season assimilation -- stream a growing season's observations
# (satellite, soil-moisture) into a running ensemble and read off a
# tightening parameter posterior after each window, rather than only at
# season end.
#
# It reuses the §1 forward-model contract (so it is parallel- and
# multi-fidelity-ready through pesto_evaluate()) and the C++
# ensemble_solution() GLM kernel (Chen & Oliver 2013).

#' Run a Sequential (Filter-Mode) Iterative Ensemble Smoother
#'
#' Assimilates observations in time-ordered **windows** against a static
#' parameter ensemble. Each window runs an IES Gauss-Levenberg-Marquardt
#' update using only that window's observation block; the updated ensemble
#' is carried forward as the prior for the next window. This is the
#' filtering analogue of [pesto_ies_callback()]: instead of one batch
#' assimilation of all observations, the posterior is refined window by
#' window and is available after each.
#'
#' @section Filter vs smoother:
#' [pesto_ies_callback()] forms parameter residuals against the *original*
#' prior mean and assimilates every observation together. Here each window
#' forms residuals against the *current* (carried-forward) ensemble mean
#' and assimilates only its own block, so information accrues sequentially.
#' With the default `use_approx = TRUE` the carried-forward ensemble itself
#' encodes the background covariance, so the sequential behaviour comes
#' from propagating the ensemble, not from an explicit prior term.
#'
#' @section Multi-fidelity:
#' When `forward_model` is a [pesto_multifidelity_model()], `fidelity_schedule`
#' selects the fidelity level evaluated at each *window* (recycled / padded
#' to the number of windows; default: highest fidelity throughout). The
#' final ensemble refresh always uses the highest fidelity.
#'
#' @param forward_model A function `function(theta) -> obs`, a
#'   [pesto_forward_model()], or a [pesto_multifidelity_model()] (then see
#'   `fidelity_schedule`). The model maps the full `nreal x npar` parameter
#'   matrix to the full `nreal x nobs` observation matrix; windows select
#'   columns of that output. A bare function is auto-wrapped via
#'   [as_forward_model()]. See [pesto_ies_callback()] for the contract.
#' @param prior_ensemble Matrix or data.table, `nreal x npar`.
#' @param obs Named numeric vector of length `nobs`. The full set of
#'   target observations; `windows` indexes into it.
#' @param obs_sd Numeric scalar or length-`nobs` vector. Observation
#'   standard deviation(s); weights are `1 / obs_sd`.
#' @param windows A list of integer vectors. Each element gives the
#'   observation indices (into `obs`, `1`-based) assimilated at that
#'   window, in assimilation order. Windows must be **disjoint** (no
#'   observation assimilated twice) and need not cover all of `obs`.
#' @param window_noptmax Integer scalar or vector. IES iterations *within*
#'   each window (default `1L`, the pure filter; `> 1` gives an iterated
#'   filter per window). A scalar is recycled; a short vector is
#'   right-padded with its last value.
#' @param lambda Numeric scalar or per-window vector. Marquardt lambda
#'   (default `1.0`; recycled / right-padded across windows).
#' @param fidelity_schedule Integer vector or `NULL`. Per-window fidelity
#'   level for a [pesto_multifidelity_model()] (see *Multi-fidelity*).
#' @param parcov Numeric vector of length `npar`, the diagonal of the
#'   prior parameter covariance. Defaults to the prior ensemble's
#'   column-wise variance (non-positive entries replaced with `1`).
#' @param eigthresh Numeric. SVD eigenvalue truncation (default `1e-6`).
#' @param use_approx Logical. If `TRUE` (default) skip the prior-scaling
#'   correction, matching the `pestpp-ies` default.
#' @param inflation A [pesto_inflation()] specification, or `NULL` (default).
#'   Applied within each window's inner update to counteract ensemble
#'   under-dispersion. See [pesto_ies_callback()].
#' @param localisation A [pesto_localisation()] specification, or `NULL`
#'   (default). Tapers the per-window Kalman gain to suppress spurious
#'   finite-sample correlations. See [pesto_ies_callback()].
#' @param on_failure Character. `"na"` (default) tolerates failed
#'   realisations; `"stop"` aborts on any failure.
#' @param verbose Logical. Print a per-window phi summary.
#'
#' @return A list of class
#'   `c("pesto_ies_filter_result", "pesto_ies_result")` with components:
#'   \describe{
#'     \item{phi}{data.table of per-realisation phi by window (on that
#'       window's observation block).}
#'     \item{par_ensemble}{Final parameter ensemble (data.table).}
#'     \item{obs_ensemble}{Final simulated-observation ensemble, full
#'       `nobs` columns (data.table).}
#'     \item{windows}{List of per-window metadata: assimilated indices,
#'       lambda, mean phi, per-parameter ensemble mean and standard
#'       deviation (the sd trace shows the posterior tightening), and
#'       failure count.}
#'     \item{runtime_seconds, n_forward_evals, failure_rate}{Run totals.}
#'     \item{fidelity}{Multi-fidelity provenance (or `NULL`), as in
#'       [pesto_ies_callback()].}
#'   }
#' @references
#' Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
#' iterative ensemble smoother for efficient history matching and
#' uncertainty quantification. *Computational Geosciences*, 17(4),
#' 689--703.
#' @seealso [pesto_ies_callback()] for the batch smoother;
#'   [pesto_multifidelity_model()] for fidelity stacks; [as_manifest()]
#'   to wrap the result in the ensemble-manifest contract.
#' @export
#' @examples
#' # Linear-Gaussian recovery, assimilated in three observation windows.
#' set.seed(1)
#' npar <- 3; nobs <- 9; nreal <- 80
#' G <- matrix(rnorm(nobs * npar), nobs, npar)
#' theta_true <- c(1.0, -0.5, 2.0)
#' y <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = 0.05)
#' f <- function(theta) theta %*% t(G)
#' prior <- matrix(rnorm(nreal * npar), nreal, npar,
#'                 dimnames = list(NULL, paste0("p", 1:npar)))
#' fit <- pesto_ies_filter(
#'   forward_model  = f, prior_ensemble = prior,
#'   obs = setNames(y, paste0("o", 1:nobs)), obs_sd = 0.05,
#'   windows = list(1:3, 4:6, 7:9), verbose = FALSE
#' )
#' # Posterior sd should shrink window over window:
#' vapply(fit$windows, function(w) mean(w$par_sd), numeric(1))
pesto_ies_filter <- function(forward_model,
                             prior_ensemble,
                             obs,
                             obs_sd,
                             windows,
                             window_noptmax = 1L,
                             lambda = 1.0,
                             fidelity_schedule = NULL,
                             parcov = NULL,
                             eigthresh = 1e-6,
                             use_approx = TRUE,
                             inflation = NULL,
                             localisation = NULL,
                             on_failure = c("na", "stop"),
                             verbose = TRUE) {

  # Validate inputs -----------------------------------------------------
  on_failure <- match.arg(on_failure)
  .check_pesto_ies_callback_inputs(forward_model, 1L, eigthresh)
  .check_inflation(inflation)
  .check_localisation(localisation)

  # Coerce prior ensemble -----------------------------------------------
  if (data.table::is.data.table(prior_ensemble) ||
      is.data.frame(prior_ensemble)) {
    par_names_local <- setdiff(names(prior_ensemble), "real_name")
    par_mat <- as.matrix(
      data.table::as.data.table(prior_ensemble)[, par_names_local,
                                                with = FALSE]
    )
  } else {
    par_mat <- as.matrix(prior_ensemble)
    par_names_local <- colnames(par_mat)
    if (is.null(par_names_local)) {
      par_names_local <- paste0("par", seq_len(ncol(par_mat)))
      colnames(par_mat) <- par_names_local
    }
  }
  storage.mode(par_mat) <- "double"
  nreal <- nrow(par_mat)
  npar  <- ncol(par_mat)
  if (nreal < 2L) {
    stop("`prior_ensemble` must contain at least 2 realisations.",
      call. = FALSE
    )
  }

  # Coerce observations + weights ---------------------------------------
  obs_vec <- as.numeric(obs)
  obs_names_local <- names(obs)
  if (is.null(obs_names_local)) {
    obs_names_local <- paste0("obs", seq_along(obs_vec))
  }
  nobs <- length(obs_vec)

  obs_sd_vec <- as.numeric(obs_sd)
  if (length(obs_sd_vec) == 1L) obs_sd_vec <- rep(obs_sd_vec, nobs)
  if (length(obs_sd_vec) != nobs || any(obs_sd_vec <= 0)) {
    stop("`obs_sd` must be a positive scalar or length-nobs vector.",
      call. = FALSE
    )
  }
  weights <- 1.0 / obs_sd_vec

  # Validate windows ----------------------------------------------------
  .check_windows(windows, nobs)
  nwindows <- length(windows)

  # Prior covariance diagonal -------------------------------------------
  if (is.null(parcov)) {
    parcov_diag <- apply(par_mat, 2L, stats::var)
    parcov_diag[parcov_diag <= 0 | !is.finite(parcov_diag)] <- 1.0
  } else {
    parcov_diag <- as.numeric(parcov)
    if (length(parcov_diag) != npar || any(parcov_diag <= 0)) {
      stop(sprintf("`parcov` must be a positive length-%d vector.", npar),
        call. = FALSE
      )
    }
  }
  parcov_inv <- 1.0 / parcov_diag

  # Per-window lambda + inner-iteration schedules -----------------------
  lambda_seq    <- .recycle_to(as.numeric(lambda), nwindows)
  inner_seq     <- .recycle_to(as.integer(window_noptmax), nwindows)
  if (any(inner_seq < 1L)) {
    stop("`window_noptmax` entries must be positive integers (>= 1).",
      call. = FALSE
    )
  }

  # Forward-model evaluation contract (mirrors pesto_ies_callback) ------
  is_mf <- S7::S7_inherits(forward_model, pesto_multifidelity_model)
  if (is_mf) {
    n_levels   <- length(forward_model@levels)
    fid_sched  <- .resolve_fidelity_schedule(fidelity_schedule,
                                             n_levels, nwindows)
    top_level  <- n_levels - 1L
    eval_model <- forward_model
  } else {
    fid_sched  <- NULL
    top_level  <- 0L
    eval_model <- .fm_with_nobs(
      as_forward_model(forward_model, on_failure = on_failure), nobs
    )
  }
  lvl_for <- function(k) if (is_mf) fid_sched[k] else 0L
  eval_at <- function(pm, level) {
    if (is_mf) pesto_evaluate(eval_model, pm, level = level)
    else       pesto_evaluate(eval_model, pm)
  }

  # Sequential window loop ----------------------------------------------
  phi_history    <- vector("list", nwindows)
  win_meta       <- vector("list", nwindows)
  total_evals    <- 0L
  total_failures <- 0L

  t0 <- proc.time()["elapsed"]

  obs_mat <- eval_at(par_mat, lvl_for(1L))
  total_evals    <- total_evals + nreal
  total_failures <- total_failures + attr(obs_mat, "n_failures")

  for (k in seq_len(nwindows)) {
    idx          <- as.integer(windows[[k]])
    obs_blk_vec  <- obs_vec[idx]
    weights_blk  <- weights[idx]
    # Filter prior: residuals are taken against the ensemble carried in
    # to this window (recomputed once, held across its inner iterations).
    prior_mean_k <- colMeans(par_mat)

    last_phi  <- NULL
    last_ok   <- NULL
    last_diag <- NULL
    for (j in seq_len(inner_seq[k])) {
      blk <- .ies_glm_block(
        par_mat       = par_mat,
        obs_block     = obs_mat[, idx, drop = FALSE],
        obs_vec_block = obs_blk_vec,
        weights_block = weights_blk,
        parcov_inv    = parcov_inv,
        prior_mean    = prior_mean_k,
        lambda        = lambda_seq[k],
        eigthresh     = eigthresh,
        use_approx    = use_approx,
        window        = k,
        inflation     = inflation,
        localisation  = localisation
      )
      par_mat   <- blk$par_mat
      last_phi  <- blk$phi
      last_ok   <- blk$ok
      last_diag <- blk$diag
      if (j < inner_seq[k]) {
        obs_mat <- eval_at(par_mat, lvl_for(k))
        total_evals    <- total_evals + nreal
        total_failures <- total_failures + attr(obs_mat, "n_failures")
      }
    }

    phi_history[[k]] <- data.table::data.table(
      window      = k,
      realisation = last_ok,
      phi         = last_phi
    )
    par_mean_k <- colMeans(par_mat)
    par_sd_k   <- apply(par_mat, 2L, stats::sd)
    win_meta[[k]] <- c(
      list(
        window      = k,
        obs_indices = idx,
        lambda      = lambda_seq[k],
        mean_phi    = mean(last_phi),
        par_mean    = stats::setNames(par_mean_k, par_names_local),
        par_sd      = stats::setNames(par_sd_k, par_names_local),
        n_failures  = nreal - length(last_ok)
      ),
      if (is.null(last_diag)) list() else last_diag
    )
    if (verbose) {
      message(sprintf(
        paste0("[pesto_ies_filter] window %d/%d: %d obs, phi mean=%.4g, ",
               "mean par sd=%.4g (lambda=%.3g)"),
        k, nwindows, length(idx), mean(last_phi), mean(par_sd_k),
        lambda_seq[k]
      ))
    }

    # Carry the updated ensemble forward to the next window.
    if (k < nwindows) {
      obs_mat <- eval_at(par_mat, lvl_for(k + 1L))
      total_evals    <- total_evals + nreal
      total_failures <- total_failures + attr(obs_mat, "n_failures")
    }
  }

  # Final refresh at the highest fidelity -------------------------------
  obs_mat_final <- eval_at(par_mat, top_level)
  total_evals    <- total_evals + nreal
  total_failures <- total_failures + attr(obs_mat_final, "n_failures")

  runtime <- as.numeric(proc.time()["elapsed"] - t0)

  # Assemble result -----------------------------------------------------
  par_dt <- data.table::as.data.table(par_mat)
  data.table::setnames(par_dt, par_names_local)
  par_dt[, real_name := paste0("real_", seq_len(nreal))]
  data.table::setcolorder(par_dt, c("real_name", par_names_local))

  obs_dt <- data.table::as.data.table(obs_mat_final)
  data.table::setnames(obs_dt, obs_names_local)
  obs_dt[, real_name := paste0("real_", seq_len(nreal))]
  data.table::setcolorder(obs_dt, c("real_name", obs_names_local))

  fidelity_record <- if (is_mf) {
    list(
      type        = "multifidelity",
      schedule    = as.integer(fid_sched),
      final_level = as.integer(top_level),
      n_levels    = as.integer(n_levels),
      costs       = as.numeric(eval_model@costs)
    )
  } else {
    NULL
  }

  output <- list(
    phi             = data.table::rbindlist(phi_history),
    par_ensemble    = par_dt,
    obs_ensemble    = obs_dt,
    windows         = win_meta,
    runtime_seconds = runtime,
    n_forward_evals = total_evals,
    failure_rate    = total_failures / total_evals,
    fidelity        = fidelity_record,
    # Assimilation inputs preserved for as_manifest(). `iterations` is the
    # manifest-facing per-step record; for the filter, one entry per window,
    # carrying the dispersion / countermeasure diagnostics when present.
    iterations      = lapply(win_meta, function(w) {
      w[intersect(
        c("lambda", "spread_ess", "spread_ess_ratio", "inflation_method",
          "inflation_factor", "retention", "localisation", "loc_threshold",
          "loc_frac_active"),
        names(w)
      )]
    }),
    obs_target      = stats::setNames(obs_vec, obs_names_local),
    obs_sd          = stats::setNames(obs_sd_vec, obs_names_local),
    weights         = stats::setNames(weights, obs_names_local)
  )
  class(output) <- c("pesto_ies_filter_result", "pesto_ies_result")
  output
}


# Single-window IES GLM update on an observation block. Delegates the update
# math (anomalies, localised / standard solve, inflation, spread-ESS) to the
# shared `.ies_apply_update()` core, restricted to this window's column subset
# and against a window-local prior mean. Returns the updated ensemble, block
# phi, the complete-case indices, and the step diagnostics.
.ies_glm_block <- function(par_mat, obs_block, obs_vec_block,
                           weights_block, parcov_inv, prior_mean,
                           lambda, eigthresh, use_approx, window,
                           inflation = NULL, localisation = NULL) {
  ok <- stats::complete.cases(obs_block)
  if (sum(ok) < 2L) {
    stop(
      sprintf(
        paste0(
          "Window %d: fewer than 2 successful realisations. ",
          "Cannot continue."
        ),
        window
      ),
      call. = FALSE
    )
  }
  par_ok <- par_mat[ok, , drop = FALSE]
  obs_ok <- obs_block[ok, , drop = FALSE]

  step <- .ies_apply_update(
    par_ok       = par_ok,
    obs_ok       = obs_ok,
    obs_vec      = obs_vec_block,
    weights      = weights_block,
    parcov_inv   = parcov_inv,
    prior_mean   = prior_mean,
    lambda       = lambda,
    eigthresh    = eigthresh,
    use_approx   = use_approx,
    iter         = window,
    inflation    = inflation,
    localisation = localisation
  )
  par_mat[ok, ] <- step$par_ok_new
  list(par_mat = par_mat, phi = step$phi, ok = which(ok),
       nreal_ok = nrow(par_ok), diag = step$diag)
}


# Recycle / right-pad a vector to length n (last value repeated).
.recycle_to <- function(x, n) {
  if (length(x) >= n) return(x[seq_len(n)])
  c(x, rep(x[length(x)], n - length(x)))
}


# Validate the `windows` list: non-empty list of disjoint, in-range,
# 1-based integer index vectors.
.check_windows <- function(windows, nobs) {
  if (!is.list(windows) || length(windows) == 0L) {
    stop("`windows` must be a non-empty list of integer index vectors.",
      call. = FALSE
    )
  }
  seen <- integer(0)
  for (k in seq_along(windows)) {
    wi <- suppressWarnings(as.integer(windows[[k]]))
    if (length(wi) == 0L || anyNA(wi) || any(wi < 1L) || any(wi > nobs)) {
      stop(
        sprintf(
          "`windows[[%d]]` must be a non-empty vector of integers in [1, %d].",
          k, nobs
        ),
        call. = FALSE
      )
    }
    if (anyDuplicated(wi)) {
      stop(sprintf("`windows[[%d]]` has duplicate indices.", k),
        call. = FALSE
      )
    }
    if (length(intersect(wi, seen)) > 0L) {
      stop(
        sprintf(
          "`windows[[%d]]` overlaps an earlier window; windows must be disjoint.",
          k
        ),
        call. = FALSE
      )
    }
    seen <- c(seen, wi)
  }
  invisible(TRUE)
}
