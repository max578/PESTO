#' Surrogate-Accelerated IES Iteration
#'
#' Performs a single IES iteration using a Gaussian Process surrogate
#' to reduce the number of expensive full-model evaluations.
#'
#' **How it works:**
#' 1. A GP surrogate is trained on existing parameter-observation pairs
#' 2. The surrogate predicts model outputs for all ensemble members
#' 3. Only members with high prediction uncertainty trigger full model runs
#' 4. Control-variate bias correction blends surrogate and model results
#' 5. Standard IES update is computed on the blended ensemble
#'
#' This typically saves 50-90% of model evaluations per iteration.
#'
#' @param par_ensemble data.table or matrix. Current parameter ensemble
#'   (rows = realisations, columns = parameters).
#' @param obs_ensemble data.table or matrix. Current observation ensemble
#'   from model evaluations.
#' @param obs_target Named numeric vector. Target observation values.
#' @param weights Numeric vector. Observation weights.
#' @param parcov_inv Numeric vector. Diagonal of inverse parameter covariance.
#' @param lambda Numeric. Marquardt lambda (default 1.0).
#' @param uncertainty_threshold Numeric. Threshold for surrogate/model switching.
#'   Fraction of signal variance (default 0.1 = 10%).
#' @param eigthresh Numeric. SVD eigenvalue threshold.
#' @return A list containing:
#'   \describe{
#'     \item{upgrade}{Matrix of parameter upgrades}
#'     \item{n_model_runs}{Number of full model evaluations needed}
#'     \item{n_surrogate_runs}{Number of surrogate-only evaluations}
#'     \item{savings_pct}{Percentage of model runs saved}
#'     \item{gp_diagnostics}{GP training diagnostics}
#'   }
#' @references
#' Rasmussen, C.E. & Williams, C.K.I. (2006). Gaussian Processes for
#' Machine Learning. MIT Press.
#'
#' Liu, F. & Guillas, S. (2017). Dimension reduction for Gaussian process
#' emulation. *Statistics and Computing*, 27(3), 785-802.
#' @examples
#' \donttest{
#' set.seed(7L)
#' n_real <- 15L; n_par <- 5L; n_obs <- 8L
#' par_ens <- matrix(rnorm(n_real * n_par), n_real, n_par,
#'                   dimnames = list(NULL, paste0("k", 1:n_par)))
#' obs_ens <- matrix(rnorm(n_real * n_obs), n_real, n_obs,
#'                   dimnames = list(NULL, paste0("h", 1:n_obs)))
#' obs_target <- rnorm(n_obs)
#' weights    <- rep(1.0, n_obs)
#' parcov_inv <- rep(1.0, n_par)
#' res <- pesto_surrogate_ies(
#'   par_ensemble = par_ens,
#'   obs_ensemble = obs_ens,
#'   obs_target   = obs_target,
#'   weights      = weights,
#'   parcov_inv   = parcov_inv,
#'   lambda       = 1.0
#' )
#' res$savings_pct
#' }
#' @export
pesto_surrogate_ies <- function(par_ensemble,
                                obs_ensemble,
                                obs_target,
                                weights,
                                parcov_inv,
                                lambda = 1.0,
                                uncertainty_threshold = 0.1,
                                eigthresh = 1e-6) {

  if (data.table::is.data.table(par_ensemble)) {
    par_mat <- as.matrix(par_ensemble[, .SD, .SDcols = is.numeric])
  } else {
    par_mat <- as.matrix(par_ensemble)
  }

  if (data.table::is.data.table(obs_ensemble)) {
    obs_mat <- as.matrix(obs_ensemble[, .SD, .SDcols = is.numeric])
  } else {
    obs_mat <- as.matrix(obs_ensemble)
  }

  result <- surrogate_ensemble_update(
    par_ensemble = par_mat,
    obs_ensemble = obs_mat,
    obs_target   = as.numeric(obs_target),
    weights      = as.numeric(weights),
    parcov_inv   = as.numeric(parcov_inv),
    cur_lam      = lambda,
    uncertainty_threshold = uncertainty_threshold,
    eigthresh    = eigthresh
  )

  result
}


#' Plot Surrogate Diagnostics
#'
#' Visualises the surrogate-accelerated IES performance including
#' model savings, uncertainty distribution, and GP quality metrics.
#'
#' @param results List of surrogate update results from multiple iterations.
#' @param title Character. Plot title.
#' @return A ggplot2 object.
#' @examples
#' iter1 <- list(n_model_runs = 12L, n_surrogate_runs = 38L,
#'               savings_pct = 76.0, mean_uncertainty = 0.18)
#' iter2 <- list(n_model_runs = 8L,  n_surrogate_runs = 42L,
#'               savings_pct = 84.0, mean_uncertainty = 0.11)
#' iter3 <- list(n_model_runs = 5L,  n_surrogate_runs = 45L,
#'               savings_pct = 90.0, mean_uncertainty = 0.07)
#' p <- plot_surrogate_diagnostics(list(iter1, iter2, iter3))
#' inherits(p, "ggplot")
#' @export
plot_surrogate_diagnostics <- function(results,
                                       title = "Surrogate IES Diagnostics") {
  if (!is.list(results)) {
    stop("`results` must be a list of surrogate update results.",
      call. = FALSE
    )
  }

  # Accept either a single-iteration result or a list of them ------------
  if (!is.null(results$savings_pct)) {
    results <- list(results)
  }

  iter_data <- data.table::data.table(
    iteration    = seq_along(results),
    model_runs   = vapply(results, function(r) as.integer(r$n_model_runs), integer(1)),
    surr_runs    = vapply(results, function(r) as.integer(r$n_surrogate_runs), integer(1)),
    savings_pct  = vapply(results, function(r) r$savings_pct, numeric(1)),
    mean_uncert  = vapply(results, function(r) r$mean_uncertainty, numeric(1))
  )

  p <- ggplot2::ggplot(iter_data, ggplot2::aes(x = iteration)) +
    ggplot2::geom_col(
      ggplot2::aes(y = model_runs, fill = "Full model"),
      alpha = 0.8
    ) +
    ggplot2::geom_col(
      ggplot2::aes(y = surr_runs, fill = "Surrogate"),
      alpha = 0.8,
      position = "stack"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = model_runs + surr_runs,
                   label = sprintf("%.0f%% saved", savings_pct)),
      vjust = -0.5, size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      values = c("Full model" = "#D55E00", "Surrogate" = "#009E73")
    ) +
    ggplot2::labs(
      title = title,
      x = "Iteration",
      y = "Number of Evaluations",
      fill = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  p
}


#' Train a GP Surrogate with Maximum-Likelihood (Anisotropic) Length Scales
#'
#' [train_gp_surrogate()] defaults to a single median-heuristic length scale:
#' fast and robust, but on a strongly anisotropic or multi-scale response it can
#' be several times less accurate than length scales tuned to the data. This
#' helper keeps the same fast C++ GP but selects the length scale(s) by
#' **maximising the GP's own log marginal likelihood** -- the criterion a
#' maximum-likelihood Gaussian-process library (for example `DiceKriging`)
#' optimises.
#'
#' By default the fit is **anisotropic**: one length scale is estimated per input
#' dimension by pre-scaling each coordinate (the isotropic C++ kernel applied to
#' coordinates divided by per-axis length scales is an anisotropic kernel on the
#' originals), optimised with [stats::optim()] from several starts. On a strongly
#' anisotropic response this recovers most of the accuracy a single length scale
#' leaves on the table: on the Branin function it cuts held-out error roughly
#' three-fold and brings the surrogate to within a small factor of an anisotropic
#' MLE oracle, where a single length scale sits about seven-fold worse. With
#' `anisotropic = FALSE`, or a single input dimension, one length scale is tuned
#' by a log-spaced grid plus a [stats::optimize()] refinement.
#'
#' The response is centred before fitting (the C++ GP is zero-mean), and the
#' marginal variance is left at the GP's automatic value unless `signal_var` is
#' supplied.
#'
#' Because an anisotropic fit stores the GP on pre-scaled coordinates, **predict
#' with [predict_gp_surrogate_tuned()]**, which reapplies the pre-scaling and
#' adds the centred mean back. Calling [predict_gp_surrogate()] directly on a
#' tuned surrogate is correct only in the isotropic case.
#'
#' @param X_train Numeric matrix of training inputs (rows = points,
#'   columns = parameters).
#' @param Y_train Numeric matrix (or vector) of training outputs.
#' @param anisotropic Logical. If `TRUE` (default) and there are at least two
#'   input dimensions, estimate one length scale per dimension; otherwise a
#'   single length scale.
#' @param signal_var Numeric or `NULL` (default, the GP's automatic mean
#'   per-output response variance).
#' @param noise_var Numeric observation-noise variance / nugget. Default `1e-4`.
#' @param n_restarts Integer. Random restarts for the anisotropic optimisation;
#'   the best marginal likelihood is kept. Default `5`.
#' @param n_grid Integer. Length scales in the isotropic grid search. Default
#'   `40`.
#' @param length_scale_bounds Numeric `c(lower, upper)` for the isotropic search,
#'   or `NULL` (default) to derive it from the pairwise distances of `X_train`.
#'
#' @return The list from [train_gp_surrogate()] at the maximum-likelihood length
#'   scale(s), trained on centred (and, when anisotropic, pre-scaled) data, with
#'   an added `tuning` element: `anisotropic`, `length_scale` (per dimension when
#'   anisotropic, scalar otherwise), `input_scale` (the per-axis divisor applied
#'   before training; all ones when isotropic), `y_mean` (the per-output centring
#'   offset), `length_scale_median`, `log_marginal_likelihood` (at the optimum)
#'   and `log_marginal_likelihood_median` (at the single median heuristic, the
#'   baseline this improves on).
#'
#' @seealso [predict_gp_surrogate_tuned()] to predict from the result;
#'   [train_gp_surrogate()] for the default single-heuristic fit.
#'
#' @examples
#' set.seed(1L)
#' X <- matrix(runif(40L * 2L), 40L, 2L)
#' y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2
#' gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
#' gp$tuning$length_scale
#' pred <- predict_gp_surrogate_tuned(gp, X)
#' @export
train_gp_surrogate_tuned <- function(X_train, Y_train,
                                     anisotropic = TRUE,
                                     signal_var = NULL,
                                     noise_var = 1e-4,
                                     n_restarts = 5L,
                                     n_grid = 40L,
                                     length_scale_bounds = NULL) {
  X_train <- as.matrix(X_train)
  Y_train <- as.matrix(Y_train)
  if (nrow(X_train) != nrow(Y_train)) {
    stop("`X_train` and `Y_train` must have the same number of rows.",
      call. = FALSE
    )
  }
  if (nrow(X_train) < 3L) {
    stop("at least 3 training points are needed to tune a length scale.",
      call. = FALSE
    )
  }
  npar <- ncol(X_train)
  sv <- if (is.null(signal_var)) 0.0 else signal_var
  # The C++ GP is zero-mean; centre the response so it models the residual.
  y_mean <- colMeans(Y_train)
  Yc <- sweep(Y_train, 2L, y_mean, `-`)

  dists <- as.numeric(stats::dist(X_train))
  dists <- dists[dists > 0]
  med <- if (length(dists)) stats::median(dists) else 1
  if (!is.finite(med) || med <= 0) {
    med <- 1
  }
  lml_median <- train_gp_surrogate(X_train, Yc, length_scale = med,
    signal_var = sv, noise_var = noise_var)$log_marginal_likelihood

  if (isTRUE(anisotropic) && npar >= 2L) {
    # Per-axis median distance as the starting per-dimension length scales.
    axis_med <- vapply(seq_len(npar), function(j) {
      dj <- as.numeric(stats::dist(X_train[, j, drop = FALSE]))
      dj <- dj[dj > 0]
      m <- if (length(dj)) stats::median(dj) else med
      if (!is.finite(m) || m <= 0) med else m
    }, numeric(1L))
    neg_lml <- function(log_scale) {
      xs <- sweep(X_train, 2L, exp(log_scale), `/`)
      -train_gp_surrogate(xs, Yc, length_scale = 1, signal_var = sv,
        noise_var = noise_var)$log_marginal_likelihood
    }
    starts <- c(
      list(log(axis_med)),
      lapply(seq_len(max(0L, as.integer(n_restarts) - 1L)),
        function(i) log(axis_med) + stats::runif(npar, -1, 1))
    )
    best <- NULL
    for (s0 in starts) {
      op <- tryCatch(
        stats::optim(s0, neg_lml, method = "Nelder-Mead",
          control = list(maxit = 300L, reltol = 1e-7)),
        error = function(e) NULL
      )
      if (!is.null(op) && (is.null(best) || op$value < best$value)) {
        best <- op
      }
    }
    if (is.null(best)) {
      stop("anisotropic length-scale optimisation did not converge.",
        call. = FALSE
      )
    }
    input_scale <- exp(best$par)
    xs <- sweep(X_train, 2L, input_scale, `/`)
    gp <- train_gp_surrogate(xs, Yc, length_scale = 1, signal_var = sv,
      noise_var = noise_var)
    length_scale <- input_scale
  } else {
    if (is.null(length_scale_bounds)) {
      length_scale_bounds <- c(med / 50, med * 5)
    }
    grid <- exp(seq(log(length_scale_bounds[1L]), log(length_scale_bounds[2L]),
      length.out = as.integer(n_grid)))
    lml_at <- function(ls) {
      train_gp_surrogate(X_train, Yc, length_scale = ls, signal_var = sv,
        noise_var = noise_var)$log_marginal_likelihood
    }
    lml_grid <- vapply(grid, lml_at, numeric(1L))
    bestls <- grid[which.max(lml_grid)]
    step <- grid[2L] / grid[1L]
    refine <- stats::optimize(function(log_ls) -lml_at(exp(log_ls)),
      interval = log(c(bestls / step, bestls * step)))
    ls_mle <- if (-refine$objective > max(lml_grid)) {
      exp(refine$minimum)
    } else {
      bestls
    }
    gp <- train_gp_surrogate(X_train, Yc, length_scale = ls_mle,
      signal_var = sv, noise_var = noise_var)
    input_scale <- rep(1, npar)
    length_scale <- ls_mle
  }

  gp$tuning <- list(
    anisotropic = isTRUE(anisotropic) && npar >= 2L,
    length_scale = length_scale,
    input_scale = input_scale,
    y_mean = y_mean,
    length_scale_median = med,
    log_marginal_likelihood = gp$log_marginal_likelihood,
    log_marginal_likelihood_median = lml_median
  )
  gp
}


#' Predict from an MLE-Tuned GP Surrogate
#'
#' Companion predictor for [train_gp_surrogate_tuned()]. It reapplies the
#' per-axis pre-scaling the tuner stored (so an anisotropic surrogate predicts on
#' the same geometry it was fitted on) and adds the centred response mean back,
#' returning predictions on the original response scale.
#'
#' @param gp A surrogate from [train_gp_surrogate_tuned()].
#' @param X_new Numeric matrix of inputs to predict at.
#'
#' @return The list [predict_gp_surrogate()] returns (`mean`, `variance`,
#'   `uncertainty`), with `mean` on the original response scale.
#'
#' @seealso [train_gp_surrogate_tuned()].
#'
#' @examples
#' set.seed(1L)
#' X <- matrix(runif(40L * 2L), 40L, 2L)
#' y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2
#' gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
#' pred <- predict_gp_surrogate_tuned(gp, X)
#' @export
predict_gp_surrogate_tuned <- function(gp, X_new) {
  if (is.null(gp$tuning)) {
    stop("`gp` is not a tuned surrogate (no `tuning`); use ",
      "`predict_gp_surrogate()` for a plain surrogate.", call. = FALSE)
  }
  X_new <- as.matrix(X_new)
  xs <- sweep(X_new, 2L, gp$tuning$input_scale, `/`)
  pred <- predict_gp_surrogate(gp, xs)
  pred$mean <- sweep(pred$mean, 2L, gp$tuning$y_mean, `+`)
  pred
}
