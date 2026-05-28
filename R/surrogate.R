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
