# Default print / plot methods for the IES result objects shared by the
# callback and sequential-filter drivers (class `pesto_ies_result`).

#' Print and plot methods for PESTO ensemble-smoother results
#'
#' `print()` gives a one-screen summary of a fitted iterative ensemble-smoother
#' run (driver, problem dimensions, the phi-convergence trace, the spread-ESS
#' dispersion diagnostic and the failure rate). `plot()` draws the
#' objective-function (phi) convergence trace by delegating to [plot_phi()]. For
#' posterior parameter distributions (prior vs posterior) use [plot_ensemble()].
#'
#' @param x A `pesto_ies_result`, as returned by [pesto_ies_callback()] or
#'   [pesto_ies_filter()].
#' @param ... Further arguments. For `plot()` these are passed to [plot_phi()]
#'   (e.g. `log_scale`, `title`); `print()` ignores them.
#' @return `print()` returns `x` invisibly; `plot()` returns a `ggplot2` object.
#'
#' @srrstats {BS6.0} Default `print` method for the return object, summarising
#'   the driver, problem dimensions, the per-step phi trace, the spread-ESS
#'   dispersion diagnostic and the failure rate.
#' @srrstats {BS6.1} Default `plot` method for the return object: the phi
#'   convergence trace (delegates to [plot_phi()]).
#' @srrstats {BS6.3} Posterior distributional estimates (prior vs posterior
#'   parameter distributions) are plotted by [plot_ensemble()], demonstrated in
#'   the getting-started vignette.
#' @name pesto_ies_result-methods
#' @examples
#' set.seed(1)
#' G <- matrix(rnorm(18L), 6L, 3L)
#' prior <- matrix(rnorm(150L), 50L, 3L, dimnames = list(NULL, paste0("p", 1:3)))
#' obs <- stats::setNames(as.numeric(G %*% c(1, -0.5, 2)) + rnorm(6L, sd = 0.05),
#'                        paste0("o", 1:6))
#' fit <- pesto_ies_callback(function(t) t %*% t(G), prior, obs,
#'                           obs_sd = 0.05, noptmax = 4L, verbose = FALSE)
#' print(fit)
#' p <- plot(fit)        # phi convergence (a ggplot2 object)
NULL

#' @rdname pesto_ies_result-methods
#' @export
print.pesto_ies_result <- function(x, ...) {
  npar      <- ncol(x$par_ensemble) - 1L
  nobs      <- ncol(x$obs_ensemble) - 1L
  nreal     <- nrow(x$par_ensemble)
  is_filter <- inherits(x, "pesto_ies_filter_result")
  n_steps   <- length(x$iterations)
  unit      <- if (is_filter) "window(s)" else "iteration(s)"

  cat(sprintf("<pesto_ies_result>  driver: %s\n",
              if (is_filter) "sequential filter" else "callback (in-process)"))
  cat(sprintf("  %d parameter(s), %d observation(s), %d realisation(s)\n",
              npar, nobs, nreal))
  cat(sprintf("  %d %s, %d forward evaluation(s), %.3gs runtime\n",
              n_steps, unit, x$n_forward_evals, x$runtime_seconds))

  if (!is.null(x$phi) && nrow(x$phi) > 0L) {
    mphi <- tapply(x$phi$phi, x$phi$iteration, mean, na.rm = TRUE)
    cat(sprintf("  mean phi: %.4g -> %.4g over %d step(s)\n",
                mphi[1L], mphi[length(mphi)], length(mphi)))
  }

  ess <- vapply(x$iterations, function(it) {
    v <- it[["spread_ess_ratio"]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1L))
  fin <- ess[length(ess)]
  if (length(fin) == 1L && is.finite(fin)) {
    cat(sprintf("  final spread-ESS ratio: %.3f\n", fin))
  }
  cat(sprintf("  failure rate: %.1f%%\n", 100 * x$failure_rate))
  if (!is.null(x$fidelity)) {
    cat(sprintf("  multi-fidelity: %d level(s), final level %d\n",
                x$fidelity$n_levels, x$fidelity$final_level))
  }
  invisible(x)
}

#' @rdname pesto_ies_result-methods
#' @export
plot.pesto_ies_result <- function(x, ...) {
  if (is.null(x$phi) || nrow(x$phi) == 0L) {
    stop("`x` has no `$phi` trace to plot.", call. = FALSE)
  }
  # Summarise the long-format phi trace (iteration / realisation / phi) into the
  # per-iteration band that plot_phi() expects, using base aggregation so no
  # data.table non-standard evaluation is introduced.
  ph  <- as.data.frame(x$phi)
  spl <- split(ph$phi, ph$iteration)
  summ <- data.frame(
    iteration = as.integer(names(spl)),
    mean      = vapply(spl, mean, numeric(1L)),
    min       = vapply(spl, min, numeric(1L)),
    max       = vapply(spl, max, numeric(1L)),
    median    = vapply(spl, stats::median, numeric(1L))
  )
  plot_phi(summ, ...)
}
