#' Check whether a surrogate-IES regime is favourable
#'
#' Issues a warning when the ratio of training points to parameters falls
#' below an empirical threshold, where the Gaussian-process surrogate
#' inside [pesto_surrogate_ies()] and [surrogate_ensemble_update()] is
#' unlikely to repay its training cost. The check is a soft guardrail —
#' it does not modify the run, only flags an unfavourable regime so the
#' caller can decide whether to fall back to pure IES.
#'
#' The default threshold of `5` corresponds to the soft floor
#' `n_train >= 5 * n_params` documented in
#' `vignette("surrogate-ies", package = "PESTO")`. Below that floor the
#' GP posterior variance typically stays above the uncertainty-driven
#' switching threshold and surrogate savings collapse to near zero.
#'
#' This is exposed as a stand-alone helper so users can call it
#' explicitly before scheduling an expensive ensemble. It is **not**
#' invoked automatically by [pesto_surrogate_ies()] in the current
#' release; that wiring is tracked as a v0.2 enhancement candidate.
#'
#' @param n_params Integer. Number of estimated parameters.
#' @param n_train Integer. Number of training samples available to the
#'   surrogate (typically the ensemble size).
#' @param threshold Numeric. Minimum acceptable `n_train / n_params`
#'   ratio. Default `5`, the empirical soft floor from the surrogate-IES
#'   vignette.
#'
#' @return Invisibly returns `TRUE` when the regime is favourable
#'   (`n_train >= threshold * n_params`) and `FALSE` otherwise.
#'   Called for the warning side-effect.
#'
#' @seealso [pesto_surrogate_ies()], [surrogate_ensemble_update()],
#'   `vignette("surrogate-ies", package = "PESTO")`
#'
#' @references
#' Rasmussen, C. E. & Williams, C. K. I. (2006). *Gaussian Processes for
#' Machine Learning*. MIT Press.
#'
#' @examples
#' # Favourable regime: 100 training points for 10 parameters.
#' check_surrogate_regime(n_params = 10L, n_train = 100L)
#'
#' # Unfavourable regime: 30 training points for 30 parameters
#' # (the curse-of-dimensionality case from Scenario C of the
#' # comparison-and-simulation vignette). Emits a warning.
#' suppressWarnings(
#'   check_surrogate_regime(n_params = 30L, n_train = 30L)
#' )
#'
#' # Custom threshold for users with a smoother forward model.
#' check_surrogate_regime(n_params = 20L, n_train = 60L, threshold = 3)
#'
#' @export
check_surrogate_regime <- function(n_params, n_train, threshold = 5) {

  if (length(n_params) != 1L || !is.numeric(n_params) ||
      !is.finite(n_params) || n_params <= 0) {
    stop("`n_params` must be a single positive number.", call. = FALSE)
  }
  if (length(n_train) != 1L || !is.numeric(n_train) ||
      !is.finite(n_train) || n_train <= 0) {
    stop("`n_train` must be a single positive number.", call. = FALSE)
  }
  if (length(threshold) != 1L || !is.numeric(threshold) ||
      !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive number.", call. = FALSE)
  }

  ratio_required <- threshold * n_params
  favourable <- n_train >= ratio_required

  if (!favourable) {
    warning(
      sprintf(
        paste0(
          "Surrogate regime is unfavourable: n_train = %g for n_params = %g ",
          "(ratio %.2f < threshold %g). The Gaussian-process surrogate is ",
          "likely to yield near-zero model-call savings in this regime; ",
          "consider running pure IES instead. See ?pesto_surrogate_ies and ",
          "vignette(\"surrogate-ies\", package = \"PESTO\") for guidance."
        ),
        n_train, n_params, n_train / n_params, threshold
      ),
      call. = FALSE
    )
  }

  invisible(favourable)
}
