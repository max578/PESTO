# multifidelity.R -- First-class multi-fidelity forward-model bridge.
#
# Bridge-invariant #3 of the APSIM roadmap: "the bridge layer accepts a
# fidelity vector (cheap, expensive) so future surrogate cascades plug
# in without refactor." This file makes that vector first-class. A
# `pesto_multifidelity_model` is an ordered stack of
# `pesto_forward_model` levels (cheapest first) plus their relative
# costs. The IES driver consumes it through a `fidelity_schedule`
# (which level each iteration runs at); the affine control-variate
# combiner [mf_control_variate()] is the statistical primitive a
# surrogate cascade plugs into to debias a cheap level against a sparse
# expensive sample.

#' Multi-Fidelity Forward Model (S7 class)
#'
#' Bundles an ordered stack of [pesto_forward_model()] levels --
#' cheapest (`level = 0`) to most expensive (`level = n - 1`) -- with
#' their relative per-evaluation costs. This is the first-class form of
#' the bridge's `(cheap, expensive)` fidelity vector: the IES driver
#' [pesto_ies_callback()] selects a level per iteration via
#' `fidelity_schedule`, and [mf_control_variate()] debiases a cheap
#' level against a sparse expensive sample for surrogate cascades.
#'
#' Each element of `levels` may be a bare `function(theta) -> obs` or a
#' fully-specified [pesto_forward_model()]; bare functions are coerced
#' via [as_forward_model()]. Levels must be ordered by ascending
#' fidelity (cheapest first) -- the convention the `level` index and the
#' `costs` vector both follow.
#'
#' @param levels List of [pesto_forward_model()] objects (or bare
#'   functions, which are coerced), ordered cheapest first.
#' @param costs Numeric vector of relative per-evaluation costs, one per
#'   level, ascending by convention. Defaults to `seq_along(levels)`.
#'   Carried for cost-aware allocation; not yet used to schedule
#'   automatically (that is the documented extension point).
#' @param label Character. Optional human label.
#' @return A `pesto_multifidelity_model` S7 object.
#' @seealso [pesto_forward_model()], [pesto_evaluate()],
#'   [mf_control_variate()], [pesto_ies_callback()].
#' @examples
#' cheap     <- function(theta) theta %*% c(1, 1)        # fast, biased
#' expensive <- function(theta) theta %*% c(1, 1) + 0.5  # slow, truth
#' mf <- pesto_multifidelity_model(
#'   levels = list(
#'     pesto_forward_model(fn = cheap,     n_obs = 1L, fidelity = 0L),
#'     pesto_forward_model(fn = expensive, n_obs = 1L, fidelity = 1L)
#'   ),
#'   costs = c(1, 25)
#' )
#' theta <- matrix(c(1, 0, 0, 1), nrow = 2L, byrow = TRUE)
#' pesto_evaluate(mf, theta, level = 0L)  # cheap
#' pesto_evaluate(mf, theta, level = 1L)  # expensive
#' @export
pesto_multifidelity_model <- S7::new_class(
  "pesto_multifidelity_model",
  package = "PESTO",
  properties = list(
    levels = S7::class_list,
    costs  = S7::class_numeric,
    label  = S7::new_property(S7::class_character, default = NA_character_)
  ),
  constructor = function(levels, costs = NULL, label = NA_character_) {
    levels <- lapply(levels, as_forward_model)
    if (is.null(costs)) costs <- as.numeric(seq_along(levels))
    S7::new_object(
      S7::S7_object(),
      levels = levels,
      costs  = as.numeric(costs),
      label  = as.character(label)
    )
  },
  validator = function(self) {
    errs <- character(0)
    if (length(self@levels) < 1L) {
      errs <- c(errs, "`levels` must contain at least one fidelity level")
    }
    ok_levels <- vapply(
      self@levels,
      function(x) S7::S7_inherits(x, pesto_forward_model),
      logical(1L)
    )
    if (!all(ok_levels)) {
      errs <- c(errs,
                "every element of `levels` must be a `pesto_forward_model`")
    }
    if (length(self@costs) != length(self@levels)) {
      errs <- c(errs, "`costs` must have one entry per fidelity level")
    } else if (any(!is.finite(self@costs)) || any(self@costs <= 0)) {
      errs <- c(errs, "`costs` must be positive and finite")
    }
    if (length(errs) == 0L) NULL else paste(errs, collapse = "; ")
  }
)


S7::method(pesto_evaluate, pesto_multifidelity_model) <-
  function(model, theta, ..., level = NULL) {
    n_levels <- length(model@levels)
    if (is.null(level)) level <- n_levels - 1L
    .assert_fidelity_level(level, n_levels)
    pesto_evaluate(model@levels[[as.integer(level) + 1L]], theta)
  }


#' Affine control-variate bias correction across fidelities
#'
#' Debiases a cheap-fidelity output ensemble against a sparse expensive
#' sample, per observation dimension. For each observation `j` it fits
#' the first-order autoregressive control variate (the linear term of
#' the Kennedy-O'Hagan multi-fidelity model)
#' `high_j ~ a_j + b_j * low_j` on the paired subset, with
#' `b_j = cov(high_j, low_j) / var(low_j)` the variance-minimising
#' coefficient, then predicts the corrected high-fidelity output for
#' every realisation as `a_j + b_j * low_all_j`.
#'
#' The estimator degrades gracefully: where the cheap output has zero
#' variance on the subset it falls back to the expensive subset mean
#' (`b_j = 0`), and where the two fidelities are weakly correlated the
#' correction shrinks toward that mean rather than amplifying noise.
#'
#' This is the plug-in primitive for surrogate cascades: a cascade runs
#' the cheap level over the full ensemble, the expensive level over a
#' chosen subset, and calls this to lift the cheap ensemble toward the
#' expensive one at a fraction of the cost.
#'
#' @param low_all Numeric matrix, `nreal x nobs`. Cheap-fidelity output
#'   for every realisation.
#' @param high_sub Numeric matrix, `nsub x nobs`. Expensive-fidelity
#'   output for the paired subset.
#' @param low_sub Numeric matrix, `nsub x nobs`. Cheap-fidelity output
#'   for the same subset, row-aligned with `high_sub`.
#' @return A `nreal x nobs` matrix of bias-corrected outputs, with
#'   attributes `"intercept"` (`a_j`), `"slope"` (`b_j`), and
#'   `"subset_cor"` (per-dimension subset correlation; `NA` for a
#'   degenerate dimension).
#' @references
#' Kennedy, M. C. & O'Hagan, A. (2000). Predicting the output from a
#' complex computer code when fast approximations are available.
#' *Biometrika*, 87(1), 1--13.
#' @seealso [pesto_multifidelity_model()].
#' @examples
#' set.seed(1L)
#' low_all  <- matrix(rnorm(40L), ncol = 2L)
#' sub      <- 1:5
#' low_sub  <- low_all[sub, , drop = FALSE]
#' high_sub <- 0.3 + 1.2 * low_sub + matrix(rnorm(10L, sd = 0.01), ncol = 2L)
#' corrected <- mf_control_variate(low_all, high_sub, low_sub)
#' attr(corrected, "slope")
#' @export
mf_control_variate <- function(low_all, high_sub, low_sub) {
  .assert_matrix(low_all, "low_all")
  .assert_matrix(high_sub, "high_sub")
  .assert_matrix(low_sub, "low_sub")
  .assert_same_ncol(low_all, high_sub, "low_all", "high_sub")
  .assert_same_nrow(high_sub, low_sub, "high_sub", "low_sub")

  nobs <- ncol(low_all)
  intercept  <- numeric(nobs)
  slope      <- numeric(nobs)
  subset_cor <- numeric(nobs)
  corrected  <- matrix(NA_real_, nrow = nrow(low_all), ncol = nobs)

  for (j in seq_len(nobs)) {
    lj <- low_sub[, j]
    hj <- high_sub[, j]
    var_l <- stats::var(lj)
    if (!is.finite(var_l) || var_l <= 0) {
      slope[j]      <- 0
      intercept[j]  <- mean(hj)
      subset_cor[j] <- NA_real_
    } else {
      slope[j]      <- stats::cov(hj, lj) / var_l
      intercept[j]  <- mean(hj) - slope[j] * mean(lj)
      subset_cor[j] <- stats::cor(hj, lj)
    }
    corrected[, j] <- intercept[j] + slope[j] * low_all[, j]
  }

  attr(corrected, "intercept")  <- intercept
  attr(corrected, "slope")      <- slope
  attr(corrected, "subset_cor") <- subset_cor
  corrected
}


# Resolve a fidelity schedule into a length-noptmax integer vector of
# 0-based level indices, padding the tail with the last entry. NULL
# means "highest fidelity every iteration".
.resolve_fidelity_schedule <- function(schedule, n_levels, noptmax) {
  top <- n_levels - 1L
  if (is.null(schedule)) return(rep(top, noptmax))
  s <- as.integer(schedule)
  if (length(s) == 0L || any(is.na(s)) || any(s < 0L) || any(s > top)) {
    stop(
      sprintf(
        "`fidelity_schedule` entries must be integers in [0, %d].", top
      ),
      call. = FALSE
    )
  }
  if (length(s) < noptmax) {
    s <- c(s, rep(s[length(s)], noptmax - length(s)))
  }
  s[seq_len(noptmax)]
}


# Validate a single 0-based fidelity level against the level count.
.assert_fidelity_level <- function(level, n_levels) {
  lv <- suppressWarnings(as.integer(level))
  if (length(lv) != 1L || is.na(lv) || lv < 0L || lv > n_levels - 1L) {
    stop(
      sprintf(
        "`level` must be a single integer in [0, %d].", n_levels - 1L
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
