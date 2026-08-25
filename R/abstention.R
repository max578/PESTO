# abstention.R -- the typed "I cannot answer" result.
#
# PESTO's reliability gates (a degenerate ensemble, a structurally
# over-determined problem, a surrogate run outside its favourable training
# regime) previously either raised a bare `stop()` or -- for the surrogate
# regime -- only warned and proceeded. Neither shape is visible to the
# orchestra's cross-member refusal contract, which recognises a decline by
# class-name convention (`is_orchestra_decline()` in
# `ORCHESTRA_dev/integration/refusal_contract.R`: any class ending in
# `_refusal` or `_abstention`). `pesto_abstention()` gives PESTO's declines
# that shape without changing which conditions are declined.

#' Construct a typed PESTO abstention
#'
#' Returns a small classed object signalling that a PESTO reliability gate
#' declined to produce a result, with a machine-readable `reason` a caller
#' (or the orchestra's `is_orchestra_decline()`) can branch on, instead of an
#' uncaught error or a half-populated result.
#'
#' @param reason Character scalar reason code naming the declined
#'   condition, e.g. `"degenerate_ensemble"` (fewer than two realisations
#'   survived a step and the update cannot be formed, raised by
#'   [pesto_ies_callback()] / [pesto_ies_filter()]) or
#'   `"surrogate_off_design"` (the GP surrogate's
#'   training-point-to-parameter ratio falls below
#'   [check_surrogate_regime()]'s favourable floor, raised by
#'   [pesto_surrogate_ies()]). Not a closed enum: callers may mint further
#'   reason codes (e.g. `"over_determined"` for a future gate) as long as
#'   the class stays `pesto_abstention`.
#' @param detail Character scalar human-readable detail, or `NA`.
#' @param diagnostics Named list of scalar diagnostics supporting the
#'   decline (e.g. `list(iteration = 3L, n_ok = 1L, nreal = 40L)`). Default
#'   an empty list.
#' @param scope Character scalar naming what was refused (default `"call"`).
#'
#' @returns An object of class `"pesto_abstention"`: a list with `reason`,
#'   `detail`, `diagnostics`, `scope` and `abstained = TRUE`.
#'
#' @examples
#' a <- pesto_abstention("degenerate_ensemble",
#'                       detail = "1 of 40 realisations succeeded",
#'                       diagnostics = list(iteration = 2L, n_ok = 1L))
#' is_pesto_abstention(a)
#'
#' @export
pesto_abstention <- function(reason, detail = NA_character_,
                             diagnostics = list(), scope = "call") {
  stopifnot(is.character(reason), length(reason) == 1L, nzchar(reason))
  stopifnot(is.list(diagnostics))
  structure(
    list(reason = reason, detail = as.character(detail),
         diagnostics = diagnostics, scope = as.character(scope),
         abstained = TRUE),
    class = "pesto_abstention"
  )
}

#' Is an object a PESTO abstention?
#'
#' @param x Any object.
#' @returns A single logical: `TRUE` when `x` is a `pesto_abstention`.
#' @examples
#' is_pesto_abstention(pesto_abstention("over_determined"))
#' is_pesto_abstention(42)
#' @export
is_pesto_abstention <- function(x) {
  inherits(x, "pesto_abstention")
}

#' @export
print.pesto_abstention <- function(x, ...) {
  cat("<pesto_abstention>\n")
  cat(sprintf("  scope:  %s\n", x$scope))
  cat(sprintf("  reason: %s\n", x$reason))
  if (!is.na(x$detail) && nzchar(x$detail)) {
    cat(sprintf("  detail: %s\n", x$detail))
  }
  if (length(x$diagnostics) > 0L) {
    for (nm in names(x$diagnostics)) {
      cat(sprintf("  %s: %s\n", nm, format(x$diagnostics[[nm]])))
    }
  }
  invisible(x)
}
