# Test helpers for the rOpenSci statistical-software (srr) testing battery.
#
# `skip_if_not_extended()` gates the long-running members of the battery
# (many-seed parameter recovery, finer data-size scaling) behind the
# `PESTO_EXTENDED_TESTS` environment variable, so the regular suite stays fast
# on CRAN while the full statistical battery runs under one flag in CI.
# The conditions for running the extended tests are documented in
# tests/README.md.

#' @srrstats {G5.10} Extended tests run under the same testthat framework as the
#'   regular suite but are switched on by the `PESTO_EXTENDED_TESTS` environment
#'   variable via `skip_if_not_extended()`. CI enables them by setting
#'   `PESTO_EXTENDED_TESTS=true` in the workflow `env` block.
#' @srrstats {G5.12} The conditions necessary to run the extended tests --- the
#'   `PESTO_EXTENDED_TESTS` flag, optional-package requirements (apsimx,
#'   deSolve, yaml), the multicore requirement, and approximate runtimes --- are
#'   documented in tests/README.md.
#' @noRd
NULL

.pesto_extended_tests <- function() {
  flag <- Sys.getenv("PESTO_EXTENDED_TESTS", unset = "false")
  isTRUE(tolower(flag) %in% c("true", "1", "yes", "on"))
}

skip_if_not_extended <- function() {
  if (!.pesto_extended_tests()) {
    testthat::skip("extended srr test; set PESTO_EXTENDED_TESTS=true to run")
  }
  invisible(TRUE)
}
