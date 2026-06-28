# Deprecated functions -------------------------------------------------------

#' Deprecated: ensemble solution with adaptive SVD backend
#'
#' `ensemble_solution_gpu()` was renamed to [ensemble_solution_adaptive()].
#' The original name implied GPU computation; the function performs CPU
#' adaptive-SVD backend selection (randomised SVD versus a dense LAPACK /
#' Accelerate decomposition). The alias is retained for backward compatibility
#' and will be removed in a future release.
#'
#' @param ... Arguments passed on to [ensemble_solution_adaptive()].
#' @return The list returned by [ensemble_solution_adaptive()].
#' @seealso [ensemble_solution_adaptive()]
#' @examples
#' set.seed(1L)
#' np <- 4L; no <- 30L; nr <- 20L
#' pd <- matrix(rnorm(np * nr), np, nr)
#' od <- matrix(rnorm(no * nr), no, nr)
#' or_ <- matrix(rnorm(no * nr, sd = 0.5), no, nr)
#' Am <- matrix(0, 0, 0)
#' suppressWarnings(
#'   res <- ensemble_solution_gpu(pd, od, or_, pd, rep(1, no), rep(1, np),
#'                                Am, cur_lam = 1.0)
#' )
#' dim(res$upgrade)
#' @export
ensemble_solution_gpu <- function(...) {
  .Deprecated("ensemble_solution_adaptive", package = "PESTO")
  ensemble_solution_adaptive(...)
}
