# internal_validation.R -- Shared argument validators.
#
# Internal helpers used by exported functions to delegate validation
# out of function bodies. Each helper signals failure via
# `stop(call. = FALSE, ...)` with a backticked argument name so the
# caller sees the offending argument exactly as it was named at the
# call site. None of these are exported; they are loaded into every
# function body via `.assert_*()` and `.check_*()` calls.

#' Assert that a value is a single finite positive number
#'
#' @param x Value to check.
#' @param name Character. Argument name used in the error message.
#'
#' @return Invisibly returns `TRUE` when the check passes.
#'
#' @noRd
#' @keywords internal
.assert_positive_scalar <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive number.", name),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a value is a single non-negative finite number
#'
#' @noRd
#' @keywords internal
.assert_nonneg_scalar <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be a single non-negative number.", name),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a value is a single non-empty character string
#'
#' @noRd
#' @keywords internal
.assert_character_scalar <- function(x, name) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    stop(
      sprintf("`%s` must be a single non-empty character string.", name),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a value is a single logical (TRUE / FALSE; not NA)
#'
#' @noRd
#' @keywords internal
.assert_logical_scalar <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a path exists on disk
#'
#' @noRd
#' @keywords internal
.assert_path_exists <- function(path, name) {
  .assert_character_scalar(path, name)
  if (!file.exists(path)) {
    stop(sprintf("`%s` does not exist: %s", name, path), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a value is a numeric matrix
#'
#' @noRd
#' @keywords internal
.assert_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a value is a non-empty numeric vector
#'
#' @noRd
#' @keywords internal
.assert_numeric_vector <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop(sprintf("`%s` must be a non-empty numeric vector.", name),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a value is a function
#'
#' @noRd
#' @keywords internal
.assert_function <- function(x, name) {
  if (!is.function(x)) {
    stop(sprintf("`%s` must be a function.", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert that a value is a data.frame or data.table
#'
#' @noRd
#' @keywords internal
.assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame or data.table.", name),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a value is one of an allowed set
#'
#' @noRd
#' @keywords internal
.assert_choice <- function(x, name, choices) {
  if (length(x) != 1L || !(x %in% choices)) {
    stop(
      sprintf(
        "`%s` must be one of %s; got %s.",
        name,
        paste0("\"", choices, "\"", collapse = ", "),
        deparse(x)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that two matrices have the same number of columns
#'
#' @noRd
#' @keywords internal
.assert_same_ncol <- function(a, b, name_a, name_b) {
  if (ncol(a) != ncol(b)) {
    stop(
      sprintf(
        "`%s` and `%s` must have the same number of columns; got %d and %d.",
        name_a, name_b, ncol(a), ncol(b)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that two matrices have the same number of rows
#'
#' @noRd
#' @keywords internal
.assert_same_nrow <- function(a, b, name_a, name_b) {
  if (nrow(a) != nrow(b)) {
    stop(
      sprintf(
        "`%s` and `%s` must have the same number of rows; got %d and %d.",
        name_a, name_b, nrow(a), nrow(b)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assert that a data.frame contains every required column
#'
#' @noRd
#' @keywords internal
.assert_required_cols <- function(x, cols, name) {
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "`%s` is missing required column(s): %s.",
        name, paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
