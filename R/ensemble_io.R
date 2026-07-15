#' Read an Ensemble File
#'
#' Reads PEST++ ensemble files in CSV or binary (.jcb/.jco) format.
#'
#' @param file Character. Path to ensemble file.
#' @param format Character. One of "csv" (default) or "binary".
#' @return A data.table with realisations as rows and parameters/observations as columns.
#'   The first column `real_name` contains realisation names.
#' @examples
#' ens <- data.table::data.table(
#'   real_name = sprintf("real_%02d", 1:10),
#'   k1 = rnorm(10, mean = 1.0, sd = 0.2),
#'   k2 = rnorm(10, mean = 0.5, sd = 0.1),
#'   k3 = rnorm(10, mean = 2.0, sd = 0.3)
#' )
#' tf <- tempfile(fileext = ".csv")
#' on.exit(unlink(tf), add = TRUE)
#' write_ensemble(ens, tf)
#' ens_back <- read_ensemble(tf, format = "csv")
#' identical(names(ens_back), names(ens))
#' nrow(ens_back) == nrow(ens)
#' @seealso [write_ensemble()]
#' @export
read_ensemble <- function(file, format = c("csv", "binary")) {
  format <- match.arg(format)
  .assert_path_exists(file, "file")

  if (format == "csv") {
    dt <- data.table::fread(file, header = TRUE)
    # First column is typically the realisation name
    if (ncol(dt) > 1 && is.character(dt[[1]])) {
      data.table::setnames(dt, 1, "real_name")
    }
    return(dt)
  } else {
    # Binary format: use C++ reader
    return(.read_ensemble_binary(file))
  }
}

#' Write an Ensemble File
#'
#' Writes an ensemble data.table to CSV format compatible with PEST++.
#'
#' @param ensemble A data.table with realisation data.
#' @param file Character. Output file path.
#' @param format Character. Only `"csv"` is supported; any other value is an
#'   error. Note the asymmetry with [read_ensemble()], which also reads PEST++
#'   binary: writing that format is not implemented.
#' @return Invisible `NULL`.
#' @examples
#' ens <- data.table::data.table(
#'   real_name = sprintf("real_%02d", 1:5),
#'   k1 = runif(5, 0.1, 10),
#'   k2 = runif(5, 0.01, 1)
#' )
#' tf <- tempfile(fileext = ".csv")
#' on.exit(unlink(tf), add = TRUE)
#' write_ensemble(ens, tf)
#' file.exists(tf)
#' @seealso [read_ensemble()]
#' @export
write_ensemble <- function(ensemble, file, format = "csv") {
  # `csv` is the only implemented writer, so refuse anything else rather than
  # writing CSV under the name the caller asked for. `read_ensemble()` gates
  # its own `format` the same way.
  format <- match.arg(format, "csv")

  if (!data.table::is.data.table(ensemble)) {
    ensemble <- data.table::as.data.table(ensemble)
  }
  data.table::fwrite(ensemble, file)
  invisible(NULL)
}

#' Internal binary ensemble reader
#' @param file Path to binary file
#' @return data.table
#' @keywords internal
.read_ensemble_binary <- function(file) {
  # Read PEST++ binary format (JCB/JCO)
  con <- file(file, "rb")
  on.exit(close(con))

  # Read header
  n_neg <- readBin(con, "integer", n = 1, size = 4)
  nrow_val <- abs(n_neg)
  ncol_val <- readBin(con, "integer", n = 1, size = 4)

  # Read data matrix (column-major, doubles)
  mat <- matrix(
    readBin(con, "double", n = nrow_val * ncol_val, size = 8),
    nrow = nrow_val, ncol = ncol_val, byrow = FALSE
  )

  # Read row names
  row_names <- character(nrow_val)
  for (i in seq_len(nrow_val)) {
    nm_len <- readBin(con, "integer", n = 1, size = 4)
    row_names[i] <- trimws(readChar(con, nm_len))
  }

  # Read column names
  col_names <- character(ncol_val)
  for (i in seq_len(ncol_val)) {
    nm_len <- readBin(con, "integer", n = 1, size = 4)
    col_names[i] <- trimws(readChar(con, nm_len))
  }

  # A file written without column names reads back with blank labels, and
  # `set()` keys on the name: every blank column would land on the same ""
  # column, so all but the last would be silently discarded -- fewer columns
  # than the header declared, with no error. Name the blanks positionally so
  # each column survives. `plot_identifiability()` substitutes real labels from
  # a `pst` when one is supplied.
  blank <- !nzchar(col_names)
  if (any(blank)) {
    col_names[blank] <- paste0("p", seq_len(ncol_val)[blank])
  }

  dt <- data.table::data.table(real_name = row_names)
  for (j in seq_len(ncol_val)) {
    data.table::set(dt, j = col_names[j], value = mat[, j])
  }

  # Record which labels are placeholders rather than names the file carried, so
  # a caller holding the authoritative names can substitute them. Without this
  # the substitution is unavailable: `p1` is indistinguishable from a parameter
  # genuinely called `p1`.
  data.table::setattr(dt, "unnamed_columns", blank)

  dt
}
