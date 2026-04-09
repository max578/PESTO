#' Read an Ensemble File
#'
#' Reads PEST++ ensemble files in CSV or binary (.jcb/.jco) format.
#'
#' @param file Character. Path to ensemble file.
#' @param format Character. One of "csv" (default) or "binary".
#' @return A data.table with realisations as rows and parameters/observations as columns.
#'   The first column `real_name` contains realisation names.
#' @export
#' @examples
#' # ens <- read_ensemble("prior.csv")
#' # print(ens)
read_ensemble <- function(file, format = c("csv", "binary")) {

  format <- match.arg(format)
  if (!file.exists(file)) {
    stop("Ensemble file not found: ", file, call. = FALSE)
  }

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
#' @param format Character. Currently only "csv" is supported.
#' @return Invisible `NULL`.
#' @export
write_ensemble <- function(ensemble, file, format = "csv") {
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

  dt <- data.table::data.table(real_name = row_names)
  for (j in seq_len(ncol_val)) {
    data.table::set(dt, j = col_names[j], value = mat[, j])
  }

  dt
}
