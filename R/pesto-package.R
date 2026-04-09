#' @keywords internal
"_PACKAGE"

#' @useDynLib PESTO, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 .data
#' @importFrom data.table data.table as.data.table fread fwrite setDT
#'   setnames is.data.table rbindlist melt dcast set `:=` .SD .I .N
#' @importFrom stats var median
NULL

# Suppress R CMD check notes for data.table / ggplot2 non-standard evaluation
utils::globalVariables(c(
  ".", ".I", ".SD", ".N", ".data",
  "iteration", "phi", "realisation", "parameter", "value", "source",
  "identifiability", "model_runs", "surr_runs", "savings_pct",
  "offset", "dercom", "scale", "real_name"
))
