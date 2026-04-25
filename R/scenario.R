#' Create a PEST Scenario Programmatically
#'
#' Builds a `pesto_pst` object from R data structures, without
#' requiring an existing .pst file.
#'
#' @param parameters data.table. Parameter definitions with columns:
#'   `parnme`, `partrans`, `parchglim`, `parval1`, `parlbnd`, `parubnd`, `pargp`.
#' @param observations data.table. Observation definitions with columns:
#'   `obsnme`, `obsval`, `weight`, `obgnme`.
#' @param model_command Character. Model command line(s).
#' @param template_files data.table. Template/model input file pairs
#'   with columns `template_file`, `model_file`.
#' @param instruction_files data.table. Instruction/model output file pairs
#'   with columns `instruction_file`, `model_file`.
#' @param pestpp_options Named list. PEST++ options.
#' @return A `pesto_pst` object.
#' @examples
#' pars <- data.table::data.table(
#'   parnme = c("k1", "k2"),
#'   partrans = "log",
#'   parchglim = "factor",
#'   parval1 = c(1.0, 0.5),
#'   parlbnd = c(0.01, 0.001),
#'   parubnd = c(100, 50),
#'   pargp = "hydraulic"
#' )
#' obs <- data.table::data.table(
#'   obsnme = c("h1", "h2"),
#'   obsval = c(1.0, 2.0),
#'   weight = c(1.0, 1.0),
#'   obgnme = "head"
#' )
#' pst <- create_pest_scenario(
#'   parameters    = pars,
#'   observations  = obs,
#'   model_command = "python model.py"
#' )
#' inherits(pst, "pesto_pst")
#' pst$control_data$npar
#' @export
create_pest_scenario <- function(parameters,
                                  observations,
                                  model_command,
                                  template_files = NULL,
                                  instruction_files = NULL,
                                  pestpp_options = list()) {

  if (!data.table::is.data.table(parameters)) {
    parameters <- data.table::as.data.table(parameters)
  } else {
    parameters <- data.table::copy(parameters)
  }
  if (!data.table::is.data.table(observations)) {
    observations <- data.table::as.data.table(observations)
  } else {
    observations <- data.table::copy(observations)
  }

  # Validate required columns
  req_par <- c("parnme", "partrans", "parchglim", "parval1", "parlbnd", "parubnd", "pargp")
  missing_par <- setdiff(req_par, names(parameters))
  if (length(missing_par) > 0) {
    stop("Missing parameter columns: ", paste(missing_par, collapse = ", "), call. = FALSE)
  }

  req_obs <- c("obsnme", "obsval", "weight", "obgnme")
  missing_obs <- setdiff(req_obs, names(observations))
  if (length(missing_obs) > 0) {
    stop("Missing observation columns: ", paste(missing_obs, collapse = ", "), call. = FALSE)
  }

  # Add defaults
  if (!"scale" %in% names(parameters)) parameters[, scale := 1.0]
  if (!"offset" %in% names(parameters)) parameters[, offset := 0.0]
  if (!"dercom" %in% names(parameters)) parameters[, dercom := 1L]

  npar <- nrow(parameters)
  nobs <- nrow(observations)
  npargp <- length(unique(parameters$pargp))
  nobsgp <- length(unique(observations$obgnme))
  ntpl <- if (!is.null(template_files)) nrow(template_files) else 0L
  nins <- if (!is.null(instruction_files)) nrow(instruction_files) else 0L

  pst <- list(
    control_data = list(
      rstfle   = "restart",
      pestmode = "estimation",
      npar     = npar,
      nobs     = nobs,
      npargp   = npargp,
      nprior   = 0L,
      nobsgp   = nobsgp,
      ntplfle  = ntpl,
      ninsfle  = nins
    ),
    parameters = parameters,
    observations = observations,
    model_command = model_command,
    io_files = if (!is.null(template_files)) template_files else data.table::data.table(),
    pestpp_options = pestpp_options,
    file = NA_character_
  )

  class(pst) <- "pesto_pst"
  pst
}
