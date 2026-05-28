#' Read a PEST Control File (.pst)
#'
#' Parses a PEST/PEST++ control file and returns a structured list
#' containing all sections: control data, parameter groups, parameter
#' data, observation groups, observation data, model command, and
#' template/instruction file pairs.
#'
#' @param file Character. Path to the .pst file.
#' @return A list of class `pesto_pst` containing:
#'   \describe{
#'     \item{control_data}{Control section parameters (NPAR, NOBS, etc.)}
#'     \item{parameter_groups}{data.table of parameter group definitions}
#'     \item{parameters}{data.table of parameter data}
#'     \item{observation_groups}{data.table of observation group definitions}
#'     \item{observations}{data.table of observation data}
#'     \item{model_command}{Character vector of model command lines}
#'     \item{template_files}{data.table of template/model input file pairs}
#'     \item{instruction_files}{data.table of instruction/model output file pairs}
#'     \item{prior_information}{data.table of prior information (if present)}
#'     \item{pestpp_options}{Named list of ++ options}
#'   }
#' @examples
#' pars <- data.table::data.table(
#'   parnme = c("k1", "k2", "k3"),
#'   partrans = c("log", "log", "none"),
#'   parchglim = "factor",
#'   parval1 = c(1.0, 0.5, 0.1),
#'   parlbnd = c(0.01, 0.001, 0.0),
#'   parubnd = c(100, 50, 1.0),
#'   pargp = c("hydraulic", "hydraulic", "storage")
#' )
#' obs <- data.table::data.table(
#'   obsnme = c("h1", "h2", "h3"),
#'   obsval = c(1.0, 2.0, 1.5),
#'   weight = c(1.0, 1.0, 0.5),
#'   obgnme = "head"
#' )
#' pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#' tf <- tempfile(fileext = ".pst")
#' on.exit(unlink(tf), add = TRUE)
#' write_pst(pst, tf)
#' pst_back <- read_pst(tf)
#' pst_back$control_data$npar
#' pst_back$control_data$nobs
#' @seealso [write_pst()], [create_pest_scenario()]
#' @export
read_pst <- function(file) {
  .assert_path_exists(file, "file")

  # Read and locate section markers --------------------------------------
  lines <- readLines(file, warn = FALSE)
  lines <- trimws(lines)

  sections <- c(
    "pcf", "* control data", "* parameter groups",
    "* parameter data", "* observation groups",
    "* observation data", "* model command line",
    "* model input/output", "* prior information",
    "* regularisation", "* regularization"
  )

  section_idx <- list()
  for (i in seq_along(lines)) {
    ll <- tolower(lines[i])
    for (s in sections) {
      if (grepl(paste0("^\\*?\\s*", gsub("\\*", "\\\\*", s)), ll)) {
        section_idx[[s]] <- i
      }
    }
  }

  result <- list()

  # Parse control data ---------------------------------------------------
  cd_start <- section_idx[["* control data"]]
  if (!is.null(cd_start)) {
    cd_lines <- lines[(cd_start + 1):(cd_start + 4)]
    # Line 1: RSTFLE PESTMODE
    tokens1 <- strsplit(cd_lines[1], "\\s+")[[1]]
    # Line 2: NPAR NOBS NPARGP NPRIOR NOBSGP [MAXCOMPDIM]
    tokens2 <- as.integer(strsplit(cd_lines[2], "\\s+")[[1]])
    # Line 3: NTPLFLE NINSFLE PRECIS DPOINT [NUMCOM JACFILE MESSFILE] [OBSREREF]
    tokens3 <- strsplit(cd_lines[3], "\\s+")[[1]]
    # Line 4: RLAMBDA1 RLAMFAC PHIRATSUF PHIREDLAM NUMLAM [JACUPDATE] [LAMFORGIVE] [DERFORGIVE]
    tokens4 <- strsplit(cd_lines[4], "\\s+")[[1]]

    result$control_data <- list(
      rstfle   = if (length(tokens1) >= 1) tokens1[1] else NA_character_,
      pestmode = if (length(tokens1) >= 2) tokens1[2] else NA_character_,
      npar     = if (length(tokens2) >= 1) tokens2[1] else NA_integer_,
      nobs     = if (length(tokens2) >= 2) tokens2[2] else NA_integer_,
      npargp   = if (length(tokens2) >= 3) tokens2[3] else NA_integer_,
      nprior   = if (length(tokens2) >= 4) tokens2[4] else NA_integer_,
      nobsgp   = if (length(tokens2) >= 5) tokens2[5] else NA_integer_,
      ntplfle  = if (length(tokens3) >= 1) as.integer(tokens3[1]) else NA_integer_,
      ninsfle  = if (length(tokens3) >= 2) as.integer(tokens3[2]) else NA_integer_
    )
  }

  # Parse parameter data -------------------------------------------------
  pd_start <- section_idx[["* parameter data"]]
  og_start <- section_idx[["* observation groups"]]
  if (!is.null(pd_start) && !is.null(og_start)) {
    pd_lines <- lines[(pd_start + 1):(og_start - 1)]
    pd_lines <- pd_lines[nchar(pd_lines) > 0]

    if (length(pd_lines) > 0) {
      par_list <- lapply(pd_lines, function(l) {
        tokens <- strsplit(l, "\\s+")[[1]]
        if (length(tokens) >= 7) {
          data.table::data.table(
            parnme    = tokens[1],
            partrans  = tokens[2],
            parchglim = tokens[3],
            parval1   = as.numeric(tokens[4]),
            parlbnd   = as.numeric(tokens[5]),
            parubnd   = as.numeric(tokens[6]),
            pargp     = tokens[7],
            scale     = if (length(tokens) >= 8) as.numeric(tokens[8]) else 1.0,
            offset    = if (length(tokens) >= 9) as.numeric(tokens[9]) else 0.0,
            dercom    = if (length(tokens) >= 10) as.integer(tokens[10]) else 1L
          )
        }
      })
      par_list <- Filter(Negate(is.null), par_list)
      if (length(par_list) > 0) {
        result$parameters <- data.table::rbindlist(par_list)
      }
    }
  }

  # Parse observation data -----------------------------------------------
  od_start <- section_idx[["* observation data"]]
  mc_start <- section_idx[["* model command line"]]
  if (!is.null(od_start) && !is.null(mc_start)) {
    od_lines <- lines[(od_start + 1):(mc_start - 1)]
    od_lines <- od_lines[nchar(od_lines) > 0]

    if (length(od_lines) > 0) {
      obs_list <- lapply(od_lines, function(l) {
        tokens <- strsplit(l, "\\s+")[[1]]
        if (length(tokens) >= 3) {
          data.table::data.table(
            obsnme  = tokens[1],
            obsval  = as.numeric(tokens[2]),
            weight  = as.numeric(tokens[3]),
            obgnme  = if (length(tokens) >= 4) tokens[4] else NA_character_
          )
        }
      })
      obs_list <- Filter(Negate(is.null), obs_list)
      if (length(obs_list) > 0) {
        result$observations <- data.table::rbindlist(obs_list)
      }
    }
  }

  # Parse model command --------------------------------------------------
  if (!is.null(mc_start)) {
    mio_start <- section_idx[["* model input/output"]]
    if (!is.null(mio_start)) {
      mc_lines <- lines[(mc_start + 1):(mio_start - 1)]
      mc_lines <- mc_lines[nchar(mc_lines) > 0]
      result$model_command <- mc_lines
    }
  }

  # Parse template/instruction files -------------------------------------
  mio_start <- section_idx[["* model input/output"]]
  if (!is.null(mio_start)) {
    # Find next section or end of file
    next_section <- length(lines)
    for (s in names(section_idx)) {
      if (section_idx[[s]] > mio_start) {
        next_section <- min(next_section, section_idx[[s]] - 1)
      }
    }
    mio_lines <- lines[(mio_start + 1):next_section]
    mio_lines <- mio_lines[nchar(mio_lines) > 0 & !grepl("^\\+\\+", mio_lines)]

    if (length(mio_lines) > 0) {
      tpl_list <- lapply(mio_lines, function(l) {
        tokens <- strsplit(l, "\\s+")[[1]]
        if (length(tokens) >= 2) {
          data.table::data.table(
            template_file = tokens[1],
            model_file    = tokens[2]
          )
        }
      })
      tpl_list <- Filter(Negate(is.null), tpl_list)
      if (length(tpl_list) > 0) {
        result$io_files <- data.table::rbindlist(tpl_list)
      }
    }
  }

  # Parse PEST++ ++ options ----------------------------------------------
  pp_lines <- lines[grepl("^\\+\\+", lines)]
  if (length(pp_lines) > 0) {
    pp_opts <- list()
    for (l in pp_lines) {
      l <- sub("^\\+\\+\\s*", "", l)
      parts <- strsplit(l, "\\s*[=()]\\s*")[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        val <- trimws(parts[2])
        # Try numeric
        num_val <- suppressWarnings(as.numeric(val))
        if (!is.na(num_val)) {
          pp_opts[[key]] <- num_val
        } else if (tolower(val) %in% c("true", "false")) {
          pp_opts[[key]] <- tolower(val) == "true"
        } else {
          pp_opts[[key]] <- val
        }
      }
    }
    result$pestpp_options <- pp_opts
  }

  result$file <- normalizePath(file, mustWork = FALSE)
  class(result) <- "pesto_pst"
  result
}

#' Write a PEST Control File (.pst)
#'
#' Writes a `pesto_pst` object to a PEST-format control file.
#'
#' @param pst A `pesto_pst` object (as returned by [read_pst()]).
#' @param file Character. Output file path.
#' @return Invisible `NULL`. File is written as a side effect.
#' @examples
#' pars <- data.table::data.table(
#'   parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'   parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'   parubnd = c(100, 50), pargp = "hydraulic"
#' )
#' obs <- data.table::data.table(
#'   obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'   weight = c(1.0, 1.0), obgnme = "head"
#' )
#' pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#' tf <- tempfile(fileext = ".pst")
#' on.exit(unlink(tf), add = TRUE)
#' write_pst(pst, tf)
#' file.exists(tf)
#' @seealso [read_pst()]
#' @export
write_pst <- function(pst, file) {
  # Validate inputs ------------------------------------------------------
  if (!inherits(pst, "pesto_pst")) {
    stop("`pst` must be a `pesto_pst` object.", call. = FALSE)
  }
  .assert_character_scalar(file, "file")

  lines <- character()
  lines <- c(lines, "pcf")

  # Control data ---------------------------------------------------------
  cd <- pst$control_data
  lines <- c(lines, "* control data")
  lines <- c(lines, paste(cd$rstfle, cd$pestmode))
  lines <- c(lines, paste(cd$npar, cd$nobs, cd$npargp, cd$nprior, cd$nobsgp))
  lines <- c(lines, paste(cd$ntplfle, cd$ninsfle, "double point"))
  lines <- c(lines, "5.0 2.0 0.3 0.03 10")
  lines <- c(lines, "5.0 5.0 0.001")
  lines <- c(lines, "0.1")
  lines <- c(lines, "30 0.01 4 4 0.01 4")
  lines <- c(lines, "1 1 1")

  # Parameter groups -----------------------------------------------------
  if (!is.null(pst$parameter_groups)) {
    lines <- c(lines, "* parameter groups")
    for (i in seq_len(nrow(pst$parameter_groups))) {
      row <- pst$parameter_groups[i, ]
      lines <- c(lines, paste(row, collapse = " "))
    }
  }

  # Parameter data -------------------------------------------------------
  if (!is.null(pst$parameters)) {
    lines <- c(lines, "* parameter data")
    for (i in seq_len(nrow(pst$parameters))) {
      p <- pst$parameters[i, ]
      lines <- c(lines, paste(
        p$parnme, p$partrans, p$parchglim,
        format(p$parval1, scientific = FALSE),
        format(p$parlbnd, scientific = FALSE),
        format(p$parubnd, scientific = FALSE),
        p$pargp, p$scale, p$offset, p$dercom
      ))
    }
  }

  # Observation groups ---------------------------------------------------
  if (!is.null(pst$observation_groups)) {
    lines <- c(lines, "* observation groups")
    for (i in seq_len(nrow(pst$observation_groups))) {
      lines <- c(lines, paste(pst$observation_groups[i, ], collapse = " "))
    }
  }

  # Observation data -----------------------------------------------------
  if (!is.null(pst$observations)) {
    lines <- c(lines, "* observation data")
    for (i in seq_len(nrow(pst$observations))) {
      o <- pst$observations[i, ]
      lines <- c(lines, paste(o$obsnme, o$obsval, o$weight, o$obgnme))
    }
  }

  # Model command --------------------------------------------------------
  lines <- c(lines, "* model command line")
  lines <- c(lines, pst$model_command)

  # Template / instruction file pairs ------------------------------------
  lines <- c(lines, "* model input/output")
  if (!is.null(pst$io_files)) {
    for (i in seq_len(nrow(pst$io_files))) {
      lines <- c(lines, paste(pst$io_files[i, ], collapse = " "))
    }
  }

  # PEST++ ++ options ----------------------------------------------------
  if (!is.null(pst$pestpp_options)) {
    for (nm in names(pst$pestpp_options)) {
      lines <- c(lines, paste0("++", nm, "(", pst$pestpp_options[[nm]], ")"))
    }
  }

  writeLines(lines, file)
  invisible(NULL)
}

#' Print method for pesto_pst objects
#' @param x A `pesto_pst` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`. Called for the side effect of printing.
#' @examples
#' pars <- data.table::data.table(
#'   parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'   parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'   parubnd = c(100, 50), pargp = "hydraulic"
#' )
#' obs <- data.table::data.table(
#'   obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'   weight = c(1.0, 1.0), obgnme = "head"
#' )
#' pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#' print(pst)
#' @export
print.pesto_pst <- function(x, ...) {
  cat("PESTO PST Control File\n")
  cat("---------------------\n")
  if (!is.null(x$control_data)) {
    cat(sprintf("Parameters:   %d\n", x$control_data$npar))
    cat(sprintf("Observations: %d\n", x$control_data$nobs))
    cat(sprintf("Par groups:   %d\n", x$control_data$npargp))
    cat(sprintf("Obs groups:   %d\n", x$control_data$nobsgp))
  }
  if (!is.null(x$pestpp_options)) {
    cat(sprintf("PEST++ options: %d\n", length(x$pestpp_options)))
  }
  cat(sprintf("Source: %s\n", x$file))
  invisible(x)
}
