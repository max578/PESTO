#' Run PEST++ IES (Iterative Ensemble Smoother)
#'
#' Executes the pestpp-ies algorithm on a PEST control file.
#' This is the primary method for ensemble-based parameter estimation
#' and uncertainty quantification.
#'
#' @param pst_file Character. Path to the .pst control file.
#' @param exe Character. Path to pestpp-ies executable. If NULL,
#'   uses the bundled binary.
#' @param num_reals Integer. Number of ensemble realisations (overrides
#'   the value in the .pst file).
#' @param noptmax Integer. Maximum number of iterations.
#' @param lambda_scale_fac Numeric vector. Lambda scaling factors.
#' @param ies_par_en Character. Path to existing parameter ensemble file.
#' @param extra_args Named list. Additional PEST++ arguments.
#' @param working_dir Character. Working directory for the run.
#'   Defaults to the directory containing the .pst file.
#' @param verbose Logical. Print stdout/stderr from pestpp-ies.
#' @return A list of class `pesto_ies_result` containing:
#'   \describe{
#'     \item{phi}{data.table of objective function values per iteration}
#'     \item{par_ensemble}{Final parameter ensemble (data.table)}
#'     \item{obs_ensemble}{Final observation ensemble (data.table)}
#'     \item{exit_code}{Integer exit code from pestpp-ies}
#'     \item{runtime_seconds}{Total wall-clock runtime}
#'   }
#' @export
#' @examples
#' # result <- pesto_ies("model.pst", num_reals = 50, noptmax = 3)
#' # plot_phi(result)
pesto_ies <- function(pst_file,
                     exe = NULL,
                     num_reals = 50,
                     noptmax = 4,
                     lambda_scale_fac = c(0.1, 0.5, 1.0),
                     ies_par_en = NULL,
                     extra_args = list(),
                     working_dir = NULL,
                     verbose = TRUE) {

  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) {
    working_dir <- dirname(pst_file)
  }

  exe <- .find_pestpp_exe("pestpp-ies", exe)

  # Build command with ++ args
  args <- c(pst_file)

  # Construct overrides
  overrides <- list(
    ies_num_reals = num_reals,
    noptmax = noptmax,
    ies_lambda_mults = paste(lambda_scale_fac, collapse = ",")
  )
  if (!is.null(ies_par_en)) {
    overrides$ies_par_en <- ies_par_en
  }
  overrides <- c(overrides, extra_args)

  # Write overrides to a temporary ++ section or pass via command line
  override_args <- vapply(names(overrides), function(nm) {
    sprintf("/h :%s=%s", nm, as.character(overrides[[nm]]))
  }, character(1))

  t0 <- proc.time()["elapsed"]

  # Build full args: pst file + override arguments
  all_args <- c(basename(pst_file), override_args)

  result <- system2(
    command = exe,
    args = all_args,
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    env = paste0("PATH=", dirname(exe), ":", Sys.getenv("PATH")),
    wait = TRUE
  )

  runtime <- proc.time()["elapsed"] - t0

  # Parse outputs
  base_name <- tools::file_path_sans_ext(basename(pst_file))
  phi_file <- file.path(working_dir, paste0(base_name, ".phi.actual.csv"))
  par_file <- file.path(working_dir, paste0(base_name, ".0.par.csv"))
  obs_file <- file.path(working_dir, paste0(base_name, ".0.obs.csv"))

  output <- list(
    exit_code = result,
    runtime_seconds = as.numeric(runtime)
  )

  if (file.exists(phi_file)) {
    output$phi <- data.table::fread(phi_file)
  }

  # Find the last iteration's ensemble files
  par_files <- sort(list.files(working_dir, paste0(base_name, "\\.[0-9]+\\.par\\.csv"), full.names = TRUE))
  obs_files <- sort(list.files(working_dir, paste0(base_name, "\\.[0-9]+\\.obs\\.csv"), full.names = TRUE))

  if (length(par_files) > 0) {
    output$par_ensemble <- data.table::fread(par_files[length(par_files)])
  }
  if (length(obs_files) > 0) {
    output$obs_ensemble <- data.table::fread(obs_files[length(obs_files)])
  }

  class(output) <- "pesto_ies_result"
  output
}

#' Run PEST++ GLM (Gauss-Levenberg-Marquardt)
#'
#' Executes the pestpp-glm algorithm for deterministic parameter estimation.
#'
#' @param pst_file Character. Path to the .pst control file.
#' @param exe Character. Path to pestpp-glm executable.
#' @param noptmax Integer. Maximum number of iterations.
#' @param extra_args Named list. Additional PEST++ arguments.
#' @param working_dir Character. Working directory.
#' @param verbose Logical. Print output.
#' @return A list of class `pesto_glm_result`.
#' @export
pesto_glm <- function(pst_file,
                     exe = NULL,
                     noptmax = 20,
                     extra_args = list(),
                     working_dir = NULL,
                     verbose = TRUE) {

  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-glm", exe)

  t0 <- proc.time()["elapsed"]

  result <- system2(
    command = exe,
    args = c(basename(pst_file)),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    wait = TRUE
  )

  runtime <- proc.time()["elapsed"] - t0

  base_name <- tools::file_path_sans_ext(basename(pst_file))

  output <- list(
    exit_code = result,
    runtime_seconds = as.numeric(runtime)
  )

  # Parse .iobj file (iteration objective function)
  iobj_file <- file.path(working_dir, paste0(base_name, ".iobj"))
  if (file.exists(iobj_file)) {
    output$iterations <- data.table::fread(iobj_file)
  }

  # Parse .par file
  par_file <- file.path(working_dir, paste0(base_name, ".par"))
  if (file.exists(par_file)) {
    par_lines <- readLines(par_file)[-1]  # skip header "single point"
    par_dt <- data.table::fread(text = par_lines, header = FALSE,
                                 col.names = c("parnme", "parval", "scale", "offset"))
    output$parameters <- par_dt
  }

  # Parse .rei file (residuals)
  rei_file <- file.path(working_dir, paste0(base_name, ".rei"))
  if (file.exists(rei_file)) {
    output$residuals_file <- rei_file
  }

  # Parse .jco file (Jacobian)
  jco_file <- file.path(working_dir, paste0(base_name, ".jco"))
  if (file.exists(jco_file)) {
    output$jacobian_file <- jco_file
  }

  class(output) <- "pesto_glm_result"
  output
}

#' Run PEST++ SWP (Parametric Sweep)
#'
#' Executes pestpp-swp for embarrassingly parallel model runs across
#' a parameter ensemble.
#'
#' @param pst_file Character. Path to the .pst control file.
#' @param par_ensemble data.table or path. Parameter ensemble.
#' @param exe Character. Path to pestpp-swp executable.
#' @param working_dir Character. Working directory.
#' @param verbose Logical. Print output.
#' @return A list containing observation outputs for each realisation.
#' @export
pesto_sweep <- function(pst_file,
                       par_ensemble,
                       exe = NULL,
                       working_dir = NULL,
                       verbose = TRUE) {

  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-swp", exe)

  # Write parameter ensemble if data.table

  if (data.table::is.data.table(par_ensemble) || is.data.frame(par_ensemble)) {
    sweep_in <- file.path(working_dir, "sweep_in.csv")
    data.table::fwrite(par_ensemble, sweep_in)
  }

  t0 <- proc.time()["elapsed"]
  result <- system2(
    command = exe,
    args = c(basename(pst_file)),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    wait = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  output <- list(
    exit_code = result,
    runtime_seconds = as.numeric(runtime)
  )

  sweep_out <- file.path(working_dir, "sweep_out.csv")
  if (file.exists(sweep_out)) {
    output$results <- data.table::fread(sweep_out)
  }

  class(output) <- "pesto_sweep_result"
  output
}

#' Run PEST++ SEN (Global Sensitivity Analysis)
#'
#' Executes pestpp-sen for Morris or Sobol sensitivity analysis.
#'
#' @param pst_file Character. Path to the .pst control file.
#' @param method Character. "morris" or "sobol".
#' @param exe Character. Path to pestpp-sen executable.
#' @param extra_args Named list. Additional options.
#' @param working_dir Character. Working directory.
#' @param verbose Logical. Print output.
#' @return A list of class `pesto_sen_result`.
#' @export
pesto_sensitivity <- function(pst_file,
                             method = c("morris", "sobol"),
                             exe = NULL,
                             extra_args = list(),
                             working_dir = NULL,
                             verbose = TRUE) {

  method <- match.arg(method)
  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-sen", exe)

  t0 <- proc.time()["elapsed"]
  result <- system2(
    command = exe,
    args = c(basename(pst_file)),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    wait = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  output <- list(
    exit_code = result,
    runtime_seconds = as.numeric(runtime),
    method = method
  )

  base_name <- tools::file_path_sans_ext(basename(pst_file))
  msn_file <- file.path(working_dir, paste0(base_name, ".msn"))
  if (file.exists(msn_file)) {
    output$sensitivity <- data.table::fread(msn_file)
  }

  class(output) <- "pesto_sen_result"
  output
}


#' Find PEST++ executable
#' @param name Name of the executable (e.g., "pestpp-ies")
#' @param user_path User-specified path (overrides search)
#' @return Path to executable
#' @keywords internal
.find_pestpp_exe <- function(name, user_path = NULL) {
  if (!is.null(user_path)) {
    if (!file.exists(user_path)) {
      stop("Executable not found: ", user_path, call. = FALSE)
    }
    return(normalizePath(user_path))
  }

  # Check bundled binaries
  pkg_bin <- system.file("bin", name, package = "PESTO")
  if (nchar(pkg_bin) > 0 && file.exists(pkg_bin)) {
    return(pkg_bin)
  }

  # Check PATH
  sys_exe <- Sys.which(name)
  if (nchar(sys_exe) > 0) {
    return(sys_exe)
  }

  stop(
    "Cannot find ", name, " executable. Either:\n",
    "  1. Install PEST++ and ensure it is on your PATH, or\n",
    "  2. Specify the path via the 'exe' argument.",
    call. = FALSE
  )
}

#' Get PESTO package version information
#'
#' Returns version info for both the PESTO R package and the
#' bundled PEST++ binaries.
#'
#' @return A list with version strings.
#' @export
pesto_version <- function() {
  pkg_ver <- utils::packageVersion("PESTO")

  # Try to get PEST++ version
  pestpp_ver <- tryCatch({
    exe <- .find_pestpp_exe("pestpp-ies", NULL)
    out <- system2(exe, "--version", stdout = TRUE, stderr = TRUE)
    paste(out, collapse = " ")
  }, error = function(e) "not found")

  list(
    pesto_version = as.character(pkg_ver),
    pestpp_version = pestpp_ver,
    platform = R.version$platform,
    r_version = R.version.string
  )
}
