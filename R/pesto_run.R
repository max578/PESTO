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
#' @examples
#' \donttest{
#' if (nzchar(Sys.which("pestpp-ies"))) {
#'   pars <- data.table::data.table(
#'     parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'     parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'     parubnd = c(100, 50), pargp = "hydraulic"
#'   )
#'   obs <- data.table::data.table(
#'     obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'     weight = c(1.0, 1.0), obgnme = "head"
#'   )
#'   pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#'   tf <- tempfile(fileext = ".pst")
#'   on.exit(unlink(tf), add = TRUE)
#'   write_pst(pst, tf)
#'   res <- pesto_ies(tf, num_reals = 3, noptmax = 1, verbose = FALSE)
#'   res$exit_code
#' }
#' }
#' @export
pesto_ies <- function(pst_file,
                     exe = NULL,
                     num_reals = 50,
                     noptmax = 4,
                     lambda_scale_fac = c(0.1, 0.5, 1.0),
                     ies_par_en = NULL,
                     extra_args = list(),
                     working_dir = NULL,
                     verbose = TRUE) {

  # Resolve paths and binary ---------------------------------------------
  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) {
    working_dir <- dirname(pst_file)
  }
  exe <- .find_pestpp_exe("pestpp-ies", exe)

  # Assemble PEST++ command-line overrides --------------------------------
  overrides <- list(
    ies_num_reals    = num_reals,
    noptmax          = noptmax,
    ies_lambda_mults = paste(lambda_scale_fac, collapse = ",")
  )
  if (!is.null(ies_par_en)) {
    overrides$ies_par_en <- ies_par_en
  }
  overrides <- c(overrides, extra_args)

  override_args <- vapply(names(overrides), function(nm) {
    sprintf("/h :%s=%s", nm, as.character(overrides[[nm]]))
  }, character(1))

  # Run pestpp-ies --------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  all_args <- c(basename(pst_file), override_args)

  result <- system2(
    command = exe,
    args    = all_args,
    stdout  = if (verbose) "" else TRUE,
    stderr  = if (verbose) "" else TRUE,
    env     = paste0("PATH=", dirname(exe), ":", Sys.getenv("PATH")),
    wait    = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  # Parse outputs ---------------------------------------------------------
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
#' @examples
#' \donttest{
#' if (nzchar(Sys.which("pestpp-glm"))) {
#'   pars <- data.table::data.table(
#'     parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'     parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'     parubnd = c(100, 50), pargp = "hydraulic"
#'   )
#'   obs <- data.table::data.table(
#'     obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'     weight = c(1.0, 1.0), obgnme = "head"
#'   )
#'   pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#'   tf <- tempfile(fileext = ".pst")
#'   on.exit(unlink(tf), add = TRUE)
#'   write_pst(pst, tf)
#'   res <- pesto_glm(tf, noptmax = 1, verbose = FALSE)
#'   res$exit_code
#' }
#' }
#' @export
pesto_glm <- function(pst_file,
                     exe = NULL,
                     noptmax = 20,
                     extra_args = list(),
                     working_dir = NULL,
                     verbose = TRUE) {

  # Resolve paths and binary ---------------------------------------------
  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-glm", exe)

  # Run pestpp-glm --------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  result <- system2(
    command = exe,
    args    = c(basename(pst_file)),
    stdout  = if (verbose) "" else TRUE,
    stderr  = if (verbose) "" else TRUE,
    wait    = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  # Parse outputs ---------------------------------------------------------
  base_name <- tools::file_path_sans_ext(basename(pst_file))

  output <- list(
    exit_code       = result,
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
#' @examples
#' \donttest{
#' if (nzchar(Sys.which("pestpp-swp"))) {
#'   pars <- data.table::data.table(
#'     parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'     parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'     parubnd = c(100, 50), pargp = "hydraulic"
#'   )
#'   obs <- data.table::data.table(
#'     obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'     weight = c(1.0, 1.0), obgnme = "head"
#'   )
#'   pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#'   tf <- tempfile(fileext = ".pst")
#'   on.exit(unlink(tf), add = TRUE)
#'   write_pst(pst, tf)
#'   par_ens <- data.table::data.table(
#'     real_name = c("r1", "r2", "r3"),
#'     k1 = c(0.8, 1.0, 1.2),
#'     k2 = c(0.4, 0.5, 0.6)
#'   )
#'   res <- pesto_sweep(tf, par_ensemble = par_ens, verbose = FALSE)
#'   res$exit_code
#' }
#' }
#' @export
pesto_sweep <- function(pst_file,
                       par_ensemble,
                       exe = NULL,
                       working_dir = NULL,
                       verbose = TRUE) {

  # Resolve paths and binary ---------------------------------------------
  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-swp", exe)

  # Materialise the parameter ensemble for the binary --------------------
  if (data.table::is.data.table(par_ensemble) ||
      is.data.frame(par_ensemble)) {
    sweep_in <- file.path(working_dir, "sweep_in.csv")
    data.table::fwrite(par_ensemble, sweep_in)
  }

  # Run pestpp-swp --------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  result <- system2(
    command = exe,
    args = c(basename(pst_file)),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    wait = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  # Parse outputs ---------------------------------------------------------
  output <- list(
    exit_code       = result,
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
#' @examples
#' \donttest{
#' if (nzchar(Sys.which("pestpp-sen"))) {
#'   pars <- data.table::data.table(
#'     parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
#'     parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
#'     parubnd = c(100, 50), pargp = "hydraulic"
#'   )
#'   obs <- data.table::data.table(
#'     obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
#'     weight = c(1.0, 1.0), obgnme = "head"
#'   )
#'   pst <- create_pest_scenario(pars, obs, model_command = "echo run")
#'   tf <- tempfile(fileext = ".pst")
#'   on.exit(unlink(tf), add = TRUE)
#'   write_pst(pst, tf)
#'   res <- pesto_sensitivity(tf, method = "morris", verbose = FALSE)
#'   res$method
#' }
#' }
#' @export
pesto_sensitivity <- function(pst_file,
                             method = c("morris", "sobol"),
                             exe = NULL,
                             extra_args = list(),
                             working_dir = NULL,
                             verbose = TRUE) {

  # Resolve paths and binary ---------------------------------------------
  method <- match.arg(method)
  pst_file <- normalizePath(pst_file, mustWork = TRUE)
  if (is.null(working_dir)) working_dir <- dirname(pst_file)
  exe <- .find_pestpp_exe("pestpp-sen", exe)

  # Run pestpp-sen --------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  result <- system2(
    command = exe,
    args = c(basename(pst_file)),
    stdout = if (verbose) "" else TRUE,
    stderr = if (verbose) "" else TRUE,
    wait = TRUE
  )
  runtime <- proc.time()["elapsed"] - t0

  # Parse outputs ---------------------------------------------------------
  output <- list(
    exit_code       = result,
    runtime_seconds = as.numeric(runtime),
    method          = method
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
      stop(
        sprintf("Executable not found: %s", user_path),
        call. = FALSE
      )
    }
    return(normalizePath(user_path))
  }

  # Bundled binary (preferred)  -----------------------------------------
  pkg_bin <- system.file("bin", name, package = "PESTO")
  if (nchar(pkg_bin) > 0 && file.exists(pkg_bin)) {
    return(pkg_bin)
  }

  # Fall back to PATH  --------------------------------------------------
  sys_exe <- Sys.which(name)
  if (nchar(sys_exe) > 0) {
    return(sys_exe)
  }

  stop(
    sprintf(
      paste0(
        "Cannot find %s executable. Either:\n",
        "  1. Install PEST++ and ensure it is on your PATH, or\n",
        "  2. Specify the path via the `exe` argument."
      ),
      name
    ),
    call. = FALSE
  )
}

#' Run IES with an In-Process R Callback Forward Model
#'
#' Drives an Iterative Ensemble Smoother entirely in R, using a user-supplied
#' forward model callable instead of the PEST++ `.pst`-file invocation cycle.
#' Each iteration:
#'
#' 1. Evaluates `forward_model(par_ensemble)` to obtain simulated observations.
#' 2. Forms parameter / observation anomalies and residuals.
#' 3. Calls the C++ kernel [ensemble_solution()] (GLM form, Chen & Oliver 2013)
#'    to compute an `nreal x npar` upgrade.
#' 4. Adds the upgrade to the current ensemble.
#'
#' The classic `.pst`-file path remains available via [pesto_ies()] for full
#' PEST++ compatibility. Use this callback driver when the forward model is
#' itself an R function (e.g. an `apsimx` wrapper from [apsim_callback()],
#' a Python bridge, or a synthetic test problem) and the per-realisation
#' file-I/O overhead of the `.pst` path is the bottleneck.
#'
#' Phase-1 behaviour: single lambda per iteration (or a user-supplied schedule).
#' A line-search over `ies_lambda_mults` matching `pestpp-ies` is a planned
#' Phase-2 enhancement; for the common case of a well-behaved forward model
#' with `lambda = 1`, the GLM update reduces phi reliably (see vignette
#' `apsim-callback`).
#'
#' @section Multi-fidelity:
#' When `forward_model` is a [pesto_multifidelity_model()], `fidelity_schedule`
#' selects which fidelity level each iteration evaluates -- a recycled / padded
#' integer vector, one entry per iteration (default: the highest fidelity every
#' iteration, i.e. exactly the single-fidelity behaviour). This supports
#' fidelity ramping (cheap early iterations, expensive late ones); the final
#' ensemble refresh always uses the highest fidelity so the returned posterior
#' is at full resolution. The control-variate combiner [mf_control_variate()]
#' is the plug-in point for bias-corrected surrogate cascades.
#'
#' @param forward_model One of: a function with signature
#'   `function(theta) -> obs` (where `theta` is an `nreal x npar` numeric
#'   matrix and `obs` an `nreal x nobs` numeric matrix); a
#'   [pesto_forward_model()] (the typed contract object, e.g. one carrying a
#'   `parallel = "multicore"` evaluation strategy); or a
#'   [pesto_multifidelity_model()] (then see `fidelity_schedule`). A bare
#'   function is auto-wrapped via [as_forward_model()] with this call's
#'   `on_failure`. Failed realisations may return rows of `NA`; the driver
#'   tolerates them (see `on_failure`). To run realisations in parallel, pass a
#'   `pesto_forward_model` built with `parallel = "multicore"` rather than a
#'   bare function.
#' @param prior_ensemble Matrix or data.table, `nreal x npar`. Columns are
#'   parameters; an optional `real_name` column is preserved if present.
#'   Column names supply parameter names.
#' @param obs Named numeric vector. Target observations.
#' @param obs_sd Numeric scalar or vector of length `nobs`. Observation
#'   standard deviation(s); the IES weights are `1/obs_sd`.
#' @param noptmax Integer. Number of IES iterations (default 4).
#' @param lambda Numeric scalar or vector. Marquardt lambda per iteration.
#'   A scalar is recycled; a vector shorter than `noptmax` is right-padded
#'   with its last value (default 1.0).
#' @param fidelity_schedule Integer vector or `NULL`. Only consulted when
#'   `forward_model` is a [pesto_multifidelity_model()]: the fidelity level to
#'   evaluate at each iteration (recycled / right-padded to `noptmax`). `NULL`
#'   (default) uses the highest fidelity every iteration. Ignored otherwise.
#' @param parcov Numeric vector of length `npar`, the diagonal of the prior
#'   parameter covariance. Defaults to the column-wise variance of
#'   `prior_ensemble`; zero or negative entries are replaced with 1.0.
#' @param eigthresh Numeric. SVD eigenvalue truncation threshold (default 1e-6).
#' @param use_approx Logical. If TRUE (default), skip the prior-scaling
#'   correction (upgrade_2); matches the typical `pestpp-ies` default.
#' @param on_failure Character. `"na"` (default) carries failed realisations
#'   forward unchanged and proceeds; `"stop"` aborts on any failure.
#' @param verbose Logical. Print per-iteration phi summaries.
#' @return A list of class `c("pesto_ies_callback_result", "pesto_ies_result")`
#'   with components:
#'   \describe{
#'     \item{phi}{data.table of per-realisation phi by iteration.}
#'     \item{par_ensemble}{Final parameter ensemble (data.table).}
#'     \item{obs_ensemble}{Final simulated-observation ensemble (data.table).}
#'     \item{iterations}{List of per-iteration metadata (lambda, mean phi,
#'       failure count).}
#'     \item{runtime_seconds}{Total wall-clock runtime.}
#'     \item{n_forward_evals}{Total number of realisation-level forward
#'       evaluations across all iterations (including the final refresh).}
#'     \item{failure_rate}{Fraction of forward evaluations that returned NA.}
#'     \item{fidelity}{For a [pesto_multifidelity_model()] run, a provenance
#'       list `list(type, schedule, final_level, n_levels, costs)` recording the
#'       realised per-iteration fidelity schedule; `NULL` for a single-fidelity
#'       run. Consumed by [as_manifest()] to populate the manifest contract.}
#'   }
#' @references
#' Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
#' iterative ensemble smoother for efficient history matching and
#' uncertainty quantification. *Computational Geosciences*, 17(4), 689--703.
#' @seealso [pesto_ies()] for the `.pst`-file path; [apsim_callback()]
#'   for the apsimx adapter.
#' @export
#' @examples
#' # Linear-Gaussian recovery toy
#' set.seed(1)
#' npar <- 3; nobs <- 6; nreal <- 80
#' G <- matrix(rnorm(nobs * npar), nobs, npar)
#' theta_true <- c(1.0, -0.5, 2.0)
#' y <- as.numeric(G %*% theta_true) + rnorm(nobs, sd = 0.05)
#' f <- function(theta) theta %*% t(G)
#' prior <- matrix(rnorm(nreal * npar), nreal, npar,
#'                 dimnames = list(NULL, paste0("p", 1:npar)))
#' fit <- pesto_ies_callback(
#'   forward_model = f, prior_ensemble = prior,
#'   obs = setNames(y, paste0("o", 1:nobs)), obs_sd = 0.05,
#'   noptmax = 5, verbose = FALSE
#' )
#' colMeans(as.matrix(fit$par_ensemble[, -1]))  # should approach theta_true
pesto_ies_callback <- function(forward_model,
                               prior_ensemble,
                               obs,
                               obs_sd,
                               noptmax = 4L,
                               lambda = 1.0,
                               fidelity_schedule = NULL,
                               parcov = NULL,
                               eigthresh = 1e-6,
                               use_approx = TRUE,
                               on_failure = c("na", "stop"),
                               verbose = TRUE) {

  # Validate inputs -----------------------------------------------------
  on_failure <- match.arg(on_failure)
  .check_pesto_ies_callback_inputs(forward_model, noptmax, eigthresh)
  noptmax <- as.integer(noptmax)

  # Coerce prior ensemble -----------------------------------------------
  if (data.table::is.data.table(prior_ensemble) ||
      is.data.frame(prior_ensemble)) {
    par_names_local <- setdiff(names(prior_ensemble), "real_name")
    par_mat <- as.matrix(
      data.table::as.data.table(prior_ensemble)[, par_names_local, with = FALSE]
    )
  } else {
    par_mat <- as.matrix(prior_ensemble)
    par_names_local <- colnames(par_mat)
    if (is.null(par_names_local)) {
      par_names_local <- paste0("par", seq_len(ncol(par_mat)))
      colnames(par_mat) <- par_names_local
    }
  }
  storage.mode(par_mat) <- "double"
  nreal <- nrow(par_mat)
  npar  <- ncol(par_mat)
  if (nreal < 2L) {
    stop(
      "`prior_ensemble` must contain at least 2 realisations.",
      call. = FALSE
    )
  }

  # Coerce observations -------------------------------------------------
  obs_vec <- as.numeric(obs)
  obs_names_local <- names(obs)
  if (is.null(obs_names_local)) {
    obs_names_local <- paste0("obs", seq_along(obs_vec))
  }
  nobs <- length(obs_vec)

  obs_sd_vec <- as.numeric(obs_sd)
  if (length(obs_sd_vec) == 1L) obs_sd_vec <- rep(obs_sd_vec, nobs)
  if (length(obs_sd_vec) != nobs || any(obs_sd_vec <= 0)) {
    stop(
      "`obs_sd` must be a positive scalar or length-nobs vector.",
      call. = FALSE
    )
  }
  weights <- 1.0 / obs_sd_vec

  # Prior covariance diagonal -------------------------------------------
  if (is.null(parcov)) {
    parcov_diag <- apply(par_mat, 2L, stats::var)
    parcov_diag[parcov_diag <= 0 | !is.finite(parcov_diag)] <- 1.0
  } else {
    parcov_diag <- as.numeric(parcov)
    if (length(parcov_diag) != npar || any(parcov_diag <= 0)) {
      stop(
        sprintf(
          "`parcov` must be a positive length-%d vector.", npar
        ),
        call. = FALSE
      )
    }
  }
  parcov_inv <- 1.0 / parcov_diag

  # Lambda schedule -----------------------------------------------------
  lambda_seq <- as.numeric(lambda)
  if (length(lambda_seq) < noptmax) {
    lambda_seq <- c(
      lambda_seq,
      rep(lambda_seq[length(lambda_seq)], noptmax - length(lambda_seq))
    )
  }
  lambda_seq <- lambda_seq[seq_len(noptmax)]

  # Forward-model evaluation contract -----------------------------------
  # Accept a bare function, a typed pesto_forward_model, or a
  # multi-fidelity container. The schedule is only meaningful for the
  # latter; eval_at()/lvl_for() hide the dispatch from the loop below.
  is_mf <- S7::S7_inherits(forward_model, pesto_multifidelity_model)
  if (is_mf) {
    n_levels   <- length(forward_model@levels)
    fid_sched  <- .resolve_fidelity_schedule(fidelity_schedule,
                                             n_levels, noptmax)
    top_level  <- n_levels - 1L
    eval_model <- forward_model
  } else {
    fid_sched  <- NULL
    top_level  <- 0L
    eval_model <- .fm_with_nobs(
      as_forward_model(forward_model, on_failure = on_failure), nobs
    )
  }
  lvl_for <- function(k) if (is_mf) fid_sched[k] else 0L
  eval_at <- function(pm, level) {
    if (is_mf) {
      pesto_evaluate(eval_model, pm, level = level)
    } else {
      pesto_evaluate(eval_model, pm)
    }
  }

  # IES iteration loop --------------------------------------------------
  par_prior_mean <- colMeans(par_mat)

  phi_history <- vector("list", noptmax)
  iter_meta   <- vector("list", noptmax)
  total_evals <- 0L
  total_failures <- 0L

  t0 <- proc.time()["elapsed"]

  obs_mat <- eval_at(par_mat, lvl_for(1L))
  total_evals    <- total_evals + nreal
  total_failures <- total_failures + attr(obs_mat, "n_failures")

  for (k in seq_len(noptmax)) {
    ok <- stats::complete.cases(obs_mat)
    if (sum(ok) < 2L) {
      stop(
        sprintf(
          paste0(
            "Iteration %d: fewer than 2 successful realisations. ",
            "Cannot continue."
          ),
          k
        ),
        call. = FALSE
      )
    }
    par_ok <- par_mat[ok, , drop = FALSE]
    obs_ok <- obs_mat[ok, , drop = FALSE]
    nreal_ok <- nrow(par_ok)

    par_mean_iter <- colMeans(par_ok)
    obs_mean_iter <- colMeans(obs_ok)
    par_diff  <- t(sweep(par_ok, 2L, par_mean_iter, "-"))       # npar x nreal_ok
    obs_diff  <- t(sweep(obs_ok, 2L, obs_mean_iter, "-"))       # nobs x nreal_ok
    obs_resid <- matrix(obs_vec, nrow = nobs, ncol = nreal_ok) - t(obs_ok)
    par_resid <- t(sweep(par_ok, 2L, par_prior_mean, "-"))      # npar x nreal_ok

    phi_vec <- compute_phi(obs_resid, weights)
    phi_history[[k]] <- data.table::data.table(
      iteration   = k,
      realisation = which(ok),
      phi         = phi_vec
    )
    if (verbose) {
      message(sprintf(
        "[pesto_ies_callback] iter %d/%d: phi mean=%.4g min=%.4g max=%.4g (lambda=%.3g, nreal_ok=%d/%d)",
        k, noptmax, mean(phi_vec), min(phi_vec), max(phi_vec),
        lambda_seq[k], nreal_ok, nreal
      ))
    }

    Am_k <- matrix(0.0, nrow = npar, ncol = max(nreal_ok - 1L, 1L))

    upgrade <- ensemble_solution(
      par_diff          = par_diff,
      obs_diff          = obs_diff,
      obs_resid         = obs_resid,
      par_resid         = par_resid,
      weights           = weights,
      parcov_inv        = parcov_inv,
      Am                = Am_k,
      cur_lam           = lambda_seq[k],
      eigthresh         = eigthresh,
      use_approx        = use_approx,
      use_prior_scaling = FALSE,
      iter              = as.integer(k),
      reg_factor        = -1.0
    )
    # upgrade: nreal_ok x npar. The C++ kernel returns the *negative-direction*
    # step (the Chen-Oliver 2013 GLM update formula carries an explicit leading
    # minus sign); apply by subtraction so phi descends.
    par_ok_new <- par_ok - upgrade
    par_mat[ok, ] <- par_ok_new

    iter_meta[[k]] <- list(
      lambda     = lambda_seq[k],
      mean_phi   = mean(phi_vec),
      n_failures = nreal - nreal_ok
    )

    if (k < noptmax) {
      obs_mat <- eval_at(par_mat, lvl_for(k + 1L))
      total_evals    <- total_evals + nreal
      total_failures <- total_failures + attr(obs_mat, "n_failures")
    }
  }

  # Final refresh and assemble result -----------------------------------
  # Always refresh at the highest fidelity so the returned ensemble is
  # at full resolution, regardless of the iteration schedule.
  obs_mat_final <- eval_at(par_mat, top_level)
  total_evals    <- total_evals + nreal
  total_failures <- total_failures + attr(obs_mat_final, "n_failures")

  runtime <- as.numeric(proc.time()["elapsed"] - t0)

  par_dt <- data.table::as.data.table(par_mat)
  data.table::setnames(par_dt, par_names_local)
  par_dt[, real_name := paste0("real_", seq_len(nreal))]
  data.table::setcolorder(par_dt, c("real_name", par_names_local))

  obs_dt <- data.table::as.data.table(obs_mat_final)
  data.table::setnames(obs_dt, obs_names_local)
  obs_dt[, real_name := paste0("real_", seq_len(nreal))]
  data.table::setcolorder(obs_dt, c("real_name", obs_names_local))

  # Fidelity provenance: a structured record for a multi-fidelity run
  # (so `as_manifest()` can close the C2 lineage), NULL for a plain
  # single-fidelity run (manifests stay byte-identical to before).
  fidelity_record <- if (is_mf) {
    list(
      type        = "multifidelity",
      schedule    = as.integer(fid_sched),
      final_level = as.integer(top_level),
      n_levels    = as.integer(n_levels),
      costs       = as.numeric(eval_model@costs)
    )
  } else {
    NULL
  }

  output <- list(
    phi             = data.table::rbindlist(phi_history),
    par_ensemble    = par_dt,
    obs_ensemble    = obs_dt,
    iterations      = iter_meta,
    runtime_seconds = runtime,
    n_forward_evals = total_evals,
    failure_rate    = total_failures / total_evals,
    fidelity        = fidelity_record,
    # Assimilation inputs preserved so downstream code (e.g. A5
    # manifest emitter) can reconstruct the full run context.
    obs_target      = stats::setNames(obs_vec, obs_names_local),
    obs_sd          = stats::setNames(obs_sd_vec, obs_names_local),
    weights         = stats::setNames(weights, obs_names_local)
  )
  class(output) <- c("pesto_ies_callback_result", "pesto_ies_result")
  output
}

#' Get PESTO package version information
#'
#' Returns version info for both the PESTO R package and the
#' bundled PEST++ binaries.
#'
#' @return A list with version strings.
#' @examples
#' v <- pesto_version()
#' v$pesto_version
#' v$platform
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


# Internal helpers ------------------------------------------------------

#' Validate the pre-coercion inputs to `pesto_ies_callback()`
#'
#' Catches the cheap, type-and-shape-level errors before the function
#' starts coercing the prior ensemble. Post-coercion checks (`nreal`,
#' `npar`, `nobs` agreement) stay inline in the caller, where the
#' derived shapes are visible.
#'
#' @noRd
#' @keywords internal
.check_pesto_ies_callback_inputs <- function(forward_model, noptmax,
                                             eigthresh) {
  ok_model <- is.function(forward_model) ||
    S7::S7_inherits(forward_model, pesto_forward_model) ||
    S7::S7_inherits(forward_model, pesto_multifidelity_model)
  if (!ok_model) {
    stop(
      paste0(
        "`forward_model` must be a function, a `pesto_forward_model`, ",
        "or a `pesto_multifidelity_model`."
      ),
      call. = FALSE
    )
  }
  noptmax_int <- suppressWarnings(as.integer(noptmax))
  if (length(noptmax_int) != 1L || is.na(noptmax_int) ||
      noptmax_int < 1L) {
    stop(
      "`noptmax` must be a positive integer scalar (>= 1).",
      call. = FALSE
    )
  }
  .assert_positive_scalar(eigthresh, "eigthresh")
  invisible(TRUE)
}
