#' apsimx Forward-Model Adapter for PESTO IES
#'
#' Builds an in-process forward-model closure for [pesto_ies_callback()]
#' that drives APSIM (Next Gen `.apsimx` or Classic `.apsim`) through the
#' `apsimx` package without going via the `.pst`-file path. Used as
#' Year-1 §D4 of the UQ ag-stack roadmap (`uq_ag_stack_roadmap_v0.md`).
#'
#' The returned closure has signature `function(theta) -> obs`, where
#' `theta` is an `nreal x npar` matrix with column names matching
#' `names(param_map)`, and `obs` is an `nreal x nobs` matrix. Each row is
#' produced by:
#'
#' 1. Copying `template` into a fresh per-realisation file under `workdir`.
#' 2. For each `(par_name, node_path)` in `param_map`, calling the
#'    appropriate `apsimx::edit_*` function to write `theta[i, par_name]`
#'    to that node.
#' 3. Calling `apsimx::apsimx()` (or `apsim()` for Classic) and passing
#'    the returned simulation object to `output_extractor()`, which must
#'    return a length-`nobs` numeric vector.
#'
#' Per-realisation failures (APSIM crash, missing output, extractor
#' error) populate an `NA` row; [pesto_ies_callback()] then carries that
#' realisation forward unchanged or aborts, depending on its
#' `on_failure` setting.
#'
#' @section APSIM version compatibility:
#' Tested against the `apsimx` package API as of CRAN 2.7.x. The exact
#' editor function differs between APSIM Next Gen (`edit_apsimx*`) and
#' APSIM Classic (`edit_apsim`); selection is by file extension of
#' `template`. If your installed `apsimx` version exposes a different
#' editor signature, supply `param_writer` to override the default
#' per-parameter writer.
#'
#' @section Concurrency:
#' The returned closure evaluates whatever block of realisations it is
#' handed and writes each to a **unique** per-realisation file under
#' `workdir`, so it is safe to drive in parallel. To run an ensemble
#' concurrently, wrap the closure in a [pesto_forward_model()] with
#' `parallel = "multicore"` (or a custom `map_fn`); the PESTO evaluation
#' engine then dispatches realisations across forked workers, each
#' invoking APSIM on its own input file. Called directly (the bulk
#' path), realisations run serially. `apsimx`'s own thread-safety under
#' heavy ensemble load has not been independently verified, so start
#' with a modest `n_cores`.
#'
#' @param template Character. Path to a working `.apsimx` (Next Gen) or
#'   `.apsim` (Classic) template file. Per-realisation copies are made
#'   into `workdir`; the template itself is never modified.
#' @param param_map Named list. Names are the parameter columns expected
#'   in `theta`; values are character node paths understood by the
#'   appropriate `apsimx::edit_*` function (e.g.
#'   `"Manager.SowingRule.Rule.Population"` for Next Gen).
#' @param output_extractor Function. Takes the object returned by
#'   `apsimx::apsimx()` / `apsim()` (typically a data.frame of report
#'   variables) and returns a length-`nobs` numeric vector. The first
#'   successful realisation defines `nobs`.
#' @param workdir Character. Per-run working directory. Created if it
#'   does not exist; not cleaned automatically (so failures are
#'   inspectable). Default: a fresh `tempfile("apsim_cb_")`.
#' @param param_writer Optional function with signature
#'   `function(file, src.dir, node, value)`. Overrides the default
#'   apsimx editor dispatch. Useful for unusual node paths or
#'   non-standard apsimx versions.
#' @param simulation_runner Optional function with signature
#'   `function(file, src.dir)` returning the simulation object. Overrides
#'   the default `apsimx::apsimx()` / `apsim()` dispatch.
#' @param verbose Logical. Forward verbose flag into apsimx calls and
#'   print per-realisation status (default `FALSE`).
#' @return A closure of signature `function(theta) -> obs` suitable for
#'   the `forward_model =` argument of [pesto_ies_callback()]. The closure
#'   carries an `"apsim_version"` attribute recording the in-use APSIM
#'   binary version (read from `Models --version` on the configured engine,
#'   `NA_character_` if it cannot be determined), so a calibrated run can be
#'   grounded to the exact simulator that produced it. Thread it into the
#'   manifest with
#'   `as_manifest(fit, apsim_version = attr(fm, "apsim_version"))`, so a
#'   downstream consumer can refuse to compare two manifests built against
#'   incompatible APSIM major versions.
#' @seealso [pesto_ies_callback()] for the IES driver; [pesto_ies()]
#'   for the classic `.pst`-file path.
#' @export
#' @examples
#' \dontrun{
#' # Requires apsimx and a working APSIM installation
#' fm <- apsim_callback(
#'   template  = "wheat_wagga.apsimx",
#'   param_map = list(
#'     RUE       = "Wheat.Leaf.Photosynthesis.RUE.FixedValue",
#'     CN2       = "Soil.SoilWater.CN2Bare"
#'   ),
#'   output_extractor = function(sim) {
#'     # sim is a data.frame; extract end-of-season yield trajectory
#'     as.numeric(sim$Wheat.Grain.Total.Wt)
#'   }
#' )
#' prior <- matrix(c(runif(40, 1.0, 2.0), runif(40, 60, 90)),
#'                 ncol = 2, dimnames = list(NULL, c("RUE", "CN2")))
#' fit <- pesto_ies_callback(
#'   forward_model  = fm,
#'   prior_ensemble = prior,
#'   obs            = c(y1 = 4500, y2 = 5200),
#'   obs_sd         = 200,
#'   noptmax        = 4
#' )
#' }
apsim_callback <- function(template,
                           param_map,
                           output_extractor,
                           workdir = tempfile("apsim_cb_"),
                           param_writer = NULL,
                           simulation_runner = NULL,
                           verbose = FALSE) {
  .check_apsim_callback_inputs(template, param_map, output_extractor)

  template <- normalizePath(template, mustWork = TRUE)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  workdir <- normalizePath(workdir)

  ext <- tolower(tools::file_ext(template))
  is_apsimx <- ext == "apsimx"
  if (!is_apsimx && ext != "apsim") {
    stop(
      sprintf(
        "`template` extension must be .apsimx or .apsim; got '.%s'.", ext
      ),
      call. = FALSE
    )
  }

  par_names_expected <- names(param_map)

  writer <- param_writer %||% .default_apsim_writer(is_apsimx, verbose)
  runner <- simulation_runner %||% .default_apsim_runner(is_apsimx, verbose)

  # Ground the run to the simulator that produced it: capture the in-use
  # APSIM binary version once, at construction. Stamped on the returned
  # closure as an attribute so the caller can thread it into the manifest.
  apsim_ver <- .capture_apsim_version()

  fm <- function(theta) {
    # Coerce theta and validate parameter columns ------------------------
    if (!is.matrix(theta)) theta <- as.matrix(theta)
    nreal <- nrow(theta)

    cols <- colnames(theta)
    if (is.null(cols)) {
      if (ncol(theta) != length(par_names_expected)) {
        stop(
          sprintf(
            "`theta` has %d unnamed columns; `param_map` expects %d.",
            ncol(theta), length(par_names_expected)
          ),
          call. = FALSE
        )
      }
      cols <- par_names_expected
      colnames(theta) <- cols
    }
    missing_pars <- setdiff(par_names_expected, cols)
    if (length(missing_pars) > 0L) {
      stop(
        sprintf(
          "`theta` is missing parameters required by `param_map`: %s.",
          paste(missing_pars, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Per-realisation evaluation loop ------------------------------------
    obs_list <- vector("list", nreal)
    for (i in seq_len(nreal)) {
      # Unique per-realisation file: keeps a human-readable `real_<i>_`
      # prefix but appends a tempfile token so concurrent workers (which
      # each see a single-row block, i.e. i == 1) never collide.
      run_file <- tempfile(
        pattern = sprintf("real_%05d_", i),
        tmpdir  = workdir,
        fileext = paste0(".", ext)
      )
      file.copy(template, run_file, overwrite = TRUE)
      src_dir  <- dirname(run_file)
      run_base <- basename(run_file)

      ok_edit <- tryCatch({
        for (par_name in par_names_expected) {
          writer(file  = run_base,
                 src.dir = src_dir,
                 node    = param_map[[par_name]],
                 value   = as.numeric(theta[i, par_name]))
        }
        TRUE
      }, error = function(e) {
        if (verbose) message("[apsim_callback] real ", i,
                             " edit failed: ", conditionMessage(e))
        FALSE
      })

      if (!ok_edit) { obs_list[[i]] <- NULL; next }

      sim <- tryCatch(
        runner(file = run_base, src.dir = src_dir),
        error = function(e) {
          if (verbose) message("[apsim_callback] real ", i,
                               " run failed: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(sim)) { obs_list[[i]] <- NULL; next }

      obs_list[[i]] <- tryCatch(
        as.numeric(output_extractor(sim)),
        error = function(e) {
          if (verbose) message("[apsim_callback] real ", i,
                               " extractor failed: ", conditionMessage(e))
          NULL
        }
      )
    }

    # Assemble nreal x nobs result matrix --------------------------------
    nobs <- max(c(0L, vapply(obs_list,
                             function(v) if (is.null(v)) 0L else length(v),
                             integer(1L))))
    if (nobs == 0L) {
      stop(
        sprintf(
          paste0(
            "apsim_callback: no realisation returned a usable output ",
            "vector. Inspect `%s` for failed runs."
          ),
          workdir
        ),
        call. = FALSE
      )
    }

    out_mat <- matrix(NA_real_, nrow = nreal, ncol = nobs)
    for (i in seq_len(nreal)) {
      v <- obs_list[[i]]
      if (!is.null(v) && length(v) == nobs && all(is.finite(v))) {
        out_mat[i, ] <- v
      }
    }
    out_mat
  }

  attr(fm, "apsim_version") <- apsim_ver
  fm
}

# Read the APSIM engine path apsimx is currently configured to use, without
# disturbing it. apsimx exposes no public getter for `exe.path`, and calling
# `apsimx::apsimx_options()` with no arguments *resets* every option to its
# default -- so its option store (an unexported environment) is read directly.
# Guarded: any change to apsimx's internals degrades to NA, never an error.
.apsimx_exe_path <- function() {
  # 1. What the caller configured explicitly always wins.
  opt_env <- tryCatch(
    get("apsimx.options", envir = asNamespace("apsimx")),
    error = function(e) NULL
  )
  if (is.environment(opt_env)) {
    exe <- tryCatch(
      get0("exe.path", envir = opt_env, inherits = FALSE, ifnotfound = NA),
      error = function(e) NA
    )
    if (length(exe) == 1L && !is.na(exe) && nzchar(exe)) return(exe)
  }

  # 2. APSIM_EXE_PATH. APSIM Next Gen has no discoverable install location on
  # macOS -- apsimx::apsim_version() scans `/Applications` folder names and so
  # sees neither a `~/Applications` install nor a source build, and never
  # consults exe.path at all. Pointing at the engine by environment variable is
  # how it is actually reached from R here. Without this, a machine with a
  # working APSIM reported "no APSIM" and the manifest recorded the provenance
  # of the run as unknown.
  env <- Sys.getenv("APSIM_EXE_PATH", unset = "")
  if (nzchar(env) && file.exists(env)) return(env)

  NA_character_
}

# Capture the in-use APSIM binary version for run provenance, read straight
# from the configured engine via `Models --version`. This deliberately avoids
# `apsimx::apsim_version()`, which on macOS derives the version from install-
# folder names under `/Applications` and never inspects the binary in use --
# so a source build, or an install under `~/Applications`, is invisible to it.
# Defensive: returns NA_character_ when apsimx, the configured path, or a
# working .NET runtime is absent, so a non-APSIM or CI context never errors.
# The captured string lets a downstream consumer refuse a manifest pair built
# against incompatible APSIM major versions.
.capture_apsim_version <- function() {
  if (!requireNamespace("apsimx", quietly = TRUE)) return(NA_character_)

  # apsimx stores either the Models binary itself or the directory holding it.
  exe <- .apsimx_exe_path()
  if (is.na(exe)) return(NA_character_)
  if (dir.exists(exe)) {
    hit <- list.files(exe, pattern = "^Models(\\.exe)?$", full.names = TRUE)
    if (length(hit) == 0L) return(NA_character_)
    exe <- hit[[1L]]
  }
  if (!file.exists(exe)) return(NA_character_)

  # `Models --version` prints e.g. "APSIM 2026.5.8046.0". A missing .NET
  # runtime makes it exit non-zero with empty stdout; that "couldn't
  # determine -> NA" outcome is handled here so it never leaks to the caller.
  out <- tryCatch(
    suppressWarnings(system2(exe, "--version", stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0L)
  )
  out <- out[nzchar(out)]
  if (length(out) == 0L) return(NA_character_)
  sub("^APSIM[[:space:]]+", "", out[[length(out)]])
}

# Default editor dispatch. Returns a function(file, src.dir, node, value).
.default_apsim_writer <- function(is_apsimx, verbose) {
  if (is_apsimx) {
    function(file, src.dir, node, value) {
      apsimx::edit_apsimx(
        file     = file,
        src.dir  = src.dir,
        wrt.dir  = src.dir,
        node     = "Other",
        parm.path = node,
        value    = value,
        overwrite = TRUE,
        verbose   = verbose
      )
    }
  } else {
    function(file, src.dir, node, value) {
      apsimx::edit_apsim(
        file      = file,
        src.dir   = src.dir,
        wrt.dir   = src.dir,
        node      = "Other",
        parm.path = node,
        value     = value,
        overwrite = TRUE,
        verbose   = verbose
      )
    }
  }
}

# Default simulation runner.
.default_apsim_runner <- function(is_apsimx, verbose) {
  if (is_apsimx) {
    function(file, src.dir) {
      apsimx::apsimx(file = file, src.dir = src.dir, silent = !verbose,
                     value = "report")
    }
  } else {
    function(file, src.dir) {
      apsimx::apsim(file = file, src.dir = src.dir, silent = !verbose,
                    value = "report")
    }
  }
}

# Internal helper: a tiny `%||%` so we don't take a rlang dep.
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Validate inputs to `apsim_callback()`
#'
#' Checks the optional `apsimx` Suggests dependency, template path,
#' `param_map` shape (named non-empty list of character node paths),
#' and that `output_extractor` is callable.
#'
#' @noRd
#' @keywords internal
.check_apsim_callback_inputs <- function(template, param_map,
                                         output_extractor) {
  if (!requireNamespace("apsimx", quietly = TRUE)) {
    stop(
      paste0(
        "`apsim_callback()` requires the 'apsimx' package (>= 2.7.0). ",
        "Install with: install.packages('apsimx')."
      ),
      call. = FALSE
    )
  }
  .assert_path_exists(template, "template")
  if (!is.list(param_map) || length(param_map) == 0L ||
    is.null(names(param_map)) || any(names(param_map) == "") ||
    !all(vapply(param_map, is.character, logical(1L)))) {
    stop(
      "`param_map` must be a non-empty named list of character node paths.",
      call. = FALSE
    )
  }
  .assert_function(output_extractor, "output_extractor")
  invisible(TRUE)
}
