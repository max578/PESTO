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
#' Phase-1 D4 runs realisations **serially**. Parallel execution via
#' `future` or `mirai` is a planned follow-up; `apsimx`'s thread-safety
#' under ensemble load has not been verified.
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
#'   the `forward_model =` argument of [pesto_ies_callback()].
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

  function(theta) {
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
      run_file <- file.path(workdir, sprintf("real_%05d.%s", i, ext))
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
