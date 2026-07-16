# pestpp_invocation.R -- how PESTO hands control variables to a PEST++ binary.
#
# PEST++ takes NO control variables on the command line. Its parser accepts
#
#     <control_file> [/r|/j] [/e|/g|/h <host:port>]
#
# and nothing else; `/h` selects the PANTHER parallel run manager and its
# argument is a host:port, not a setting. Anything else is a fatal command-line
# error. Grounded 2026-07-15 against the authority -- usgs/pestpp,
# `src/libs/common/utilities.cpp`, `CmdLine::CmdLine` -- not against PESTO's
# own belief about PEST++, which is what got this wrong in the first place.
#
# Control variables therefore reach PEST++ only through the control file, and
# they live in two different places:
#
#   * `++key(value)` lines carry PestppOptions keys. The parser requires the
#     parentheses (`pest_utils::parse_plusplus_line` throws "incorrect format
#     for '++' line (missing'(')" without them), and rejects a duplicate key.
#     An unrecognised key is an error unless `forgive_unknown_args` is set.
#   * NOPTMAX is NOT a PestppOptions key. It belongs to `ControlInfo`
#     (`pest_data_structs.cpp`, `ControlInfo::assign_value_by_key`), which the
#     `++` path never falls back to, so `++noptmax(4)` is REJECTED. It sits in
#     `* control data`, line 7, field 1.
#
# The injection below is textual on purpose. PESTO's own read_pst()/write_pst()
# round-trip is lossy -- read_pst() reads 4 of the 8 `* control data` lines and
# write_pst() re-emits the rest from hard-coded literals -- so routing a user's
# real control file through it would silently replace their settings with
# PESTO's defaults. Editing the lines in place preserves every byte PESTO does
# not deliberately change.


# NOPTMAX is field 1 of line 7 of `* control data`:
#   RSTFLE PESTMODE / NPAR NOBS ... / NTPLFLE NINSFLE ... / RLAMBDA1 ... /
#   RELPARMAX ... / PHIREDSWH ... / NOPTMAX PHIREDSTP ... / ICOV ICOR IEIG
.PESTO_NOPTMAX_LINE <- 7L


#' Assemble the `* model input/output` table
#'
#' PEST++ reads that section positionally: the first NTPLFLE lines are
#' template/model-input pairs and every line after them is an
#' instruction/model-output pair (`Pest.cpp`, `MODEL INPUT/OUTPUT`, which keys
#' off `i_tpl_ins < num_tpl_file`). So the two tables are bound by position and
#' templates must come first. Column names differ between them (`tpl`/`inp` vs
#' `ins`/`out`) and carry no meaning in the file, so they are bound by position
#' rather than by name.
#'
#' @param template_files data.table or `NULL`. Template/model-input pairs.
#' @param instruction_files data.table or `NULL`. Instruction/model-output
#'   pairs.
#'
#' @returns A two-column data.table, templates first.
#' @noRd
.pst_io_files <- function(template_files, instruction_files) {
  parts <- Filter(Negate(is.null), list(template_files, instruction_files))
  if (length(parts) == 0L) {
    return(data.table::data.table())
  }

  data.table::rbindlist(parts, use.names = FALSE)
}


#' Index of a control-file section header
#'
#' @param lines Character vector. The control file.
#' @param section Character. Section name without the leading `*`.
#'
#' @returns Integer index of the header line, or `NA_integer_` if absent.
#' @noRd
.pst_section_index <- function(lines, section) {
  hit <- grep(
    paste0("^[[:space:]]*\\*[[:space:]]*", section, "[[:space:]]*$"),
    lines,
    ignore.case = TRUE
  )
  if (length(hit) == 0L) NA_integer_ else hit[1L]
}


#' Data lines belonging to a control-file section
#'
#' Returns the indices of the lines after a section header and before the next
#' one, ignoring blanks so positional line numbers count data, not whitespace.
#'
#' @param lines Character vector. The control file.
#' @param start Integer. Index of the section header.
#'
#' @returns Integer vector of line indices, in file order.
#' @noRd
.pst_section_lines <- function(lines, start) {
  if (is.na(start) || start >= length(lines)) {
    return(integer(0))
  }

  rest <- seq.int(start + 1L, length(lines))
  next_header <- rest[grepl("^[[:space:]]*\\*", lines[rest])]
  last <- if (length(next_header) > 0L) next_header[1L] - 1L else length(lines)
  if (last < start + 1L) {
    return(integer(0))
  }

  idx <- seq.int(start + 1L, last)
  idx[nzchar(trimws(lines[idx]))]
}


#' Read NOPTMAX from a control file's `* control data` section
#'
#' @param lines Character vector. The control file.
#' @param start Integer. Index of the `* control data` header.
#'
#' @returns Integer NOPTMAX, or `NA_integer_` if the section is too short to
#'   carry one.
#' @noRd
.pst_read_noptmax <- function(lines, start) {
  idx <- .pst_section_lines(lines, start)
  if (length(idx) < .PESTO_NOPTMAX_LINE) {
    return(NA_integer_)
  }

  tokens <- strsplit(
    trimws(lines[idx[.PESTO_NOPTMAX_LINE]]), "[[:space:]]+"
  )[[1L]]
  suppressWarnings(as.integer(tokens[1L]))
}


#' Set NOPTMAX in a control file's `* control data` section
#'
#' @param lines Character vector. The control file.
#' @param noptmax Integer. Iteration cap to write.
#'
#' @returns `lines`, with line 7 of `* control data` re-emitted.
#' @noRd
.pst_set_noptmax <- function(lines, noptmax) {
  start <- .pst_section_index(lines, "control data")
  if (is.na(start)) {
    stop(
      "control file has no `* control data` section, so `noptmax` cannot be ",
      "set.",
      call. = FALSE
    )
  }

  idx <- .pst_section_lines(lines, start)
  if (length(idx) < .PESTO_NOPTMAX_LINE) {
    stop(
      "`* control data` has ", length(idx), " lines; NOPTMAX is on line ",
      .PESTO_NOPTMAX_LINE, ", so this control file is not one PEST++ accepts.",
      call. = FALSE
    )
  }

  target <- idx[.PESTO_NOPTMAX_LINE]
  tokens <- strsplit(trimws(lines[target]), "[[:space:]]+")[[1L]]
  tokens[1L] <- as.character(as.integer(noptmax))
  lines[target] <- paste(tokens, collapse = " ")

  lines
}


#' Merge `++` options into a control file
#'
#' Existing `++` lines for the same key are dropped before the new ones are
#' appended: PEST++ treats a duplicated `++` key as a control-file error, so
#' appending blindly would break a file that already sets the option.
#'
#' @param lines Character vector. The control file.
#' @param options Named list. PestppOptions keys and values.
#'
#' @returns `lines`, with the requested options set.
#' @noRd
.pst_set_plusplus <- function(lines, options) {
  if (length(options) == 0L) {
    return(lines)
  }
  if (is.null(names(options)) || any(!nzchar(names(options)))) {
    stop("every element of `extra_args` must be named.", call. = FALSE)
  }

  keys <- tolower(names(options))
  keep <- vapply(
    lines,
    function(ln) {
      if (!grepl("^[[:space:]]*\\+\\+", ln)) {
        return(TRUE)
      }
      key <- tolower(trimws(sub("^[[:space:]]*\\+\\+([^(]*)\\(.*$", "\\1", ln)))
      !(key %in% keys)
    },
    logical(1),
    USE.NAMES = FALSE
  )
  lines <- lines[keep]

  emitted <- vapply(
    seq_along(options),
    function(i) {
      value <- options[[i]]
      if (length(value) != 1L) {
        stop(
          "`", names(options)[i], "` must be a single value; PEST++ reads a ",
          "`++` option as one `key(value)` pair.",
          call. = FALSE
        )
      }
      paste0("++", names(options)[i], "(", as.character(value), ")")
    },
    character(1)
  )

  c(lines, emitted)
}


#' Write the control file a PEST++ run will actually read
#'
#' Injects the caller's control variables into a copy of `pst_file` and writes
#' it beside the original, so relative template, instruction, and model-command
#' paths still resolve. The copy is left in place deliberately: it is the exact
#' input PEST++ was given, and the run's outputs are named after it.
#'
#' @param pst_file Character. Path to the user's control file.
#' @param working_dir Character. Directory to write the run file into.
#' @param noptmax Integer or `NULL`. `NULL` leaves the file's own value.
#' @param pestpp_options Named list. `++` options to set.
#'
#' @returns Character. Path to the control file to run -- `pst_file` itself when
#'   there is nothing to inject.
#' @noRd
.pesto_run_control_file <- function(pst_file,
                                    working_dir,
                                    noptmax = NULL,
                                    pestpp_options = list()) {
  # Nothing to change means nothing to write: run the caller's own file, so a
  # default call leaves no artefact behind and names its outputs after the file
  # the caller passed.
  if (is.null(noptmax) && length(pestpp_options) == 0L) {
    return(pst_file)
  }

  lines <- readLines(pst_file, warn = FALSE)

  # `* control data keyword` is the alternative to the positional section, and
  # PEST++ refuses to read a file that mixes it with `++` args. Both of this
  # function's injections would be wrong there, so refuse rather than write a
  # file the binary will reject for reasons the caller cannot see.
  if (any(grepl("^[[:space:]]*\\*[[:space:]]*control data keyword",
                lines, ignore.case = TRUE))) {
    stop(
      "control file uses `* control data keyword`, which PEST++ does not ",
      "allow alongside the `++` options PESTO sets. Convert it to the ",
      "positional `* control data` section, or set the options in the file ",
      "and leave PESTO's arguments at their defaults.",
      call. = FALSE
    )
  }

  if (!is.null(noptmax)) {
    lines <- .pst_set_noptmax(lines, noptmax)
  }
  lines <- .pst_set_plusplus(lines, pestpp_options)

  base <- tools::file_path_sans_ext(basename(pst_file))
  run_file <- file.path(working_dir, paste0(base, "_pesto.pst"))
  writeLines(lines, run_file)

  run_file
}


#' Read a PEST++ binary's version from its startup banner
#'
#' PEST++ has no `--version` flag. `pestpp-ies --version` is read as a *control
#' file name*: the binary prints its banner, fails on the missing
#' `--version.pst`, exits 1, and leaves `--version.log`, `--version.rec` and
#' `--version.rst` in the caller's working directory. Asking for the version
#' that way returned the whole failed run -- `std::exception` and all -- as the
#' version string, and looked right only because the banner is printed before
#' the error.
#'
#' Invoked with no arguments the binary prints the same banner, exits 1, and
#' writes nothing. That is the probe. Run from a temporary directory anyway, so
#' a future PEST++ that does emit files cannot litter the caller's.
#'
#' @param exe Character. Path to the binary.
#'
#' @returns Character. The version string, or `"unknown"` if the banner carries
#'   no `version:` line.
#' @noRd
.pestpp_binary_version <- function(exe) {
  probe_dir <- tempfile("pestpp_version_")
  dir.create(probe_dir)
  on.exit(unlink(probe_dir, recursive = TRUE), add = TRUE)

  old <- setwd(probe_dir)
  on.exit(setwd(old), add = TRUE, after = FALSE)

  # Exit status is 1 even on the clean banner path, so it carries no signal
  # here and is deliberately ignored; the banner is the answer.
  out <- suppressWarnings(system2(
    exe, args = character(0), stdout = TRUE, stderr = TRUE
  ))

  hit <- grep("^[[:space:]]*version:", out, value = TRUE)
  if (length(hit) == 0L) {
    return("unknown")
  }

  trimws(sub("^[[:space:]]*version:[[:space:]]*", "", hit[1L]))
}


#' Run a PEST++ binary on a control file
#'
#' Runs in `working_dir`, because the control file is passed by basename and
#' PEST++ resolves it -- and every relative path inside it -- against the
#' process's working directory.
#'
#' @param exe Character. Path to the binary.
#' @param pst_file Character. Control file to run; passed by basename.
#' @param working_dir Character. Directory to run in.
#' @param verbose Logical. Stream PEST++ output.
#'
#' @returns Integer exit code, in both the streaming and the capturing case.
#' @noRd
.pesto_run_pestpp <- function(exe, pst_file, working_dir, verbose = TRUE) {
  old <- setwd(working_dir)
  on.exit(setwd(old), add = TRUE)

  # The only argument PEST++ accepts here is the control file: no run-manager
  # switch (serial), and no control variables -- those are in the file.
  res <- suppressWarnings(system2(
    command = exe,
    args    = basename(pst_file),
    stdout  = if (verbose) "" else TRUE,
    stderr  = if (verbose) "" else TRUE,
    wait    = TRUE
  ))

  # system2() returns the exit code when output is streamed, but the captured
  # OUTPUT when it is not -- with the status hidden in an attribute that is
  # absent on success. Returning that verbatim put PEST++'s entire log in
  # `$exit_code` whenever `verbose = FALSE`, which is what the examples use.
  # Normalise here so callers get an exit code either way.
  if (is.character(res)) {
    status <- attr(res, "status")
    return(if (is.null(status)) 0L else as.integer(status))
  }

  as.integer(res)
}
