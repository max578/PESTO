# test-pestpp-invocation.R -- independent-oracle tests for the PEST++ interface
# (recipe 52 / the Independent Oracle Principle).
#
# PESTO asserts facts about an external authority it did not author: how PEST++
# parses its command line, and where PEST++ reads control variables from. Every
# test in this package that touched the invocation path checked PESTO against
# PESTO, so a fabricated interface stayed green -- and did. `pesto_ies()` built
# `/h :noptmax=4` switches that PEST++ has never accepted; the binary is not a
# test dependency, so nothing ever contradicted it.
#
# The oracle here is `.pestpp_cmdline_verdict()`: a transcription of the REAL
# parser, from
#
#   usgs/pestpp, src/libs/common/utilities.cpp, CmdLine::CmdLine (read
#   2026-07-15)
#
# It is deliberately written from the C++ control flow rather than from what
# PESTO believes, so an assertion failing here is a finding about PESTO. It is
# a stand-in for the binary, not a substitute: it settles grammar, not
# behaviour. Only a real pestpp run can do that, and none is available here --
# which is exactly the gap that let the fabrication live.


# ---- fixture ---------------------------------------------------------------

# A minimal scenario, only so write_pst() has something to emit. The tests are
# about where control variables land, not about the model.
.make_test_scenario <- function() {
  create_pest_scenario(
    parameters = data.table::data.table(
      parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
      parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
      parubnd = c(100, 50), pargp = "hydraulic"
    ),
    observations = data.table::data.table(
      obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
      weight = c(1.0, 1.0), obgnme = "head"
    ),
    model_command = "echo run"
  )
}


# ---- the oracle ------------------------------------------------------------

# Transcribed from CmdLine::CmdLine. Returns "ok" or the error PEST++ raises.
#
# `argv` is the FULL vector the process receives, program name first, so the
# positions line up with the C++ (whose argv[0] is the binary). Passing only
# the arguments would shift every index by one and quietly change the verdict.
#
#   argv[0]                      -> binary
#   argv[1]                      -> control file (required)
#   "/r", "/j"                   -> stripped as flags
#   argc 3, or > 4               -> error unless argv[2] == "/e"
#   argc 2                       -> serial run manager
#   argv[2] in /e /g /h          -> run manager; anything else is an error
#   /h -> argv[3] must contain ":"; leading ":" = master, port must parse int
.pestpp_cmdline_verdict <- function(argv) {
  n <- length(argv)
  if (n < 2L) return("too few args, no control file name found")

  lower <- tolower(argv)
  keep <- !(lower %in% c("/r", "/j"))
  argv <- argv[keep]
  lower <- lower[keep]
  n <- length(argv)

  if ((n == 3L) || (n > 4L)) {
    if (lower[3L] != "/e") {
      return("wrong number of args, expecting 2 or 4")
    }
  }
  if (n == 2L) return("ok")

  third <- lower[3L]
  if (!third %in% c("/e", "/g", "/h")) {
    return(paste0("unrecognized commandline arg '", third, "'"))
  }
  if (third != "/h") return("ok")

  forth <- argv[4L]
  if (!grepl(":", forth, fixed = TRUE)) {
    return(paste0("panther master/worker arg '", forth, "' doesn't have a ':'"))
  }
  if (startsWith(forth, ":")) {
    port <- substring(forth, 2L)
    if (is.na(suppressWarnings(as.integer(port)))) {
      return(paste0("error casting master port number '", port, "' to int"))
    }
    return("ok")
  }
  tokens <- strsplit(forth, ":", fixed = TRUE)[[1L]]
  if (length(tokens) != 2L) return("wrong number of colon-delimited tokens")
  if (is.na(suppressWarnings(as.integer(tokens[2L])))) {
    return(paste0("error casting master port number '", tokens[2L], "' to int"))
  }
  "ok"
}


test_that("the oracle rejects the syntax PESTO used to emit", {
  # Guards the oracle itself: if this passes, the transcription is not simply
  # accepting everything. These are the exact argv the old override code built.
  expect_equal(
    .pestpp_cmdline_verdict(c("pestpp-ies", "model.pst", "/h", ":noptmax=4")),
    "error casting master port number 'noptmax=4' to int"
  )
  expect_equal(
    .pestpp_cmdline_verdict(c(
      "pestpp-ies", "model.pst",
      "/h", ":ies_num_reals=50", "/h", ":noptmax=4"
    )),
    "wrong number of args, expecting 2 or 4"
  )
  # And that it accepts what the parser accepts.
  expect_equal(.pestpp_cmdline_verdict(c("pestpp-ies", "model.pst")), "ok")
  expect_equal(
    .pestpp_cmdline_verdict(c("pestpp-ies", "model.pst", "/h", "localhost:4004")),
    "ok"
  )
})


# ---- the emitted command line ----------------------------------------------

test_that("every run function emits a command line PEST++ accepts", {
  # The regression that matters: whatever PESTO sends must parse. Captured at
  # the system2() boundary, so this asserts what the binary would receive.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "oracle.pst")
  on.exit(unlink(list.files(tempdir(), "^oracle", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  seen <- new.env(parent = emptyenv())
  seen$argv <- list()

  local_mocked_bindings(
    .find_pestpp_exe = function(name, exe = NULL) name,
    .pesto_run_pestpp = function(exe, pst_file, working_dir, verbose = TRUE) {
      # Mirrors what .pesto_run_pestpp() passes to system2(): the binary, then
      # the control file's basename and nothing else.
      seen$argv <- c(seen$argv, list(c(exe, basename(pst_file))))
      0L
    }
  )

  pesto_ies(tf, num_reals = 3, noptmax = 1, verbose = FALSE)
  pesto_glm(tf, noptmax = 2, verbose = FALSE)
  pesto_sensitivity(tf, method = "sobol", verbose = FALSE)

  expect_length(seen$argv, 3L)
  for (argv in seen$argv) {
    expect_equal(.pestpp_cmdline_verdict(argv), "ok")
    # No control variable may ride on the command line: PEST++ has no such
    # mechanism, so anything beyond the control file is a fabrication.
    expect_length(argv, 2L)
    expect_false(any(grepl("=", argv, fixed = TRUE)))
    expect_false(any(startsWith(argv[-1L], "/")))
  }
})


# ---- control variables reach the control file ------------------------------

test_that("noptmax is written to * control data, not as a ++ option", {
  # `++noptmax()` is rejected by PEST++: NOPTMAX belongs to ControlInfo, and
  # the ++ path never falls back to it. Asserting BOTH halves -- present in the
  # section, absent from the ++ lines -- because either alone would pass while
  # the value went to the wrong place.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "noptmax.pst")
  on.exit(unlink(list.files(tempdir(), "^noptmax", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  run_file <- .pesto_run_control_file(tf, tempdir(), noptmax = 7L)
  lines <- readLines(run_file)

  start <- .pst_section_index(lines, "control data")
  cd <- lines[.pst_section_lines(lines, start)]
  expect_equal(strsplit(trimws(cd[7L]), "[[:space:]]+")[[1L]][1L], "7")

  expect_false(any(grepl("^\\+\\+[[:space:]]*noptmax", lines, ignore.case = TRUE)))
})

test_that("noptmax changes the control file PEST++ reads", {
  # The parameter must change the output; at a fixed value an honoured and an
  # ignored noptmax are indistinguishable, which is how it stayed dead.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "vary.pst")
  on.exit(unlink(list.files(tempdir(), "^vary", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  read_noptmax <- function(n) {
    f <- .pesto_run_control_file(tf, tempdir(), noptmax = n)
    lines <- readLines(f)
    .pst_read_noptmax(lines, .pst_section_index(lines, "control data"))
  }

  expect_equal(read_noptmax(1L), 1L)
  expect_equal(read_noptmax(99L), 99L)
})

test_that("noptmax = NULL leaves the file's own iteration cap alone", {
  # The default. A caller who set NOPTMAX in their own control file must get
  # that value, not one of PESTO's -- the old defaults (4 / 20) would have
  # silently overridden it on every default call.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "keepcap.pst")
  on.exit(unlink(list.files(tempdir(), "^keepcap", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  lines <- .pst_set_noptmax(readLines(tf), 13L)
  writeLines(lines, tf)

  # Force a run file to be written anyway, so this tests NULL rather than the
  # no-injection short-circuit below.
  run_file <- .pesto_run_control_file(
    tf, tempdir(),
    noptmax = NULL, pestpp_options = list(ies_num_reals = 6L)
  )
  out <- readLines(run_file)

  expect_equal(
    .pst_read_noptmax(out, .pst_section_index(out, "control data")), 13L
  )
  expect_true(any(grepl("^\\+\\+ies_num_reals\\(6\\)$", out)))
})

test_that("pesto_ies and pesto_glm default to the file's own iteration cap", {
  # Grades the DEFAULT, not just the NULL path: the old defaults (4 for
  # pesto_ies, 20 for pesto_glm) were documented as overrides but had never
  # once been applied, so honouring them now would silently retune every
  # caller's control file. Tested through the exported functions because a
  # regression would live in the signature, which .pesto_run_control_file
  # tests cannot see.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "defcap.pst")
  on.exit(unlink(list.files(tempdir(), "^defcap", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)
  writeLines(.pst_set_noptmax(readLines(tf), 13L), tf)

  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .find_pestpp_exe = function(name, exe = NULL) name,
    .pesto_run_pestpp = function(exe, pst_file, working_dir, verbose = TRUE) {
      lines <- readLines(pst_file)
      seen$noptmax <- .pst_read_noptmax(
        lines, .pst_section_index(lines, "control data")
      )
      0L
    }
  )

  pesto_ies(tf, verbose = FALSE)
  expect_equal(seen$noptmax, 13L)

  pesto_glm(tf, verbose = FALSE)
  expect_equal(seen$noptmax, 13L)

  # And an explicit value still overrides, or the parameter would be dead again.
  pesto_glm(tf, noptmax = 2L, verbose = FALSE)
  expect_equal(seen$noptmax, 2L)
})

test_that("nothing to inject writes no file and runs the caller's own", {
  # With no noptmax and no options there is nothing to change, so PESTO must
  # not manufacture a `_pesto.pst` the caller never asked for.
  scenario <- .make_test_scenario()
  wd <- file.path(tempdir(), "noinject")
  dir.create(wd, showWarnings = FALSE)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  tf <- file.path(wd, "plain.pst")
  write_pst(scenario, tf)

  run_file <- .pesto_run_control_file(tf, wd)

  expect_identical(run_file, tf)
  expect_false(file.exists(file.path(wd, "plain_pesto.pst")))
})

test_that("extra_args become ++key(value) lines", {
  # The parentheses are not cosmetic: parse_plusplus_line throws "incorrect
  # format for '++' line (missing'(')" without them.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "extra.pst")
  on.exit(unlink(list.files(tempdir(), "^extra", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  run_file <- .pesto_run_control_file(
    tf, tempdir(),
    pestpp_options = list(ies_num_reals = 12, ies_bad_phi = 1e10)
  )
  lines <- readLines(run_file)

  expect_true(any(lines == "++ies_num_reals(12)"))
  expect_true(any(grepl("^\\+\\+ies_bad_phi\\(", lines)))
})

test_that("an option the file already sets is replaced, not duplicated", {
  # PEST++ treats a duplicated ++ key as a control-file error, so appending
  # blindly would break a file that already carries the option.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "dup.pst")
  on.exit(unlink(list.files(tempdir(), "^dup", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)
  writeLines(c(readLines(tf), "++ies_num_reals(999)"), tf)

  run_file <- .pesto_run_control_file(
    tf, tempdir(), pestpp_options = list(ies_num_reals = 5)
  )
  lines <- readLines(run_file)

  expect_length(grep("^\\+\\+ies_num_reals\\(", lines), 1L)
  expect_true(any(lines == "++ies_num_reals(5)"))
})

test_that("pesto_sensitivity routes method through GSA_METHOD", {
  # `method` was read (match.arg) and stored in the result, but never reached
  # the binary: pestpp-sen defaults to Morris, so `method = "sobol"` labelled
  # Morris output "sobol". Both branches asserted -- equality would mean the
  # parameter is inert again.
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "gsa.pst")
  on.exit(unlink(list.files(tempdir(), "^gsa", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .find_pestpp_exe = function(name, exe = NULL) name,
    .pesto_run_pestpp = function(exe, pst_file, working_dir, verbose = TRUE) {
      seen$lines <- readLines(pst_file)
      0L
    }
  )

  pesto_sensitivity(tf, method = "sobol", verbose = FALSE)
  sobol <- grep("^\\+\\+gsa_method\\(", seen$lines, value = TRUE)

  pesto_sensitivity(tf, method = "morris", verbose = FALSE)
  morris <- grep("^\\+\\+gsa_method\\(", seen$lines, value = TRUE)

  expect_equal(sobol, "++gsa_method(SOBOL)")
  expect_equal(morris, "++gsa_method(MORRIS)")
  expect_false(identical(sobol, morris))
})


# ---- refusals --------------------------------------------------------------

test_that("a keyword-style control file is refused, not silently mangled", {
  # `* control data keyword` cannot carry ++ args (PEST++ rejects the mix) and
  # its NOPTMAX is a keyword, not a positional field, so both injections would
  # be wrong. Refusing beats writing a file the binary rejects opaquely.
  tf <- tempfile(fileext = ".pst")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("pcf", "* control data keyword", "noptmax 5"), tf)

  expect_error(
    .pesto_run_control_file(tf, tempdir(), noptmax = 3L),
    "control data keyword"
  )
})

test_that("an unnamed extra_args element is refused", {
  scenario <- .make_test_scenario()
  tf <- file.path(tempdir(), "unnamed.pst")
  on.exit(unlink(list.files(tempdir(), "^unnamed", full.names = TRUE)), add = TRUE)
  write_pst(scenario, tf)

  expect_error(
    .pesto_run_control_file(tf, tempdir(), pestpp_options = list(5)),
    "must be named"
  )
})


# ---- the run happens where the control file lives --------------------------

test_that("the binary runs in the control file's directory, with a bare argv", {
  # PEST++ is handed a basename, and resolves it -- and every relative template,
  # instruction, and model-command path inside it -- against its own working
  # directory. Running from R's cwd meant it could not find the file at all.
  #
  # This is the only test that reaches the real system2() boundary: the ones
  # above mock .pesto_run_pestpp(), so they cannot see what it actually sends.
  # Both halves of the contract are asserted right here at the hand-off.
  run_dir <- file.path(tempdir(), "rundir")
  dir.create(run_dir, showWarnings = FALSE)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  probe <- new.env(parent = emptyenv())
  fake <- function(command, args, ...) {
    probe$cwd <- normalizePath(getwd())
    probe$argv <- c(command, args)
    0L
  }
  local_mocked_bindings(system2 = fake, .package = "base")

  before <- normalizePath(getwd())
  .pesto_run_pestpp("pestpp-ies", file.path(run_dir, "x.pst"), run_dir,
                    verbose = FALSE)

  expect_equal(probe$cwd, normalizePath(run_dir))
  expect_equal(normalizePath(getwd()), before)  # and restored afterwards

  # The argv PEST++ would really receive: binary, control file, nothing else.
  expect_equal(probe$argv, c("pestpp-ies", "x.pst"))
  expect_equal(.pestpp_cmdline_verdict(probe$argv), "ok")
})
