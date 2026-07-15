# Structural tests for apsim_callback(). Real APSIM integration is gated on
# `skip_if_not_installed("apsimx")` AND availability of a working APSIM
# binary; absent either, the tests still exercise input validation and
# the closure-construction path with user-supplied stubs.

test_that("apsim_callback requires apsimx package", {
  if (requireNamespace("apsimx", quietly = TRUE)) {
    skip("apsimx is installed; this branch only runs when it is missing")
  }
  tmpf <- tempfile(fileext = ".apsimx"); file.create(tmpf); on.exit(unlink(tmpf))
  expect_error(
    apsim_callback(
      template         = tmpf,
      param_map        = list(a = "Some.Path"),
      output_extractor = identity
    ),
    "apsimx"
  )
})

test_that("apsim_callback validates inputs", {
  skip_if_not_installed("apsimx")
  tmpf <- tempfile(fileext = ".apsimx"); file.create(tmpf); on.exit(unlink(tmpf))

  expect_error(
    apsim_callback("nonexistent.apsimx",
                   list(a = "X"), identity),
    "`template` does not exist"
  )
  expect_error(
    apsim_callback(tmpf, list(), identity),
    "param_map"
  )
  expect_error(
    apsim_callback(tmpf, list(a = 123), identity),  # non-character node path
    "param_map"
  )
  expect_error(
    apsim_callback(tmpf, list(a = "X"), output_extractor = "not_a_fn"),
    "output_extractor"
  )

  tmpf_bad <- tempfile(fileext = ".txt"); file.create(tmpf_bad)
  on.exit(unlink(tmpf_bad), add = TRUE)
  expect_error(
    apsim_callback(tmpf_bad, list(a = "X"), identity),
    "extension"
  )
})

test_that("apsim_callback closure runs with stub writer + runner", {
  skip_if_not_installed("apsimx")
  # Use param_writer / simulation_runner overrides so the test does not
  # require a real APSIM installation. This validates the per-realisation
  # plumbing of the closure (file-copy, edit, run, extract, NA-handling).
  tmpf <- tempfile(fileext = ".apsimx")
  writeLines("stub-apsimx-template", tmpf)
  on.exit(unlink(tmpf))

  workdir <- tempfile("apsim_cb_test_")
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  edits <- new.env(parent = emptyenv())
  edits$values <- list()

  stub_writer <- function(file, src.dir, node, value) {
    edits$values <- c(edits$values, list(list(file = file, node = node,
                                              value = value)))
    invisible(TRUE)
  }
  # Stub runner reads the realisation file path back from its working
  # directory and returns a deterministic 2-vector that depends on the
  # most recently written parameter values.
  stub_runner <- function(file, src.dir) {
    list(file = file, src.dir = src.dir)
  }
  stub_extractor <- function(sim) {
    # Pull the last 2 written values; deterministic per realisation.
    last <- utils::tail(edits$values, 2L)
    vapply(last, `[[`, numeric(1L), "value")
  }

  fm <- apsim_callback(
    template          = tmpf,
    param_map         = list(a = "Node.A", b = "Node.B"),
    output_extractor  = stub_extractor,
    workdir           = workdir,
    param_writer      = stub_writer,
    simulation_runner = stub_runner
  )

  theta <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2L, byrow = TRUE,
                  dimnames = list(NULL, c("a", "b")))
  out <- fm(theta)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 2L))
  expect_true(all(is.finite(out)))
  # One unique per-realisation file per row (parallel-safe naming).
  run_files <- list.files(workdir, pattern = "^real_\\d{5}_.*\\.apsimx$")
  expect_length(run_files, 2L)
})

test_that("apsim_callback rejects theta missing required parameters", {
  skip_if_not_installed("apsimx")
  tmpf <- tempfile(fileext = ".apsimx"); writeLines("stub", tmpf)
  on.exit(unlink(tmpf))
  workdir <- tempfile("apsim_cb_"); on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  fm <- apsim_callback(
    template          = tmpf,
    param_map         = list(a = "A", b = "B"),
    output_extractor  = function(x) c(0, 0),
    workdir           = workdir,
    param_writer      = function(...) invisible(NULL),
    simulation_runner = function(...) list()
  )
  bad_theta <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(NULL, "a"))
  expect_error(fm(bad_theta), "missing parameters")
})

test_that(".capture_apsim_version reads the version from `Models --version`", {
  skip_if_not_installed("apsimx")
  skip_on_os("windows")

  # A stand-in `Models` that prints an APSIM-style version banner, so the
  # capture path is exercised deterministically without a real APSIM install.
  stub_dir <- tempfile("apsim_stub_")
  dir.create(stub_dir)
  on.exit(unlink(stub_dir, recursive = TRUE), add = TRUE)
  stub <- file.path(stub_dir, "Models")
  writeLines(c("#!/bin/sh", "echo 'APSIM 9.9.9-test+deadbeef'"), stub)
  Sys.chmod(stub, mode = "0755")

  old_exe <- .apsimx_exe_path()
  on.exit(suppressWarnings(
    apsimx::apsimx_options(exe.path = old_exe, warn.find.apsimx = FALSE)
  ), add = TRUE)
  suppressWarnings(
    apsimx::apsimx_options(exe.path = stub, warn.find.apsimx = FALSE)
  )

  expect_identical(.capture_apsim_version(), "9.9.9-test+deadbeef")
})

test_that(".capture_apsim_version returns NA when the engine path is unusable", {
  skip_if_not_installed("apsimx")

  old_exe <- .apsimx_exe_path()
  on.exit(suppressWarnings(
    apsimx::apsimx_options(exe.path = old_exe, warn.find.apsimx = FALSE)
  ), add = TRUE)
  suppressWarnings(
    apsimx::apsimx_options(
      exe.path = file.path(tempdir(), "no_such_Models"),
      warn.find.apsimx = FALSE
    )
  )

  expect_true(is.na(.capture_apsim_version()))
})
