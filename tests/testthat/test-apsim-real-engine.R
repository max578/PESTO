# test-apsim-real-engine.R -- the oracle no stub can be: APSIM itself.
#
# test-apsim-callback.R exercises the closure with user-supplied stubs, so it
# checks PESTO against PESTO. That is the same shape that let PESTO ship a
# PEST++ command line PEST++ has never accepted for three months: the authority
# was not a test dependency, so nothing could contradict the package's belief
# about it.
#
# Opt-in: needs a real APSIM Next Gen engine, found via APSIM_EXE_PATH (the
# directory holding `Models`, as ~/.Renviron sets it) or an explicit
# apsimx_options(exe.path=). Skipped otherwise -- including on CRAN. The skip
# is honest; a stubbed "pass" would not be, since stubbing the engine is
# precisely the gap these tests exist to close.
#
#   APSIM_EXE_PATH=/path/to/net8.0 R -e 'devtools::test(filter="apsim-real")'

.apsim_models_exe <- function() {
  dir <- Sys.getenv("APSIM_EXE_PATH", unset = "")
  if (!nzchar(dir)) return("")
  exe <- if (dir.exists(dir)) file.path(dir, "Models") else dir
  if (file.exists(exe)) exe else ""
}

test_that("APSIM_EXE_PATH alone is enough to find the engine", {
  # The engine is unreachable by any automatic scan on macOS: apsimx's own
  # apsim_version() reads `/Applications` folder names, so it sees neither a
  # ~/Applications install nor a source build, and never consults exe.path.
  # PESTO consulted only exe.path, so a machine with APSIM installed and
  # APSIM_EXE_PATH exported still reported no engine, and the manifest recorded
  # the run's simulator provenance as unknown.
  exe <- .apsim_models_exe()
  skip_if(!nzchar(exe), "no APSIM engine (set APSIM_EXE_PATH)")
  skip_if_not_installed("apsimx")

  # Clear anything a previous test left configured, so this proves the env var
  # is doing the work rather than riding on apsimx's option.
  opt_env <- tryCatch(get("apsimx.options", envir = asNamespace("apsimx")),
                      error = function(e) NULL)
  skip_if(!is.environment(opt_env), "apsimx internals not as expected")
  old <- tryCatch(get0("exe.path", envir = opt_env, inherits = FALSE,
                       ifnotfound = NA), error = function(e) NA)
  assign("exe.path", NA, envir = opt_env)
  on.exit(assign("exe.path", old, envir = opt_env), add = TRUE)

  found <- PESTO:::.apsimx_exe_path()

  expect_false(is.na(found))
  expect_true(file.exists(if (dir.exists(found)) file.path(found, "Models") else found))
})

test_that("the captured version comes from the running binary", {
  # Not from apsimx::apsim_version(), which cannot see either engine on this
  # platform, and not from a version-bearing folder name, which a source build
  # does not have. Only the binary knows what it is.
  exe <- .apsim_models_exe()
  skip_if(!nzchar(exe), "no APSIM engine (set APSIM_EXE_PATH)")
  skip_if_not_installed("apsimx")

  v <- PESTO:::.capture_apsim_version()
  expect_false(is.na(v))

  # The independent oracle: ask the binary ourselves and require agreement.
  raw <- suppressWarnings(
    system2(exe, "--version", stdout = TRUE, stderr = FALSE)
  )
  raw <- raw[nzchar(raw)]
  expect_gt(length(raw), 0L)
  expect_equal(v, sub("^APSIM[[:space:]]+", "", raw[[length(raw)]]))

  # And it must be a version, not an error transcript -- the failure mode the
  # PEST++ probe actually had.
  expect_false(grepl("error|exception|not found", v, ignore.case = TRUE))
})

test_that("apsim_callback stamps the real engine version on the closure", {
  # The provenance claim: a manifest built from this closure should be able to
  # say which simulator produced the numbers. NA here means the run is
  # unattributable.
  exe <- .apsim_models_exe()
  skip_if(!nzchar(exe), "no APSIM engine (set APSIM_EXE_PATH)")
  skip_if_not_installed("apsimx")
  apsimx::apsimx_options(exe.path = exe)

  tmpl <- file.path(tempdir(), "probe.apsimx")
  writeLines("{}", tmpl)
  on.exit(unlink(tmpl), add = TRUE)

  fm <- apsim_callback(
    template         = tmpl,
    param_map        = list(k = "x"),
    output_extractor = function(out) 1,
    # Stubs here are correct: this test is about the version stamp, not the
    # run. The engine is still the one that answered.
    param_writer      = function(...) invisible(NULL),
    simulation_runner = function(...) data.frame(y = 1)
  )

  ver <- attr(fm, "apsim_version")
  expect_false(is.na(ver))
  expect_match(ver, "^[0-9]{4}\\.")
})
