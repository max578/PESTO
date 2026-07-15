# test-pestpp-real-binary.R -- the oracle no mock can be: PEST++ itself.
#
# Every other test in this package checks PESTO against PESTO. That is exactly
# how PESTO shipped a PEST++ command line PEST++ has never accepted: the binary
# was not a test dependency, so nothing could contradict the package's belief
# about it. These tests run the real thing.
#
# Opt-in: they need a pestpp binary on PATH (or PESTPP_BIN set), and are
# skipped otherwise -- including on CRAN. A skip is honest here; a mocked
# "pass" would not be, since mocking the authority is precisely the failure
# this file exists to prevent.
#
#   PATH=/path/to/pestpp/bin:$PATH R -e 'devtools::test(filter="real-binary")'
#
# The problem is linear with an exact answer, so "did it work" is a number:
#
#   h1 = 2*k1 + 3*k2 ;  h2 = k1 + k2 ;  obs (13, 6)  ->  k1 = 5, k2 = 1
#
# started at (1, 1), well away from the solution, so an honoured `noptmax` and
# an ignored one give visibly different answers.

.pestpp_exe <- function(name) {
  dir <- Sys.getenv("PESTPP_BIN")
  if (nzchar(dir) && file.exists(file.path(dir, name))) {
    return(file.path(dir, name))
  }
  found <- Sys.which(name)
  if (nzchar(found)) found else ""
}

skip_on_cran()
skip_if(!nzchar(.pestpp_exe("pestpp-glm")),
        "no pestpp binary: put one on PATH or set PESTPP_BIN to run the real-binary oracle")

# A fresh directory per run. Base R rather than withr, to keep the test suite's
# dependency surface as it is.
.real_dir <- function(tag) {
  wd <- file.path(tempfile("pestpp_"), tag)
  dir.create(wd, recursive = TRUE)
  wd
}

# Build a complete, runnable PEST problem: a real model executable, a real
# template, a real instruction file, and a control file written by PESTO.
.real_problem <- function(wd) {
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "#!/bin/sh",
    "k1=$(sed -n '1p' model.in)",
    "k2=$(sed -n '2p' model.in)",
    "awk -v a=\"$k1\" -v b=\"$k2\" 'BEGIN{printf \"%.6f\\n%.6f\\n\", 2*a+3*b, a+b}' > model.out"
  ), file.path(wd, "model.sh"))
  Sys.chmod(file.path(wd, "model.sh"), "0755")

  writeLines(c("ptf $", "$k1            $", "$k2            $"),
             file.path(wd, "model.tpl"))
  # `l1` advances one line; it is not "go to line 1".
  writeLines(c("pif $", "l1 !h1!", "l1 !h2!"), file.path(wd, "model.ins"))

  pst <- create_pest_scenario(
    parameters = data.table::data.table(
      parnme = c("k1", "k2"), partrans = "none", parchglim = "relative",
      parval1 = c(1.0, 1.0), parlbnd = c(-20, -20), parubnd = c(20, 20),
      pargp = "g", scale = 1.0, offset = 0.0, dercom = 1L
    ),
    observations = data.table::data.table(
      obsnme = c("h1", "h2"), obsval = c(13.0, 6.0),
      weight = c(1.0, 1.0), obgnme = "og"
    ),
    model_command = "./model.sh",
    template_files = data.table::data.table(tpl = "model.tpl", inp = "model.in"),
    instruction_files = data.table::data.table(ins = "model.ins", out = "model.out")
  )

  f <- file.path(wd, "m.pst")
  write_pst(pst, f)
  f
}


test_that("PESTO writes a control file PEST++ accepts", {
  # PESTO declared NINSFLE but wrote no instruction lines, so PEST++ refused
  # every control file it produced: "number of instruction files = 0". Nothing
  # local could see it -- the control file is only ever read by the binary.
  wd <- .real_dir("accept")
  pst_file <- .real_problem(wd)

  # Run IN `wd`: the control file is passed by basename, so from anywhere else
  # PEST++ simply would not find it -- and would then report a *different*
  # error, quietly satisfying the assertions below without ever reading the
  # file. Assert it got far enough to matter, not merely that it complained
  # about something else.
  #
  # The directory is restored inside local(), before anything else runs: an
  # on.exit() restore registered alongside a cleanup would fire in registration
  # order, and deleting the directory while it is still the working directory
  # leaves the session pointing at nothing -- which sprayed PEST++ output into
  # tests/testthat/ when a later test resolved a relative path.
  out <- local({
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)
    suppressWarnings(system2(
      .pestpp_exe("pestpp-glm"), basename(pst_file),
      stdout = TRUE, stderr = TRUE
    ))
  })

  expect_true(any(grepl("checking model IO files...done", out, fixed = TRUE)))
  expect_false(any(grepl("control file parsing error", out, fixed = TRUE)))
  expect_false(any(grepl("COMMAND LINE ERROR", out, fixed = TRUE)))
})


test_that("pesto_glm() recovers the true parameters, and noptmax changes them", {
  # The RB25 defect, settled by the authority. `noptmax` was documented,
  # exported, and read by nothing; a unit test can show it now reaches the
  # control file, but only a real run shows it changes the answer.
  one <- pesto_glm(
    .real_problem(.real_dir("n1")),
    exe = .pestpp_exe("pestpp-glm"), noptmax = 1L, verbose = FALSE
  )
  many <- pesto_glm(
    .real_problem(.real_dir("n20")),
    exe = .pestpp_exe("pestpp-glm"), noptmax = 20L, verbose = FALSE
  )

  expect_equal(one$exit_code, 0L)
  expect_equal(many$exit_code, 0L)

  k <- function(r, nm) r$parameters$parval[r$parameters$parnme == nm]

  # 20 iterations reach the exact solution; 1 iteration cannot.
  expect_equal(k(many, "k1"), 5.0, tolerance = 1e-2)
  expect_equal(k(many, "k2"), 1.0, tolerance = 1e-2)
  expect_false(isTRUE(all.equal(k(one, "k1"), 5.0, tolerance = 1e-2)))

  # And the two runs genuinely differ: if `noptmax` were ignored again, both
  # would run the control file's own value and agree.
  expect_false(isTRUE(all.equal(k(one, "k1"), k(many, "k1"))))
})


test_that("pesto_ies() runs, and num_reals sizes the ensemble", {
  # This function could never have worked: every call built at least one
  # `/h :name=value` switch, and PEST++ exited on the command line.
  res <- pesto_ies(
    .real_problem(.real_dir("ies")),
    exe = .pestpp_exe("pestpp-ies"), num_reals = 8L, noptmax = 2L,
    verbose = FALSE
  )

  expect_equal(res$exit_code, 0L)
  expect_false(is.null(res$par_ensemble))
  expect_equal(nrow(res$par_ensemble), 8L)   # ies_num_reals reached the binary
  expect_true(mean(res$par_ensemble$k1) > 3.5)   # moved from 1 toward 5
})


test_that("pesto_sensitivity(method) selects the algorithm PEST++ runs", {
  # `method` was read and stored on the result but never sent, so pestpp-sen
  # ran its Morris default and the output was labelled whatever was asked for.
  # PEST++'s own record file is the oracle for what it actually ran.
  gsa_from_rec <- function(m) {
    wd <- .real_dir(m)
    pst_file <- .real_problem(wd)
    res <- pesto_sensitivity(pst_file, method = m,
                             exe = .pestpp_exe("pestpp-sen"), verbose = FALSE)
    expect_equal(res$exit_code, 0L)
    rec <- file.path(wd, "m_pesto.rec")
    expect_true(file.exists(rec))
    trimws(sub(".*gsa_method:", "",
               grep("gsa_method", readLines(rec), value = TRUE,
                    ignore.case = TRUE)[1]))
  }

  expect_equal(gsa_from_rec("morris"), "MORRIS")
  expect_equal(gsa_from_rec("sobol"), "SOBOL")
})
