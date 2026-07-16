# RB25 -- formals that were documented, exported, and read by nothing.
#
# Each test below asserts that a parameter CHANGES THE OUTPUT. That is the
# point: a dead formal is invisible to any test that never varies it, so the
# suite passed for as long as the parameter was ignored. Asserting the default
# behaviour would re-admit the defect, because at the default an ignored
# parameter and an honoured one agree.
#
# Found by the RB25 dead-formal gate:
#   Rscript ~/.claude/skills/rpkg/scripts/rpkg_rb25_gate.R <pkg-root>


# ---- helper ----------------------------------------------------------------

# Write a minimal PEST++ binary matrix in the layout .read_ensemble_binary()
# reads: n_neg, ncol, column-major doubles, then (len, chars) per row name and
# per column name. `col_names = NULL` writes zero-length names, which is the
# "column names are absent" case `pst` is documented to cover.
.write_test_jco <- function(mat, row_names, col_names = NULL) {
  path <- tempfile(fileext = ".jco")
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)

  writeBin(-nrow(mat), con, size = 4L)
  writeBin(ncol(mat), con, size = 4L)
  writeBin(as.double(as.vector(mat)), con, size = 8L)

  for (nm in row_names) {
    writeBin(nchar(nm), con, size = 4L)
    writeChar(nm, con, nchars = nchar(nm), eos = NULL)
  }
  for (j in seq_len(ncol(mat))) {
    nm <- if (is.null(col_names)) "" else col_names[j]
    writeBin(nchar(nm), con, size = 4L)
    if (nchar(nm) > 0L) writeChar(nm, con, nchars = nchar(nm), eos = NULL)
  }

  path
}


# ---- write_ensemble(format) ------------------------------------------------

test_that("write_ensemble() refuses a format it cannot write", {
  # Was: `format` was never read, so "binary" wrote a CSV and returned quietly.
  # The caller got a CSV named as a binary -- a wrong file, not an error.
  ens <- data.table::data.table(real_name = c("r1", "r2"), k1 = c(1.0, 2.0))
  tf <- tempfile(fileext = ".jcb")
  on.exit(unlink(tf), add = TRUE)

  expect_error(write_ensemble(ens, tf, format = "binary"))
  expect_false(file.exists(tf))
})

test_that("write_ensemble() still writes the format it does support", {
  ens <- data.table::data.table(real_name = c("r1", "r2"), k1 = c(1.25, 2.5))
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)

  expect_silent(write_ensemble(ens, tf, format = "csv"))
  expect_identical(read_ensemble(tf, format = "csv")$k1, c(1.25, 2.5))
})


# ---- .read_ensemble_binary(): unnamed columns ------------------------------

test_that("an unnamed .jco keeps every column", {
  # Was: `set(j = "")` keys on the name, so every blank-named column landed on
  # the same "" column and all but the last were discarded -- silently, with a
  # narrower table than the header declared. Values are distinct per column so
  # a collapse cannot pass by coincidence.
  mat <- matrix(c(1, 2, 3, 10, 20, 30, 100, 200, 300), nrow = 3L, ncol = 3L)
  jco <- .write_test_jco(mat, row_names = c("r1", "r2", "r3"))
  on.exit(unlink(jco), add = TRUE)

  dt <- PESTO:::.read_ensemble_binary(jco)

  expect_identical(ncol(dt), 4L)                      # real_name + 3 columns
  expect_identical(names(dt)[-1], c("p1", "p2", "p3"))
  expect_identical(dt$p1, c(1, 2, 3))
  expect_identical(dt$p3, c(100, 200, 300))
  expect_identical(attr(dt, "unnamed_columns"), rep(TRUE, 3L))
})

test_that("a named .jco keeps the names the file carries", {
  mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3L, ncol = 2L)
  jco <- .write_test_jco(mat, c("r1", "r2", "r3"), col_names = c("kx", "ky"))
  on.exit(unlink(jco), add = TRUE)

  dt <- PESTO:::.read_ensemble_binary(jco)

  expect_identical(names(dt)[-1], c("kx", "ky"))
  expect_identical(attr(dt, "unnamed_columns"), c(FALSE, FALSE))
})


# ---- plot_identifiability(pst) ---------------------------------------------

test_that("pst labels an unnamed .jco, and changes the plot", {
  # Was: `pst` was never read, so the plot showed placeholders whether or not a
  # pst was supplied. Both branches are asserted here -- equality of the two
  # would mean `pst` is inert again.
  set.seed(42L)
  mat <- matrix(rnorm(12L), nrow = 4L, ncol = 3L)
  jco <- .write_test_jco(mat, row_names = paste0("r", 1:4))
  on.exit(unlink(jco), add = TRUE)

  pst <- create_pest_scenario(
    parameters = data.table::data.table(
      parnme = c("alpha", "beta", "gamma"), partrans = "log",
      parchglim = "factor", parval1 = c(1, 1, 1), parlbnd = c(0.1, 0.1, 0.1),
      parubnd = c(10, 10, 10), pargp = "g"
    ),
    observations = data.table::data.table(
      obsnme = c("o1", "o2"), obsval = c(1, 2), weight = 1, obgnme = "og"
    ),
    model_command = "echo run"
  )

  without <- plot_identifiability(jco_file = jco)
  with_pst <- plot_identifiability(jco_file = jco, pst = pst)

  labs_without <- levels(without$data$parameter)
  labs_with <- levels(with_pst$data$parameter)

  expect_setequal(labs_without, c("p1", "p2", "p3"))
  expect_setequal(labs_with, c("alpha", "beta", "gamma"))
  expect_false(identical(labs_without, labs_with))
})

test_that("a pst that does not describe the Jacobian is refused", {
  # Positional labelling is only defensible when the counts agree; otherwise
  # the plot would carry confident wrong names on real numbers.
  set.seed(42L)
  mat <- matrix(rnorm(12L), nrow = 4L, ncol = 3L)
  jco <- .write_test_jco(mat, row_names = paste0("r", 1:4))
  on.exit(unlink(jco), add = TRUE)

  pst <- create_pest_scenario(
    parameters = data.table::data.table(
      parnme = c("alpha", "beta"), partrans = "log", parchglim = "factor",
      parval1 = c(1, 1), parlbnd = c(0.1, 0.1), parubnd = c(10, 10),
      pargp = "g"
    ),
    observations = data.table::data.table(
      obsnme = c("o1", "o2"), obsval = c(1, 2), weight = 1, obgnme = "og"
    ),
    model_command = "echo run"
  )

  expect_error(
    plot_identifiability(jco_file = jco, pst = pst),
    "cannot be matched"
  )
})

test_that("pst does not override the names a .jco carries", {
  # The documented contract is a fallback, not an override: a file that names
  # its own columns is the authority on their order.
  set.seed(42L)
  mat <- matrix(rnorm(12L), nrow = 4L, ncol = 3L)
  jco <- .write_test_jco(mat, paste0("r", 1:4), col_names = c("kx", "ky", "kz"))
  on.exit(unlink(jco), add = TRUE)

  pst <- create_pest_scenario(
    parameters = data.table::data.table(
      parnme = c("alpha", "beta", "gamma"), partrans = "log",
      parchglim = "factor", parval1 = c(1, 1, 1), parlbnd = c(0.1, 0.1, 0.1),
      parubnd = c(10, 10, 10), pargp = "g"
    ),
    observations = data.table::data.table(
      obsnme = c("o1", "o2"), obsval = c(1, 2), weight = 1, obgnme = "og"
    ),
    model_command = "echo run"
  )

  p <- plot_identifiability(jco_file = jco, pst = pst)

  expect_setequal(levels(p$data$parameter), c("kx", "ky", "kz"))
})
