# Tests for the pesto_ensemble_manifest S7 class and its YAML+CSV
# round-trip. Year-1 §A5 of the UQ ag-stack roadmap.

make_fit <- function(seed = 1L, nreal = 30L, npar = 3L, nobs = 5L,
                     noptmax = 3L) {
  set.seed(seed)
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  theta_true <- stats::rnorm(npar)
  y <- as.numeric(G %*% theta_true) + stats::rnorm(nobs, sd = 0.05)
  names(y) <- paste0("o", seq_len(nobs))
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, paste0("p", seq_len(npar))))
  pesto_ies_callback(
    forward_model  = function(t) t %*% t(G),
    prior_ensemble = prior,
    obs            = y,
    obs_sd         = 0.05,
    noptmax        = noptmax,
    verbose        = FALSE
  )
}

test_that("as_manifest constructs a valid S7 object with consistent hash", {
  fit <- make_fit()
  m <- as_manifest(fit, seed = 1L)

  expect_s7_class <- function(object, class_def) {
    expect_true(S7::S7_inherits(object, class_def))
  }
  expect_s7_class(m, pesto_ensemble_manifest)
  expect_true(grepl("^sha256:", m@data_hash))
  expect_equal(m@method, "ies_callback")
  expect_equal(m@noptmax, 3L)
  expect_equal(nrow(m@params), nrow(fit$par_ensemble))
  expect_equal(nrow(m@outputs), nrow(fit$obs_ensemble))

  v <- verify_manifest(m)
  expect_true(v$ok)
  expect_identical(v$stored, v$recomputed)
})

test_that("write_manifest + read_manifest round-trip preserves hash", {
  fit <- make_fit()
  m1  <- as_manifest(fit, seed = 1L, run_id = "test_roundtrip_001")

  dir <- tempfile("pesto_manifest_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  paths <- write_manifest(m1, file.path(dir, "run.yaml"))
  expect_length(paths, 4L)
  expect_true(all(file.exists(paths)))
  expect_true(file.exists(file.path(dir, "run_assim.rds")))

  m2 <- read_manifest(file.path(dir, "run.yaml"))

  # All slots round-trip
  expect_equal(m2@run_id,          m1@run_id)
  expect_equal(m2@method,          m1@method)
  expect_equal(m2@noptmax,         m1@noptmax)
  expect_equal(m2@seed,            m1@seed)
  expect_equal(m2@pesto_version,   m1@pesto_version)
  expect_equal(m2@data_hash,       m1@data_hash)
  expect_equal(as.numeric(m2@weights),    as.numeric(m1@weights))
  expect_equal(as.numeric(m2@obs_target), as.numeric(m1@obs_target))
  expect_equal(as.numeric(m2@lambda_schedule),
               as.numeric(m1@lambda_schedule))

  # Hash recomputes correctly on the reloaded object
  expect_true(verify_manifest(m2)$ok)
})

test_that("verify_manifest detects post-write tampering with the CSV", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "test_tamper_001")

  dir <- tempfile("pesto_manifest_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  yaml_path <- file.path(dir, "run.yaml")
  write_manifest(m, yaml_path)

  # Tamper with the outputs RDS sidecar by perturbing one value
  outputs_rds <- file.path(dir, "run_outputs.rds")
  df <- readRDS(outputs_rds)
  num_cols <- setdiff(names(df), "real_name")
  df[1L, num_cols[1L]] <- df[1L, num_cols[1L]] + 1e-3
  saveRDS(df, outputs_rds, version = 3L)

  m2 <- read_manifest(yaml_path)
  v  <- verify_manifest(m2)
  expect_false(v$ok)
  expect_false(identical(v$stored, v$recomputed))
})

test_that("constructor validator rejects malformed manifests", {
  expect_error(
    pesto_ensemble_manifest(
      run_id          = "",
      params          = data.frame(real_name = "r1", p1 = 0),
      outputs         = data.frame(real_name = "r1", o1 = 0),
      weights         = 1, obs_target = 0,
      data_hash       = "sha256:x", pesto_version = "0.0",
      timestamp       = Sys.time(),
      method          = "ies_callback",
      noptmax         = 1L, lambda_schedule = 1
    ),
    "run_id"
  )
  expect_error(
    pesto_ensemble_manifest(
      run_id          = "ok",
      params          = data.frame(real_name = c("r1", "r2"),
                                   p1 = 0:1),
      outputs         = data.frame(real_name = "r1", o1 = 0),
      weights         = 1, obs_target = 0,
      data_hash       = "sha256:x", pesto_version = "0.0",
      timestamp       = Sys.time(),
      method          = "ies_callback",
      noptmax         = 1L, lambda_schedule = 1
    ),
    "equal row counts"
  )
  expect_error(
    pesto_ensemble_manifest(
      run_id          = "ok",
      params          = data.frame(real_name = "r1", p1 = 0),
      outputs         = data.frame(real_name = "r1", o1 = 0),
      weights         = 1, obs_target = 0,
      data_hash       = "sha256:x", pesto_version = "0.0",
      timestamp       = Sys.time(),
      method          = "BOGUS",
      noptmax         = 1L, lambda_schedule = 1
    ),
    "method"
  )
})

test_that("print method emits the expected header lines", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "print_test_001")
  out <- utils::capture.output(print(m))
  expect_true(any(grepl("pesto_ensemble_manifest", out)))
  expect_true(any(grepl("print_test_001", out)))
  expect_true(any(grepl("ies_callback", out)))
  expect_true(any(grepl("data hash", out)))
})

test_that("write_manifest errors when directory does not exist", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L)
  expect_error(
    write_manifest(m, "/nonexistent_dir_xyz/run.yaml"),
    "Directory does not exist"
  )
})

test_that("read_manifest errors when file is missing", {
  expect_error(read_manifest("/nonexistent_file.yaml"), "not found")
})

# ---- format = "rds" / "both" / "csv" -----------------------------------

test_that("write_manifest format='rds' is the bit-exact default", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "fmt_rds_001")
  dir <- tempfile("pesto_manifest_rds_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  paths <- write_manifest(m, file.path(dir, "run.yaml"))  # default
  # YAML + 3 RDS, no CSV
  expect_length(paths, 4L)
  expect_true(all(file.exists(paths)))
  expect_false(any(grepl("\\.csv$", paths)))

  m2 <- read_manifest(file.path(dir, "run.yaml"))
  expect_equal(m2@format, "rds")
  v <- verify_manifest(m2)
  expect_true(isTRUE(v$ok))
  expect_null(v$message)
})

test_that("write_manifest format='both' adds inspection CSVs alongside RDS", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "fmt_both_001")
  dir <- tempfile("pesto_manifest_both_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  paths <- write_manifest(m, file.path(dir, "run.yaml"), format = "both")
  # YAML + 3 RDS + 3 inspection CSV
  expect_length(paths, 7L)
  expect_true(all(file.exists(paths)))
  expect_equal(sum(grepl("_inspection\\.csv$", paths)), 3L)
  expect_equal(sum(grepl("\\.rds$", paths)),            3L)

  # YAML records the inspection CSVs in their own block
  yl <- yaml::read_yaml(file.path(dir, "run.yaml"))
  expect_equal(yl$format, "both")
  expect_named(yl$inspection_csv, c("params", "outputs", "assim"),
               ignore.order = TRUE)

  # read_manifest reads from RDS (the hash-bearing files), so the
  # integrity contract holds bit-exactly.
  m2 <- read_manifest(file.path(dir, "run.yaml"))
  expect_equal(m2@format, "both")
  v <- verify_manifest(m2)
  expect_true(isTRUE(v$ok))
  expect_identical(v$stored, v$recomputed)
})

test_that("format='csv_unverified' writes CSV-only sidecars with integrity flag", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "fmt_csv_001")
  dir <- tempfile("pesto_manifest_csv_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  paths <- write_manifest(m, file.path(dir, "run.yaml"),
                          format = "csv_unverified")
  # YAML + 3 CSV, no RDS
  expect_length(paths, 4L)
  expect_true(all(file.exists(paths)))
  expect_false(any(grepl("\\.rds$", paths)))
  expect_equal(sum(grepl("\\.csv$", paths)), 3L)

  yl <- yaml::read_yaml(file.path(dir, "run.yaml"))
  expect_equal(yl$format, "csv_unverified")
  expect_equal(yl$integrity, "not_verifiable")
  expect_match(yl$artefacts$params,  "\\.csv$")
  expect_match(yl$artefacts$outputs, "\\.csv$")
  expect_match(yl$artefacts$assim,   "\\.csv$")

  m2 <- read_manifest(file.path(dir, "run.yaml"))
  expect_equal(m2@format, "csv_unverified")

  v <- verify_manifest(m2)
  expect_true(is.na(v$ok))
  expect_match(v$message, "not_verifiable")
})

test_that("legacy format='csv' is accepted with a deprecation warning", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L, run_id = "fmt_legacy_csv_001")
  dir <- tempfile("pesto_manifest_legacy_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  expect_warning(
    paths <- write_manifest(m, file.path(dir, "run.yaml"),
                            format = "csv"),
    "csv_unverified"
  )
  # Persisted form uses the new spelling.
  yl <- yaml::read_yaml(file.path(dir, "run.yaml"))
  expect_equal(yl$format, "csv_unverified")
  expect_equal(yl$integrity, "not_verifiable")
})

test_that("verifiable modes record integrity=verifiable in YAML", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L)
  dir <- tempfile("pesto_manifest_integ_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  write_manifest(m, file.path(dir, "rds.yaml"))
  write_manifest(m, file.path(dir, "both.yaml"), format = "both")
  expect_equal(
    yaml::read_yaml(file.path(dir, "rds.yaml"))$integrity,
    "verifiable"
  )
  expect_equal(
    yaml::read_yaml(file.path(dir, "both.yaml"))$integrity,
    "verifiable"
  )
})

test_that("write_manifest rejects unknown format", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L)
  dir <- tempfile("pesto_manifest_bad_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  expect_error(
    suppressWarnings(
      write_manifest(m, file.path(dir, "run.yaml"), format = "parquet")
    ),
    "should be one of"
  )
})

test_that("manifest constructor rejects invalid format slot value", {
  fit <- make_fit()
  m   <- as_manifest(fit, seed = 1L)
  expect_error(
    {
      m_bad <- m
      m_bad@format <- "parquet"
    },
    "format"
  )
})
