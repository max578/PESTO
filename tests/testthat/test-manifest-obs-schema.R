# Tests for the obs_schema grounded semantic descriptor (schema 1.1.0).
# Phase-0 FX-1 of the orchestra hallucination-leakage fix: column meaning
# becomes a machine-checkable field instead of an out-of-band convention.

.fixture_manifest <- function(obs_schema = NULL) {
  set.seed(1)
  npar <- 2L; nobs <- 2L; nreal <- 20L
  G <- matrix(stats::rnorm(nobs * npar), nobs, npar)
  y <- as.numeric(G %*% c(0.5, -1.0)) + stats::rnorm(nobs, sd = 0.05)
  names(y) <- c("yield", "grain_n")
  prior <- matrix(stats::rnorm(nreal * npar), nreal, npar,
                  dimnames = list(NULL, c("rue", "cn2")))
  fit <- pesto_ies_callback(function(t) t %*% t(G), prior, y, 0.05,
                            noptmax = 2, verbose = FALSE)
  as_manifest(fit, seed = 1L, obs_schema = obs_schema)
}

test_that("default schema_version is 1.1.0", {
  m <- .fixture_manifest()
  expect_identical(m@schema_version, "1.1.0")
  expect_null(m@obs_schema)
})

test_that("pesto_obs_schema() builds a valid descriptor and defaults provenance", {
  s <- pesto_obs_schema(
    outputs = data.frame(name = c("yield", "grain_n"),
                         quantity = c("grain_yield", "grain_nitrogen"),
                         unit = c("t/ha", "kg/ha"), stringsAsFactors = FALSE),
    params  = data.frame(name = c("rue", "cn2"),
                         apsim_node = c("[Leaf].Photosynthesis.RUE", "CN2Bare"),
                         unit = c("g/MJ", "unitless"), stringsAsFactors = FALSE)
  )
  expect_named(s, c("outputs", "params"))
  expect_true(all(is.na(s$outputs$verified_on)))      # unverified by default
  expect_s3_class(s$outputs$verified_on, "Date")
  expect_true("oracle_kind" %in% names(s$outputs))
})

test_that("pesto_obs_schema() errors on a missing required column", {
  expect_error(
    pesto_obs_schema(outputs = data.frame(name = "yield", unit = "t/ha")),
    "missing required column"
  )
})

test_that("a valid obs_schema constructs and is carried on the manifest", {
  s <- pesto_obs_schema(
    outputs = data.frame(name = c("yield", "grain_n"),
                         quantity = c("grain_yield", "grain_nitrogen"),
                         unit = c("t/ha", "kg/ha"), stringsAsFactors = FALSE))
  m <- .fixture_manifest(obs_schema = s)
  expect_equal(m@obs_schema$outputs$unit, c("t/ha", "kg/ha"))
})

test_that("a schema naming a non-existent column is rejected (hallucinated column)", {
  bad <- pesto_obs_schema(
    outputs = data.frame(name = c("yield", "evapotranspiration"),
                         quantity = c("grain_yield", "et"),
                         unit = c("t/ha", "mm"), stringsAsFactors = FALSE))
  expect_error(.fixture_manifest(obs_schema = bad),
               "absent from the manifest")
})

test_that("obs_schema is provenance metadata — not folded into data_hash", {
  s <- pesto_obs_schema(
    outputs = data.frame(name = "yield", quantity = "grain_yield",
                         unit = "t/ha", stringsAsFactors = FALSE))
  m0 <- .fixture_manifest()
  m1 <- .fixture_manifest(obs_schema = s)
  # Same data, schema present vs absent -> identical integrity hash.
  expect_identical(m0@data_hash, m1@data_hash)
  expect_true(verify_manifest(m1)$ok)
})

test_that("obs_schema round-trips through YAML write/read", {
  s <- pesto_obs_schema(
    outputs = data.frame(name = c("yield", "grain_n"),
                         quantity = c("grain_yield", "grain_nitrogen"),
                         unit = c("t/ha", "kg/ha"),
                         verified_on = as.Date(c("2026-06-09", NA)),
                         stringsAsFactors = FALSE),
    params  = data.frame(name = c("rue", "cn2"),
                         apsim_node = c("[Leaf].Photosynthesis.RUE", "CN2Bare"),
                         unit = c("g/MJ", "unitless"), stringsAsFactors = FALSE))
  m <- .fixture_manifest(obs_schema = s)
  tmp <- file.path(tempdir(), "obs_schema_rt.yaml")
  write_manifest(m, tmp, format = "rds")
  back <- read_manifest(tmp)
  expect_equal(back@obs_schema$outputs$unit, c("t/ha", "kg/ha"))
  expect_equal(back@obs_schema$outputs$quantity,
               c("grain_yield", "grain_nitrogen"))
  expect_equal(back@obs_schema$outputs$verified_on,
               as.Date(c("2026-06-09", NA)))
  expect_equal(back@obs_schema$params$apsim_node,
               c("[Leaf].Photosynthesis.RUE", "CN2Bare"))
})
