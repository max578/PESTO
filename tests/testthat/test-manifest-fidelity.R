# Fidelity provenance: pesto_ies_callback() records the realised
# fidelity schedule for a multi-fidelity run, as_manifest() inherits it,
# and write/read round-trips the structured record faithfully (closing
# the C2 manifest lineage).

mf_fit <- function(noptmax = 4L, schedule = c(0L, 0L, 1L, 1L)) {
  set.seed(31L)
  G <- matrix(stats::rnorm(8L * 3L), 8L, 3L)
  truth  <- function(theta) theta %*% t(G)
  biased <- function(theta) theta %*% t(G) + 0.4
  mf <- pesto_multifidelity_model(
    levels = list(
      pesto_forward_model(fn = biased, n_obs = 8L, fidelity = 0L),
      pesto_forward_model(fn = truth,  n_obs = 8L, fidelity = 1L)
    ),
    costs = c(1, 20)
  )
  prior <- matrix(stats::rnorm(30L * 3L), 30L, 3L,
                  dimnames = list(NULL, paste0("p", 1:3)))
  y <- as.numeric(truth(matrix(c(1, -0.5, 2), nrow = 1L)))
  pesto_ies_callback(
    mf, prior, stats::setNames(y, paste0("o", 1:8)), obs_sd = 0.1,
    noptmax = noptmax, fidelity_schedule = schedule, verbose = FALSE
  )
}

single_fit <- function() {
  set.seed(31L)
  G <- matrix(stats::rnorm(8L * 3L), 8L, 3L)
  f <- function(theta) theta %*% t(G)
  prior <- matrix(stats::rnorm(30L * 3L), 30L, 3L,
                  dimnames = list(NULL, paste0("p", 1:3)))
  y <- as.numeric(f(matrix(c(1, -0.5, 2), nrow = 1L)))
  pesto_ies_callback(f, prior, stats::setNames(y, paste0("o", 1:8)),
                     obs_sd = 0.1, noptmax = 4L, verbose = FALSE)
}

test_that("a multi-fidelity run records its realised schedule", {
  fit <- mf_fit()
  fd <- fit$fidelity
  expect_type(fd, "list")
  expect_identical(fd$type, "multifidelity")
  expect_identical(fd$schedule, c(0L, 0L, 1L, 1L))
  expect_identical(fd$final_level, 1L)
  expect_identical(fd$n_levels, 2L)
  expect_equal(fd$costs, c(1, 20))
})

test_that("a single-fidelity run records no fidelity provenance", {
  expect_null(single_fit()$fidelity)
})

test_that("as_manifest inherits the recorded fidelity provenance", {
  m <- as_manifest(mf_fit(), run_id = "mf_run")
  expect_identical(m@fidelity$type, "multifidelity")
  expect_identical(m@fidelity$schedule, c(0L, 0L, 1L, 1L))
  expect_null(as_manifest(single_fit(), run_id = "single_run")@fidelity)
})

test_that("an explicit fidelity argument overrides the recorded one", {
  m <- as_manifest(mf_fit(), run_id = "mf_run",
                   fidelity = c(custom = 1))
  expect_equal(m@fidelity, c(custom = 1))
})

test_that("write/read round-trips the structured fidelity record", {
  m <- as_manifest(mf_fit(), run_id = "mf_rt")
  dir <- tempfile("mf_manifest_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  yaml_path <- file.path(dir, "manifest.yaml")
  write_manifest(m, yaml_path)               # default rds sidecars
  back <- read_manifest(yaml_path)

  expect_identical(back@fidelity$type, "multifidelity")
  expect_identical(back@fidelity$schedule, c(0L, 0L, 1L, 1L))
  expect_identical(back@fidelity$final_level, 1L)
  expect_identical(back@fidelity$n_levels, 2L)
  expect_equal(back@fidelity$costs, c(1, 20))
  # Integrity (params/outputs/weights/obs_target/seed) is unaffected.
  expect_true(isTRUE(verify_manifest(back)$ok))
})
