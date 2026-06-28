# Coverage for exports that the "comparison and simulation" vignette used to be
# the sole exerciser of (its Scenario C, a coverage-driving Monte-Carlo study,
# was removed when the vignette was reduced to a publication artefact). These
# tests move that coverage into the suite, where it belongs.

test_that("ensemble CSV round-trips through write_ensemble / read_ensemble", {
  ens <- data.table::data.table(
    real_name = sprintf("r%02d", seq_len(8)),
    k1 = rnorm(8), k2 = rnorm(8)
  )
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)

  expect_null(write_ensemble(ens, tf))
  back <- read_ensemble(tf)
  expect_equal(nrow(back), 8L)
  expect_true(all(c("k1", "k2") %in% names(back)))
  expect_equal(back$k1, ens$k1, tolerance = 1e-8)
})

test_that("SVD backends agree on the leading singular value", {
  set.seed(1L)
  A <- matrix(rnorm(120 * 40), 120, 40)

  auto <- adaptive_svd(A, k = 10L, method = "auto")
  expect_true(nzchar(auto$method_used))
  expect_gte(length(auto$d), 1L)

  acc <- accelerate_svd(A, thin = TRUE)
  expect_gte(length(acc$d), 1L)

  rs <- rsvd(A, k = 10L)
  expect_equal(length(rs$d), 10L)
  # Randomised and dense decompositions must agree on the top singular value.
  expect_lt(abs(rs$d[1] - acc$d[1]) / acc$d[1], 0.1)
})

test_that("ensemble_solution_adaptive returns an upgrade and timing diagnostics", {
  set.seed(2L)
  np <- 20L; no <- 40L; nr <- 30L
  pd <- matrix(rnorm(np * nr), np, nr)
  od <- matrix(rnorm(no * nr), no, nr)
  or_ <- matrix(rnorm(no * nr), no, nr)
  Am <- matrix(rnorm(np * (nr - 1L)), np, nr - 1L)

  g <- ensemble_solution_adaptive(pd, od, or_, pd, rep(1, no), rep(1, np), Am,
                                  cur_lam = 1.0, svd_method = "auto")
  expect_equal(sort(dim(g$upgrade)), sort(c(np, nr)))
  expect_true(is.numeric(g$total_time_ms))
  expect_true(nzchar(g$svd_method))

  # The former name is a deprecated alias that warns and forwards.
  expect_warning(
    g2 <- ensemble_solution_gpu(pd, od, or_, pd, rep(1, no), rep(1, np), Am,
                                cur_lam = 1.0, svd_method = "auto"),
    "deprecated"
  )
  expect_equal(sort(dim(g2$upgrade)), sort(c(np, nr)))
})

test_that("surrogate_ensemble_update and pesto_surrogate_ies report savings", {
  set.seed(3L)
  np <- 4L; no <- 12L; nr <- 80L
  G <- matrix(rnorm(no * np, sd = 1 / sqrt(np)), no, np)
  fwd <- function(t) as.numeric(G %*% t)
  pe <- matrix(rnorm(nr * np), nr, np)
  oe <- t(apply(pe, 1L, function(p) fwd(p) + rnorm(no, sd = 0.05)))
  yt <- fwd(rnorm(np)) + rnorm(no, sd = 0.05)

  r <- surrogate_ensemble_update(pe, oe, yt, rep(1 / 0.05, no), rep(1, np),
                                 uncertainty_threshold = 0.2)
  expect_gte(r$savings_pct, 0)
  expect_lte(r$savings_pct, 100)
  expect_equal(r$n_model_runs + r$n_surrogate_runs, r$n_total)
  expect_equal(sort(dim(r$upgrade)), sort(c(np, nr)))

  r2 <- pesto_surrogate_ies(pe, oe, yt, rep(1 / 0.05, no), rep(1, np),
                            uncertainty_threshold = 0.2)
  expect_false(is.null(r2$savings_pct))
})

test_that("GP and RFF surrogates predict at the training inputs", {
  set.seed(4L)
  nr <- 60L; np <- 8L; no <- 15L
  X <- matrix(rnorm(nr * np), nr, np)
  Y <- matrix(rnorm(nr * no), nr, no)

  gp <- train_gp_surrogate(X, Y)
  gpp <- predict_gp_surrogate(gp, X)
  expect_equal(dim(gpp$mean), c(nr, no))

  rff <- train_rff_surrogate(X, Y, n_features = 100L)
  rp <- predict_rff_surrogate(rff, X)
  expect_equal(dim(rp$mean), c(nr, no))
})

test_that("adaptive_ensemble_size returns a bounded recommendation", {
  set.seed(5L)
  phi <- rlnorm(50, 3, 0.5)
  sz <- adaptive_ensemble_size(phi, current_size = 50L)
  expect_gte(sz$recommended_size, 20L)
  expect_lte(sz$recommended_size, 500L)
  expect_true(is.finite(sz$ess))
})

test_that("check_surrogate_regime flags an unfavourable regime", {
  expect_true(check_surrogate_regime(n_params = 4, n_train = 100))
  expect_warning(check_surrogate_regime(n_params = 30, n_train = 30))
  expect_false(suppressWarnings(check_surrogate_regime(n_params = 30, n_train = 30)))
})

test_that("pesto_reference_ies returns an upgrade of the parameter shape", {
  set.seed(6L)
  np <- 6L; no <- 12L; nr <- 40L
  G <- matrix(rnorm(no * np), no, np)
  yt <- as.numeric(G %*% rnorm(np)) + rnorm(no, sd = 0.05)
  pe <- matrix(rnorm(np * nr), np, nr)               # n_par x n_real
  oe <- apply(pe, 2L, function(p) as.numeric(G %*% p))

  upgrade <- pesto_reference_ies(pe, oe, yt, rep(1 / 0.05, no))
  expect_equal(dim(upgrade), c(np, nr))
  expect_true(all(is.finite(upgrade)))
})

test_that("plot helpers return ggplot objects, including the high-dim path", {
  skip_if_not_installed("ggplot2")

  phi_dt <- data.table::data.table(
    iteration = 0:3, mean = c(10, 5, 2, 1),
    min = c(8, 4, 1, 0.5), max = c(12, 6, 3, 1.5)
  )
  expect_s3_class(plot_phi(phi_dt), "ggplot")

  ens <- data.table::data.table(k1 = rnorm(30), k2 = rnorm(30), k3 = rnorm(30))
  expect_s3_class(plot_ensemble(ens), "ggplot")

  J <- matrix(rnorm(30 * 8), 30, 8)
  colnames(J) <- paste0("k", seq_len(8))
  expect_s3_class(plot_identifiability(jacobian = J), "ggplot")
  expect_s3_class(plot_identifiability(jacobian = J, top_n = 3L), "ggplot")

  # High-dimensional Jacobian: the default top_n cap must keep it legible.
  J_big <- matrix(rnorm(120 * 60), 120, 60)
  expect_s3_class(plot_identifiability(jacobian = J_big), "ggplot")

  set.seed(8L)
  Gp <- matrix(rnorm(3 * 4, sd = 0.5), 3, 4)
  pe <- matrix(rnorm(60 * 4), 60, 4)
  oe <- t(apply(pe, 1L, function(p) as.numeric(Gp %*% p) + rnorm(3, sd = 0.05)))
  surr <- surrogate_ensemble_update(
    pe, oe, as.numeric(Gp %*% rnorm(4)), rep(1, 3), rep(1, 4),
    uncertainty_threshold = 0.2
  )
  expect_s3_class(plot_surrogate_diagnostics(surr), "ggplot")
})

test_that("external-binary wrappers are exported and callable", {
  # These shell out to PEST++ executables, so a live run is only possible
  # where the binaries are installed; the function objects are always checked.
  expect_true(is.function(pesto_ies))
  expect_true(is.function(pesto_glm))
  expect_true(is.function(pesto_sweep))
  expect_true(is.function(pesto_sensitivity))
})
