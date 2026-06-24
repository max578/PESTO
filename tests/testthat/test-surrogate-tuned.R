# Tests for the MLE-tuned (optionally anisotropic) GP surrogate.

test_that("train_gp_surrogate_tuned returns a tuned, predictable surrogate", {
  set.seed(1L)
  X <- matrix(runif(40L * 2L), 40L, 2L)
  y <- sin(3 * X[, 1]) + 0.5 * X[, 2]^2
  gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))

  expect_true(is.list(gp$tuning))
  expect_true(gp$tuning$anisotropic)
  expect_length(gp$tuning$length_scale, 2L)      # one per dimension
  expect_length(gp$tuning$input_scale, 2L)
  expect_true(all(gp$tuning$input_scale > 0))

  # Predicts the training response back to high accuracy (smooth function).
  pred <- predict_gp_surrogate_tuned(gp, X)
  expect_equal(as.numeric(pred$mean), y, tolerance = 0.05)
  # The tuner improved on the single median-heuristic baseline.
  expect_gt(gp$tuning$log_marginal_likelihood,
            gp$tuning$log_marginal_likelihood_median)
})

test_that("anisotropic tuning beats a single length scale on an anisotropic function", {
  set.seed(2L)
  # Strongly different scales per axis: fast in x1, slow in x2.
  X <- matrix(runif(50L * 2L), 50L, 2L)
  f <- function(M) sin(8 * M[, 1]) + 0.3 * M[, 2]
  y <- f(X)
  Xte <- matrix(runif(200L * 2L), 200L, 2L)
  yte <- f(Xte)
  rmse <- function(a, b) sqrt(mean((a - b)^2))

  gp_iso <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L), anisotropic = FALSE)
  gp_ani <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
  r_iso <- rmse(predict_gp_surrogate_tuned(gp_iso, Xte)$mean[, 1], yte)
  r_ani <- rmse(predict_gp_surrogate_tuned(gp_ani, Xte)$mean[, 1], yte)

  expect_false(gp_iso$tuning$anisotropic)
  expect_lt(r_ani, r_iso)            # anisotropy helps on an anisotropic target
})

test_that("the isotropic path centres the response and predicts on the original scale", {
  set.seed(3L)
  X <- matrix(runif(30L), 30L, 1L)             # 1-D -> isotropic by construction
  y <- 100 + 5 * X[, 1]                         # large non-zero mean
  gp <- train_gp_surrogate_tuned(X, matrix(y, ncol = 1L))
  expect_false(gp$tuning$anisotropic)           # single dimension
  expect_equal(gp$tuning$y_mean, mean(y), tolerance = 1e-8)
  pred <- predict_gp_surrogate_tuned(gp, X)
  expect_equal(as.numeric(pred$mean), y, tolerance = 0.5)   # mean added back
})

test_that("input validation and the non-tuned guard fire", {
  expect_error(train_gp_surrogate_tuned(matrix(1:4, 2L, 2L), matrix(1:2, 1L)),
               "same number of rows")
  expect_error(train_gp_surrogate_tuned(matrix(1:2, 2L, 1L), matrix(1:2, 2L)),
               "at least 3")
  plain <- train_gp_surrogate(matrix(runif(20L), 10L, 2L),
                              matrix(runif(10L), ncol = 1L))
  expect_error(predict_gp_surrogate_tuned(plain, matrix(runif(4L), 2L, 2L)),
               "not a tuned surrogate")
})
