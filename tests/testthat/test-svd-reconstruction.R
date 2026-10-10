test_that("accelerate_svd reconstructs tall, wide and square matrices", {
  set.seed(11L)
  for (dims in list(c(8L, 5L), c(5L, 8L), c(6L, 6L))) {
    A <- matrix(rnorm(prod(dims)), dims[1], dims[2])
    k <- min(dims)
    thin <- accelerate_svd(A, thin = TRUE)
    expect_equal(dim(thin$u), c(dims[1], k))
    expect_equal(dim(thin$v), c(dims[2], k))
    expect_equal(thin$u %*% diag(thin$d) %*% t(thin$v), A, tolerance = 1e-10)

    full <- accelerate_svd(A, thin = FALSE)
    expect_equal(dim(full$u), c(dims[1], dims[1]))
    expect_equal(dim(full$v), c(dims[2], dims[2]))
    expect_equal(crossprod(full$u), diag(dims[1]), tolerance = 1e-10)
    expect_equal(crossprod(full$v), diag(dims[2]), tolerance = 1e-10)
    expect_equal(full$u[, seq_len(k)] %*% diag(full$d) %*% t(full$v[, seq_len(k)]),
                 A, tolerance = 1e-10)
  }
})

test_that("the LAPACK and Eigen paths of the adaptive solver agree", {
  set.seed(12L)
  np <- 20L; no <- 40L; nr <- 30L
  pd <- matrix(rnorm(np * nr), np, nr)
  od <- matrix(rnorm(no * nr), no, nr)
  or_ <- matrix(rnorm(no * nr), no, nr)
  Am <- matrix(rnorm(np * (nr - 1L)), np, nr - 1L)

  sv_lapack <- adaptive_svd(od, method = "accelerate")
  expect_equal(sv_lapack$u %*% diag(sv_lapack$d) %*% t(sv_lapack$v), od,
               tolerance = 1e-10)

  solve_with <- function(method) {
    ensemble_solution_adaptive(pd, od, or_, pd, rep(1, no), rep(1, np), Am,
                               cur_lam = 1.0, svd_method = method)$upgrade
  }
  expect_equal(solve_with("accelerate"), solve_with("eigen"), tolerance = 1e-8)
})
