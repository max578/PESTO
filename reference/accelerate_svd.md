# Hardware-Accelerated SVD via LAPACK

Uses R's linked LAPACK (which on macOS is Apple Accelerate/AMX, on Linux
is typically OpenBLAS or MKL) for hardware-optimised SVD computation.

## Usage

``` r
accelerate_svd(A, thin = TRUE)
```

## Arguments

- A:

  Matrix (m x n). Input matrix.

- thin:

  Logical. If TRUE (default), compute thin SVD.

## Value

A list with components `u`, `d` and `v`, so that
`A = u %*% diag(d) %*% t(v)`. With `thin = TRUE`, `u` is m x k and `v`
is n x k, where k = min(m, n). With `thin = FALSE`, `u` is m x m and `v`
is n x n, and only their first k columns pair with `d`.

## Examples

``` r
set.seed(1L)
A <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
res <- accelerate_svd(A, thin = TRUE)
length(res$d)
#> [1] 5
all.equal(res$d, svd(A)$d)
#> [1] TRUE
max(abs(A - res$u %*% diag(res$d) %*% t(res$v)))
#> [1] 7.827072e-15
```
