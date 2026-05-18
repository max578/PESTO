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

A list with components U, d, V.

## Examples

``` r
set.seed(1L)
A <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
res <- accelerate_svd(A, thin = TRUE)
length(res$d)
#> [1] 5
all.equal(sort(res$d, decreasing = TRUE), svd(A)$d)
#> [1] TRUE
```
