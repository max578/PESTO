# Adaptive SVD with Automatic Backend Selection

Automatically selects the optimal SVD algorithm based on matrix
dimensions and available hardware:

## Usage

``` r
adaptive_svd(A, k = 0L, method = "auto")
```

## Arguments

- A:

  Matrix (m x n). Input matrix.

- k:

  Integer. Target rank. If 0 (default), computes full SVD.

- method:

  Character. Force a specific method: "auto" (default), "rsvd",
  "accelerate", "eigen", "cuda".

## Value

A list with components U (m x k), d (k), V (n x k), plus `method_used`
and `time_ms`.

## Details

- **Randomised SVD** (Halko et al., 2011): when target rank k \<\<
  min(m,n)

- **Apple Accelerate**: on macOS, for full SVD leveraging AMX
  coprocessor

- **Eigen BDCSVD**: cross-platform fallback with divide-and-conquer

- **CUDA cuSOLVER**: on GPU-equipped systems (when compiled with
  PESTO_USE_CUDA)

## Examples

``` r
set.seed(1L)
A <- matrix(rnorm(20 * 12), nrow = 20, ncol = 12)
res <- adaptive_svd(A, k = 5L, method = "auto")
length(res$d)
#> [1] 5
is.character(res$method_used)
#> [1] TRUE
```
