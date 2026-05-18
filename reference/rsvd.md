# Randomised SVD (Halko-Martinsson-Tropp Algorithm)

Computes a rank-k approximation to the SVD using randomised projections.
This is asymptotically faster than full SVD for problems where the
target rank k is much smaller than min(m,n).

## Usage

``` r
rsvd(A, k, p = 10L, q = 2L)
```

## Arguments

- A:

  Matrix (m x n). Input matrix.

- k:

  Integer. Target rank (number of singular values to compute).

- p:

  Integer. Oversampling parameter (default 10). Higher = more accurate.

- q:

  Integer. Number of power iterations (default 2). Higher = better for
  matrices with slowly decaying singular values.

## Value

A list with components U (m x k), d (k), V (n x k).

## Details

Complexity: O(mn*k) vs O(mn*min(m,n)) for full SVD.

## References

Halko, N., Martinsson, P.G., & Tropp, J.A. (2011). Finding structure
with randomness: Probabilistic algorithms for constructing approximate
matrix decompositions. SIAM Review, 53(2), 217-288.

## Examples

``` r
set.seed(1L)
A <- matrix(rnorm(10 * 6), nrow = 10, ncol = 6)
res <- rsvd(A, k = 3L)
length(res$d)
#> [1] 3
A_hat <- res$u %*% diag(res$d) %*% t(res$v)
mean((A - A_hat)^2)
#> [1] 0.1065446
```
