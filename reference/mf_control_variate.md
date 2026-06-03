# Affine control-variate bias correction across fidelities

Debiases a cheap-fidelity output ensemble against a sparse expensive
sample, per observation dimension. For each observation `j` it fits the
first-order autoregressive control variate (the linear term of the
Kennedy-O'Hagan multi-fidelity model) `high_j ~ a_j + b_j * low_j` on
the paired subset, with `b_j = cov(high_j, low_j) / var(low_j)` the
variance-minimising coefficient, then predicts the corrected
high-fidelity output for every realisation as `a_j + b_j * low_all_j`.

## Usage

``` r
mf_control_variate(low_all, high_sub, low_sub)
```

## Arguments

- low_all:

  Numeric matrix, `nreal x nobs`. Cheap-fidelity output for every
  realisation.

- high_sub:

  Numeric matrix, `nsub x nobs`. Expensive-fidelity output for the
  paired subset.

- low_sub:

  Numeric matrix, `nsub x nobs`. Cheap-fidelity output for the same
  subset, row-aligned with `high_sub`.

## Value

A `nreal x nobs` matrix of bias-corrected outputs, with attributes
`"intercept"` (`a_j`), `"slope"` (`b_j`), and `"subset_cor"`
(per-dimension subset correlation; `NA` for a degenerate dimension).

## Details

The estimator degrades gracefully: where the cheap output has zero
variance on the subset it falls back to the expensive subset mean
(`b_j = 0`), and where the two fidelities are weakly correlated the
correction shrinks toward that mean rather than amplifying noise.

This is the plug-in primitive for surrogate cascades: a cascade runs the
cheap level over the full ensemble, the expensive level over a chosen
subset, and calls this to lift the cheap ensemble toward the expensive
one at a fraction of the cost.

## References

Kennedy, M. C. & O'Hagan, A. (2000). Predicting the output from a
complex computer code when fast approximations are available.
*Biometrika*, 87(1), 1–13.

## See also

[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md).

## Examples

``` r
set.seed(1L)
low_all  <- matrix(rnorm(40L), ncol = 2L)
sub      <- 1:5
low_sub  <- low_all[sub, , drop = FALSE]
high_sub <- 0.3 + 1.2 * low_sub + matrix(rnorm(10L, sd = 0.01), ncol = 2L)
corrected <- mf_control_variate(low_all, high_sub, low_sub)
attr(corrected, "slope")
#> [1] 1.200315 1.200632
```
