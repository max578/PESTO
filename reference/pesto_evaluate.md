# Evaluate a PESTO forward model

Runs the forward model on a parameter matrix under its own failure
policy and concurrency strategy, returning a shape-guaranteed
`nreal x nobs` observation matrix. Failed realisations populate `NA`
rows (when `on_failure = "na"`). The returned matrix carries two
attributes: `"n_failures"` (integer count of `NA` rows) and `"fail_idx"`
(integer realisation indices that failed).

## Usage

``` r
pesto_evaluate(model, theta, ...)
```

## Arguments

- model:

  A `pesto_forward_model` (or, for the multi-fidelity method, a
  `pesto_multifidelity_model`).

- theta:

  Numeric matrix, `nreal x npar`. Column names, when present, are
  checked against `model@param_names`.

- ...:

  Method-specific arguments. The multi-fidelity method accepts `level`
  (integer fidelity level to evaluate).

## Value

An `nreal x nobs` numeric matrix with attributes `"n_failures"` and
`"fail_idx"`.

## See also

[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md),
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md).

## Examples

``` r
fm <- pesto_forward_model(fn = function(theta) theta[, 1, drop = FALSE],
                          n_obs = 1L)
pesto_evaluate(fm, matrix(c(1, 2, 3), ncol = 1L))
#>      [,1]
#> [1,]    1
#> [2,]    2
#> [3,]    3
#> attr(,"n_failures")
#> [1] 0
#> attr(,"fail_idx")
#> integer(0)
```
