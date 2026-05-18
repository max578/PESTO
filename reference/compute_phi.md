# Compute Phi (Objective Function) for Ensemble

Calculates the weighted sum of squared residuals for each realisation in
the ensemble.

## Usage

``` r
compute_phi(residuals, weights)
```

## Arguments

- residuals:

  Matrix (nobs x nreal). Observation residuals.

- weights:

  Numeric vector (nobs). Observation weights.

## Value

Numeric vector (nreal). Phi value per realisation.

## Examples

``` r
set.seed(1L)
residuals <- matrix(rnorm(5 * 4), 5, 4)
weights   <- rep(1, 5)
phi <- compute_phi(residuals, weights)
length(phi)
#> [1] 4
phi
#> [1] 3.777941 1.880665 8.993765 1.920231
```
