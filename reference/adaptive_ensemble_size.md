# Adaptive Ensemble Sizing

Dynamically determines the optimal ensemble size based on convergence
diagnostics and information-theoretic criteria.

## Usage

``` r
adaptive_ensemble_size(
  phi_values,
  current_size,
  min_size = 20L,
  max_size = 500L,
  cv_target = 0.3
)
```

## Arguments

- phi_values:

  Numeric vector. Current phi values per realisation.

- current_size:

  Integer. Current ensemble size.

- min_size:

  Integer. Minimum ensemble size (default 20).

- max_size:

  Integer. Maximum ensemble size (default 500).

- cv_target:

  Numeric. Target coefficient of variation for phi (default 0.3).

## Value

A list with recommended_size, reasoning, and diagnostics.

## Details

Uses the effective sample size (ESS) and coefficient of variation of phi
to determine whether the ensemble is too large (wasting compute) or too
small (poor UQ coverage).

## Examples

``` r
set.seed(1L)
phi_values <- rnorm(50L, mean = 100, sd = 20)^2
res <- adaptive_ensemble_size(
  phi_values   = phi_values,
  current_size = 50L
)
res$recommended_size
#> [1] 65
res$cv_phi
#> [1] 0.3029322
res$ess_ratio
#> [1] 0.02
```
