# Surrogate-Accelerated Ensemble Update

Performs an IES update step using a GP surrogate for cheap
pre-screening, with adaptive switching to the full model based on
prediction uncertainty.

## Usage

``` r
surrogate_ensemble_update(
  par_ensemble,
  obs_ensemble,
  obs_target,
  weights,
  parcov_inv,
  cur_lam = 1,
  uncertainty_threshold = 0.1,
  eigthresh = 1e-06
)
```

## Arguments

- par_ensemble:

  Matrix (nreal x npar). Current parameter ensemble.

- obs_ensemble:

  Matrix (nreal x nobs). Current observation ensemble (from model).

- obs_target:

  Numeric vector (nobs). Target observations.

- weights:

  Numeric vector (nobs). Observation weights.

- parcov_inv:

  Numeric vector (npar). Inverse parameter covariance diagonal.

- cur_lam:

  Numeric. Marquardt lambda.

- uncertainty_threshold:

  Numeric. Threshold for surrogate/model switching. Realisations with GP
  uncertainty above this are re-evaluated with full model. Default 0.1
  (10% of signal variance).

- eigthresh:

  Numeric. SVD eigenvalue threshold.

## Value

A list with:

- upgrade:

  Matrix. Parameter upgrades.

- n_model_runs:

  Integer. Number of full model evaluations needed.

- n_surrogate_runs:

  Integer. Number of surrogate evaluations.

- savings_pct:

  Numeric. Percentage of model runs saved.

- gp_diagnostics:

  List. GP training diagnostics.

## Details

**Algorithm:**

1.  Train GP surrogate from current ensemble evaluations

2.  Generate candidate upgrades using surrogate predictions

3.  Evaluate uncertainty of surrogate predictions

4.  Run full model only for realisations where uncertainty exceeds
    threshold

5.  Blend surrogate and model results using control-variate correction

This typically reduces full model evaluations by 50-90%.

## Examples

``` r
# \donttest{
set.seed(1L)
npar  <- 5L
nreal <- 15L
nobs  <- 8L
par_ensemble <- matrix(rnorm(nreal * npar), nreal, npar)
obs_ensemble <- matrix(rnorm(nreal * nobs), nreal, nobs)
obs_target   <- rnorm(nobs)
weights      <- rep(1, nobs)
parcov_inv   <- rep(1, npar)
res <- surrogate_ensemble_update(
  par_ensemble = par_ensemble,
  obs_ensemble = obs_ensemble,
  obs_target   = obs_target,
  weights      = weights,
  parcov_inv   = parcov_inv,
  cur_lam      = 1.0
)
dim(res$upgrade)
#> [1] 15  5
res$savings_pct
#> [1] 100
# }
```
