# Ensemble Solution Kernel (MDA / Evensen form)

Implements the Multiple Data Assimilation (MDA) update from Evensen
(2018), Section 14.3.2. Uses low-rank representation of the error
covariance.

## Usage

``` r
ensemble_solution_mda(
  par_diff,
  obs_diff,
  obs_resid,
  obs_err,
  cur_lam = 1,
  eigthresh = 1e-06
)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies.

- obs_diff:

  Matrix (nobs x nreal). Observation anomalies.

- obs_resid:

  Matrix (nobs x nreal). Observation residuals.

- obs_err:

  Matrix (nobs x nreal). Observation error realisations.

- cur_lam:

  Numeric. Inflation factor.

- eigthresh:

  Numeric. Eigenvalue truncation threshold.

## Value

Matrix (nreal x npar). Parameter upgrade vectors.

## Examples

``` r
set.seed(1L)
npar  <- 4L
nreal <- 20L
nobs  <- 30L
par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
obs_err   <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
upgrade <- ensemble_solution_mda(
  par_diff  = par_diff,
  obs_diff  = obs_diff,
  obs_resid = obs_resid,
  obs_err   = obs_err,
  cur_lam   = 1.0
)
dim(upgrade)
#> [1] 20  4
```
