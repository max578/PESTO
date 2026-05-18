# Surrogate-Accelerated IES Iteration

Performs a single IES iteration using a Gaussian Process surrogate to
reduce the number of expensive full-model evaluations.

## Usage

``` r
pesto_surrogate_ies(
  par_ensemble,
  obs_ensemble,
  obs_target,
  weights,
  parcov_inv,
  lambda = 1,
  uncertainty_threshold = 0.1,
  eigthresh = 1e-06
)
```

## Arguments

- par_ensemble:

  data.table or matrix. Current parameter ensemble (rows = realisations,
  columns = parameters).

- obs_ensemble:

  data.table or matrix. Current observation ensemble from model
  evaluations.

- obs_target:

  Named numeric vector. Target observation values.

- weights:

  Numeric vector. Observation weights.

- parcov_inv:

  Numeric vector. Diagonal of inverse parameter covariance.

- lambda:

  Numeric. Marquardt lambda (default 1.0).

- uncertainty_threshold:

  Numeric. Threshold for surrogate/model switching. Fraction of signal
  variance (default 0.1 = 10%).

- eigthresh:

  Numeric. SVD eigenvalue threshold.

## Value

A list containing:

- upgrade:

  Matrix of parameter upgrades

- n_model_runs:

  Number of full model evaluations needed

- n_surrogate_runs:

  Number of surrogate-only evaluations

- savings_pct:

  Percentage of model runs saved

- gp_diagnostics:

  GP training diagnostics

## Details

**How it works:**

1.  A GP surrogate is trained on existing parameter-observation pairs

2.  The surrogate predicts model outputs for all ensemble members

3.  Only members with high prediction uncertainty trigger full model
    runs

4.  Control-variate bias correction blends surrogate and model results

5.  Standard IES update is computed on the blended ensemble

This typically saves 50-90% of model evaluations per iteration.

## References

Rasmussen, C.E. & Williams, C.K.I. (2006). Gaussian Processes for
Machine Learning. MIT Press.

Liu, F. & Guillas, S. (2017). Dimension reduction for Gaussian process
emulation. *Statistics and Computing*, 27(3), 785-802.

## Examples

``` r
# \donttest{
set.seed(7L)
n_real <- 15L; n_par <- 5L; n_obs <- 8L
par_ens <- matrix(rnorm(n_real * n_par), n_real, n_par,
                  dimnames = list(NULL, paste0("k", 1:n_par)))
obs_ens <- matrix(rnorm(n_real * n_obs), n_real, n_obs,
                  dimnames = list(NULL, paste0("h", 1:n_obs)))
obs_target <- rnorm(n_obs)
weights    <- rep(1.0, n_obs)
parcov_inv <- rep(1.0, n_par)
res <- pesto_surrogate_ies(
  par_ensemble = par_ens,
  obs_ensemble = obs_ens,
  obs_target   = obs_target,
  weights      = weights,
  parcov_inv   = parcov_inv,
  lambda       = 1.0
)
res$savings_pct
#> [1] 100
# }
```
