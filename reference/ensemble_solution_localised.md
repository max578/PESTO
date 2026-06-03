# Localised Ensemble Solution Kernel (explicit-gain GLM form)

Computes the IES Gauss-Levenberg-Marquardt update with state-space
covariance localisation applied as a Schur (elementwise) product on the
explicit Kalman gain. The standard SVD kernel
[`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
works in the reduced observation-anomaly subspace and never forms the
`npar x nobs` gain, so it cannot host localisation; this kernel
reconstructs the gain \$\$K = \Delta\Theta\\ V\\ \mathrm{diag}(s)\\
\mathrm{diag}((s^2 + (\lambda+1))^{-1})\\ U^{T},\$\$ (with \\U s V^{T}\\
the thin SVD of the weight-scaled observation anomalies) tapers it as
\\K \circ \rho\\, and applies it to the weighted residuals.

## Usage

``` r
ensemble_solution_localised(
  par_diff,
  obs_diff,
  obs_resid,
  weights,
  rho,
  cur_lam,
  eigthresh = 1e-06
)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies.

- obs_diff:

  Matrix (nobs x nreal). Observation anomalies.

- obs_resid:

  Matrix (nobs x nreal). Observation residuals (sim - obs); see
  [`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
  for the sign rationale.

- weights:

  Numeric vector (nobs). Observation weights (1 / sqrt(variance)).

- rho:

  Matrix (npar x nobs). Localisation taper in \\\[0, 1\]\\, e.g. from
  [`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md)
  or
  [`gaspari_cohn()`](https://max578.github.io/PESTO/reference/gaspari_cohn.md).

- cur_lam:

  Numeric. Current Marquardt lambda.

- eigthresh:

  Numeric. Eigenvalue truncation threshold (0-1).

## Value

Matrix (nreal x npar). Negative-direction parameter upgrade, applied by
subtraction.

## Details

When \\\rho \equiv 1\\ the result is identical (to truncation tolerance)
to
[`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
with `use_approx = TRUE`; the prior-scaling null-space correction
(`upgrade_2`) is not part of the localised path. The returned matrix
follows the same sign convention as
[`ensemble_solution()`](https://max578.github.io/PESTO/reference/ensemble_solution.md)
— it is the negative-direction step, applied to the ensemble by
subtraction (`par_new = par_old - upgrade`).

## References

Chen, Y. & Oliver, D.S. (2013). Levenberg-Marquardt forms of the
iterative ensemble smoother for efficient history matching and
uncertainty quantification. *Computational Geosciences*, 17(4), 689–703.

## Examples

``` r
set.seed(1L)
npar <- 4L; nreal <- 20L; nobs <- 6L
par_diff  <- matrix(rnorm(npar * nreal), npar, nreal)
obs_diff  <- matrix(rnorm(nobs * nreal), nobs, nreal)
obs_resid <- matrix(rnorm(nobs * nreal, sd = 0.5), nobs, nreal)
weights   <- rep(1, nobs)
rho       <- matrix(1, npar, nobs)          # no localisation
upg <- ensemble_solution_localised(
  par_diff, obs_diff, obs_resid, weights, rho, cur_lam = 1.0
)
dim(upg)
#> [1] 20  4
```
