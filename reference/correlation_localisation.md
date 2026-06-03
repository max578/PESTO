# Correlation-Based Automatic Localisation Taper

Builds an `npar x nobs` localisation taper directly from the ensemble,
with no parameter / observation coordinates required. This is the
iterative-ensemble-smoother-native localisation of Luo & Bhakta (2020):
spurious sample correlations between a parameter and an observation (an
artefact of finite ensemble size) are damped, while genuine correlations
that stand above an estimated noise floor are retained.

## Usage

``` r
correlation_localisation(
  par_diff,
  obs_diff,
  threshold = -1,
  taper = "hard",
  n_shuffle = 1L,
  quantile = 0.95
)
```

## Arguments

- par_diff:

  Matrix (npar x nreal). Parameter anomalies.

- obs_diff:

  Matrix (nobs x nreal). Observation anomalies.

- threshold:

  Numeric. Noise floor on \\\|\rho\|\\. Negative (default `-1`) triggers
  automatic estimation by permutation.

- taper:

  Character. `"hard"` (indicator) or `"soft"` (linear ramp above the
  floor).

- n_shuffle:

  Integer. Number of permutation replicates for the automatic floor
  (default 1; each replicate yields `npar * nobs` spurious samples).
  Ignored when `threshold >= 0`.

- quantile:

  Numeric in (0, 1). Quantile of the spurious-correlation distribution
  used as the floor (default 0.95). Ignored when `threshold >= 0`.

## Value

A list with `rho` (the npar x nobs taper), `threshold` (the floor used),
`n_active` (count of entries with non-zero weight), and `frac_active`
(that count over `npar * nobs`).

## Details

The sample correlation \\\rho\_{ij}\\ between parameter-anomaly row
\\i\\ and observation-anomaly row \\j\\ is compared against a noise
floor \\\theta\\. When `threshold < 0` the floor is estimated by
destroying the parameter-observation link — the realisation order of the
observation anomalies is randomly permuted, independently per replicate,
and the floor is taken as a high quantile (default 0.95) of the
resulting spurious \\\|\rho\|\\ values. The permutation uses R's RNG, so
the estimate is reproducible under
[`set.seed()`](https://rdrr.io/r/base/Random.html).

Two tapers are offered. `"hard"` keeps correlations above the floor
unchanged (weight 1) and zeroes the rest. `"soft"` applies a smooth,
monotone ramp \\w\_{ij} = \mathrm{clip}((\|\rho\_{ij}\| - \theta) / (1 -
\theta), 0, 1)\\, which downweights near-floor correlations rather than
thresholding them sharply.

## References

Luo, X. & Bhakta, T. (2020). Automatic and adaptive localization for
ensemble-based history matching. *Journal of Petroleum Science and
Engineering*, 184, 106559.

## Examples

``` r
set.seed(1L)
npar <- 8L; nobs <- 5L; nreal <- 40L
pd <- matrix(rnorm(npar * nreal), npar, nreal)
# Make parameter 1 genuinely correlated with observation 1.
od <- matrix(rnorm(nobs * nreal), nobs, nreal)
od[1L, ] <- od[1L, ] + 2 * pd[1L, ]
loc <- correlation_localisation(pd, od)
loc$threshold
#> [1] 0.2559031
loc$rho[1L, 1L]
#> [1] 1
```
