# Covariance Localisation Specification for IES

Builds a control object that tells
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
and
[`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
how to taper the ensemble Kalman gain, suppressing the spurious
long-range parameter-observation correlations that a finite ensemble
manufactures. Localisation is applied as a Schur (elementwise) product
on the explicit gain inside
[`ensemble_solution_localised()`](https://max578.github.io/PESTO/reference/ensemble_solution_localised.md);
the default `method = "none"` leaves the standard SVD update untouched.

## Usage

``` r
pesto_localisation(
  method = c("none", "correlation", "distance"),
  taper = c("hard", "soft"),
  threshold = -1,
  n_shuffle = 1L,
  quantile = 0.95,
  distances = NULL,
  par_coords = NULL,
  obs_coords = NULL,
  radius = NULL
)
```

## Arguments

- method:

  Character. One of `"none"` (default), `"correlation"`, `"distance"`.

- taper:

  Character. `"hard"` (default) or `"soft"`; passed to
  [`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md)
  for `method = "correlation"`.

- threshold:

  Numeric. Correlation noise floor; negative (default -1) triggers
  automatic per-iteration estimation. `method = "correlation"`.

- n_shuffle:

  Integer \\\ge\\ 1. Permutation replicates for the automatic floor
  (default 1). `method = "correlation"`.

- quantile:

  Numeric in (0, 1). Quantile of the spurious-correlation distribution
  used as the floor (default 0.95). `method = "correlation"`.

- distances:

  Matrix (npar x nobs) or `NULL`. Precomputed parameter-to-observation
  distances for `method = "distance"`.

- par_coords, obs_coords:

  Matrices (npar x d, nobs x d) or `NULL`. Parameter / observation
  coordinates; Euclidean distances are derived when `distances` is
  `NULL`. `method = "distance"`.

- radius:

  Numeric (\> 0) or `NULL`. Gaspari-Cohn localisation radius; required
  for `method = "distance"`.

## Value

An object of class `"pesto_localisation"`.

## Details

`"correlation"` is the iterative-ensemble-smoother-native automatic
localisation of Luo & Bhakta (2020): it needs no parameter or
observation coordinates, estimating a noise floor from the ensemble
itself and damping sample correlations that fall below it (see
[`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md)).
This is the recommended default for parameter-estimation problems whose
parameters carry no spatial metric. `"distance"` is classical
distance-based localisation: a Gaspari-Cohn taper
([`gaspari_cohn()`](https://max578.github.io/PESTO/reference/gaspari_cohn.md))
of a parameter-to-observation distance matrix, for problems where such a
metric exists — supply either `distances` directly or `par_coords` +
`obs_coords` (Euclidean distances are then computed), together with
`radius`.

## References

Luo, X. & Bhakta, T. (2020). Automatic and adaptive localization for
ensemble-based history matching. *Journal of Petroleum Science and
Engineering*, 184, 106559.

## See also

[`pesto_inflation()`](https://max578.github.io/PESTO/reference/pesto_inflation.md),
[`correlation_localisation()`](https://max578.github.io/PESTO/reference/correlation_localisation.md),
[`gaspari_cohn()`](https://max578.github.io/PESTO/reference/gaspari_cohn.md).

## Examples

``` r
loc <- pesto_localisation("correlation", taper = "soft")
loc
#> $method
#> [1] "correlation"
#> 
#> $taper
#> [1] "soft"
#> 
#> $threshold
#> [1] -1
#> 
#> $n_shuffle
#> [1] 1
#> 
#> $quantile
#> [1] 0.95
#> 
#> $distances
#> NULL
#> 
#> $par_coords
#> NULL
#> 
#> $obs_coords
#> NULL
#> 
#> $radius
#> NULL
#> 
#> attr(,"class")
#> [1] "pesto_localisation"
```
