# Covariance Inflation Specification for IES

Builds a control object that tells
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
and
[`pesto_ies_filter()`](https://max578.github.io/PESTO/reference/pesto_ies_filter.md)
how to counteract ensemble under-dispersion — the progressive collapse
of posterior spread that an iterative ensemble smoother suffers with a
finite ensemble. Inflation re-expands the post-update parameter spread
each iteration; the default `method = "none"` leaves the update
byte-identical to the un-inflated smoother.

## Usage

``` r
pesto_inflation(
  method = c("none", "rtps", "adaptive", "multiplicative"),
  alpha = 0.5,
  factor = 1,
  retention_floor = 0.5,
  max_factor = 5
)
```

## Arguments

- method:

  Character. One of `"none"` (default), `"rtps"`, `"adaptive"`,
  `"multiplicative"`.

- alpha:

  Numeric in \[0, 1\]. RTPS relaxation coefficient (default 0.5); used
  only when `method = "rtps"`.

- factor:

  Numeric \\\ge\\ 1. Fixed inflation factor for
  `method = "multiplicative"` (default 1, i.e. no inflation).

- retention_floor:

  Numeric in (0, 1\]. Target floor on the mean spread-retention ratio
  for `method = "adaptive"` (default 0.5).

- max_factor:

  Numeric \\\ge\\ 1. Upper bound on any single-iteration inflation
  factor for `"rtps"` and `"adaptive"` (default 5).

## Value

An object of class `"pesto_inflation"`.

## Details

Four methods are offered. `"rtps"` is relaxation-to-prior-spread
(Whitaker & Hamill 2012): each parameter's posterior anomalies are
rescaled by \\\alpha(\sigma^{b} - \sigma^{a})/\sigma^{a} + 1\\, where
\\\sigma^{b}\\ and \\\sigma^{a}\\ are the background (pre-update) and
analysis (post-update) standard deviations. Being per-parameter, it
re-inflates the directions that collapsed hardest, so it is the
spectrally-aware workhorse. `"adaptive"` is a global,
magnitude-targeting scheme: it measures the mean spread-retention ratio
\\q = \mathrm{mean}\_j(\sigma^{a}\_j/\sigma^{b}\_j)\\ and, when `q`
falls below `retention_floor`, applies a single multiplicative factor
\\\min(\texttt{max\\factor}, \texttt{retention\\floor}/q)\\ to restore
the lost variance magnitude. `"multiplicative"` applies a fixed `factor`
every iteration. `"none"` disables inflation.

The companion *diagnostic* is the spectral spread-ESS
([`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md)),
recorded each iteration regardless of method: it reports the effective
number of variance-carrying directions and is what detects directional
collapse. Because that participation ratio is invariant to a global
rescaling, a global (`"multiplicative"` / `"adaptive"`) inflation
restores variance *magnitude* but not the spectral *shape*; `"rtps"` is
the method that reshapes the spectrum. The two compose well.

## References

Whitaker, J.S. & Hamill, T.M. (2012). Evaluating methods to account for
system errors in ensemble data assimilation. *Monthly Weather Review*,
140(9), 3078–3089.

## See also

[`pesto_localisation()`](https://max578.github.io/PESTO/reference/pesto_localisation.md),
[`ensemble_spread_ess()`](https://max578.github.io/PESTO/reference/ensemble_spread_ess.md),
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).

## Examples

``` r
inf <- pesto_inflation("rtps", alpha = 0.5)
inf
#> $method
#> [1] "rtps"
#> 
#> $alpha
#> [1] 0.5
#> 
#> $factor
#> [1] 1
#> 
#> $retention_floor
#> [1] 0.5
#> 
#> $max_factor
#> [1] 5
#> 
#> attr(,"class")
#> [1] "pesto_inflation"
```
