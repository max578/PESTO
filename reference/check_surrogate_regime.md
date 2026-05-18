# Check whether a surrogate-IES regime is favourable

Issues a warning when the ratio of training points to parameters falls
below an empirical threshold, where the Gaussian-process surrogate
inside
[`pesto_surrogate_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_surrogate_ies.md)
and
[`surrogate_ensemble_update()`](https://AAGI-AUS.github.io/PESTO/reference/surrogate_ensemble_update.md)
is unlikely to repay its training cost. The check is a soft guardrail —
it does not modify the run, only flags an unfavourable regime so the
caller can decide whether to fall back to pure IES.

## Usage

``` r
check_surrogate_regime(n_params, n_train, threshold = 5)
```

## Arguments

- n_params:

  Integer. Number of estimated parameters.

- n_train:

  Integer. Number of training samples available to the surrogate
  (typically the ensemble size).

- threshold:

  Numeric. Minimum acceptable `n_train / n_params` ratio. Default `5`,
  the empirical soft floor from the surrogate-IES vignette.

## Value

Invisibly returns `TRUE` when the regime is favourable
(`n_train >= threshold * n_params`) and `FALSE` otherwise. Called for
the warning side-effect.

## Details

The default threshold of `5` corresponds to the soft floor
`n_train >= 5 * n_params` documented in
[`vignette("surrogate-ies", package = "PESTO")`](https://AAGI-AUS.github.io/PESTO/articles/surrogate-ies.md).
Below that floor the GP posterior variance typically stays above the
uncertainty-driven switching threshold and surrogate savings collapse to
near zero.

This is exposed as a stand-alone helper so users can call it explicitly
before scheduling an expensive ensemble. It is **not** invoked
automatically by
[`pesto_surrogate_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_surrogate_ies.md)
in the current release; that wiring is tracked as a v0.2 enhancement
candidate.

## References

Rasmussen, C. E. & Williams, C. K. I. (2006). *Gaussian Processes for
Machine Learning*. MIT Press.

## See also

[`pesto_surrogate_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_surrogate_ies.md),
[`surrogate_ensemble_update()`](https://AAGI-AUS.github.io/PESTO/reference/surrogate_ensemble_update.md),
[`vignette("surrogate-ies", package = "PESTO")`](https://AAGI-AUS.github.io/PESTO/articles/surrogate-ies.md)

## Examples

``` r
# Favourable regime: 100 training points for 10 parameters.
check_surrogate_regime(n_params = 10L, n_train = 100L)

# Unfavourable regime: 30 training points for 30 parameters
# (the curse-of-dimensionality case from Scenario C of the
# comparison-and-simulation vignette). Emits a warning.
suppressWarnings(
  check_surrogate_regime(n_params = 30L, n_train = 30L)
)

# Custom threshold for users with a smoother forward model.
check_surrogate_regime(n_params = 20L, n_train = 60L, threshold = 3)
```
