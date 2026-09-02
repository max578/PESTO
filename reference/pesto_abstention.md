# Construct a typed PESTO abstention

Returns a small classed object signalling that a PESTO reliability gate
declined to produce a result, with a machine-readable `reason` a caller
(or the orchestra's `is_orchestra_decline()`) can branch on, instead of
an uncaught error or a half-populated result.

## Usage

``` r
pesto_abstention(
  reason,
  detail = NA_character_,
  diagnostics = list(),
  scope = "call"
)
```

## Arguments

- reason:

  Character scalar reason code. One of `"degenerate_ensemble"` (fewer
  than two realisations survived a step and the update cannot be
  formed), `"over_determined"` (the ensemble is too small relative to
  the number of parameters to support the update – `nreal <= npar`), or
  `"surrogate_off_design"` (the GP surrogate's
  training-point-to-parameter ratio falls below
  [`check_surrogate_regime()`](https://max578.github.io/PESTO/reference/check_surrogate_regime.md)'s
  favourable floor).

- detail:

  Character scalar human-readable detail, or `NA`.

- diagnostics:

  Named list of scalar diagnostics supporting the decline (e.g.
  `list(iteration = 3L, n_ok = 1L, nreal = 40L)`). Default an empty
  list.

- scope:

  Character scalar naming what was refused (default `"call"`).

## Value

An object of class `"pesto_abstention"`: a list with `reason`, `detail`,
`diagnostics`, `scope` and `abstained = TRUE`.

## Examples

``` r
a <- pesto_abstention("degenerate_ensemble",
                      detail = "1 of 40 realisations succeeded",
                      diagnostics = list(iteration = 2L, n_ok = 1L))
is_pesto_abstention(a)
#> [1] TRUE
```
