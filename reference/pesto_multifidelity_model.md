# Multi-Fidelity Forward Model (S7 class)

Bundles an ordered stack of
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
levels – cheapest (`level = 0`) to most expensive (`level = n - 1`) –
with their relative per-evaluation costs. This is the first-class form
of the bridge's `(cheap, expensive)` fidelity vector: the IES driver
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
selects a level per iteration via `fidelity_schedule`, and
[`mf_control_variate()`](https://max578.github.io/PESTO/reference/mf_control_variate.md)
debiases a cheap level against a sparse expensive sample for surrogate
cascades.

## Usage

``` r
pesto_multifidelity_model(levels, costs = NULL, label = NA_character_)
```

## Arguments

- levels:

  List of
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  objects (or bare functions, which are coerced), ordered cheapest
  first.

- costs:

  Numeric vector of relative per-evaluation costs, one per level,
  ascending by convention. Defaults to `seq_along(levels)`. Carried for
  cost-aware allocation; not yet used to schedule automatically (that is
  the documented extension point).

- label:

  Character. Optional human label.

## Value

A `pesto_multifidelity_model` S7 object.

## Details

Each element of `levels` may be a bare `function(theta) -> obs` or a
fully-specified
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md);
bare functions are coerced via
[`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md).
Levels must be ordered by ascending fidelity (cheapest first) – the
convention the `level` index and the `costs` vector both follow.

## See also

[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md),
[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md),
[`mf_control_variate()`](https://max578.github.io/PESTO/reference/mf_control_variate.md),
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).

## Examples

``` r
cheap     <- function(theta) theta %*% c(1, 1)        # fast, biased
expensive <- function(theta) theta %*% c(1, 1) + 0.5  # slow, truth
mf <- pesto_multifidelity_model(
  levels = list(
    pesto_forward_model(fn = cheap,     n_obs = 1L, fidelity = 0L),
    pesto_forward_model(fn = expensive, n_obs = 1L, fidelity = 1L)
  ),
  costs = c(1, 25)
)
theta <- matrix(c(1, 0, 0, 1), nrow = 2L, byrow = TRUE)
pesto_evaluate(mf, theta, level = 0L)  # cheap
#>      [,1]
#> [1,]    1
#> [2,]    1
#> attr(,"n_failures")
#> [1] 0
#> attr(,"fail_idx")
#> integer(0)
pesto_evaluate(mf, theta, level = 1L)  # expensive
#>      [,1]
#> [1,]  1.5
#> [2,]  1.5
#> attr(,"n_failures")
#> [1] 0
#> attr(,"fail_idx")
#> integer(0)
```
