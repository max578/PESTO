# Convert a PESTO ensemble result into a `pesto_ensemble_manifest`

Wraps the plain-list result returned by
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
(and, eventually,
[`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md))
in the S7 manifest contract object so downstream packages can consume it
via S7 dispatch without reaching into PESTO-specific list internals.

## Usage

``` r
as_manifest(x, ...)
```

## Arguments

- x:

  A `pesto_ies_callback_result`, a `pesto_ies_filter_result`, a
  [`pesto_abstention()`](https://max578.github.io/PESTO/reference/pesto_abstention.md)
  (any object with a method registered against this generic).

- ...:

  Method-specific arguments. For `pesto_ies_callback_result`: `run_id`
  (character, defaults to a timestamp+hash slug), `seed` (integer,
  defaults to `NA_integer_`), `fidelity` (structured provenance list or
  `NULL`; defaults to the record captured by
  [`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
  for a multi-fidelity run), `apsim_version` (character, defaults to
  `NA_character_`; pass `attr(fm, "apsim_version")` from an
  [`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
  forward model to ground the run to its simulator), and `obs_schema` (a
  grounded semantic descriptor from
  [`pesto_obs_schema()`](https://max578.github.io/PESTO/reference/pesto_obs_schema.md)
  or `NULL`). For a `pesto_abstention`: the same `run_id` / `seed` /
  `apsim_version` / `obs_schema` arguments plus `method` (character, one
  of `"ies_callback"`, `"ies_filter"`, `"ies_pst"`, `"mda"`,
  `"surrogate"`; default `"ies_callback"` – there is no ensemble to
  infer it from).

## Value

A `pesto_ensemble_manifest` S7 object. When `x` is a
[`pesto_abstention()`](https://max578.github.io/PESTO/reference/pesto_abstention.md),
`params` / `outputs` are empty (zero-row) data frames, `weights` /
`obs_target` are empty numeric vectors, `failure_rate` is `1`, and
`summary` carries
`list(abstained = TRUE, reason = x$reason, detail = x$detail)` – the
abstention's reason travels into the manifest's typed verdict home so a
downstream consumer sees the decline without inspecting PESTO-internal
objects.

## Examples

``` r
fit <- pesto_ies_callback(
  function(t) t %*% t(matrix(1, 2L, 2L)),
  matrix(rnorm(20L), 10L, 2L, dimnames = list(NULL, c("a", "b"))),
  c(o1 = 0.1, o2 = -0.2), obs_sd = 0.1, noptmax = 2L, verbose = FALSE
)
m <- as_manifest(fit, seed = 1L)
class(m)
#> [1] "PESTO::pesto_ensemble_manifest" "S7_object"                     
```
