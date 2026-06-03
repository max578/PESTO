# Forward-Model Contract (S7 class)

A `pesto_forward_model` wraps a user callable of signature
`function(theta) -> obs` – where `theta` is an `nreal x npar` numeric
matrix and `obs` is an `nreal x nobs` numeric matrix – in a typed object
that owns the evaluation contract: output dimensionality, expected
parameter names, the failure policy, the concurrency strategy, and a
fidelity tag.

## Usage

``` r
pesto_forward_model(
  fn,
  n_obs = NA_integer_,
  param_names = character(0),
  on_failure = "na",
  parallel = "serial",
  n_cores = NA_integer_,
  map_fn = NULL,
  max_fail_frac = 1,
  fidelity = 0L,
  label = NA_character_
)
```

## Arguments

- fn:

  Function. The forward model, signature `function(theta) -> obs`.

- n_obs:

  Integer or `NA`. Known observation dimensionality. If `NA` (default)
  it is inferred from the first successful evaluation.

- param_names:

  Character. Expected parameter column names. Empty (default) disables
  the column check.

- on_failure:

  Character. `"na"` (default) records failed realisations as `NA` rows
  and proceeds; `"stop"` aborts on any failure.

- parallel:

  Character. `"serial"` (default) or `"multicore"`.

- n_cores:

  Integer or `NA`. Worker count for `"multicore"`. `NA` (default)
  resolves to `parallel::detectCores() - 1L` at evaluation time.

- map_fn:

  Function or `NULL`. Optional `lapply`-shaped override
  `function(X, FUN, ...)`; when supplied it drives per-realisation
  dispatch regardless of `parallel`.

- max_fail_frac:

  Numeric in `[0, 1]`. Abort if the fraction of failed realisations in
  any single evaluation exceeds this. Default `1` (never abort on
  fraction; `on_failure` still governs the zero-success case).

- fidelity:

  Integer. Fidelity level tag (default `0L`).

- label:

  Character. Optional human label carried into diagnostics.

## Value

A `pesto_forward_model` S7 object.

## Details

This is the single contract both PESTO adapter modes honour. The native
callback driver
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
accepts either a bare function (auto-wrapped via
[`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md)
with that driver's failure policy) or a `pesto_forward_model` built
here; the apsimx adapter
[`apsim_callback()`](https://max578.github.io/PESTO/reference/apsim_callback.md)
can emit one directly. Evaluation is via
[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md),
which guarantees the returned shape and accounts for failed realisations
as `NA` rows.

## Concurrency

With `parallel = "serial"` (the default) and no `map_fn`, evaluation
attempts a single bulk call `fn(theta)` and falls back to a serial
per-realisation loop only if the bulk call errors – preserving the fast
path for vectorised forward models. With `parallel = "multicore"`
realisations are dispatched per row through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(fork-based; silently serial on Windows). A custom `map_fn` (an
`lapply`-shaped `function(X, FUN, ...)`) overrides both and lets callers
plug in `future.apply::future_lapply`, `mirai`, or a cluster backend.
For reproducible parallel runs set `RNGkind("L'Ecuyer-CMRG")` and
[`set.seed()`](https://rdrr.io/r/base/Random.html) before evaluating;
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) then
draws independent streams per realisation.

## Fidelity

`fidelity` is an integer level tag (`0L` = base / cheapest by
convention; higher = more expensive / higher resolution). A
single-fidelity model leaves it at `0L`. The tag is what
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
uses to order levels and what the manifest emitter threads into
provenance.

## See also

[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md)
to evaluate one;
[`as_forward_model()`](https://max578.github.io/PESTO/reference/as_forward_model.md)
to coerce a bare function;
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md)
to compose several across fidelity levels;
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
for the IES driver that consumes it.

## Examples

``` r
# A vectorised linear forward model wrapped as a contract object.
G  <- matrix(c(1, 0, 0, 1, 1, -1), nrow = 3L, byrow = TRUE)
fm <- pesto_forward_model(
  fn          = function(theta) theta %*% t(G),
  n_obs       = 3L,
  param_names = c("a", "b")
)
theta <- matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE,
                dimnames = list(NULL, c("a", "b")))
pesto_evaluate(fm, theta)
#>      [,1] [,2] [,3]
#> [1,]    1    2   -1
#> [2,]    3    4   -1
#> attr(,"n_failures")
#> [1] 0
#> attr(,"fail_idx")
#> integer(0)
```
