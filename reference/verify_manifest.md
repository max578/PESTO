# Verify the integrity of a manifest

Recomputes the SHA-256 hash over
`(params, outputs, weights, obs_target, seed)` and compares against the
stored `data_hash`. Use this after
[`read_manifest()`](https://max578.github.io/PESTO/reference/read_manifest.md)
to confirm the data files have not been tampered with or silently
re-saved.

## Usage

``` r
verify_manifest(manifest, ...)
```

## Arguments

- manifest:

  A `pesto_ensemble_manifest`.

- ...:

  Reserved.

## Value

A list with `ok` (`TRUE`, `FALSE`, or `NA` – see Details), `stored` (the
manifest's recorded hash), `recomputed` (the hash computed from current
data), and `message` (`NULL` for verifiable formats, otherwise an
explanation).

## Details

When the manifest's `format` slot is `"csv"`, the on-disk data has been
round-tripped through
[`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) and IEEE
754 doubles have been truncated at the formatter's precision (~15-17
digits). That precision loss is enough to flip the SHA-256 hash, so
`verify_manifest()` returns `ok = NA` with an explanatory `message`
field rather than reporting a spurious `FALSE`.

## Examples

``` r
m <- as_manifest(
  pesto_ies_callback(
    function(t) t %*% t(matrix(1, 2L, 2L)),
    matrix(rnorm(20L), 10L, 2L, dimnames = list(NULL, c("a", "b"))),
    c(o1 = 0.1, o2 = -0.2), obs_sd = 0.1, noptmax = 2L, verbose = FALSE
  ),
  seed = 1L
)
verify_manifest(m)$ok
#> [1] TRUE
```
