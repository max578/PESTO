# Contributing to PESTO

Thank you for your interest in contributing to PESTO. This document
covers the developer-side workflows that are not part of the CRAN
shipped package.

## Developer-only PEST++ benchmark

The CRAN tarball is intentionally self-contained: every example, test,
and vignette runs without the upstream `pestpp-ies` binary. The
`vignettes/pestpp-comparison-and-simulation.Rmd` vignette compares PESTO
native IES against the pure-R textbook reference
([`pesto_reference_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_reference_ies.md))
by default; this comparison is what every end-user sees.

If you have `pestpp-ies` installed and want the vignette to extend its
agreement plot to also include the upstream binary’s posterior, set the
`PESTO_PESTPP_BIN` environment variable (or have `pestpp-ies` on `PATH`)
and ensure a developer-side cache exists at
`tools/pestpp_benchmark/scenario_a_pestpp_ies.rds`.

To regenerate the cache:

``` bash
export PESTO_PESTPP_BIN=/path/to/pestpp-ies
cd /path/to/PESTO
Rscript tools/pestpp_benchmark/run_benchmark.R
```

The script regenerates the always-shipped reference cache at
`inst/extdata/pestpp_cache/scenario_a_reference.rds` deterministically
(SHA-256-pinned prior ensemble) and, when the binary is available, also
refreshes the developer-only `pestpp-ies` cache. Both files are keyed to
the same prior digest so the vignette refuses to compare mismatched
caches.

`tools/` is `.Rbuildignore`-d and never enters the CRAN tarball.

## Re-run policy

Re-run `tools/pestpp_benchmark/run_benchmark.R` after any change to:

1.  The Scenario A problem definition in
    `vignettes/pestpp-comparison-and-simulation.Rmd`.
2.  `R/pesto_reference_ies.R` or `src/ensemble_solve.cpp`.
3.  The pestpp-ies binary version (when refreshing the developer-side
    cache).

Commit the regenerated reference cache; commit the regenerated
developer-side `pestpp-ies` cache only if you also want to update the
documentation snapshot for other developers.

## R package conventions

- snake_case for functions, arguments, and exported objects.
- `data.table` is the default for tabular work; tidyverse permitted
  where it raises clarity (e.g. `dplyr` for grammar-shaped APIs).
- `cli` for user-facing messages with bare-symbol keys (`{.fn fn}`,
  `{.path /tmp}`, `{.val 0.05}`).
- Australian / British English (en-AU); ISO-8601 dates.
- Run `devtools::document()` before every commit that touches roxygen.
- Run
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  before every commit that touches `// [[Rcpp::export]]` blocks in
  `src/`.

## Quality gates before opening a PR

``` r

# House-style + spelling + URLs
lintr::lint_package()
devtools::spell_check()
urlchecker::url_check()

# Test surface
devtools::test()
covr::package_coverage()  # target >= 85%

# Full pre-CRAN pass
devtools::check(args = "--as-cran")
```

The Apple Silicon toolchain currently lacks `libgcov`, so
`covr::package_coverage()` is meaningful only on Linux CI.

## Issue and discussion

PEST++ upstream lives at <https://github.com/usgs/pestpp>; PESTO is not
a fork — it is a separate R-native re-engineering. Issues specific to
PESTO go to <https://github.com/AAGI-AUS/PESTO/issues>.
