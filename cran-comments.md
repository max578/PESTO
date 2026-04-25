# cran-comments.md

## Test environments

- Local: macOS 26.3.1 (Apple Silicon, arm64), R 4.5.2 (2025-10-31).
- (Pre-submission) Will be expanded to `devtools::check_win_devel()` and
  `rhub::rhub_check()` matrix prior to first CRAN release.

## R CMD check results

`R CMD check --as-cran` on the source tarball: **0 errors, 1 warning, 2 notes**.

### WARNING (upstream, not from PESTO source)

```
* checking whether package ‘PESTO’ can be installed ... WARNING
Found the following significant warnings:
  warning: unknown warning group '-Wfixed-enum-extension', ignored
  [-Wunknown-warning-option]
```

These warnings originate in upstream headers, not PESTO source code:

- `-Wfixed-enum-extension` is emitted by `R_ext/Boolean.h:62` (R itself)
  under macOS 26 / Apple Clang 21 — a known transitional incompatibility
  between R 4.5 and the Tahoe SDK. The warning will disappear once R
  updates the header guard.
- The remaining unused-but-set-variable warnings come from RcppEigen's
  bundled Eigen headers (`SparseLU`, `IterativeSolvers`,
  `SparseExtra/MarketIO`).

PESTO has no path to fix either source. Adding `-Wno-*` flags to
`PKG_CXXFLAGS` resolves the WARNING but introduces a non-portable-flags
WARNING (CRAN policy), so we leave them. The CRAN reference machine and
Windows builds do not surface this warning.

### NOTE 1 — new submission and development version

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Max Moldovan <max.moldovan@adelaide.edu.au>'
New submission
Version contains large components (0.1.0.9000)
```

This is a development pre-release tarball. The version will be reset to
`0.1.0` (semver-compliant) at CRAN submission and the "New submission"
note is by design for a first submission.

### NOTE 2 — local HTML Tidy version

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent
enough HTML Tidy.
```

Local toolchain note only — CRAN reference machines bundle a current
HTML Tidy, so this note will not appear there.

## Downstream dependencies

None — first submission, no reverse dependencies.

## Numerical-validation regression test

`tests/testthat/test-ensemble-solution-sign.R` (added 2026-04-25 under
investigation I1) asserts strict monotone phi descent across at least
three IES iterations on a well-conditioned linear inverse problem with
the documented `obs_resid = sim - obs` convention, and geometric
divergence under the inverted convention. This regression-tests the
IES kernel sign convention against future docstring or refactor drift.

The PESTO IES kernel is independently verified against Chen & Oliver
(2013) eq. 12 at machine precision (max element-wise iteration-1 delta
= 5.8 x 10^-15 against a pure-R textbook reference; see
`inst/scripts/i2_paired_seed_check.R` and the I2 information note).
