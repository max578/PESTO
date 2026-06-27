# Is a PEST++ Executable Available?

A non-erroring capability probe for a PEST++ family executable. It is
the documented way for examples, vignettes, and conditional tests to
skip gracefully when no binary is installed: every PESTO algorithm runs
natively in R, so the external binaries are needed only for the optional
cross-checking and `.pst`-file paths.

## Usage

``` r
pestpp_available(which = "pestpp-ies")
```

## Arguments

- which:

  Character scalar naming the executable to probe. Defaults to
  `"pestpp-ies"`; any PEST++ tool name is accepted, e.g. `"pestpp-glm"`,
  `"pestpp-swp"`, `"pestpp-sen"`.

## Value

A length-one logical: `TRUE` if the named executable is resolvable,
`FALSE` otherwise. Never errors.

## Details

The probe looks first for a copy bundled with PESTO (`inst/bin`), then
for the executable on the system `PATH`. It mirrors the resolution used
by
[`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md)
and friends but never throws.

## See also

[`pesto_ies()`](https://max578.github.io/PESTO/reference/pesto_ies.md),
[`pesto_glm()`](https://max578.github.io/PESTO/reference/pesto_glm.md)

## Examples

``` r
# FALSE on a machine without PEST++ installed -- and that is fine:
pestpp_available()
#> [1] FALSE
pestpp_available("pestpp-glm")
#> [1] FALSE
```
