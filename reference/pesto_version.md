# Get PESTO package version information

Returns the PESTO package version, plus the version of the PEST++
install PESTO resolves (see
[`pestpp_available()`](https://max578.github.io/PESTO/reference/pestpp_available.md)
for how it is found). PESTO bundles no binaries, so `pestpp_version`
reports whatever install is configured, or `"not found"` when there is
none.

## Usage

``` r
pesto_version()
```

## Value

A list with version strings: `pesto_version`, `pestpp_version`,
`platform`, and `r_version`.

## Examples

``` r
v <- pesto_version()
v$pesto_version
#> [1] "0.10.0"
v$platform
#> [1] "x86_64-pc-linux-gnu"
```
