# Get PESTO package version information

Returns version info for both the PESTO R package and the bundled PEST++
binaries.

## Usage

``` r
pesto_version()
```

## Value

A list with version strings.

## Examples

``` r
v <- pesto_version()
v$pesto_version
#> [1] "0.9.0.9000"
v$platform
#> [1] "x86_64-pc-linux-gnu"
```
