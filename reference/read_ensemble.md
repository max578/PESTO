# Read an Ensemble File

Reads PEST++ ensemble files in CSV or binary (.jcb/.jco) format.

## Usage

``` r
read_ensemble(file, format = c("csv", "binary"))
```

## Arguments

- file:

  Character. Path to ensemble file.

- format:

  Character. One of "csv" (default) or "binary".

## Value

A data.table with realisations as rows and parameters/observations as
columns. The first column `real_name` contains realisation names.

## See also

[`write_ensemble()`](https://AAGI-AUS.github.io/PESTO/reference/write_ensemble.md)

## Examples

``` r
ens <- data.table::data.table(
  real_name = sprintf("real_%02d", 1:10),
  k1 = rnorm(10, mean = 1.0, sd = 0.2),
  k2 = rnorm(10, mean = 0.5, sd = 0.1),
  k3 = rnorm(10, mean = 2.0, sd = 0.3)
)
tf <- tempfile(fileext = ".csv")
on.exit(unlink(tf), add = TRUE)
write_ensemble(ens, tf)
ens_back <- read_ensemble(tf, format = "csv")
identical(names(ens_back), names(ens))
#> [1] TRUE
nrow(ens_back) == nrow(ens)
#> [1] TRUE
```
