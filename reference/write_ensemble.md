# Write an Ensemble File

Writes an ensemble data.table to CSV format compatible with PEST++.

## Usage

``` r
write_ensemble(ensemble, file, format = "csv")
```

## Arguments

- ensemble:

  A data.table with realisation data.

- file:

  Character. Output file path.

- format:

  Character. Currently only "csv" is supported.

## Value

Invisible `NULL`.

## See also

[`read_ensemble()`](https://max578.github.io/PESTO/reference/read_ensemble.md)

## Examples

``` r
ens <- data.table::data.table(
  real_name = sprintf("real_%02d", 1:5),
  k1 = runif(5, 0.1, 10),
  k2 = runif(5, 0.01, 1)
)
tf <- tempfile(fileext = ".csv")
on.exit(unlink(tf), add = TRUE)
write_ensemble(ens, tf)
file.exists(tf)
#> [1] TRUE
```
