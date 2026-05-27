# Write a PEST Control File (.pst)

Writes a `pesto_pst` object to a PEST-format control file.

## Usage

``` r
write_pst(pst, file)
```

## Arguments

- pst:

  A `pesto_pst` object (as returned by
  [`read_pst()`](https://max578.github.io/PESTO/reference/read_pst.md)).

- file:

  Character. Output file path.

## Value

Invisible `NULL`. File is written as a side effect.

## See also

[`read_pst()`](https://max578.github.io/PESTO/reference/read_pst.md)

## Examples

``` r
pars <- data.table::data.table(
  parnme = c("k1", "k2"), partrans = "log", parchglim = "factor",
  parval1 = c(1.0, 0.5), parlbnd = c(0.01, 0.001),
  parubnd = c(100, 50), pargp = "hydraulic"
)
obs <- data.table::data.table(
  obsnme = c("h1", "h2"), obsval = c(1.0, 2.0),
  weight = c(1.0, 1.0), obgnme = "head"
)
pst <- create_pest_scenario(pars, obs, model_command = "echo run")
tf <- tempfile(fileext = ".pst")
on.exit(unlink(tf), add = TRUE)
write_pst(pst, tf)
file.exists(tf)
#> [1] TRUE
```
