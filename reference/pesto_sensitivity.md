# Run PEST++ SEN (Global Sensitivity Analysis)

Executes pestpp-sen for Morris or Sobol sensitivity analysis.

## Usage

``` r
pesto_sensitivity(
  pst_file,
  method = c("morris", "sobol"),
  exe = NULL,
  extra_args = list(),
  working_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- pst_file:

  Character. Path to the .pst control file.

- method:

  Character. `"morris"` or `"sobol"`. Selects the algorithm by setting
  the `GSA_METHOD` option in the control file; pestpp-sen defaults to
  Morris, so without this the result would carry the requested label
  whatever was actually computed.

- exe:

  Character. Path to pestpp-sen executable.

- extra_args:

  Named list. Additional options.

- working_dir:

  Character. Working directory.

- verbose:

  Logical. Print output.

## Value

A list of class `pesto_sen_result`.

## Examples

``` r
# \donttest{
if (nzchar(Sys.which("pestpp-sen"))) {
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
  res <- pesto_sensitivity(tf, method = "morris", verbose = FALSE)
  res$method
}
# }
```
