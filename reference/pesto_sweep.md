# Run PEST++ SWP (Parametric Sweep)

Executes pestpp-swp for embarrassingly parallel model runs across a
parameter ensemble.

## Usage

``` r
pesto_sweep(
  pst_file,
  par_ensemble,
  exe = NULL,
  working_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- pst_file:

  Character. Path to the .pst control file.

- par_ensemble:

  data.table or path. Parameter ensemble.

- exe:

  Character. Path to pestpp-swp executable.

- working_dir:

  Character. Working directory.

- verbose:

  Logical. Print output.

## Value

A list containing observation outputs for each realisation.

## Examples

``` r
# \donttest{
if (nzchar(Sys.which("pestpp-swp"))) {
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
  par_ens <- data.table::data.table(
    real_name = c("r1", "r2", "r3"),
    k1 = c(0.8, 1.0, 1.2),
    k2 = c(0.4, 0.5, 0.6)
  )
  res <- pesto_sweep(tf, par_ensemble = par_ens, verbose = FALSE)
  res$exit_code
}
# }
```
