# Read a PEST Control File (.pst)

Parses a PEST/PEST++ control file and returns a structured list
containing all sections: control data, parameter groups, parameter data,
observation groups, observation data, model command, and
template/instruction file pairs.

## Usage

``` r
read_pst(file)
```

## Arguments

- file:

  Character. Path to the .pst file.

## Value

A list of class `pesto_pst` containing:

- control_data:

  Control section parameters (NPAR, NOBS, etc.)

- parameter_groups:

  data.table of parameter group definitions

- parameters:

  data.table of parameter data

- observation_groups:

  data.table of observation group definitions

- observations:

  data.table of observation data

- model_command:

  Character vector of model command lines

- template_files:

  data.table of template/model input file pairs

- instruction_files:

  data.table of instruction/model output file pairs

- prior_information:

  data.table of prior information (if present)

- pestpp_options:

  Named list of ++ options

## See also

[`write_pst()`](https://AAGI-AUS.github.io/PESTO/reference/write_pst.md),
[`create_pest_scenario()`](https://AAGI-AUS.github.io/PESTO/reference/create_pest_scenario.md)

## Examples

``` r
pars <- data.table::data.table(
  parnme = c("k1", "k2", "k3"),
  partrans = c("log", "log", "none"),
  parchglim = "factor",
  parval1 = c(1.0, 0.5, 0.1),
  parlbnd = c(0.01, 0.001, 0.0),
  parubnd = c(100, 50, 1.0),
  pargp = c("hydraulic", "hydraulic", "storage")
)
obs <- data.table::data.table(
  obsnme = c("h1", "h2", "h3"),
  obsval = c(1.0, 2.0, 1.5),
  weight = c(1.0, 1.0, 0.5),
  obgnme = "head"
)
pst <- create_pest_scenario(pars, obs, model_command = "echo run")
tf <- tempfile(fileext = ".pst")
on.exit(unlink(tf), add = TRUE)
write_pst(pst, tf)
pst_back <- read_pst(tf)
pst_back$control_data$npar
#> [1] 3
pst_back$control_data$nobs
#> [1] 3
```
