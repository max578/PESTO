# Print method for pesto_pst objects

Print method for pesto_pst objects

## Usage

``` r
# S3 method for class 'pesto_pst'
print(x, ...)
```

## Arguments

- x:

  A `pesto_pst` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`. Called for the side effect of printing.

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
print(pst)
#> $control_data
#> $control_data$rstfle
#> [1] "restart"
#> 
#> $control_data$pestmode
#> [1] "estimation"
#> 
#> $control_data$npar
#> [1] 2
#> 
#> $control_data$nobs
#> [1] 2
#> 
#> $control_data$npargp
#> [1] 1
#> 
#> $control_data$nprior
#> [1] 0
#> 
#> $control_data$nobsgp
#> [1] 1
#> 
#> $control_data$ntplfle
#> [1] 0
#> 
#> $control_data$ninsfle
#> [1] 0
#> 
#> 
#> $parameters
#>    parnme partrans parchglim parval1 parlbnd parubnd     pargp scale offset
#>    <char>   <char>    <char>   <num>   <num>   <num>    <char> <num>  <num>
#> 1:     k1      log    factor     1.0   0.010     100 hydraulic     1      0
#> 2:     k2      log    factor     0.5   0.001      50 hydraulic     1      0
#>    dercom
#>     <int>
#> 1:      1
#> 2:      1
#> 
#> $observations
#>    obsnme obsval weight obgnme
#>    <char>  <num>  <num> <char>
#> 1:     h1      1      1   head
#> 2:     h2      2      1   head
#> 
#> $model_command
#> [1] "echo run"
#> 
#> $io_files
#> Null data.table (0 rows and 0 cols)
#> 
#> $pestpp_options
#> list()
#> 
#> $file
#> [1] NA
#> 
#> attr(,"class")
#> [1] "pesto_pst"
```
