# Is an object a PESTO abstention?

Is an object a PESTO abstention?

## Usage

``` r
is_pesto_abstention(x)
```

## Arguments

- x:

  Any object.

## Value

A single logical: `TRUE` when `x` is a `pesto_abstention`.

## Examples

``` r
is_pesto_abstention(pesto_abstention("over_determined"))
#> [1] TRUE
is_pesto_abstention(42)
#> [1] FALSE
```
