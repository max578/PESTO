# Coerce an object into a `pesto_forward_model`

Generic used by the IES driver so a caller can pass either a bare
`function(theta) -> obs` or a fully-specified
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md).
A bare function is wrapped with the supplied policy arguments; an
existing `pesto_forward_model` is returned unchanged (its own policy is
authoritative and the `...` are ignored).

## Usage

``` r
as_forward_model(x, ...)
```

## Arguments

- x:

  A function or a `pesto_forward_model`.

- ...:

  Policy arguments forwarded to
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  when `x` is a bare function (e.g. `on_failure`, `parallel`, `n_obs`).

## Value

A `pesto_forward_model` S7 object.

## See also

[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md),
[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md).

## Examples

``` r
fm <- as_forward_model(function(theta) theta, on_failure = "stop")
S7::S7_inherits(fm, pesto_forward_model)
#> [1] TRUE
```
