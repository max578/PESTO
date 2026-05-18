# Convert a PESTO ensemble result into a `pesto_ensemble_manifest`

Wraps the plain-list result returned by
[`pesto_ies_callback()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies_callback.md)
(and, eventually,
[`pesto_ies()`](https://AAGI-AUS.github.io/PESTO/reference/pesto_ies.md))
in the S7 manifest contract object so downstream packages can consume it
via S7 dispatch without reaching into PESTO-specific list internals.

## Usage

``` r
as_manifest(x, ...)
```

## Arguments

- x:

  A `pesto_ies_callback_result` (or any object with a method registered
  against this generic).

- ...:

  Method-specific arguments. For `pesto_ies_callback_result`: `run_id`
  (character, defaults to a timestamp+hash slug), `seed` (integer,
  defaults to `NA_integer_`), `fidelity` (named numeric or `NULL`),
  `apsim_version` (character, defaults to `NA_character_`).

## Value

A `pesto_ensemble_manifest` S7 object.
