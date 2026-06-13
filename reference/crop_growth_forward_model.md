# Crop-Growth Forward Model (Logistic / Expolinear Biomass)

A ready-to-use
[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
specialisation for the canonical single-state crop biomass-accumulation
curve. Above-ground biomass \\B\\ follows a logistic growth law
\$\$\frac{dB}{dt} = r\\B\left(1 - \frac{B}{B\_{\max}}\right),\$\$ the
standard sigmoid description of a crop's dry-matter accumulation over a
season: an early near-exponential phase at relative growth rate \\r\\,
decelerating to a canopy- and resource-limited ceiling \\B\_{\max}\\
(Goudriaan & Monteith 1990). The calibration parameters are the relative
growth rate `r`, the asymptotic biomass `b_max`, and the initial biomass
`b0`.

## Usage

``` r
crop_growth_forward_model(
  times,
  solver = c("rk4", "desolve"),
  n_steps = 10L,
  ...
)
```

## Arguments

- times:

  Numeric vector of strictly increasing observation times (for example
  thermal-time or days-after-sowing), length at least two. The first
  entry is the initial time; biomass is reported at every later entry.

- solver:

  Character. `"rk4"` (default) or `"desolve"`, as in
  [`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md).

- n_steps:

  Integer. Fixed RK4 sub-steps between observation times (default
  `10L`).

- ...:

  Further policy arguments forwarded to
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  via
  [`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
  (for example `on_failure`, `parallel`, `fidelity`).

## Value

A
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
S7 object with `param_names` `c("r", "b_max", "b0")` and
`n_obs = length(times) - 1L`.

## Details

The forward map integrates the logistic ODE across `times` and returns
the modelled biomass at every time after the first – exactly the shape a
destructive-harvest or remote-sensing biomass series takes, so the
returned object calibrates directly against an observed
biomass-over-time vector through
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).
The single-state logistic form is deliberately the simplest defensible
crop-growth template; richer multi-organ partitioning models compose
through the same
[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
entry point by supplying a vector-valued `derivs`.

## References

Goudriaan, J. & Monteith, J. L. (1990). A mathematical function for crop
growth based on light interception and leaf area expansion. *Annals of
Botany*, 66(6), 695–701.

## See also

[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
for the generic builder;
[`seir_forward_model()`](https://max578.github.io/PESTO/reference/seir_forward_model.md)
for the compartmental epidemic template;
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
for calibration.

## Examples

``` r
# Simulate a biomass series at a known parameter, then recover it.
times <- seq(0, 120, by = 15)
fm <- crop_growth_forward_model(times = times)
truth <- matrix(c(0.06, 1400, 20), nrow = 1L,
                dimnames = list(NULL, c("r", "b_max", "b0")))
biomass <- as.numeric(pesto_evaluate(fm, truth))
round(biomass)
#> [1]   48  113  248  485  793 1067 1243 1331
```
