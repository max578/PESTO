# SEIR Compartmental Forward Model

A ready-to-use
[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
specialisation for the classic Susceptible-Exposed-Infectious-Recovered
epidemic model on a closed population (Anderson & May 1991). The four
states evolve as \$\$\dot S = -\beta S I / N,\quad \dot E = \beta S I /
N - \sigma E,\$\$ \$\$\dot I = \sigma E - \gamma I,\quad \dot R = \gamma
I,\$\$ with transmission rate \\\beta\\, latency rate \\\sigma\\ (mean
incubation \\1/\sigma\\), and recovery rate \\\gamma\\ (mean infectious
period \\1/\gamma\\). The basic reproduction number is \\R_0 = \beta /
\gamma\\.

## Usage

``` r
seir_forward_model(
  times,
  n_pop = 1000,
  i0 = 1,
  solver = c("rk4", "desolve"),
  n_steps = 10L,
  ...
)
```

## Arguments

- times:

  Numeric vector of strictly increasing observation times (days), length
  at least two. The first entry is the outbreak start.

- n_pop:

  Numeric. Total (closed) population size. Default `1000`.

- i0:

  Numeric. Initial infectious count at `times[1]`. Default `1`. Must be
  positive and below `n_pop`.

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
  [`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md).

## Value

A
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
S7 object with `param_names` `c("beta", "sigma", "gamma")` and
`n_obs = length(times) - 1L`, emitting the infectious prevalence
trajectory.

## Details

The calibration parameters are `beta`, `sigma`, and `gamma`. The
population size `n_pop` and the initial infectious count `i0` are fixed
structural constants of the template (an outbreak seeded with `i0`
infectious individuals, `n_pop - i0` susceptible, and nobody exposed or
recovered). By default the forward map returns the **infectious
prevalence** \\I(t)\\ at every time after the first – the compartment a
case-count series tracks – so the object calibrates directly against an
observed epidemic curve through
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md).
Supply a custom `observe` through
[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
to read a different compartment (for example the incidence \\\sigma
E\\).

## References

Anderson, R. M. & May, R. M. (1991). *Infectious Diseases of Humans:
Dynamics and Control*. Oxford University Press.

## See also

[`ode_forward_model()`](https://max578.github.io/PESTO/reference/ode_forward_model.md)
for the generic builder;
[`crop_growth_forward_model()`](https://max578.github.io/PESTO/reference/crop_growth_forward_model.md)
for the crop template;
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
for calibration.

## Examples

``` r
# Simulate an outbreak curve at a known (beta, sigma, gamma).
times <- seq(0, 60, by = 5)
fm <- seir_forward_model(times = times, n_pop = 1000, i0 = 1)
truth <- matrix(c(0.6, 0.2, 0.1), nrow = 1L,
                dimnames = list(NULL, c("beta", "sigma", "gamma")))
round(as.numeric(pesto_evaluate(fm, truth)), 1)
#>  [1]   1.6   4.2  11.2  29.4  72.5 156.3 265.4 333.1 321.0 258.8 187.5 127.5
```
