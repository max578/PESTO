# Forward Model from a System of Ordinary Differential Equations

Each calibration parameter is a named entry of the `theta` matrix
supplied by the IES driver. `param_names` declares which columns the
model consumes and in what order. For each realisation the template
integrates \$\$\frac{d\mathbf{y}}{dt} = f(t, \mathbf{y}, \theta)\$\$
from `times[1]` with initial state `y0` – itself possibly a function of
`theta` – and applies `observe()` to the state trajectory to produce a
length-`n_obs` numeric vector.

## Usage

``` r
ode_forward_model(
  derivs,
  y0,
  times,
  param_names = character(0),
  observe = NULL,
  solver = c("rk4", "desolve"),
  n_steps = 10L,
  desolve_method = "lsoda",
  n_obs = NA_integer_,
  ...
)
```

## Arguments

- derivs:

  Function `function(t, y, theta) -> dydt`. The ODE right hand side;
  returns a numeric vector the same length as the state `y`. `theta` is
  the named parameter vector for one realisation.

- y0:

  Numeric vector, or `function(theta) -> numeric`. The initial state at
  `times[1]`. A function form lets the initial condition depend on
  calibrated parameters.

- times:

  Numeric vector of strictly increasing observation times, length at
  least two. The first entry is the initial time; the trajectory is
  recorded at every entry.

- param_names:

  Character vector of the parameter columns the model consumes from
  `theta`, in order. Empty (default) disables the column check and
  passes `theta` rows through verbatim.

- observe:

  Function `function(traj, theta) -> obs`. Maps the
  `length(times) x n_state` trajectory matrix to a length-`n_obs`
  observation vector. Default reads state variable one at every time
  after the first.

- solver:

  Character. `"rk4"` (default, self-contained) or `"desolve"` (delegates
  to [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html)).

- n_steps:

  Integer. Fixed RK4 sub-steps between successive observation times
  (default `10L`). Ignored for `solver = "desolve"`. Larger values trade
  speed for integration accuracy.

- desolve_method:

  Character. The
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) method
  when `solver = "desolve"` (default `"lsoda"`). Ignored for `"rk4"`.

- n_obs:

  Integer or `NA`. Known observation dimensionality. `NA` (default)
  infers it from the first successful realisation.

- ...:

  Further policy arguments forwarded to
  [`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
  (for example `on_failure`, `parallel`, `n_cores`, `fidelity`,
  `label`).

## Value

A
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
S7 object whose forward map integrates the ODE system. Pass it to
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
as `forward_model`, evaluate it directly with
[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md),
or bundle several across fidelity levels with
[`pesto_multifidelity_model()`](https://max578.github.io/PESTO/reference/pesto_multifidelity_model.md).

## Details

Wraps a user-supplied ODE right-hand side in a typed
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
whose forward map integrates the system over a fixed time grid and reads
the observation vector off the resulting state trajectory. This is the
generic `differential_equations` template that
[`crop_growth_forward_model()`](https://max578.github.io/PESTO/reference/crop_growth_forward_model.md)
and
[`seir_forward_model()`](https://max578.github.io/PESTO/reference/seir_forward_model.md)
specialise; use it directly for any compartmental or mechanistic ODE.

The derivative function `derivs` has signature
`function(t, y, theta) -> dydt`, returning a numeric vector the same
length as the state `y`. The initial state `y0` may be a plain numeric
vector (shared across realisations) or a function
`function(theta) -> y0` when the starting condition is itself calibrated
(for example an unknown initial inoculum). The observation map `observe`
has signature `function(traj, theta) -> obs`, where `traj` is an
`length(times) x length(y0)` matrix of states at each integration time
(column order matching `y0`); the default reads the first state variable
at every time after the first.

Two integrators are available. The default `solver = "rk4"` is a
self-contained classical fourth-order Runge-Kutta with `n_steps` fixed
sub-steps between successive observation times – no external dependency,
deterministic, and adequate for the smooth non-stiff systems these
templates target. `solver = "desolve"` delegates to
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) (an
optional `Suggests` dependency) with method `desolve_method`, which
brings adaptive step control and stiff solvers (`"lsoda"`, the default,
switches between non-stiff and stiff automatically). A realisation whose
integration fails or returns a non-finite trajectory is reported as an
`NA` row, which
[`pesto_evaluate()`](https://max578.github.io/PESTO/reference/pesto_evaluate.md)
and
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
handle under their failure policy.

## References

Soetaert, K., Petzoldt, T. & Setzer, R. W. (2010). Solving differential
equations in R: package deSolve. *Journal of Statistical Software*,
33(9), 1–25.
[doi:10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09)

## See also

[`crop_growth_forward_model()`](https://max578.github.io/PESTO/reference/crop_growth_forward_model.md)
and
[`seir_forward_model()`](https://max578.github.io/PESTO/reference/seir_forward_model.md)
for ready-made specialisations;
[`pesto_forward_model()`](https://max578.github.io/PESTO/reference/pesto_forward_model.md)
for the contract;
[`pesto_ies_callback()`](https://max578.github.io/PESTO/reference/pesto_ies_callback.md)
for the IES driver.

## Examples

``` r
# Exponential decay dy/dt = -k y, observed at five times. Calibrate k.
fm <- ode_forward_model(
  derivs      = function(t, y, theta) -theta[["k"]] * y,
  y0          = c(y = 1),
  times       = seq(0, 4, by = 1),
  param_names = "k"
)
theta <- matrix(c(0.5, 1.0), ncol = 1L, dimnames = list(NULL, "k"))
pesto_evaluate(fm, theta)
#>           [,1]      [,2]      [,3]       [,4]
#> [1,] 0.6065307 0.3678795 0.2231302 0.13533530
#> [2,] 0.3678798 0.1353355 0.0497872 0.01831571
#> attr(,"n_failures")
#> [1] 0
#> attr(,"fail_idx")
#> integer(0)
```
