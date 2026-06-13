# ode_templates.R -- Ready-to-use ODE / compartmental forward-model
# templates for the `differential_equations` model structure.
#
# PESTO's IES machinery calibrates any `function(theta) -> obs` forward
# model (see [pesto_forward_model()]). A very common forward model in the
# environmental, epidemiological, and crop sciences is a small system of
# ordinary differential equations integrated forward in time, with the
# observation vector read off the state trajectory at a set of times.
# This file supplies that machinery as drop-in templates -- the ODE
# analogue of the `apsimx` callback adapter -- so a caller writes only
# the derivative function (or picks a named template) and gets back a
# typed [pesto_forward_model()] that plugs straight into
# [pesto_ies_callback()], a [pesto_multifidelity_model()] stack, and the
# manifest emitter.
#
# Integration is by a self-contained fixed-step RK4 by default (no new
# hard dependency); a `solver = "desolve"` path delegates to the `deSolve`
# package (a Suggests-guarded optional dependency) for stiff systems and
# adaptive step control. References inline in the roxygen blocks.

#' Forward Model from a System of Ordinary Differential Equations
#'
#' Wraps a user-supplied ODE right-hand side in a typed
#' [pesto_forward_model()] whose forward map integrates the system over a
#' fixed time grid and reads the observation vector off the resulting
#' state trajectory. This is the generic `differential_equations`
#' template that [crop_growth_forward_model()] and
#' [seir_forward_model()] specialise; use it directly for any
#' compartmental or mechanistic ODE.
#'
#' @description
#' Each calibration parameter is a named entry of the `theta` matrix
#' supplied by the IES driver. `param_names` declares which columns the
#' model consumes and in what order. For each realisation the template
#' integrates
#' \deqn{\frac{d\mathbf{y}}{dt} = f(t, \mathbf{y}, \theta)}
#' from `times[1]` with initial state `y0` -- itself possibly a function
#' of `theta` -- and applies `observe()` to the state trajectory to
#' produce a length-`n_obs` numeric vector.
#'
#' @details
#' The derivative function `derivs` has signature
#' `function(t, y, theta) -> dydt`, returning a numeric vector the same
#' length as the state `y`. The initial state `y0` may be a plain numeric
#' vector (shared across realisations) or a function
#' `function(theta) -> y0` when the starting condition is itself
#' calibrated (for example an unknown initial inoculum). The observation
#' map `observe` has signature `function(traj, theta) -> obs`, where
#' `traj` is an `length(times) x length(y0)` matrix of states at each
#' integration time (column order matching `y0`); the default reads the
#' first state variable at every time after the first.
#'
#' Two integrators are available. The default `solver = "rk4"` is a
#' self-contained classical fourth-order Runge-Kutta with `n_steps`
#' fixed sub-steps between successive observation times -- no external
#' dependency, deterministic, and adequate for the smooth non-stiff
#' systems these templates target. `solver = "desolve"` delegates to
#' `deSolve::ode()` (an optional `Suggests` dependency) with method
#' `desolve_method`, which brings adaptive step control and stiff solvers
#' (`"lsoda"`, the default, switches between non-stiff and stiff
#' automatically). A realisation whose integration fails or returns a
#' non-finite trajectory is reported as an `NA` row, which
#' [pesto_evaluate()] and [pesto_ies_callback()] handle under their
#' failure policy.
#'
#' @param derivs Function `function(t, y, theta) -> dydt`. The ODE right
#'   hand side; returns a numeric vector the same length as the state
#'   `y`. `theta` is the named parameter vector for one realisation.
#' @param y0 Numeric vector, or `function(theta) -> numeric`. The initial
#'   state at `times[1]`. A function form lets the initial condition
#'   depend on calibrated parameters.
#' @param times Numeric vector of strictly increasing observation times,
#'   length at least two. The first entry is the initial time; the
#'   trajectory is recorded at every entry.
#' @param param_names Character vector of the parameter columns the model
#'   consumes from `theta`, in order. Empty (default) disables the column
#'   check and passes `theta` rows through verbatim.
#' @param observe Function `function(traj, theta) -> obs`. Maps the
#'   `length(times) x n_state` trajectory matrix to a length-`n_obs`
#'   observation vector. Default reads state variable one at every time
#'   after the first.
#' @param solver Character. `"rk4"` (default, self-contained) or
#'   `"desolve"` (delegates to `deSolve::ode()`).
#' @param n_steps Integer. Fixed RK4 sub-steps between successive
#'   observation times (default `10L`). Ignored for `solver = "desolve"`.
#'   Larger values trade speed for integration accuracy.
#' @param desolve_method Character. The `deSolve::ode()` method when
#'   `solver = "desolve"` (default `"lsoda"`). Ignored for `"rk4"`.
#' @param n_obs Integer or `NA`. Known observation dimensionality. `NA`
#'   (default) infers it from the first successful realisation.
#' @param ... Further policy arguments forwarded to
#'   [pesto_forward_model()] (for example `on_failure`, `parallel`,
#'   `n_cores`, `fidelity`, `label`).
#'
#' @return A [pesto_forward_model()] S7 object whose forward map
#'   integrates the ODE system. Pass it to [pesto_ies_callback()] as
#'   `forward_model`, evaluate it directly with [pesto_evaluate()], or
#'   bundle several across fidelity levels with
#'   [pesto_multifidelity_model()].
#' @references
#' Soetaert, K., Petzoldt, T. & Setzer, R. W. (2010). Solving
#' differential equations in R: package deSolve. *Journal of Statistical
#' Software*, 33(9), 1--25. \doi{10.18637/jss.v033.i09}
#' @seealso [crop_growth_forward_model()] and [seir_forward_model()] for
#'   ready-made specialisations; [pesto_forward_model()] for the contract;
#'   [pesto_ies_callback()] for the IES driver.
#' @examples
#' # Exponential decay dy/dt = -k y, observed at five times. Calibrate k.
#' fm <- ode_forward_model(
#'   derivs      = function(t, y, theta) -theta[["k"]] * y,
#'   y0          = c(y = 1),
#'   times       = seq(0, 4, by = 1),
#'   param_names = "k"
#' )
#' theta <- matrix(c(0.5, 1.0), ncol = 1L, dimnames = list(NULL, "k"))
#' pesto_evaluate(fm, theta)
#' @export
ode_forward_model <- function(derivs,
                              y0,
                              times,
                              param_names = character(0),
                              observe = NULL,
                              solver = c("rk4", "desolve"),
                              n_steps = 10L,
                              desolve_method = "lsoda",
                              n_obs = NA_integer_,
                              ...) {
  solver <- match.arg(solver)
  .check_ode_inputs(derivs, y0, times, observe, n_steps, solver)

  times       <- as.numeric(times)
  param_names <- as.character(param_names)
  if (is.null(observe)) observe <- .ode_default_observe

  # Build the per-realisation forward map -------------------------------
  # The closure resolves the initial state (constant or theta-dependent),
  # integrates with the chosen solver, and applies the observation map.
  # Any failure collapses to NULL, which the bulk assembler renders as an
  # NA row -- the same failure contract the apsimx adapter honours.
  fn <- function(theta) {
    if (!is.matrix(theta)) theta <- as.matrix(theta)
    nreal <- nrow(theta)
    rows  <- vector("list", nreal)
    for (i in seq_len(nreal)) {
      theta_i <- .ode_named_theta(theta[i, , drop = FALSE], param_names)
      rows[[i]] <- tryCatch({
        y0_i <- if (is.function(y0)) as.numeric(y0(theta_i)) else as.numeric(y0)
        traj <- .ode_integrate(derivs, y0_i, times, theta_i, solver,
                               as.integer(n_steps), desolve_method)
        obs  <- as.numeric(observe(traj, theta_i))
        if (length(obs) == 0L || !all(is.finite(obs))) NULL else obs
      }, error = function(e) NULL)
    }
    .ode_bind_rows(rows, nreal)
  }

  pesto_forward_model(fn = fn, n_obs = n_obs, param_names = param_names, ...)
}


#' Crop-Growth Forward Model (Logistic / Expolinear Biomass)
#'
#' A ready-to-use [ode_forward_model()] specialisation for the canonical
#' single-state crop biomass-accumulation curve. Above-ground biomass
#' \eqn{B} follows a logistic growth law
#' \deqn{\frac{dB}{dt} = r\,B\left(1 - \frac{B}{B_{\max}}\right),}
#' the standard sigmoid description of a crop's dry-matter accumulation
#' over a season: an early near-exponential phase at relative growth rate
#' \eqn{r}, decelerating to a canopy- and resource-limited ceiling
#' \eqn{B_{\max}} (Goudriaan & Monteith 1990). The calibration
#' parameters are the relative growth rate `r`, the asymptotic biomass
#' `b_max`, and the initial biomass `b0`.
#'
#' @details
#' The forward map integrates the logistic ODE across `times` and returns
#' the modelled biomass at every time after the first -- exactly the
#' shape a destructive-harvest or remote-sensing biomass series takes, so
#' the returned object calibrates directly against an observed
#' biomass-over-time vector through [pesto_ies_callback()]. The
#' single-state logistic form is deliberately the simplest defensible
#' crop-growth template; richer multi-organ partitioning models compose
#' through the same [ode_forward_model()] entry point by supplying a
#' vector-valued `derivs`.
#'
#' @param times Numeric vector of strictly increasing observation times
#'   (for example thermal-time or days-after-sowing), length at least
#'   two. The first entry is the initial time; biomass is reported at
#'   every later entry.
#' @param solver Character. `"rk4"` (default) or `"desolve"`, as in
#'   [ode_forward_model()].
#' @param n_steps Integer. Fixed RK4 sub-steps between observation times
#'   (default `10L`).
#' @param ... Further policy arguments forwarded to
#'   [pesto_forward_model()] via [ode_forward_model()] (for example
#'   `on_failure`, `parallel`, `fidelity`).
#'
#' @return A [pesto_forward_model()] S7 object with `param_names`
#'   `c("r", "b_max", "b0")` and `n_obs = length(times) - 1L`.
#' @references
#' Goudriaan, J. & Monteith, J. L. (1990). A mathematical function for
#' crop growth based on light interception and leaf area expansion.
#' *Annals of Botany*, 66(6), 695--701.
#' @seealso [ode_forward_model()] for the generic builder;
#'   [seir_forward_model()] for the compartmental epidemic template;
#'   [pesto_ies_callback()] for calibration.
#' @examples
#' # Simulate a biomass series at a known parameter, then recover it.
#' times <- seq(0, 120, by = 15)
#' fm <- crop_growth_forward_model(times = times)
#' truth <- matrix(c(0.06, 1400, 20), nrow = 1L,
#'                 dimnames = list(NULL, c("r", "b_max", "b0")))
#' biomass <- as.numeric(pesto_evaluate(fm, truth))
#' round(biomass)
#' @export
crop_growth_forward_model <- function(times,
                                      solver = c("rk4", "desolve"),
                                      n_steps = 10L,
                                      ...) {
  solver <- match.arg(solver)
  times  <- as.numeric(times)

  # Logistic dry-matter accumulation. b_max is guarded away from zero so a
  # boundary-hugging realisation cannot divide by zero mid-integration.
  derivs <- function(t, y, theta) {
    b_max <- max(theta[["b_max"]], .Machine$double.eps)
    theta[["r"]] * y * (1 - y / b_max)
  }

  ode_forward_model(
    derivs      = derivs,
    y0          = function(theta) c(B = theta[["b0"]]),
    times       = times,
    param_names = c("r", "b_max", "b0"),
    observe     = .ode_default_observe,
    solver      = solver,
    n_steps     = as.integer(n_steps),
    n_obs       = length(times) - 1L,
    ...
  )
}


#' SEIR Compartmental Forward Model
#'
#' A ready-to-use [ode_forward_model()] specialisation for the classic
#' Susceptible-Exposed-Infectious-Recovered epidemic model on a closed
#' population (Anderson & May 1991). The four states evolve as
#' \deqn{\dot S = -\beta S I / N,\quad
#'       \dot E = \beta S I / N - \sigma E,}
#' \deqn{\dot I = \sigma E - \gamma I,\quad \dot R = \gamma I,}
#' with transmission rate \eqn{\beta}, latency rate \eqn{\sigma} (mean
#' incubation \eqn{1/\sigma}), and recovery rate \eqn{\gamma} (mean
#' infectious period \eqn{1/\gamma}). The basic reproduction number is
#' \eqn{R_0 = \beta / \gamma}.
#'
#' @details
#' The calibration parameters are `beta`, `sigma`, and `gamma`. The
#' population size `n_pop` and the initial infectious count `i0` are
#' fixed structural constants of the template (an outbreak seeded with
#' `i0` infectious individuals, `n_pop - i0` susceptible, and nobody
#' exposed or recovered). By default the forward map returns the
#' **infectious prevalence** \eqn{I(t)} at every time after the first --
#' the compartment a case-count series tracks -- so the object calibrates
#' directly against an observed epidemic curve through
#' [pesto_ies_callback()]. Supply a custom `observe` through
#' [ode_forward_model()] to read a different compartment (for example the
#' incidence \eqn{\sigma E}).
#'
#' @param times Numeric vector of strictly increasing observation times
#'   (days), length at least two. The first entry is the outbreak start.
#' @param n_pop Numeric. Total (closed) population size. Default `1000`.
#' @param i0 Numeric. Initial infectious count at `times[1]`. Default
#'   `1`. Must be positive and below `n_pop`.
#' @param solver Character. `"rk4"` (default) or `"desolve"`, as in
#'   [ode_forward_model()].
#' @param n_steps Integer. Fixed RK4 sub-steps between observation times
#'   (default `10L`).
#' @param ... Further policy arguments forwarded to
#'   [pesto_forward_model()] via [ode_forward_model()].
#'
#' @return A [pesto_forward_model()] S7 object with `param_names`
#'   `c("beta", "sigma", "gamma")` and `n_obs = length(times) - 1L`,
#'   emitting the infectious prevalence trajectory.
#' @references
#' Anderson, R. M. & May, R. M. (1991). *Infectious Diseases of Humans:
#' Dynamics and Control*. Oxford University Press.
#' @seealso [ode_forward_model()] for the generic builder;
#'   [crop_growth_forward_model()] for the crop template;
#'   [pesto_ies_callback()] for calibration.
#' @examples
#' # Simulate an outbreak curve at a known (beta, sigma, gamma).
#' times <- seq(0, 60, by = 5)
#' fm <- seir_forward_model(times = times, n_pop = 1000, i0 = 1)
#' truth <- matrix(c(0.6, 0.2, 0.1), nrow = 1L,
#'                 dimnames = list(NULL, c("beta", "sigma", "gamma")))
#' round(as.numeric(pesto_evaluate(fm, truth)), 1)
#' @export
seir_forward_model <- function(times,
                               n_pop = 1000,
                               i0 = 1,
                               solver = c("rk4", "desolve"),
                               n_steps = 10L,
                               ...) {
  solver <- match.arg(solver)
  times  <- as.numeric(times)
  .assert_positive_scalar(n_pop, "n_pop")
  .assert_positive_scalar(i0, "i0")
  if (i0 >= n_pop) {
    stop("`i0` must be smaller than `n_pop`.", call. = FALSE)
  }

  # State order is (S, E, I, R); the closed-population SEIR right-hand
  # side. Compartment four (R) closes the conservation law S+E+I+R = N.
  derivs <- function(t, y, theta) {
    s <- y[[1L]]; e <- y[[2L]]; i <- y[[3L]]
    beta  <- theta[["beta"]]
    sigma <- theta[["sigma"]]
    gamma <- theta[["gamma"]]
    new_inf <- beta * s * i / n_pop
    c(
      -new_inf,
      new_inf - sigma * e,
      sigma * e - gamma * i,
      gamma * i
    )
  }

  # Read the infectious compartment (state 3) at every time after t0.
  observe_infectious <- function(traj, theta) {
    traj[-1L, 3L]
  }

  ode_forward_model(
    derivs      = derivs,
    y0          = c(S = n_pop - i0, E = 0, I = i0, R = 0),
    times       = times,
    param_names = c("beta", "sigma", "gamma"),
    observe     = observe_infectious,
    solver      = solver,
    n_steps     = as.integer(n_steps),
    n_obs       = length(times) - 1L,
    ...
  )
}


# Internal: validate the shared ODE-template inputs ---------------------

#' Validate inputs common to the ODE forward-model templates
#'
#' Checks the derivative callable, the initial state (vector or
#' `function(theta)`), the time grid (numeric, strictly increasing,
#' length >= 2), the optional observation map, the sub-step count, and --
#' for `solver = "desolve"` -- the optional `deSolve` Suggests dependency.
#'
#' @noRd
#' @keywords internal
.check_ode_inputs <- function(derivs, y0, times, observe, n_steps, solver) {
  .assert_function(derivs, "derivs")
  if (!is.function(y0) && !is.numeric(y0)) {
    stop("`y0` must be a numeric vector or a function(theta).",
         call. = FALSE)
  }
  .assert_numeric_vector(times, "times")
  if (length(times) < 2L || any(!is.finite(times)) ||
      any(diff(times) <= 0)) {
    stop("`times` must be a strictly increasing numeric vector of length >= 2.",
         call. = FALSE)
  }
  if (!is.null(observe) && !is.function(observe)) {
    stop("`observe` must be a function(traj, theta) or NULL.",
         call. = FALSE)
  }
  if (length(n_steps) != 1L || !is.finite(n_steps) ||
      as.integer(n_steps) < 1L) {
    stop("`n_steps` must be a single positive integer.", call. = FALSE)
  }
  if (solver == "desolve" && !requireNamespace("deSolve", quietly = TRUE)) {
    stop(
      paste0(
        "`solver = \"desolve\"` requires the 'deSolve' package. ",
        "Install it, or use the default solver = \"rk4\"."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# Internal: attach parameter names to a single theta row ----------------

#' Coerce one `theta` row to a named numeric vector for `derivs`
#'
#' @noRd
#' @keywords internal
.ode_named_theta <- function(theta_row, param_names) {
  v <- as.numeric(theta_row)
  nm <- colnames(theta_row)
  if (length(param_names) > 0L) {
    nm <- param_names
  } else if (is.null(nm)) {
    nm <- paste0("par", seq_along(v))
  }
  names(v) <- nm
  v
}


# Internal: integrate one realisation over the time grid ----------------

#' Integrate one ODE realisation, dispatching on the chosen solver
#'
#' Returns a `length(times) x n_state` trajectory matrix (states at each
#' observation time, including the initial time as row one).
#'
#' @noRd
#' @keywords internal
.ode_integrate <- function(derivs, y0, times, theta, solver, n_steps,
                           desolve_method) {
  if (solver == "desolve") {
    return(.ode_integrate_desolve(derivs, y0, times, theta, desolve_method))
  }
  .ode_integrate_rk4(derivs, y0, times, theta, n_steps)
}


#' Self-contained fixed-step classical RK4 integrator
#'
#' Advances the state from each observation time to the next in
#' `n_steps` equal RK4 sub-steps, recording the state at every
#' observation time. No external dependency.
#'
#' @noRd
#' @keywords internal
.ode_integrate_rk4 <- function(derivs, y0, times, theta, n_steps) {
  n_state <- length(y0)
  n_time  <- length(times)
  traj    <- matrix(NA_real_, nrow = n_time, ncol = n_state)
  traj[1L, ] <- y0
  y <- y0

  for (k in seq_len(n_time - 1L)) {
    h <- (times[k + 1L] - times[k]) / n_steps
    t <- times[k]
    for (s in seq_len(n_steps)) {
      k1 <- as.numeric(derivs(t, y, theta))
      k2 <- as.numeric(derivs(t + h / 2, y + h / 2 * k1, theta))
      k3 <- as.numeric(derivs(t + h / 2, y + h / 2 * k2, theta))
      k4 <- as.numeric(derivs(t + h, y + h * k3, theta))
      y <- y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
      t <- t + h
    }
    traj[k + 1L, ] <- y
  }
  traj
}


#' Integrate via `deSolve::ode()` (optional Suggests path)
#'
#' Delegates to `deSolve::ode()` with the requested method and returns
#' the state columns aligned to the template's trajectory convention
#' (rows = observation times, columns = states, no time column).
#'
#' @noRd
#' @keywords internal
.ode_integrate_desolve <- function(derivs, y0, times, theta, desolve_method) {
  rhs <- function(t, state, parms) {
    list(as.numeric(derivs(t, state, parms)))
  }
  out <- deSolve::ode(
    y      = y0,
    times  = times,
    func   = rhs,
    parms  = theta,
    method = desolve_method
  )
  # deSolve returns a matrix whose first column is time; drop it to match
  # the RK4 path's (n_time x n_state) state-only trajectory.
  state <- out[, -1L, drop = FALSE]
  matrix(as.numeric(state), nrow = length(times))
}


# Internal: default observation map -- state one after the initial time -

#' Default ODE observation map: state variable one at every later time
#'
#' @noRd
#' @keywords internal
.ode_default_observe <- function(traj, theta) {
  traj[-1L, 1L]
}


# Internal: assemble per-realisation outputs into an nreal x nobs matrix -

#' Bind ODE per-realisation observation vectors into a result matrix
#'
#' Mirrors the apsimx adapter's assembler: the observation width is the
#' longest successful realisation; a realisation that failed (`NULL`) or
#' returned the wrong width becomes an `NA` row.
#'
#' @noRd
#' @keywords internal
.ode_bind_rows <- function(rows, nreal) {
  widths <- vapply(rows,
                   function(v) if (is.null(v)) 0L else length(v),
                   integer(1L))
  nobs <- max(c(0L, widths))
  if (nobs == 0L) {
    stop(
      "ode_forward_model: no realisation produced a usable output vector.",
      call. = FALSE
    )
  }
  out <- matrix(NA_real_, nrow = nreal, ncol = nobs)
  for (i in seq_len(nreal)) {
    v <- rows[[i]]
    if (!is.null(v) && length(v) == nobs && all(is.finite(v))) {
      out[i, ] <- as.numeric(v)
    }
  }
  out
}
