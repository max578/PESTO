# forward_model.R -- Typed forward-model contract for PESTO's two
# adapter modes.
#
# Both the native R callback path ([pesto_ies_callback()]) and the
# classic `.pst`-file path ([pesto_ies()]) ultimately need the same
# thing: a way to turn an `nreal x npar` parameter matrix into an
# `nreal x nobs` observation matrix, tolerating per-realisation
# failures. Historically that contract lived as an implicit convention
# enforced ad hoc inside the callback driver. This file promotes it to
# a typed S7 object so the contract is enforceable, carries its own
# failure policy and evaluation strategy, and tags its fidelity level
# for the multi-fidelity bridge ([pesto_multifidelity_model()]).

#' Forward-Model Contract (S7 class)
#'
#' A `pesto_forward_model` wraps a user callable of signature
#' `function(theta) -> obs` -- where `theta` is an `nreal x npar`
#' numeric matrix and `obs` is an `nreal x nobs` numeric matrix -- in a
#' typed object that owns the evaluation contract: output dimensionality,
#' expected parameter names, the failure policy, the concurrency
#' strategy, and a fidelity tag.
#'
#' This is the single contract both PESTO adapter modes honour. The
#' native callback driver [pesto_ies_callback()] accepts either a bare
#' function (auto-wrapped via [as_forward_model()] with that driver's
#' failure policy) or a `pesto_forward_model` built here; the apsimx
#' adapter [apsim_callback()] can emit one directly. Evaluation is via
#' [pesto_evaluate()], which guarantees the returned shape and accounts
#' for failed realisations as `NA` rows.
#'
#' @section Concurrency:
#' With `parallel = "serial"` (the default) and no `map_fn`, evaluation
#' attempts a single bulk call `fn(theta)` and falls back to a serial
#' per-realisation loop only if the bulk call errors -- preserving the
#' fast path for vectorised forward models. With `parallel = "multicore"`
#' realisations are dispatched per row through
#' [parallel::mclapply()] (fork-based; silently serial on Windows). A
#' custom `map_fn` (an `lapply`-shaped `function(X, FUN, ...)`) overrides
#' both and lets callers plug in `future.apply::future_lapply`,
#' `mirai`, or a cluster backend. For reproducible parallel runs set
#' `RNGkind("L'Ecuyer-CMRG")` and `set.seed()` before evaluating;
#' [parallel::mclapply()] then draws independent streams per realisation.
#'
#' @section Fidelity:
#' `fidelity` is an integer level tag (`0L` = base / cheapest by
#' convention; higher = more expensive / higher resolution). A
#' single-fidelity model leaves it at `0L`. The tag is what
#' [pesto_multifidelity_model()] uses to order levels and what the
#' manifest emitter threads into provenance.
#'
#' @param fn Function. The forward model, signature
#'   `function(theta) -> obs`.
#' @param n_obs Integer or `NA`. Known observation dimensionality. If
#'   `NA` (default) it is inferred from the first successful evaluation.
#' @param param_names Character. Expected parameter column names. Empty
#'   (default) disables the column check.
#' @param on_failure Character. `"na"` (default) records failed
#'   realisations as `NA` rows and proceeds; `"stop"` aborts on any
#'   failure.
#' @param parallel Character. `"serial"` (default) or `"multicore"`.
#' @param n_cores Integer or `NA`. Worker count for `"multicore"`. `NA`
#'   (default) resolves to `parallel::detectCores() - 1L` at evaluation
#'   time.
#' @param map_fn Function or `NULL`. Optional `lapply`-shaped override
#'   `function(X, FUN, ...)`; when supplied it drives per-realisation
#'   dispatch regardless of `parallel`.
#' @param max_fail_frac Numeric in `[0, 1]`. Abort if the fraction of
#'   failed realisations in any single evaluation exceeds this. Default
#'   `1` (never abort on fraction; `on_failure` still governs the
#'   zero-success case).
#' @param fidelity Integer. Fidelity level tag (default `0L`).
#' @param label Character. Optional human label carried into diagnostics.
#'
#' @return A `pesto_forward_model` S7 object.
#' @seealso [pesto_evaluate()] to evaluate one; [as_forward_model()] to
#'   coerce a bare function; [pesto_multifidelity_model()] to compose
#'   several across fidelity levels; [pesto_ies_callback()] for the IES
#'   driver that consumes it.
#' @examples
#' # A vectorised linear forward model wrapped as a contract object.
#' G  <- matrix(c(1, 0, 0, 1, 1, -1), nrow = 3L, byrow = TRUE)
#' fm <- pesto_forward_model(
#'   fn          = function(theta) theta %*% t(G),
#'   n_obs       = 3L,
#'   param_names = c("a", "b")
#' )
#' theta <- matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE,
#'                 dimnames = list(NULL, c("a", "b")))
#' pesto_evaluate(fm, theta)
#' @export
pesto_forward_model <- S7::new_class(
  "pesto_forward_model",
  package = "PESTO",
  properties = list(
    fn            = S7::class_function,
    n_obs         = S7::new_property(S7::class_integer,
                                     default = NA_integer_),
    param_names   = S7::new_property(S7::class_character,
                                     default = character(0)),
    on_failure    = S7::new_property(S7::class_character,
                                     default = "na"),
    parallel      = S7::new_property(S7::class_character,
                                     default = "serial"),
    n_cores       = S7::new_property(S7::class_integer,
                                     default = NA_integer_),
    map_fn        = S7::new_property(S7::class_any, default = NULL),
    max_fail_frac = S7::new_property(S7::class_numeric, default = 1),
    fidelity      = S7::new_property(S7::class_integer, default = 0L),
    label         = S7::new_property(S7::class_character,
                                     default = NA_character_)
  ),
  constructor = function(fn, n_obs = NA_integer_,
                         param_names = character(0),
                         on_failure = "na", parallel = "serial",
                         n_cores = NA_integer_, map_fn = NULL,
                         max_fail_frac = 1, fidelity = 0L,
                         label = NA_character_) {
    # Coerce the integer-valued fields so callers can pass plain numeric
    # literals (`n_obs = 6`) without tripping S7's strict integer class.
    S7::new_object(
      S7::S7_object(),
      fn            = fn,
      n_obs         = as.integer(n_obs),
      param_names   = as.character(param_names),
      on_failure    = on_failure,
      parallel      = parallel,
      n_cores       = as.integer(n_cores),
      map_fn        = map_fn,
      max_fail_frac = as.numeric(max_fail_frac),
      fidelity      = as.integer(fidelity),
      label         = as.character(label)
    )
  },
  validator = function(self) {
    errs <- character(0)

    # Output dimensionality -------------------------------------------
    if (length(self@n_obs) != 1L ||
        (!is.na(self@n_obs) && self@n_obs < 1L)) {
      errs <- c(errs, "`n_obs` must be a single positive integer or NA")
    }

    # Failure / concurrency policy ------------------------------------
    if (length(self@on_failure) != 1L ||
        !self@on_failure %in% c("na", "stop")) {
      errs <- c(errs, "`on_failure` must be one of \"na\", \"stop\"")
    }
    if (length(self@parallel) != 1L ||
        !self@parallel %in% c("serial", "multicore")) {
      errs <- c(errs, "`parallel` must be one of \"serial\", \"multicore\"")
    }
    if (length(self@n_cores) != 1L ||
        (!is.na(self@n_cores) && self@n_cores < 1L)) {
      errs <- c(errs, "`n_cores` must be a single positive integer or NA")
    }
    if (!is.null(self@map_fn) && !is.function(self@map_fn)) {
      errs <- c(errs, "`map_fn` must be a function or NULL")
    }

    # Abort threshold + fidelity --------------------------------------
    if (length(self@max_fail_frac) != 1L ||
        !is.finite(self@max_fail_frac) ||
        self@max_fail_frac < 0 || self@max_fail_frac > 1) {
      errs <- c(errs, "`max_fail_frac` must be a single number in [0, 1]")
    }
    if (length(self@fidelity) != 1L || is.na(self@fidelity) ||
        self@fidelity < 0L) {
      errs <- c(errs, "`fidelity` must be a single non-negative integer")
    }

    if (length(errs) == 0L) NULL else paste(errs, collapse = "; ")
  }
)


#' Coerce an object into a `pesto_forward_model`
#'
#' Generic used by the IES driver so a caller can pass either a bare
#' `function(theta) -> obs` or a fully-specified [pesto_forward_model()].
#' A bare function is wrapped with the supplied policy arguments; an
#' existing `pesto_forward_model` is returned unchanged (its own policy
#' is authoritative and the `...` are ignored).
#'
#' @param x A function or a `pesto_forward_model`.
#' @param ... Policy arguments forwarded to [pesto_forward_model()] when
#'   `x` is a bare function (e.g. `on_failure`, `parallel`, `n_obs`).
#' @return A `pesto_forward_model` S7 object.
#' @seealso [pesto_forward_model()], [pesto_evaluate()].
#' @examples
#' fm <- as_forward_model(function(theta) theta, on_failure = "stop")
#' S7::S7_inherits(fm, pesto_forward_model)
#' @export
as_forward_model <- S7::new_generic("as_forward_model", "x")

S7::method(as_forward_model, S7::class_function) <-
  function(x, ...) {
    pesto_forward_model(fn = x, ...)
  }

S7::method(as_forward_model, pesto_forward_model) <-
  function(x, ...) {
    x
  }


#' Evaluate a PESTO forward model
#'
#' Runs the forward model on a parameter matrix under its own failure
#' policy and concurrency strategy, returning a shape-guaranteed
#' `nreal x nobs` observation matrix. Failed realisations populate `NA`
#' rows (when `on_failure = "na"`). The returned matrix carries two
#' attributes: `"n_failures"` (integer count of `NA` rows) and
#' `"fail_idx"` (integer realisation indices that failed).
#'
#' @param model A `pesto_forward_model` (or, for the multi-fidelity
#'   method, a `pesto_multifidelity_model`).
#' @param theta Numeric matrix, `nreal x npar`. Column names, when
#'   present, are checked against `model@param_names`.
#' @param ... Method-specific arguments. The multi-fidelity method
#'   accepts `level` (integer fidelity level to evaluate).
#' @return An `nreal x nobs` numeric matrix with attributes
#'   `"n_failures"` and `"fail_idx"`.
#' @seealso [pesto_forward_model()], [pesto_multifidelity_model()].
#' @examples
#' fm <- pesto_forward_model(fn = function(theta) theta[, 1, drop = FALSE],
#'                           n_obs = 1L)
#' pesto_evaluate(fm, matrix(c(1, 2, 3), ncol = 1L))
#' @export
pesto_evaluate <- S7::new_generic(
  "pesto_evaluate", "model",
  function(model, theta, ...) S7::S7_dispatch()
)

S7::method(pesto_evaluate, pesto_forward_model) <-
  function(model, theta, ...) {
    theta <- .fm_coerce_theta(theta, model@param_names)
    .fm_eval_matrix(model, theta)
  }


# Coerce + validate theta into a double matrix with checked columns.
.fm_coerce_theta <- function(theta, param_names) {
  if (!is.matrix(theta)) theta <- as.matrix(theta)
  storage.mode(theta) <- "double"

  if (length(param_names) > 0L) {
    cols <- colnames(theta)
    if (is.null(cols)) {
      if (ncol(theta) != length(param_names)) {
        stop(
          sprintf(
            "`theta` has %d unnamed columns; the model expects %d (%s).",
            ncol(theta), length(param_names),
            paste(param_names, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      colnames(theta) <- param_names
    } else {
      missing_pars <- setdiff(param_names, cols)
      if (length(missing_pars) > 0L) {
        stop(
          sprintf(
            "`theta` is missing parameters required by the model: %s.",
            paste(missing_pars, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      theta <- theta[, param_names, drop = FALSE]
    }
  }
  theta
}


# Core evaluation engine. Returns nreal x nobs with "n_failures" and
# "fail_idx" attributes; enforces on_failure + max_fail_frac.
.fm_eval_matrix <- function(model, theta) {
  nreal <- nrow(theta)
  use_per_row <- !is.null(model@map_fn) || model@parallel != "serial"

  if (use_per_row) {
    out <- .fm_eval_per_row(model, theta, nreal)
  } else {
    out <- .fm_eval_bulk(model, theta, nreal)
  }

  # Failure accounting + policy enforcement ---------------------------
  fail_idx <- which(!stats::complete.cases(out))
  n_fail   <- length(fail_idx)

  if (model@on_failure == "stop" && n_fail > 0L) {
    stop(
      sprintf(
        paste0(
          "forward model returned NA for %d of %d realisations ",
          "(on_failure = \"stop\")."
        ),
        n_fail, nreal
      ),
      call. = FALSE
    )
  }
  if (n_fail / nreal > model@max_fail_frac) {
    stop(
      sprintf(
        paste0(
          "forward model failed for %d of %d realisations (%.1f%%), ",
          "exceeding max_fail_frac = %.2f."
        ),
        n_fail, nreal, 100 * n_fail / nreal, model@max_fail_frac
      ),
      call. = FALSE
    )
  }

  attr(out, "n_failures") <- n_fail
  attr(out, "fail_idx")   <- fail_idx
  out
}


# Fast path: one bulk call, serial per-row fallback on bulk error.
.fm_eval_bulk <- function(model, theta, nreal) {
  bulk <- tryCatch(
    .fm_shape_bulk(model@fn(theta), model, nreal),
    error = function(e) {
      if (model@on_failure == "stop") {
        stop(
          sprintf("forward model failed: %s", conditionMessage(e)),
          call. = FALSE
        )
      }
      NULL
    }
  )
  if (!is.null(bulk)) return(bulk)

  # Bulk failed under on_failure = "na": retry row-by-row, serially.
  nobs <- .fm_resolve_nobs(model, theta)
  rows <- lapply(seq_len(nreal), function(i) {
    .fm_eval_one(model@fn, theta[i, , drop = FALSE], nobs)
  })
  .fm_bind_rows(rows, nreal, nobs)
}


# Per-realisation dispatch via the chosen concurrency strategy.
.fm_eval_per_row <- function(model, theta, nreal) {
  nobs   <- .fm_resolve_nobs(model, theta)
  mapper <- .fm_resolve_mapper(model)
  rows <- mapper(seq_len(nreal), function(i) {
    .fm_eval_one(model@fn, theta[i, , drop = FALSE], nobs)
  })
  .fm_bind_rows(rows, nreal, nobs)
}


# Resolve the lapply-shaped mapper for per-row dispatch.
.fm_resolve_mapper <- function(model) {
  if (!is.null(model@map_fn)) return(model@map_fn)
  if (model@parallel == "multicore") {
    n_cores <- model@n_cores
    if (is.na(n_cores)) {
      n_cores <- max(1L, parallel::detectCores() - 1L)
    }
    return(function(X, FUN, ...) {
      parallel::mclapply(X, FUN, ..., mc.cores = n_cores,
                         mc.set.seed = TRUE)
    })
  }
  lapply
}


# Evaluate a single realisation; return a length-nobs vector or NAs.
.fm_eval_one <- function(fn, theta_row, nobs) {
  tryCatch(
    {
      x <- fn(theta_row)
      if (!is.matrix(x)) x <- matrix(as.numeric(x), nrow = 1L)
      v <- as.numeric(x[1L, ])
      if (length(v) != nobs || !all(is.finite(v))) {
        rep(NA_real_, nobs)
      } else {
        v
      }
    },
    error = function(e) rep(NA_real_, nobs)
  )
}


# Shape + validate a bulk result against the declared/inferred nobs.
.fm_shape_bulk <- function(r, model, nreal) {
  if (!is.matrix(r)) r <- as.matrix(r)
  storage.mode(r) <- "double"
  nobs <- if (is.na(model@n_obs)) ncol(r) else model@n_obs
  if (!identical(dim(r), c(nreal, as.integer(nobs)))) {
    stop(
      sprintf(
        "forward model returned shape %dx%d; expected %dx%d.",
        nrow(r), ncol(r), nreal, nobs
      ),
      call. = FALSE
    )
  }
  r
}


# Determine nobs for the per-row path: declared, else probe row 1.
.fm_resolve_nobs <- function(model, theta) {
  if (!is.na(model@n_obs)) return(as.integer(model@n_obs))
  probe <- tryCatch(model@fn(theta[1L, , drop = FALSE]),
                    error = function(e) NULL)
  if (is.null(probe)) {
    stop(
      paste0(
        "Cannot infer `n_obs`: the first realisation failed. Supply ",
        "`n_obs` explicitly to `pesto_forward_model()`."
      ),
      call. = FALSE
    )
  }
  if (!is.matrix(probe)) probe <- matrix(as.numeric(probe), nrow = 1L)
  ncol(probe)
}


# Reconcile a model's declared n_obs against a known observation count.
# Sets it when undeclared (NA); errors on a genuine mismatch.
.fm_with_nobs <- function(model, nobs) {
  if (is.na(model@n_obs)) {
    model@n_obs <- as.integer(nobs)
  } else if (model@n_obs != as.integer(nobs)) {
    stop(
      sprintf(
        "forward model declares n_obs = %d but `obs` has length %d.",
        model@n_obs, as.integer(nobs)
      ),
      call. = FALSE
    )
  }
  model
}


# Bind a list of length-nobs vectors into an nreal x nobs matrix,
# substituting NA rows for malformed entries.
.fm_bind_rows <- function(rows, nreal, nobs) {
  out <- matrix(NA_real_, nrow = nreal, ncol = nobs)
  for (i in seq_len(nreal)) {
    v <- rows[[i]]
    if (!is.null(v) && length(v) == nobs && all(is.finite(v))) {
      out[i, ] <- as.numeric(v)
    }
  }
  out
}
