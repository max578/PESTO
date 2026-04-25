#' Plot Objective Function (Phi) Convergence
#'
#' Creates a publication-quality plot of objective function values
#' across iterations, showing mean, min, max, and individual
#' realisation traces.
#'
#' @param result A `pesto_ies_result` or `pesto_glm_result` object,
#'   or a data.table with columns `iteration` and `phi` (at minimum).
#' @param log_scale Logical. Use log10 scale for y-axis (default TRUE).
#' @param show_reals Logical. Show individual realisation traces.
#' @param title Character. Plot title.
#' @return A ggplot2 object.
#' @examples
#' phi_dt <- data.table::data.table(
#'   iteration = 0:4,
#'   total_runs = c(50L, 100L, 150L, 200L, 250L),
#'   mean = c(1200, 450, 180, 95, 72),
#'   min  = c(900, 320, 130, 70, 55),
#'   max  = c(1700, 680, 260, 140, 105),
#'   median = c(1180, 440, 175, 92, 70),
#'   std    = c(180, 80, 35, 18, 12)
#' )
#' p <- plot_phi(phi_dt, log_scale = TRUE,
#'               title = "Synthetic Phi Convergence")
#' inherits(p, "ggplot")
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon
#'   scale_y_log10 labs theme_minimal theme element_text
plot_phi <- function(result,
                     log_scale = TRUE,
                     show_reals = FALSE,
                     title = "Objective Function Convergence") {

  if (inherits(result, "pesto_ies_result")) {
    phi_dt <- result$phi
    if (is.null(phi_dt)) stop("No phi data in result", call. = FALSE)
  } else if (data.table::is.data.table(result) || is.data.frame(result)) {
    phi_dt <- data.table::copy(data.table::as.data.table(result))
  } else {
    stop("result must be a pesto_ies_result or data.table", call. = FALSE)
  }

  # Try to identify iteration and phi columns
  iter_col <- intersect(names(phi_dt), c("iteration", "iter", "i"))
  phi_cols <- setdiff(names(phi_dt), c(iter_col, "total_runs", "lambda",
                                        "min", "max", "mean", "median", "std"))

  if (length(iter_col) == 0) {
    phi_dt[, iteration := .I]
    iter_col <- "iteration"
  } else {
    iter_col <- iter_col[1]
  }

  # Summary statistics
  if ("mean" %in% names(phi_dt) && "min" %in% names(phi_dt)) {
    # Already summarised
    p <- ggplot2::ggplot(phi_dt, ggplot2::aes(x = .data[[iter_col]])) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = min, ymax = max),
        alpha = 0.2, fill = "steelblue"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = mean),
        colour = "steelblue", linewidth = 1.2
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = mean),
        colour = "steelblue", size = 2.5
      )
  } else if (length(phi_cols) > 0) {
    # Individual realisations - melt
    melt_dt <- data.table::melt(phi_dt, id.vars = iter_col,
                                 measure.vars = phi_cols,
                                 variable.name = "realisation",
                                 value.name = "phi")

    summary_dt <- melt_dt[, .(
      mean = mean(phi, na.rm = TRUE),
      min = min(phi, na.rm = TRUE),
      max = max(phi, na.rm = TRUE),
      median = stats::median(phi, na.rm = TRUE)
    ), by = iter_col]

    p <- ggplot2::ggplot()

    if (show_reals) {
      p <- p + ggplot2::geom_line(
        data = melt_dt,
        ggplot2::aes(x = .data[[iter_col]], y = phi,
                     group = realisation),
        alpha = 0.15, colour = "grey50"
      )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = summary_dt,
        ggplot2::aes(x = .data[[iter_col]], ymin = min, ymax = max),
        alpha = 0.2, fill = "steelblue"
      ) +
      ggplot2::geom_line(
        data = summary_dt,
        ggplot2::aes(x = .data[[iter_col]], y = mean),
        colour = "steelblue", linewidth = 1.2
      ) +
      ggplot2::geom_point(
        data = summary_dt,
        ggplot2::aes(x = .data[[iter_col]], y = mean),
        colour = "steelblue", size = 2.5
      )
  } else {
    stop("Cannot identify phi columns in data", call. = FALSE)
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = "Iteration",
      y = expression(Phi ~ "(Objective Function)")
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p
}


#' Plot Ensemble Parameter Distributions
#'
#' Visualises the prior and/or posterior parameter ensemble
#' distributions as violin plots or density ridges.
#'
#' @param ensemble data.table. Parameter ensemble (rows = realisations).
#' @param parameters Character vector. Parameter names to plot.
#'   If NULL, selects up to 20 parameters with highest variance.
#' @param prior_ensemble data.table. Optional prior ensemble for comparison.
#' @param max_params Integer. Maximum parameters to display.
#' @param title Character. Plot title.
#' @return A ggplot2 object.
#' @examples
#' posterior <- data.table::data.table(
#'   real_name = sprintf("real_%02d", 1:50),
#'   k1 = rnorm(50, 1.0, 0.15),
#'   k2 = rnorm(50, 0.5, 0.05),
#'   k3 = rnorm(50, 2.0, 0.40),
#'   k4 = rnorm(50, 0.1, 0.02)
#' )
#' prior <- data.table::data.table(
#'   real_name = sprintf("real_%02d", 1:50),
#'   k1 = rnorm(50, 1.0, 0.50),
#'   k2 = rnorm(50, 0.5, 0.20),
#'   k3 = rnorm(50, 2.0, 1.00),
#'   k4 = rnorm(50, 0.1, 0.08)
#' )
#' p <- plot_ensemble(posterior, prior_ensemble = prior)
#' inherits(p, "ggplot")
#' @export
plot_ensemble <- function(ensemble,
                          parameters = NULL,
                          prior_ensemble = NULL,
                          max_params = 20L,
                          title = "Parameter Ensemble Distributions") {

  if (!data.table::is.data.table(ensemble)) {
    ensemble <- data.table::as.data.table(ensemble)
  }

  # Remove non-numeric columns
  num_cols <- names(ensemble)[vapply(ensemble, is.numeric, logical(1))]

  if (is.null(parameters)) {
    # Select top-variance parameters
    vars <- vapply(num_cols, function(cn) stats::var(ensemble[[cn]], na.rm = TRUE), numeric(1))
    parameters <- names(sort(vars, decreasing = TRUE))[seq_len(min(max_params, length(vars)))]
  }

  parameters <- intersect(parameters, num_cols)
  if (length(parameters) == 0) {
    stop("No valid numeric parameter columns found", call. = FALSE)
  }

  # Melt for plotting
  plot_dt <- data.table::melt(
    ensemble[, .SD, .SDcols = parameters],
    measure.vars = parameters,
    variable.name = "parameter",
    value.name = "value"
  )
  plot_dt[, source := "Posterior"]

  if (!is.null(prior_ensemble)) {
    if (!data.table::is.data.table(prior_ensemble)) {
      prior_ensemble <- data.table::as.data.table(prior_ensemble)
    }
    prior_dt <- data.table::melt(
      prior_ensemble[, .SD, .SDcols = parameters],
      measure.vars = parameters,
      variable.name = "parameter",
      value.name = "value"
    )
    prior_dt[, source := "Prior"]
    plot_dt <- data.table::rbindlist(list(prior_dt, plot_dt))
  }

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = parameter, y = value, fill = source)) +
    ggplot2::geom_violin(alpha = 0.6, position = "identity", scale = "width") +
    ggplot2::labs(
      title = title,
      x = "Parameter",
      y = "Value",
      fill = "Ensemble"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = c("Prior" = "#E69F00", "Posterior" = "#0072B2"))

  if (is.null(prior_ensemble)) {
    p <- p + ggplot2::guides(fill = "none")
  }

  p
}


#' Plot Parameter Identifiability
#'
#' Creates a bar plot of parameter identifiability based on
#' the singular value decomposition of the Jacobian matrix. Accepts
#' either a numeric Jacobian matrix in memory or a path to a `.jco`
#' (PEST binary) file.
#'
#' @param jacobian Numeric matrix (n_obs x n_par). The Jacobian / sensitivity
#'   matrix. Column names, if present, are used as parameter labels;
#'   otherwise `p1`, `p2`, ... are generated. Either `jacobian` or
#'   `jco_file` must be supplied.
#' @param jco_file Character. Path to a `.jco` (Jacobian) binary file.
#'   Mutually exclusive with `jacobian`.
#' @param pst A `pesto_pst` object for parameter names. Optional; only
#'   used when reading from a `.jco` file and column names are absent.
#' @param n_sv Integer. Number of singular values to retain.
#' @param title Character. Plot title.
#' @return A ggplot2 object.
#' @examples
#' J <- matrix(rnorm(30 * 8), nrow = 30, ncol = 8)
#' J[, 7] <- 0.5 * J[, 1] + 0.5 * J[, 2]
#' J[, 8] <- 1e-6 * rnorm(30)
#' colnames(J) <- paste0("k", 1:8)
#' p <- plot_identifiability(jacobian = J)
#' inherits(p, "ggplot")
#' @export
plot_identifiability <- function(jacobian = NULL,
                                 jco_file = NULL,
                                 pst = NULL,
                                 n_sv = NULL,
                                 title = "Parameter Identifiability") {
  if (is.null(jacobian) && is.null(jco_file)) {
    stop("Either `jacobian` (a numeric matrix) or `jco_file` must be provided.",
         call. = FALSE)
  }
  if (!is.null(jacobian) && !is.null(jco_file)) {
    stop("Supply only one of `jacobian` or `jco_file`, not both.",
         call. = FALSE)
  }

  if (!is.null(jacobian)) {
    if (!is.matrix(jacobian) || !is.numeric(jacobian)) {
      stop("`jacobian` must be a numeric matrix.", call. = FALSE)
    }
    mat <- jacobian
    par_names <- colnames(mat)
    if (is.null(par_names)) {
      par_names <- paste0("p", seq_len(ncol(mat)))
    }
  } else {
    jco <- .read_ensemble_binary(jco_file)
    par_names <- names(jco)[-1]
    mat <- as.matrix(jco[, -1, with = FALSE])
  }

  sv <- svd(mat)

  if (is.null(n_sv)) {
    n_sv <- min(nrow(mat), ncol(mat))
  }

  # Compute identifiability: sum of squared right singular vectors
  V <- sv$v[, seq_len(n_sv), drop = FALSE]
  ident <- rowSums(V^2)

  ident_dt <- data.table::data.table(
    parameter = par_names,
    identifiability = ident
  )
  data.table::setorder(ident_dt, -identifiability)
  ident_dt[, parameter := factor(parameter, levels = parameter)]

  p <- ggplot2::ggplot(ident_dt, ggplot2::aes(x = parameter, y = identifiability)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    ggplot2::labs(
      title = title,
      x = "Parameter",
      y = "Identifiability (0-1)"
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )

  p
}
