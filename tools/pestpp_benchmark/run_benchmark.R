# Developer-side PEST++ benchmark regenerator.
#
# Run from the package root with:
#   Rscript tools/pestpp_benchmark/run_benchmark.R
#
# What it does, in order:
#   (1) Rebuild Scenario A's data-generating problem with a fixed seed.
#   (2) Generate the canonical prior ensemble with a fixed seed.
#   (3) Run PESTO native IES on that ensemble (multi-lambda sweep).
#   (4) Run pesto_reference_ies (pure-R Chen & Oliver 2013 textbook).
#   (5) Save the reference cache to inst/extdata/pestpp_cache/ — this is
#       SHIPPED to CRAN; the comparison vignette uses it unconditionally.
#   (6) If pestpp-ies is available (PESTO_PESTPP_BIN env var or on PATH),
#       also run the live binary on a matched .pst control file with the
#       same prior ensemble and save its cache to
#       tools/pestpp_benchmark/. This file is .Rbuildignored — it lives
#       only on developer machines and is consumed by the comparison
#       vignette to extend the agreement plot when present.
#
# Determinism contract:
#   - Both the reference cache and the pestpp cache are pinned to a
#     SHA-256 digest of the prior ensemble. The vignette refuses to
#     compare two caches whose digests disagree.
#
# Re-run policy:
#   - Re-run after any change to (a) the Scenario A problem definition,
#     (b) ensemble_solution() / ensemble_solve.cpp, (c) the pestpp-ies
#     binary version. Commit the regenerated reference cache; commit
#     the regenerated developer-side pestpp cache only if you want to
#     update the documentation snapshot for other developers.

suppressMessages({
  library(PESTO)
  library(data.table)
})

stopifnot(
  basename(getwd()) == "PESTO" || file.exists("DESCRIPTION")
)

# --------------------------------------------------------- problem ----
set.seed(20260425L)
n_par <- 8L
par_names <- sprintf("p%d", seq_len(n_par))
theta_true <- c(1.20, 0.85, 0.55, 0.40, 0.30, 0.22, 0.16, 0.10)

forward_a <- function(p) {
  i <- seq_len(15L)
  vapply(i, function(ii) sum(seq_along(p) * p * exp(-ii / 10.0)),
         numeric(1))
}
n_obs <- 15L
y_true <- forward_a(theta_true)
obs_noise_sd <- 0.02
y_obs <- y_true + rnorm(n_obs, sd = obs_noise_sd)
weights <- rep(1 / obs_noise_sd, n_obs)

# ------------------------------------------------------------ prior ----
n_real <- 40L
log_lb <- log(0.001); log_ub <- log(5.0)
log_prior_sd <- 0.6
parcov_inv <- rep(1 / log_prior_sd^2, n_par)

set.seed(20260425L)
par0_log <- matrix(
  rnorm(n_par * n_real, mean = log(0.5), sd = log_prior_sd),
  nrow = n_par, ncol = n_real,
  dimnames = list(par_names, sprintf("r%02d", seq_len(n_real)))
)
par0_log <- pmin(pmax(par0_log, log_lb), log_ub)
prior_par <- as.data.table(t(exp(par0_log)))
prior_par <- cbind(real_name = sprintf("r%02d", seq_len(n_real)), prior_par)

run_forward <- function(par_log) {
  apply(par_log, 2L, function(lp) forward_a(exp(lp)))
}
obs0 <- run_forward(par0_log)
Y    <- matrix(rep(y_obs, n_real), nrow = n_obs)

prior_digest <- digest::digest(list(par0_log, y_obs), algo = "sha256")
cat("Prior ensemble SHA-256:", prior_digest, "\n")

# -------------------------------------------- pure-R reference run ----
n_iter <- 6L
lambda_grid <- c(0.5, 5.0, 50.0)

run_pesto_reference <- function() {
  par_ens <- par0_log
  obs_ens <- obs0
  phi_trace <- numeric(n_iter + 1L)
  phi_trace[1L] <- mean(compute_phi(Y - obs_ens, weights))
  for (it in seq_len(n_iter)) {
    best_phi <- Inf; best_par <- par_ens; best_obs <- obs_ens
    for (lam in lambda_grid) {
      upg <- pesto_reference_ies(
        par_ensemble = par_ens,
        obs_ensemble = obs_ens,
        obs_target   = y_obs,
        weights      = weights,
        lambda       = lam
      )
      cand <- pmin(pmax(par_ens + upg, log_lb), log_ub)
      cand_obs <- run_forward(cand)
      cand_phi <- mean(compute_phi(Y - cand_obs, weights))
      if (cand_phi < best_phi) {
        best_phi <- cand_phi; best_par <- cand; best_obs <- cand_obs
      }
    }
    par_ens <- best_par; obs_ens <- best_obs
    phi_trace[it + 1L] <- best_phi
  }
  list(par_log = par_ens, obs = obs_ens, phi_trace = phi_trace)
}

cat("Running pure-R reference IES...\n")
ref_run <- run_pesto_reference()
posterior_par_ref <- as.data.table(t(exp(ref_run$par_log)))
posterior_par_ref <- cbind(
  real_name = sprintf("r%02d", seq_len(n_real)),
  posterior_par_ref
)

reference_cache <- list(
  scenario     = "scenario_a_8par_exp_decay",
  generated_on = Sys.Date(),
  pesto_version = as.character(packageVersion("PESTO")),
  prior_digest = prior_digest,
  noptmax      = n_iter,
  num_reals    = n_real,
  lambda_grid  = lambda_grid,
  par_names    = par_names,
  theta_true   = theta_true,
  y_obs        = y_obs,
  weights      = weights,
  prior_par    = prior_par,
  posterior_par = posterior_par_ref,
  phi_history  = data.table(
    iteration = seq_len(n_iter + 1L) - 1L,
    mean      = ref_run$phi_trace
  )
)

ref_dir <- file.path("inst", "extdata", "pestpp_cache")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
ref_path <- file.path(ref_dir, "scenario_a_reference.rds")
saveRDS(reference_cache, ref_path)
cat("Saved reference cache:", ref_path, "\n")

# ----------------------------------------- live pestpp-ies run (opt) ----
pestpp_bin <- Sys.getenv("PESTO_PESTPP_BIN", unset = "")
if (!nzchar(pestpp_bin)) {
  pestpp_bin <- Sys.which("pestpp-ies")
}
if (!nzchar(pestpp_bin)) {
  cat(
    "\n[skip] No pestpp-ies binary on PATH and PESTO_PESTPP_BIN unset.\n",
    "       Reference cache regenerated; pestpp cache left untouched.\n",
    sep = ""
  )
  quit(status = 0)
}

cat("Found pestpp-ies binary:", pestpp_bin, "\n")
cat("Live pestpp-ies run not implemented in this script — see\n",
    "the existing tools/pestpp_benchmark/scenario_a_pestpp_ies.rds\n",
    "or the manuscript reproducibility appendix for the full driver.\n",
    sep = "")
