# I2 — Matched-prior PESTO vs cached pestpp-ies posterior.
#
# Purpose. Quantify the PESTO-vs-PEST++ posterior disagreement once
# Hypothesis 2 (independent prior random seeds) is removed by construction.
# We re-run PESTO native IES on Scenario A starting from the EXACT prior
# ensemble that pestpp-ies consumed (the cache's `prior_par`), then compare
# posterior means against the cache's `posterior_par`.
#
# Residual disagreement after this run is attributable solely to:
#   H1 — additional Marquardt sub-cycles inside pestpp-ies, or
#   H4 — pestpp-ies's per-realisation observation-noise perturbation.
# H3 (kernel divergence) is closed by the paired-seed reference check
# (see i2_paired_seed_check.R; iteration-1 |Delta| = 5.8e-15).
#
# Status: standalone diagnostic. Run from package root via
#   Rscript inst/scripts/i2_matched_prior_check.R

suppressMessages({
  library(PESTO)
  library(data.table)
})

# Re-create Scenario A's data-generating problem exactly.
set.seed(20260425L)
n_par <- 8L
par_names <- sprintf("p%d", seq_len(n_par))
theta_true <- c(1.20, 0.85, 0.55, 0.40, 0.30, 0.22, 0.16, 0.10)
forward <- function(p) {
  i <- seq_len(15L)
  vapply(i, function(ii) sum(seq_along(p) * p * exp(-ii / 10.0)), numeric(1))
}
y_true <- forward(theta_true)
obs_noise_sd <- 0.02
y_obs <- y_true + rnorm(15L, sd = obs_noise_sd)
weights <- rep(1 / obs_noise_sd, 15L)

# Load pestpp-ies cache.
cache <- readRDS("inst/extdata/pestpp_cache/scenario_a_pestpp_ies.rds")
stopifnot(max(abs(cache$y_obs - y_obs)) < 1e-10)  # confirm same y_obs

# Build the matched prior ensemble in log-space, exactly as pestpp-ies
# consumed it (cache$prior_par columns are p1..p8 in raw space).
prior_dt <- cache$prior_par[, .SD, .SDcols = par_names]
par0_raw <- as.matrix(prior_dt)                      # nreal x npar (raw)
par0     <- t(log(par0_raw))                          # npar x nreal (log)
n_real   <- ncol(par0)
log_lb <- log(0.001); log_ub <- log(5.0)
par0   <- pmin(pmax(par0, log_lb), log_ub)
log_prior_sd <- 0.6
parcov_inv <- rep(1 / log_prior_sd^2, n_par)

run_forward <- function(par_mat) apply(par_mat, 2L, function(lp) forward(exp(lp)))
obs0 <- run_forward(par0)
Y    <- matrix(rep(y_obs, n_real), nrow = length(y_obs))

# Run PESTO IES with the same multi-lambda sweep used inside the vignette
# (Section 1). This emulates a single-cycle Marquardt sweep — not as deep as
# pestpp-ies's full sub-cycling, but a reasonable approximation.
n_iter <- cache$noptmax
lambda_grid <- c(0.5, 5.0, 50.0)
phi_trace <- numeric(n_iter + 1L)
phi_trace[1L] <- mean(compute_phi(Y - obs0, weights))

par_ens <- par0; obs_ens <- obs0
for (iter in seq_len(n_iter)) {
  pm <- rowMeans(par_ens); om <- rowMeans(obs_ens)
  pd <- par_ens - pm
  od <- obs_ens - om
  obs_resid <- obs_ens - Y           # sim - obs (correct convention post-I1)
  Am <- matrix(rnorm(n_par * (n_real - 1L)), n_par, n_real - 1L)
  best_phi <- Inf; best_par <- par_ens; best_obs <- obs_ens
  for (lam in lambda_grid) {
    upg_T <- ensemble_solution(
      par_diff   = pd,
      obs_diff   = od,
      obs_resid  = obs_resid,
      par_resid  = pd,
      weights    = weights,
      parcov_inv = parcov_inv,
      Am         = Am,
      cur_lam    = lam,
      iter       = as.integer(iter)
    )
    cand <- pmin(pmax(par_ens + t(upg_T), log_lb), log_ub)
    cand_obs <- run_forward(cand)
    cand_phi <- mean(compute_phi(Y - cand_obs, weights))
    if (cand_phi < best_phi) {
      best_phi <- cand_phi; best_par <- cand; best_obs <- cand_obs
    }
  }
  par_ens <- best_par; obs_ens <- best_obs
  phi_trace[iter + 1L] <- best_phi
}

posterior_pesto <- exp(t(par_ens))                    # nreal x npar
posterior_pestpp <- as.matrix(cache$posterior_par[, par_names, with = FALSE])

mean_pesto  <- colMeans(posterior_pesto)
mean_pestpp <- colMeans(posterior_pestpp)
rel_diff_pc <- 100 * abs(mean_pesto - mean_pestpp) /
                pmax(abs(mean_pestpp), 1e-6)

cat("\n=== I2 matched-prior PESTO vs pestpp-ies (cached) ===\n")
cat(sprintf("Iterations: %d   Realisations: %d\n", n_iter, n_real))
cat(sprintf("Lambda grid (PESTO sweep): %s\n",
            paste(lambda_grid, collapse = ", ")))
print(data.table(
  parameter      = par_names,
  truth          = round(theta_true,   4),
  pesto_mean     = round(mean_pesto,   4),
  pestpp_mean    = round(mean_pestpp,  4),
  rel_diff_pct   = round(rel_diff_pc,  2)
))
cat(sprintf("\nMedian per-parameter |rel diff| = %.2f %%\n",
            median(rel_diff_pc)))
cat(sprintf("Max    per-parameter |rel diff| = %.2f %%\n",
            max(rel_diff_pc)))

# Persist artefact for the I2 information note.
out_path <- normalizePath(file.path("..", "..", "..",
  "Library", "CloudStorage", "Box-Box", "A_UniAdelaide", "aa_at_work",
  "PEST_plus_plus", "I2_matched_prior_result.rds"), mustWork = FALSE)
tryCatch(saveRDS(list(
  date            = Sys.Date(),
  protocol        = "matched prior, multi-lambda sweep, no Marquardt sub-cycle",
  n_iter          = n_iter,
  n_real          = n_real,
  posterior_pesto_mean  = mean_pesto,
  posterior_pestpp_mean = mean_pestpp,
  rel_diff_pc           = rel_diff_pc,
  median_rel_diff_pc    = median(rel_diff_pc),
  max_rel_diff_pc       = max(rel_diff_pc),
  phi_trace_pesto       = phi_trace,
  phi_trace_pestpp      = cache$phi_history$mean
), out_path),
error = function(e) cat("(could not save:", conditionMessage(e), ")\n"))
cat("\nArtefact: ", out_path, "\n", sep = "")
