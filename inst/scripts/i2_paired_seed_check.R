# I2 — Paired-seed kernel-vs-reference check.
#
# Scope. Eliminates two of the three working hypotheses for the 15%
# PESTO-vs-PEST++ posterior disagreement reported in
# vignettes/pestpp-comparison-and-simulation.Rmd Section 1:
#
#   H1: PEST++ runs additional Marquardt sub-cycles inside each IES iteration.
#   H2: Independent prior random seeds.
#   H3: PESTO kernel sign or formula divergence (closed by I1).
#
# This script runs the Section A problem with a SINGLE shared prior ensemble
# and a SINGLE fixed lambda (no Marquardt sweep) against an independent pure-R
# reference implementation derived directly from Chen & Oliver (2013) eq. 12.
# H1 and H2 are removed by construction; the residual delta isolates kernel
# fidelity. The textbook equation is:
#
#   d_obs* = obs_resid_signed (with sign convention obs - sim)
#   K = par_diff %*% V %*% diag(s / (s^2 + (lam+1) * I)) %*% U^T
#   delta_theta_j = K %*% (W * d_obs*_j)
#
# The PESTO kernel uses obs_resid = sim - obs internally (after I1), with a
# leading minus on the upgrade. The reference uses obs - sim and no leading
# minus; the two must agree to numerical precision.
#
# Status: standalone diagnostic, NOT a unit test. Run from package root via
#   Rscript inst/scripts/i2_paired_seed_check.R
# Logs to stdout; saves a paired-seed comparison RDS to ../../<repo>/I2_paired_seed_result.rds
# in the project tracker location for traceability.

suppressMessages({
  library(PESTO)
  library(data.table)
})

stopifnot(requireNamespace("PESTO", quietly = TRUE))

# ---------------------------------------------------------------- problem ----
set.seed(20260425L)
n_par <- 8L
par_names <- sprintf("p%d", seq_len(n_par))
theta_true <- c(1.20, 0.85, 0.55, 0.40, 0.30, 0.22, 0.16, 0.10)

forward <- function(p) {
  i <- seq_len(15L)
  vapply(i, function(ii) sum(seq_along(p) * p * exp(-ii / 10.0)), numeric(1))
}
n_obs <- 15L
y_true <- forward(theta_true)
obs_noise_sd <- 0.02
y_obs <- y_true + rnorm(n_obs, sd = obs_noise_sd)
weights <- rep(1 / obs_noise_sd, n_obs)
parcov_inv <- rep(1 / 0.6^2, n_par)

n_real <- 40L
log_lb <- log(0.001); log_ub <- log(5.0)
log_prior_sd <- 0.6

# Shared prior ensemble — paired seed.
set.seed(20260425L)
par0 <- matrix(
  rnorm(n_par * n_real, mean = log(0.5), sd = log_prior_sd),
  nrow = n_par, ncol = n_real,
  dimnames = list(par_names, sprintf("r%02d", seq_len(n_real)))
)
par0 <- pmin(pmax(par0, log_lb), log_ub)
run_forward <- function(par_mat) apply(par_mat, 2L, function(lp) forward(exp(lp)))
obs0 <- run_forward(par0)
Y    <- matrix(rep(y_obs, n_real), nrow = n_obs)

# ------------------------------------------------------- reference kernel ----
# Pure-R IES upgrade derived from Chen & Oliver (2013) eq. 12.
# Convention here: obs_resid_ref = obs - sim (textbook), no leading minus.
ref_upgrade <- function(par_mat, obs_mat, y, w, lam) {
  pm <- rowMeans(par_mat); om <- rowMeans(obs_mat)
  pd <- par_mat - pm                         # npar x nreal
  od <- obs_mat - om                         # nobs x nreal
  scale <- 1 / sqrt(ncol(par_mat) - 1)
  W <- diag(w)
  Yres_signed <- matrix(rep(y, ncol(par_mat)), nrow = length(y)) - obs_mat # obs - sim
  S <- svd(scale * W %*% od)
  Sigma <- S$d
  inv <- 1 / (Sigma^2 + (lam + 1))
  X <- S$v %*% diag(Sigma * inv) %*% t(S$u) %*% (W %*% Yres_signed)
  # Upgrade in parameter space — note positive sign because Yres = obs - sim.
  upgrade <- (scale * pd) %*% X
  upgrade  # npar x nreal
}

# ------------------------------------------------------------ pesto kernel --
# Wrap the C++ kernel with the *correct* (post-I1) sign convention sim - obs.
pesto_upgrade <- function(par_mat, obs_mat, y, w, lam) {
  pm <- rowMeans(par_mat); om <- rowMeans(obs_mat)
  pd <- par_mat - pm
  od <- obs_mat - om
  obs_resid <- obs_mat - matrix(rep(y, ncol(par_mat)), nrow = length(y))  # sim - obs
  Am <- matrix(0, n_par, 0)  # use_approx = TRUE means Am ignored
  upg_T <- ensemble_solution(
    par_diff   = pd,
    obs_diff   = od,
    obs_resid  = obs_resid,
    par_resid  = pd,
    weights    = w,
    parcov_inv = parcov_inv,
    Am         = Am,
    cur_lam    = lam,
    iter       = 1L,
    use_approx = TRUE
  )
  t(upg_T)  # back to npar x nreal
}

# ------------------------------------- iteration-1 elementwise comparison ---
lam0 <- 1.0
upg_pesto <- pesto_upgrade(par0, obs0, y_obs, weights, lam0)
upg_ref   <- ref_upgrade  (par0, obs0, y_obs, weights, lam0)

elem_max <- max(abs(upg_pesto - upg_ref))
elem_mean <- mean(abs(upg_pesto - upg_ref))
elem_rel <- elem_mean / mean(abs(upg_ref))

cat("\n=== I2 paired-seed iteration-1 kernel agreement ===\n")
cat(sprintf("max |Delta|         = %.3e\n", elem_max))
cat(sprintf("mean |Delta|        = %.3e\n", elem_mean))
cat(sprintf("mean |Delta| / |ref|= %.3e\n", elem_rel))

# ----------------------------------- multi-iteration pure-step comparison ---
n_iter <- 6L
phi_pesto <- numeric(n_iter + 1L)
phi_ref   <- numeric(n_iter + 1L)
phi_pesto[1L] <- mean(compute_phi(Y - obs0, weights))
phi_ref  [1L] <- phi_pesto[1L]

p_p <- par0; o_p <- obs0
p_r <- par0; o_r <- obs0
for (it in seq_len(n_iter)) {
  # PESTO
  upg <- pesto_upgrade(p_p, o_p, y_obs, weights, lam0)
  p_p <- pmin(pmax(p_p + upg, log_lb), log_ub)
  o_p <- run_forward(p_p)
  phi_pesto[it + 1L] <- mean(compute_phi(Y - o_p, weights))
  # Reference
  upg <- ref_upgrade(p_r, o_r, y_obs, weights, lam0)
  p_r <- pmin(pmax(p_r + upg, log_lb), log_ub)
  o_r <- run_forward(p_r)
  phi_ref  [it + 1L] <- mean(compute_phi(Y - o_r, weights))
}

post_pesto <- exp(t(p_p))
post_ref   <- exp(t(p_r))
mean_pesto <- colMeans(post_pesto)
mean_ref   <- colMeans(post_ref)
rel_diff   <- 100 * abs(mean_pesto - mean_ref) / pmax(abs(mean_ref), 1e-6)

cat("\n=== I2 paired-seed posterior comparison (single lambda, no Marquardt) ===\n")
print(data.table(
  parameter   = par_names,
  truth       = round(theta_true, 4),
  pesto_mean  = round(mean_pesto, 4),
  ref_mean    = round(mean_ref,   4),
  rel_diff_pc = round(rel_diff,   3)
))

cat("\nMedian per-parameter |rel diff|: ",
    round(median(rel_diff), 4), "%\n", sep = "")
cat("Max    per-parameter |rel diff|: ",
    round(max(rel_diff),    4), "%\n", sep = "")

cat("\nPhi traces:\n")
cat("  PESTO    : ", paste(signif(phi_pesto, 4), collapse = " -> "), "\n")
cat("  reference: ", paste(signif(phi_ref,   4), collapse = " -> "), "\n")

# Save artefact for the project log.
out_path <- normalizePath(file.path("..", "..", "..",
  "Library", "CloudStorage", "Box-Box", "A_UniAdelaide", "aa_at_work",
  "PEST_plus_plus", "I2_paired_seed_result.rds"), mustWork = FALSE)
result <- list(
  date         = Sys.Date(),
  problem      = "Scenario A 8-par exponential decay",
  protocol     = "shared prior, single fixed lambda, no Marquardt sub-cycle",
  iter1_max_abs_delta  = elem_max,
  iter1_mean_abs_delta = elem_mean,
  iter1_rel_delta      = elem_rel,
  posterior_pesto_mean = mean_pesto,
  posterior_ref_mean   = mean_ref,
  posterior_rel_diff_pc = rel_diff,
  median_rel_diff_pc   = median(rel_diff),
  phi_pesto = phi_pesto,
  phi_ref   = phi_ref
)
tryCatch(
  saveRDS(result, out_path),
  error = function(e) cat("(could not save RDS:", conditionMessage(e), ")\n")
)
cat("\nArtefact: ", out_path, "\n", sep = "")
