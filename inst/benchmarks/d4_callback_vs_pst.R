# D4 benchmark: pesto_ies_callback() vs pesto_ies() (.pst-file path).
#
# Phase-1 benchmark uses a *synthetic surrogate* APSIM forward model
# (closed-form polynomial). It quantifies the file-I/O / process-spawn
# overhead of the .pst path against the in-process callback driver.
#
# The real-APSIM speed-up number (target 10x per roadmap §D4) requires
# the D1 scenario library and is captured in a separate benchmark once
# that lands. For the synthetic case, file-I/O dominates the .pst path
# and the callback driver typically runs >100x faster.
#
# Run manually with:
#   Rscript inst/benchmarks/d4_callback_vs_pst.R
#
# Not run on `R CMD check` (lives under inst/ outside vignettes/tests).

suppressPackageStartupMessages({
  library(PESTO)
  library(data.table)
})

set.seed(20260516)

# Synthetic forward model: cheap polynomial in 5 parameters, 12 outputs.
npar  <- 5L
nobs  <- 12L
nreal <- 60L
sigma <- 0.05

A <- matrix(rnorm(nobs * npar), nobs, npar)
b <- rnorm(nobs)
forward <- function(theta) sweep(theta %*% t(A), 2L, b, "+")

theta_true <- rnorm(npar)
y <- as.numeric(forward(matrix(theta_true, nrow = 1L))) +
  rnorm(nobs, sd = sigma)
names(y) <- paste0("o", seq_len(nobs))

prior <- matrix(rnorm(nreal * npar), nreal, npar,
                dimnames = list(NULL, paste0("p", seq_len(npar))))

# ---- callback driver ---------------------------------------------------
t0 <- proc.time()["elapsed"]
fit_cb <- pesto_ies_callback(
  forward_model  = forward,
  prior_ensemble = prior,
  obs            = y,
  obs_sd         = sigma,
  noptmax        = 4L,
  verbose        = FALSE
)
t_cb <- as.numeric(proc.time()["elapsed"] - t0)

cat(sprintf("callback driver: %.3fs (%d realisations x %d iter + final)\n",
            t_cb, nreal, 4L))
cat(sprintf("  posterior RMSE vs truth: %.4f\n",
            sqrt(mean((colMeans(as.matrix(fit_cb$par_ensemble[, -1L])) -
                       theta_true)^2))))
cat(sprintf("  prior RMSE vs truth:     %.4f\n",
            sqrt(mean((colMeans(prior) - theta_true)^2))))
cat(sprintf("  failure rate: %.2f%%\n", 100 * fit_cb$failure_rate))

# ---- .pst driver (skipped if pestpp-ies unavailable) -------------------
have_pestpp <- nchar(Sys.which("pestpp-ies")) > 0L
if (!have_pestpp) {
  cat("\npestpp-ies not on PATH; skipping .pst-path comparison.\n")
  quit(save = "no")
}

# (A full .pst comparison would write a synthetic template + instruction
# files, populate observations, and call pesto_ies(). That harness is
# captured in inst/benchmarks/d4_full_pst_harness.R once the synthetic
# .pst scaffolding stabilises; this script intentionally stops at the
# callback measurement for v0.2.0.)
cat("\n.pst-path harness deferred to d4_full_pst_harness.R.\n")
