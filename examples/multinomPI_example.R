###############################################################################
# Example: Model-Agnostic Predictive Intervals for Claims Reserving
#
# Demonstrates the multinomial parametric bootstrap on two contrasting
# datasets from the ChainLadder package:
#   - Taylor & Ashe (1983): stable paid losses, c_hat = 108
#   - RAA: volatile reinsurance, c_hat = 23
#
# Reference:
# Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive
# Intervals for Claims Reserving via Dirichlet-Multinomial Allocation.
###############################################################################

library(ChainLadder)
library(hgr)

set.seed(2026)

# =============================================================================
# 1. TAYLOR & ASHE: A WELL-BEHAVED TRIANGLE
# =============================================================================

cat("================================================================\n")
cat("1. Taylor & Ashe (1983) — Paid Loss Triangle\n")
cat("================================================================\n\n")

# Convert cumulative to incremental
cum_ta <- GenIns
incr_ta <- cum2incr(cum_ta)

# Diagnostic: is the Dirichlet model appropriate?
diag_ta <- diagnose_concentration(incr_ta)
cat(diag_ta$message, "\n\n")

# Bootstrap predictive intervals (using Chain-Ladder proportions)
boot_ta <- multinomial_bootstrap(incr_ta, B = 10000)

cat(sprintf("Chain-Ladder reserve:   %12.0f\n", boot_ta$reserve_cl))
cat(sprintf("Bootstrap mean:         %12.0f\n", boot_ta$reserve_mean))
cat(sprintf("Bootstrap SE:           %12.0f\n", boot_ta$reserve_se))
cat(sprintf("CV:                     %11.1f%%\n", 100 * boot_ta$cv))
cat(sprintf("95%% PI:           [%12.0f, %12.0f]\n",
            boot_ta$ci_lower, boot_ta$ci_upper))
cat(sprintf("Concentration c_hat:    %12.1f\n\n", boot_ta$c_hat))

# Per-origin breakdown
cat("Per-origin reserve estimates:\n")
print(boot_ta$by_origin[boot_ta$by_origin$reserve_cl > 0,
                          c("origin", "observed", "reserve_cl",
                            "reserve_mean", "reserve_se")],
      row.names = FALSE, digits = 0)

# Compare with Mack
mack_ta <- MackChainLadder(cum_ta)
cat(sprintf("\nMack total SE:          %12.0f\n", mack_ta$Total.Mack.S.E))
cat(sprintf("Mack CV:                %11.1f%%\n",
            100 * mack_ta$Total.Mack.S.E / mack_ta$Total.IBNR))

# Delta method approximation (fast, no bootstrap needed)
cat("\n--- Delta method approximation ---\n")
delta_ta <- delta_method_variance(incr_ta)
cat("Per-origin process SE:\n")
print(delta_ta[delta_ta$reserve > 0, c("origin", "reserve", "process_se", "cv")],
      row.names = FALSE, digits = 3)


# =============================================================================
# 2. RAA: A VOLATILE REINSURANCE TRIANGLE
# =============================================================================

cat("\n\n================================================================\n")
cat("2. RAA — Reinsurance Association of America\n")
cat("================================================================\n\n")

cum_raa <- as.triangle(ChainLadder::RAA)
incr_raa <- cum2incr(cum_raa)

# Diagnostic: c < 30 signals trouble
diag_raa <- diagnose_concentration(incr_raa)
cat(diag_raa$message, "\n\n")

# Bootstrap anyway to see the intervals
boot_raa <- multinomial_bootstrap(incr_raa, B = 10000)

cat(sprintf("Chain-Ladder reserve:   %12.0f\n", boot_raa$reserve_cl))
cat(sprintf("Bootstrap mean:         %12.0f\n", boot_raa$reserve_mean))
cat(sprintf("Bootstrap SE:           %12.0f\n", boot_raa$reserve_se))
cat(sprintf("CV:                     %11.1f%%\n", 100 * boot_raa$cv))
cat(sprintf("95%% PI:           [%12.0f, %12.0f]\n",
            boot_raa$ci_lower, boot_raa$ci_upper))
cat(sprintf("Concentration c_hat:    %12.1f\n", boot_raa$c_hat))


# =============================================================================
# 3. MODEL-AGNOSTIC: USING BORNHUETTER-FERGUSON PROPORTIONS
# =============================================================================

cat("\n\n================================================================\n")
cat("3. Model-Agnostic: BF Proportions on Taylor & Ashe\n")
cat("================================================================\n\n")

# Suppose the actuary uses Bornhuetter-Ferguson with external priors.
# The BF method still produces development proportions — we can wrap
# those in the Dirichlet framework for predictive intervals.

# Step 1: Get BF development proportions (same as CL here, but the
# point is that ANY pi vector can be plugged in)
bf_result <- estimate_development_proportions(incr_ta)
pi_bf <- bf_result$pi_hat

cat("Development proportions (same for CL and BF):\n")
cat(sprintf("  pi = (%s)\n\n",
            paste(sprintf("%.3f", pi_bf), collapse = ", ")))

# Step 2: Suppose we have external information suggesting c = 80
# (e.g., from a larger portfolio or historical data)
boot_bf <- multinomial_bootstrap(incr_ta, pi_hat = pi_bf,
                                  c_param = 80, B = 10000)

cat("With externally specified c = 80:\n")
cat(sprintf("  Bootstrap mean:       %12.0f\n", boot_bf$reserve_mean))
cat(sprintf("  Bootstrap SE:         %12.0f\n", boot_bf$reserve_se))
cat(sprintf("  95%% PI:         [%12.0f, %12.0f]\n",
            boot_bf$ci_lower, boot_bf$ci_upper))

# Step 3: Compare with data-driven c
cat(sprintf("\nWith estimated c = %.1f:\n", boot_ta$c_hat))
cat(sprintf("  Bootstrap SE:         %12.0f\n", boot_ta$reserve_se))
cat(sprintf("  95%% PI:         [%12.0f, %12.0f]\n",
            boot_ta$ci_lower, boot_ta$ci_upper))

cat("\nLower c -> wider intervals (more development variability assumed).\n")


# =============================================================================
# 4. BAYESIAN EXTENSION (for small or uncertain triangles)
# =============================================================================

cat("\n\n================================================================\n")
cat("4. Bayesian Predictive Bootstrap on RAA\n")
cat("================================================================\n\n")

# For the volatile RAA triangle, the Bayesian extension provides
# a posterior for c rather than a point estimate
bayes_raa <- bayesian_predictive_bootstrap(
  incr_raa, B = 5000,
  mu_c = log(50), sigma_c = 1,
  n_mcmc = 2000, burnin = 500
)

cat(sprintf("Bayesian reserve mean:  %12.0f\n", bayes_raa$reserve_mean))
cat(sprintf("Bayesian reserve SE:    %12.0f\n", bayes_raa$reserve_se))
cat(sprintf("95%% PI:           [%12.0f, %12.0f]\n",
            bayes_raa$ci_lower, bayes_raa$ci_upper))
cat(sprintf("Posterior c mean:       %12.1f\n", bayes_raa$c_posterior_mean))
cat(sprintf("Posterior c SD:         %12.1f\n", bayes_raa$c_posterior_sd))
cat(sprintf("MH acceptance rate:     %11.1f%%\n", 100 * bayes_raa$accept_rate))

cat("\nA wide posterior on c confirms the data are uninformative about\n")
cat("development stability — the diagnostic threshold (c < 30) is\n")
cat("the more actionable signal.\n")


# =============================================================================
# 5. BATTERY OF DATASETS: CONCENTRATION DIAGNOSTIC
# =============================================================================

cat("\n\n================================================================\n")
cat("5. Concentration Diagnostic Across Datasets\n")
cat("================================================================\n\n")

datasets <- list(
  "Taylor-Ashe"      = cum2incr(GenIns),
  "UKMotor"          = cum2incr(UKMotor),
  "RAA"              = cum2incr(as.triangle(RAA))
)

# Add MW2008 if available
if (exists("MW2008")) {
  datasets[["MW2008"]] <- cum2incr(MW2008)
}

cat(sprintf("%-20s %6s %8s %10s\n", "Dataset", "Dims", "c_hat", "Adequate?"))
cat(paste(rep("-", 48), collapse = ""), "\n")

for (nm in names(datasets)) {
  tri <- datasets[[nm]]
  diag <- diagnose_concentration(tri)
  dims <- sprintf("%dx%d", nrow(tri), ncol(tri))
  cat(sprintf("%-20s %6s %8.1f %10s\n",
              nm, dims, diag$c_hat, ifelse(diag$adequate, "Yes", "NO")))
}

cat("\nDatasets with c_hat < 30 require models with explicit\n")
cat("accident-year heterogeneity (frailty, frequency-severity\n")
cat("dependence, or micro-level approaches).\n")
