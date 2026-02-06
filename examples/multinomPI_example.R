###############################################################################
# Example: Model-Agnostic Predictive Intervals for Claims Reserving
#
# Demonstrates the Dirichlet-Multinomial parametric bootstrap on two
# contrasting datasets:
#   - Taylor & Ashe (1983): stable paid losses, c_hat ~ 108
#   - RAA: volatile reinsurance, c_hat ~ 23
#
# Reference:
# Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive
# Intervals for Claims Reserving via Dirichlet-Multinomial Allocation.
###############################################################################

library(hgr)
library(ChainLadder)

set.seed(2026)

# =============================================================================
# 1. TAYLOR & ASHE: A WELL-BEHAVED TRIANGLE
# =============================================================================

cat("================================================================\n")
cat("1. Taylor & Ashe (1983) — Paid Loss Triangle\n")
cat("================================================================\n\n")

# Convert cumulative to incremental
incr_ta <- cum2incr(GenIns)

# Step 1: Diagnostic — is the Dirichlet model appropriate?
diagnose_c(incr_ta)

# Step 2: Bootstrap predictive intervals
boot_ta <- multinomial_bootstrap(incr_ta, B = 10000)
print(boot_ta)

# Per-origin breakdown
cat("\nPer-origin reserves:\n")
print(boot_ta$by_origin[boot_ta$by_origin$reserve_cl > 0,
                          c("origin", "observed", "reserve_cl",
                            "reserve_mean", "reserve_se")],
      row.names = FALSE, digits = 0)

# Compare with Mack
mack_ta <- MackChainLadder(GenIns)
cat(sprintf("\nMack total SE:  %12.0f  (CV %.1f%%)\n",
            mack_ta$Total.Mack.S.E,
            100 * mack_ta$Total.Mack.S.E / mack_ta$Total.IBNR))

# Delta method (fast, no bootstrap)
cat("\n--- Delta method approximation ---\n")
dm <- delta_method_var(incr_ta)
print(dm[dm$reserve > 0, ], row.names = FALSE, digits = 0)


# =============================================================================
# 2. RAA: A VOLATILE REINSURANCE TRIANGLE
# =============================================================================

cat("\n\n================================================================\n")
cat("2. RAA — Reinsurance Association of America\n")
cat("================================================================\n\n")

incr_raa <- cum2incr(as.triangle(RAA))

# Diagnostic fires: c < 30
diagnose_c(incr_raa)

# Bootstrap
boot_raa <- multinomial_bootstrap(incr_raa, B = 10000)
print(boot_raa)


# =============================================================================
# 3. MODEL-AGNOSTIC: USING CUSTOM DEVELOPMENT PROPORTIONS
# =============================================================================

cat("\n\n================================================================\n")
cat("3. Model-Agnostic: Custom pi on Taylor & Ashe\n")
cat("================================================================\n\n")

# The key feature: ANY development proportions can be used.
# Here we use Chain-Ladder, but these could come from BF, Cape Cod,
# GLMs, expert judgment, or any other source.

dev <- estimate_dev_proportions(incr_ta)
cat("Chain-Ladder pi:\n")
cat(sprintf("  (%s)\n\n", paste(sprintf("%.4f", dev$pi_hat), collapse = ", ")))

# Suppose external information suggests c = 80
boot_ext <- multinomial_bootstrap(incr_ta, pi_hat = dev$pi_hat,
                                   c_param = 80, B = 10000)

cat("With externally specified c = 80:\n")
cat(sprintf("  SE:   %12.0f  (vs %.0f with estimated c = %.0f)\n",
            boot_ext$reserve_se, boot_ta$reserve_se, boot_ta$c_hat))
cat(sprintf("  95%% PI: [%.0f, %.0f]\n", boot_ext$ci_lower, boot_ext$ci_upper))
cat("\nLower c -> wider intervals (more development variability).\n")


# =============================================================================
# 4. BAYESIAN EXTENSION
# =============================================================================

cat("\n\n================================================================\n")
cat("4. Bayesian Predictive Bootstrap on RAA\n")
cat("================================================================\n\n")

bayes_raa <- bayesian_bootstrap(incr_raa, B = 5000,
                                 n_mcmc = 2000, burnin = 500)

cat(sprintf("Bayesian mean:     %12.0f\n", bayes_raa$reserve_mean))
cat(sprintf("Bayesian SE:       %12.0f\n", bayes_raa$reserve_se))
cat(sprintf("95%% PI:       [%12.0f, %12.0f]\n",
            bayes_raa$ci_lower, bayes_raa$ci_upper))
cat(sprintf("Posterior c:  %.1f (SD %.1f)\n",
            bayes_raa$c_posterior_mean, bayes_raa$c_posterior_sd))
cat(sprintf("MH accept:   %.0f%%\n", 100 * bayes_raa$accept_rate))


# =============================================================================
# 5. DIAGNOSTIC BATTERY
# =============================================================================

cat("\n\n================================================================\n")
cat("5. Concentration Diagnostic Across Datasets\n")
cat("================================================================\n\n")

datasets <- list(
  "Taylor-Ashe" = cum2incr(GenIns),
  "UKMotor"     = cum2incr(UKMotor),
  "RAA"         = cum2incr(as.triangle(RAA))
)

cat(sprintf("%-15s %6s %8s %10s\n", "Dataset", "Dims", "c_hat", "Adequate?"))
cat(paste(rep("-", 43), collapse = ""), "\n")

for (nm in names(datasets)) {
  tri <- datasets[[nm]]
  d <- diagnose_c(tri)
  cat(sprintf("%-15s %3dx%-2d %8.1f %10s\n",
              nm, nrow(tri), ncol(tri), d$c_hat,
              ifelse(d$adequate, "Yes", "NO")))
}

cat("\nDatasets with c < 30 need models with accident-year heterogeneity.\n")
