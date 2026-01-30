# =============================================================================
# Example usage of the vog package
# =============================================================================

library(vog)

# -----------------------------------------------------------------------------
# Example triangle (Taylor-Ashe, incremental)
# -----------------------------------------------------------------------------

taylor_ashe <- matrix(c(
  357848, 766940, 610542, 482940, 527326, 574398, 146342, 139950, 227229, 67948,
  352118, 884021, 933894, 1183289, 445745, 320996, 527804, 266172, 425046, NA,
  290507, 1001799, 926219, 1016654, 750816, 146923, 495992, 280405, NA, NA,
  310608, 1108250, 776189, 1562400, 272482, 352053, 206286, NA, NA, NA,
  443160, 693190, 991983, 769488, 504851, 470639, NA, NA, NA, NA,
  396132, 937085, 847498, 805037, 705960, NA, NA, NA, NA, NA,
  440832, 847631, 1131398, 1063269, NA, NA, NA, NA, NA, NA,
  359480, 1061648, 1443370, NA, NA, NA, NA, NA, NA, NA,
  376686, 986608, NA, NA, NA, NA, NA, NA, NA, NA,
  344014, NA, NA, NA, NA, NA, NA, NA, NA, NA
), nrow = 10, ncol = 10, byrow = TRUE)

# -----------------------------------------------------------------------------
# Fit NB-CL model
# -----------------------------------------------------------------------------

fit <- fit_nbcl(taylor_ashe)
print(fit)

# Summary with coefficients
summary(fit)

# -----------------------------------------------------------------------------
# Reserve estimates
# -----------------------------------------------------------------------------

reserves <- reserve_nbcl(fit)
cat("\nTotal reserve:", format(round(reserves$total), big.mark = ","), "\n")
print(reserves$by_accident_year)

# -----------------------------------------------------------------------------
# Bootstrap prediction intervals
# -----------------------------------------------------------------------------

# With bias correction (recommended)
boot_corr <- bootstrap_nbcl(fit, B = 2000, correct_kappa = TRUE, seed = 42)
print(boot_corr)

# Extract 95% prediction interval
pi_95 <- predict_interval(boot_corr, level = 0.95)
print(pi_95)

# Without bias correction (for comparison)
boot_mle <- bootstrap_nbcl(fit, B = 2000, correct_kappa = FALSE, seed = 42)
print(boot_mle)

# -----------------------------------------------------------------------------
# Model comparison
# -----------------------------------------------------------------------------

comparison <- compare_models(taylor_ashe)
print(comparison)

# Test for overdispersion
lr_test_overdispersion(fit)

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

# All diagnostic plots
plot_diagnostics(fit, which = 1:4)

# Profile likelihood for kappa
profile <- profile_kappa(fit)
plot(profile)
cat("\nKappa 95% CI: [", round(profile$ci_lower, 1), ", ", 
    round(profile$ci_upper, 1), "]\n", sep = "")

# -----------------------------------------------------------------------------
# Compare with ChainLadder package (if available)
# -----------------------------------------------------------------------------

if (requireNamespace("ChainLadder", quietly = TRUE)) {
  library(ChainLadder)
  
  # Convert to cumulative for ChainLadder
  taylor_cum <- taylor_ashe
  for (i in 1:10) {
    taylor_cum[i, ] <- cumsum(replace(taylor_ashe[i, ], is.na(taylor_ashe[i, ]), 0))
    taylor_cum[i, is.na(taylor_ashe[i, ])] <- NA
  }
  class(taylor_cum) <- c("triangle", "matrix")
  
  # ODP bootstrap from ChainLadder
  odp_boot <- BootChainLadder(taylor_cum, R = 2000, process.distr = "od.pois")
  
  cat("\nComparison of methods:\n")
  cat("ODP bootstrap mean:", format(round(mean(odp_boot$IBNR.Totals)), big.mark = ","), "\n")
  cat("NB-CL (MLE) mean:", format(round(mean(boot_mle$reserves_total)), big.mark = ","), "\n")
  cat("NB-CL (corrected) mean:", format(round(mean(boot_corr$reserves_total)), big.mark = ","), "\n")
}
