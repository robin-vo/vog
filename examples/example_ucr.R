# =============================================================================
# Example: Unified Credibility Reserving (UCR)
# =============================================================================
# 
# This example demonstrates UCR on the Taylor-Ashe triangle and compares
# it to classical methods (Chain-Ladder, Cape Cod, Bornhuetter-Ferguson).
#
# UCR is a credibility-based reserving method that adaptively estimates
# the weight on individual accident year experience vs. pooled information.

library(vog)

# =============================================================================
# 1. TAYLOR-ASHE DATA
# =============================================================================

# Taylor-Ashe triangle (cumulative paid claims)
taylor_ashe <- matrix(c(
  357848, 1124788, 1735330, 2218270, 2745596, 3319994, 3466336, 3606286, 3833515, 3901463,
  352118, 1236139, 2170033, 3353322, 3799067, 4120063, 4647867, 4914039, 5339085, NA,
  290507, 1292306, 2218525, 3235179, 3985995, 4132918, 4628910, 4909315, NA, NA,
  310608, 1418858, 2195047, 3757447, 4029929, 4381982, 4588268, NA, NA, NA,
  443160, 1136350, 2128333, 2897821, 3402672, 3873311, NA, NA, NA, NA,
  396132, 1333217, 2180715, 2985752, 3691712, NA, NA, NA, NA, NA,
  440832, 1288463, 2419861, 3483130, NA, NA, NA, NA, NA, NA,
  359480, 1421128, 2864498, NA, NA, NA, NA, NA, NA, NA,
  376686, 1363294, NA, NA, NA, NA, NA, NA, NA, NA,
  344014, NA, NA, NA, NA, NA, NA, NA, NA, NA
), nrow = 10, ncol = 10, byrow = TRUE)

rownames(taylor_ashe) <- paste0("AY", 1:10)
colnames(taylor_ashe) <- paste0("Dev", 1:10)

print(taylor_ashe)

# =============================================================================
# 2. FIT UCR
# =============================================================================

# Fit UCR (uses first period claims as exposure by default)
ucr_fit <- ucr(taylor_ashe)
print(ucr_fit)

# Detailed summary
summary(ucr_fit)

# =============================================================================
# 3. COMPARE METHODS
# =============================================================================

# Compare all methods
comparison <- compare_reserves(taylor_ashe)
print(comparison)

# =============================================================================
# 4. VISUALIZE
# =============================================================================

# Diagnostic plots
par(mfrow = c(1, 3))
plot(ucr_fit)

# =============================================================================
# 5. MACK MODEL
# =============================================================================

# Fit Mack's model for variance estimation
mack_fit <- fit_mack(taylor_ashe)
print(mack_fit)

# Prediction intervals
predict_intervals(mack_fit, level = 0.95)

# Diagnostics
par(mfrow = c(1, 3))
plot(mack_fit)

# =============================================================================
# 6. SIMULATION STUDY
# =============================================================================

# Demonstrate UCR's adaptive behavior
cat("\n=== Simulation Study ===\n\n")

# Generate triangles with different heterogeneity levels
set.seed(123)

# Homogeneous (tau^2 = 0) - UCR should favor Cape Cod
sim_homo <- simulate_triangle(I = 10, J = 10, mu = 0.1, tau_sq = 0)
ucr_homo <- ucr(sim_homo$triangle, sim_homo$exposure)
cat("Homogeneous (tau^2 = 0):\n")
cat("  tau^2_hat:", round(ucr_homo$tau_sq_hat, 6), "\n")
cat("  Mean Z1:", round(mean(ucr_homo$Z1), 3), "(low = Cape Cod)\n\n")

# Heterogeneous (tau^2 = 0.01) - UCR should favor Chain-Ladder
sim_hetero <- simulate_triangle(I = 10, J = 10, mu = 0.1, tau_sq = 0.01)
ucr_hetero <- ucr(sim_hetero$triangle, sim_hetero$exposure)
cat("Heterogeneous (tau^2 = 0.01):\n")
cat("  tau^2_hat:", round(ucr_hetero$tau_sq_hat, 6), "\n")
cat("  Mean Z1:", round(mean(ucr_hetero$Z1), 3), "(high = Chain-Ladder)\n")

# =============================================================================
# 7. CREDIBILITY THEORY EXAMPLES
# =============================================================================

cat("\n=== Credibility Theory ===\n\n")

# Bühlmann credibility
x <- c(8, 12, 9, 11)  # Individual observations
mu <- 10              # Collective mean
sigma_sq <- 4         # Within-risk variance
tau_sq <- 1           # Between-risk variance

buhl <- buhlmann(x, mu, sigma_sq, tau_sq)
print(buhl)

# Bühlmann-Straub with varying weights
x_rates <- c(0.08, 0.12, 0.095, 0.11)
w <- c(100, 200, 150, 50)  # Exposure weights

bs_params <- estimate_bs_params(x_rates, w)
cat("\nBühlmann-Straub parameter estimates:\n")
cat("  mu:", round(bs_params$mu, 4), "\n")
cat("  tau^2:", round(bs_params$tau_sq, 6), "\n")

bs_fit <- buhlmann_straub(x_rates, w, bs_params$mu, bs_params$sigma_sq, bs_params$tau_sq)
print(bs_fit)

# =============================================================================
# 8. CREDIBILITY-BAYES EQUIVALENCE
# =============================================================================

cat("\n=== Credibility = Bayesian Posterior Mean ===\n\n")

# Show that credibility estimator equals posterior mean
equiv <- credibility_bayes_equivalence(n_obs = 50, exposure = 1000, mu = 0.05, tau_sq = 0.001)
cat("Bayesian posterior mean:", round(equiv$bayesian_mean, 6), "\n")
cat("Credibility estimate:   ", round(equiv$credibility_estimate, 6), "\n")
cat("Difference:             ", round(equiv$difference, 10), "(should be ~0)\n")
