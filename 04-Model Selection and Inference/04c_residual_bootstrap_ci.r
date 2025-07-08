################################################################################
# 04c_residual_bootstrap_ci.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script applies the Residual Bootstrap method to construct a 95%
# confidence interval for the regression coefficient of the 'rm' variable
# (average number of rooms per dwelling) from the best model (M1) found
# using the AIC criterion in script 04a.
#
# The script:
# 1. Defines and fits the best model (M1) to the Boston dataset.
# 2. Extracts the original OLS estimate for the 'rm' coefficient and the
#    model's residuals.
# 3. Centers the residuals to ensure they have a mean of zero.
# 4. Performs B bootstrap iterations:
#    a. Resamples the centered residuals.
#    b. Creates a new bootstrap response variable y*.
#    c. Refits the model M1 on the bootstrap data.
#    d. Stores the new estimate of the 'rm' coefficient.
# 5. Constructs a 95% percentile confidence interval from the distribution
#    of the bootstrap estimates.
# 6. Visualizes the bootstrap distribution with the original estimate and CI.
#
# This script corresponds to Exercise 4(c) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)
library(MASS)
data(Boston)


# --- 2. FIT ORIGINAL MODEL (M1) & EXTRACT COMPONENTS ---

# Define the formula for Model M1 (best model from AIC)
formula_M1_str <- "medv ~ crim + zn + chas + nox + rm + dis + rad + tax + ptratio + black + lstat"
formula_M1 <- as.formula(formula_M1_str)

# Fit the original Model M1
model_M1_fit <- lm(formula_M1, data = Boston)

# Extract the original OLS estimate for the coefficient of 'rm'
original_beta_rm_hat <- coef(model_M1_fit)["rm"]

# Get residuals and fitted values from the original model
residuals_M1 <- residuals(model_M1_fit)
fitted_values_M1 <- fitted(model_M1_fit)

# Center the residuals (important for bootstrap validity)
centered_residuals_M1 <- residuals_M1 - mean(residuals_M1)
n_obs_boston <- length(centered_residuals_M1)


# --- 3. RESIDUAL BOOTSTRAP LOOP ---

# Parameters for Bootstrap
B_bootstrap <- 2000
bootstrap_beta_rm_estimates <- numeric(B_bootstrap) # To store estimates

cat("Starting Residual Bootstrap...\n")
set.seed(4567) # For reproducibility

for (b in 1:B_bootstrap) {
  # Step 1: Generate a bootstrap sample of residuals
  bootstrap_residuals_sample <- sample(centered_residuals_M1, size = n_obs_boston, replace = TRUE)
  
  # Step 2: Create a new bootstrap response variable y*
  y_star_bootstrap <- fitted_values_M1 + bootstrap_residuals_sample
  
  # Step 3: Create a temporary bootstrap data frame
  bootstrap_data <- Boston
  bootstrap_data$medv <- y_star_bootstrap
  
  # Step 4: Refit Model M1 using the bootstrap data
  refitted_model_M1_bootstrap <- lm(formula_M1, data = bootstrap_data)
  
  # Step 5: Store the coefficient for 'rm' from this bootstrap iteration
  bootstrap_beta_rm_estimates[b] <- coef(refitted_model_M1_bootstrap)["rm"]
  
  if (b %% 200 == 0) cat(sprintf("... Completed %d / %d bootstrap iterations.\n", b, B_bootstrap))
}
cat("Bootstrap finished.\n\n")


# --- 4. CONSTRUCT CONFIDENCE INTERVAL ---

# Construct the 95% percentile bootstrap confidence interval
alpha_CI <- 0.05
percentile_CI_beta_rm <- quantile(
  bootstrap_beta_rm_estimates,
  probs = c(alpha_CI / 2, 1 - (alpha_CI / 2)),
  na.rm = TRUE
)

cat("--- Residual Bootstrap Results for beta_rm (Coefficient of 'rm') ---\n")
cat(sprintf("Original OLS estimate for beta_rm (from Model M1): %.4f\n", original_beta_rm_hat))
cat(sprintf("Number of Bootstrap Samples (B): %d\n", B_bootstrap))
cat(sprintf(
  "95%% Percentile Bootstrap Confidence Interval for beta_rm: [%.4f, %.4f]\n\n",
  percentile_CI_beta_rm[1], percentile_CI_beta_rm[2]
))


# --- 5. VISUALIZATION ---

# Create a data frame for the histogram
hist_df_beta_rm <- data.frame(beta_rm_star = bootstrap_beta_rm_estimates)

# Create the plot
plot_hist_beta_rm <- ggplot(hist_df_beta_rm, aes(x = beta_rm_star)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, fill = "#8E3C2E", color = "black", alpha = 0.7) +
  # Add vertical line for the original OLS estimate
  geom_vline(xintercept = original_beta_rm_hat, color = "#F084C1", linetype = "dashed", linewidth = 1.5) +
  # Add vertical lines for the 95% Confidence Interval
  geom_vline(xintercept = percentile_CI_beta_rm[1], color = "#4F8E38", linetype = "dotted", linewidth = 1.2) +
  geom_vline(xintercept = percentile_CI_beta_rm[2], color = "#4F8E38", linetype = "dotted", linewidth = 1.2) +
  labs(
    title = "Bootstrap Distribution of beta_rm Coefficient",
    subtitle = "From Residual Bootstrap on Model M1",
    x = expression(paste("Bootstrap Estimate (", beta[rm]^"*", ")")),
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the plot
ggsave("plots/09_residual_bootstrap_rm_ci.png", plot = plot_hist_beta_rm, width = 8, height = 6, dpi = 300)