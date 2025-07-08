################################################################################
# 02b_monte_carlo_vs_importance_sampling.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script estimates the value of the integral E[phi(X)] where
# phi(X) = 4 / (1 + X^2) and X ~ U(0,1). The estimation is performed
# using two different stochastic methods:
#
# 1. Classic Monte Carlo Integration: Samples are drawn directly from the
#    uniform distribution U(0,1).
# 2. Importance Sampling: Samples are drawn from a carefully chosen proposal
#    distribution g(x) to reduce the variance of the estimator.
#
# The script calculates the standard error for both estimators based on
# multiple simulations and compares their efficiency, demonstrating the
# variance reduction achieved by Importance Sampling.
#
# This script corresponds to Exercise 2(b) of the project.
#
################################################################################


# --- 1. SETUP & COMMON DEFINITIONS ---

# Define the function to be integrated
phi_x <- function(x) {
  4 / (1 + x^2)
}

# Parameters for the simulation (common for both methods)
n_sample_size <- 200  # Sample size for each estimate
N_simulations <- 1000 # Number of estimates to generate for SE calculation


# --- 2. CLASSIC MONTE CARLO ESTIMATION ---

cat("--- Running Classic Monte Carlo Simulation ---\n")
# Vector to store the N_simulations estimates of theta_1
theta_1_estimates <- numeric(N_simulations)
set.seed(2024) # For reproducibility

# Loop to generate N_simulations estimates
for (k in 1:N_simulations) {
  # Step 1: Generate a sample from the target distribution f(x) = U(0,1)
  X_sample <- runif(n_sample_size, min = 0, max = 1)
  
  # Step 2: Apply the function phi to the sample
  phi_X_values <- phi_x(X_sample)
  
  # Step 3: The MC estimate is the mean
  theta_1_estimates[k] <- mean(phi_X_values)
}

# Calculate the standard error from the simulation results
se_theta_1_hat <- sd(theta_1_estimates)

cat(sprintf("Mean of Classic MC estimates: %.6f\n", mean(theta_1_estimates)))
cat(sprintf("Estimated Standard Error (Classic MC): %.6f\n", se_theta_1_hat))
cat(sprintf("True value (pi): %.6f\n\n", pi))


# --- 3. IMPORTANCE SAMPLING ESTIMATION ---

cat("--- Running Importance Sampling Simulation ---\n")
# Importance sampling density g(x)
g_importance <- function(x) {
  ifelse(x >= 0 & x <= 1, (1/3) * (4 - 2*x), 0)
}

# Function to generate samples from g_importance(x) using inversion method
generate_from_g_importance <- function() {
  u <- runif(1)
  # This is the solution to 2 - sqrt(4 - 3*u) that lies in [0, 1]
  return(2 - sqrt(4 - 3*u))
}

# Importance weight function w(x) = f(x)/g(x), where f(x)=1 for x in [0,1]
importance_weight <- function(x) {
  1 / g_importance(x)
}

# The term to average in Importance Sampling: psi(x) = phi(x) * w(x)
psi_x <- function(x) {
  phi_x(x) * importance_weight(x)
}

# Vector to store the N_simulations estimates of theta_2
theta_2_estimates <- numeric(N_simulations)
set.seed(29) # For reproducibility

# Loop to generate N_simulations estimates
for (k in 1:N_simulations) {
  # Step 1: Generate a sample from the proposal distribution g(x)
  X_sample_g <- replicate(n_sample_size, generate_from_g_importance())
  
  # Step 2: Calculate psi(X_i) for each observation
  psi_X_values <- psi_x(X_sample_g)
  
  # Step 3: The IS estimate is the mean
  theta_2_estimates[k] <- mean(psi_X_values)
}

# Calculate the standard error from the simulation results
se_theta_2_hat <- sd(theta_2_estimates)

cat(sprintf("Mean of Importance Sampling estimates: %.6f\n", mean(theta_2_estimates)))
cat(sprintf("Estimated Standard Error (Importance Sampling): %.6f\n\n", se_theta_2_hat))


# --- 4. COMPARISON OF METHODS ---

cat("--- Comparison of Standard Errors ---\n")
if (se_theta_2_hat < se_theta_1_hat) {
  reduction_percentage <- (1 - (se_theta_2_hat / se_theta_1_hat)) * 100
  cat(sprintf("Importance Sampling resulted in a variance reduction.\n"))
  cat(sprintf("Standard Error was reduced by approximately %.2f%% compared to Classic MC.\n", reduction_percentage))
} else {
  cat("Importance Sampling did NOT result in a reduction of Standard Error for this g(x).\n")
}       