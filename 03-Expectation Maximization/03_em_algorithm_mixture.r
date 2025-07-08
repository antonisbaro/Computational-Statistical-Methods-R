################################################################################
# 03_em_algorithm_mixture.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script implements the Expectation-Maximization (EM) algorithm to
# estimate the mixing probability 'p' of a mixture of two exponential
# distributions: p * Exp(1) + (1-p) * Exp(5).
#
# The script:
# 1. Loads the dataset containing observations from the mixture distribution.
# 2. Initializes the algorithm with a starting guess for p.
# 3. Iteratively performs the E-step (calculating posterior probabilities,
#    or responsibilities, gamma_i) and the M-step (updating the estimate
#    of p) until convergence.
# 4. Visualizes the convergence path of the parameter p over the iterations.
#
# This script corresponds to Exercise 3 of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)

# Load the dataset
tryCatch({
  data_em <- readRDS("data/data3em.rds")
}, error = function(e) {
  stop("Error loading data. Make sure 'data/data3em.rds' exists and the working directory is correct.")
})

# For convenience, assign data to X_obs and get n
X_obs <- data_em
n_em <- length(X_obs)


# --- 2. DEFINE DENSITY FUNCTIONS ---

# Density function for the first component: Exp(rate = 1)
f1_x <- function(x) {
  dexp(x, rate = 1)
}

# Density function for the second component: Exp(rate = 5)
f0_x <- function(x) {
  dexp(x, rate = 5)
}


# --- 3. EM ALGORITHM IMPLEMENTATION ---

# Initialization
p_current <- 0.5            # Initial guess for p
max_iterations <- 100000    # Safety limit to prevent infinite loops
tolerance <- 1e-10          # Convergence criterion: |p_new - p_old|
iterations <- 0
converged <- FALSE
p_history <- numeric(max_iterations + 1) # Vector to store the history of p
p_history[1] <- p_current

cat("Starting EM algorithm...\n")
# Pre-calculate densities as they don't change within the loop
f1_X_obs <- f1_x(X_obs)
f0_X_obs <- f0_x(X_obs)

# EM Iteration loop
for (r in 1:max_iterations) {
  iterations <- r
  p_old <- p_current
  
  # --- E-Step: Calculate responsibilities (gamma_i) ---
  # This is the expected value of the latent variable Z_i given X_obs and p_old
  numerator_gamma <- p_old * f1_X_obs
  denominator_gamma <- numerator_gamma + (1 - p_old) * f0_X_obs
  
  # Handle potential division by zero if a denominator is 0
  # (unlikely with exponential densities but good practice)
  gamma_i <- ifelse(denominator_gamma == 0, 0, numerator_gamma / denominator_gamma)
  
  # --- M-Step: Update mixing probability p ---
  # The new p is the average of the responsibilities
  p_current <- mean(gamma_i)
  p_history[r + 1] <- p_current
  
  # Check for convergence
  if (abs(p_current - p_old) < tolerance) {
    converged <- TRUE
    break # Exit loop
  }
}

# Trim history to the actual number of iterations
p_history <- p_history[1:(iterations + 1)]


# --- 4. OUTPUT RESULTS ---

cat("EM algorithm finished.\n\n")
if (converged) {
  cat(sprintf("EM algorithm converged in %d iterations.\n", iterations))
  cat(sprintf("Final estimate for p (p_hat_EM): %.12f\n\n", p_current))
} else {
  cat(sprintf("EM algorithm did NOT converge after %d iterations.\n", max_iterations))
  cat(sprintf("Last estimate for p: %.12f\n\n", p_current))
}


# --- 5. VISUALIZE CONVERGENCE ---

# Create a data frame for plotting
convergence_df <- data.frame(Iteration = 0:iterations, P_Estimate = p_history)

# Create the convergence plot
convergence_plot_em <- ggplot(convergence_df, aes(x = Iteration, y = P_Estimate)) +
  geom_line(color = "#8E3C2E", linewidth = 1.5) +
  geom_point(color = "#F084C1", size = 2.8, shape = 16) +
  geom_hline(yintercept = p_current, linetype = "dashed", color = "#4F8E38", linewidth = 1.2) +
  annotate(
    "text", x = iterations / 2, y = p_current + 0.02, # Position annotation nicely
    label = paste("Converged p =", round(p_current, 6)),
    color = "#4F8E38", size = 3.5, fontface = "bold"
  ) +
  labs(
    title = "Convergence of Parameter p in EM Algorithm",
    x = "Iteration (r)",
    y = expression(hat(p)^{(r)}) # Using hat(p) for estimate
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the plot
ggsave("plots/07_em_convergence.png", plot = convergence_plot_em, width = 8, height = 6)