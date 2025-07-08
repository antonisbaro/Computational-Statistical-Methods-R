################################################################################
# 01d_bootstrap_parametric_student_t.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script implements a parametric bootstrap to estimate the sampling
# distribution of the statistic T = min(X1, ..., Xn), assuming the data
# follows a location-scale Student's t-distribution.
#
# The script:
# 1. Defines custom functions (density, cdf, quantile, random generation)
#    for the location-scale Student's t-distribution.
# 2. Fits the t-distribution to the data using Maximum Likelihood Estimation
#    (MLE) via the `fitdistrplus` package to estimate the parameters
#    (df, mu, sigma).
# 3. Generates B bootstrap samples from this fitted distribution.
# 4. Calculates the minimum for each sample and visualizes the resulting
#    distribution, comparing it to the non-parametric approach.
#
# This script corresponds to Exercise 1(b)(ii) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)
library(fitdistrplus)

# Load the dataset for question 1.b
tryCatch({
  data_1b <- readRDS("data/data1b.rds")
}, error = function(e) {
  stop("Error loading data. Make sure 'data/data1b.rds' exists and the working directory is correct.")
})


# --- 2. DEFINE & FIT PARAMETRIC MODEL (LOCATION-SCALE T-DISTRIBUTION) ---

# Custom functions for the location-scale t-distribution for `fitdist`
dlst <- function(x, df, mu, sigma) {
  if (any(sigma <= 0) || any(df <= 0)) return(rep(NaN, length(x)))
  return(1/sigma * dt((x - mu)/sigma, df))
}
plst <- function(q, df, mu, sigma) {
  if (any(sigma <= 0) || any(df <= 0)) return(rep(NaN, length(q)))
  return(pt((q - mu)/sigma, df))
}
qlst <- function(p, df, mu, sigma) {
  if (any(sigma <= 0) || any(df <= 0)) return(rep(NaN, length(p)))
  return(mu + sigma * qt(p, df))
}
rlst <- function(n, df, mu, sigma) {
  if (any(sigma <= 0) || any(df <= 0)) return(rep(NaN, n))
  return(mu + sigma * rt(n, df))
}

# Provide reasonable starting values for the MLE fitting procedure
mu_start <- mean(data_1b)
sigma_start <- sd(data_1b)
df_start <- 5 # A common starting point for degrees of freedom

# Fit the distribution using MLE, with error handling
cat("Fitting Student's t-distribution to the data...\n")
student_fit <- NULL
tryCatch({
  student_fit <- fitdist(
    data_1b, "lst",
    start = list(df = df_start, mu = mu_start, sigma = sigma_start),
    lower = c(0.001, -Inf, 0.001) # Constraints: df > 0 and sigma > 0
  )
}, error = function(e) {
  stop(paste("Error during fitting:", e$message))
})

# Extract estimated parameters
est_params_student <- coef(student_fit)
nu_hat <- est_params_student["df"]
mu_hat <- est_params_student["mu"]
sigma_hat <- est_params_student["sigma"]

cat("Fitting complete.\n\n--- Fitted Parameters ---\n")
print(summary(student_fit))
cat("\n")


# --- 3. PARAMETRIC BOOTSTRAP ---

# Parameters for Bootstrap
B <- 2000
n_obs <- length(data_1b)

# Vector to store the T* values from parametric bootstrap
T_star_parametric_values <- numeric(B)

cat("Starting parametric bootstrap...\n")
set.seed(456) # For reproducibility
# Parametric Bootstrap loop
for (b in 1:B) {
  # Step 1: Generate a sample from the FITTED parametric distribution
  current_parametric_bootstrap_sample <- rlst(n = n_obs, df = nu_hat, mu = mu_hat, sigma = sigma_hat)
  
  # Step 2: Calculate the statistic T* for the current sample
  T_star_parametric_values[b] <- min(current_parametric_bootstrap_sample)
}
cat("Bootstrap finished.\n\n")


# --- 4. VISUALIZATION ---

# Create a data frame for ggplot
T_star_parametric_df <- data.frame(T_star_param = T_star_parametric_values)

# Get the minimum of the original sample for reference
original_sample_min <- min(data_1b)

# Create the histogram
histogram_T_star_parametric <- ggplot(T_star_parametric_df, aes(x = T_star_param)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.1,
                 fill = "#8E3C2E", color = "black", alpha = 0.7) +
  geom_density(alpha = 0.2, fill = "#F084C1", color = "#FF00B9", linewidth = 1) +
  geom_vline(aes(xintercept = original_sample_min, linetype = "Original Sample Min"),
             color = "#4F8E38", linewidth = 1.5) +
  scale_linetype_manual(name = "", values = c("Original Sample Min" = "dashed")) +
  labs(
    title = "Histogram of Parametric Bootstrap Estimates for T = min(X_i)",
    subtitle = paste0(
      "Assuming Student-t distribution; B = ", B, " samples; n = ", n_obs,
      "\nFitted params: df=", round(nu_hat, 2), ", mu=", round(mu_hat, 2), ", sigma=", round(sigma_hat, 2)
    ),
    x = "T* (Minimum of Parametric Bootstrap Sample)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, lineheight = 1.1),
    legend.position = "top",
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the histogram
ggsave("plots/05_bootstrap_parametric_min.png", plot = histogram_T_star_parametric, width = 8, height = 6)