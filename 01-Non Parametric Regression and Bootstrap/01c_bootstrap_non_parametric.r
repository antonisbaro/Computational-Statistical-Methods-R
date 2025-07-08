################################################################################
# 01c_bootstrap_non_parametric.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script implements a non-parametric bootstrap to estimate the sampling
# distribution of the statistic T = min(X1, ..., Xn).
# It generates B bootstrap samples, calculates the minimum for each,
# and visualizes the resulting distribution with a histogram.
# The script also analyzes the key limitation of non-parametric bootstrap
# for extreme value statistics, showing the high frequency of the original
# sample minimum in the bootstrap estimates.
#
# This script corresponds to Exercise 1(b)(i) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)

# Load the dataset for question 1.b
tryCatch({
  data_1b <- readRDS("data/data1b.rds")
}, error = function(e) {
  stop("Error loading data. Make sure 'data/data1b.rds' exists and the working directory is correct.")
})


# --- 2. NON-PARAMETRIC BOOTSTRAP ---

# Parameters for Bootstrap
B <- 2000 # Number of Bootstrap samples
n_obs <- length(data_1b)

# Vector to store the T* values (minimum of each Bootstrap sample)
T_star_values <- numeric(B)

cat("Starting non-parametric bootstrap...\n")
# Bootstrap loop
set.seed(123) # For reproducibility
for (b in 1:B) {
  # Step 1: Generate a Bootstrap sample by sampling with replacement
  current_bootstrap_sample <- sample(data_1b, size = n_obs, replace = TRUE)
  
  # Step 2: Calculate the statistic T* for the current Bootstrap sample
  T_star_values[b] <- min(current_bootstrap_sample)
}
cat("Bootstrap finished.\n\n")


# --- 3. ANALYSIS OF BOOTSTRAP RESULTS ---

# Calculate the minimum of the original sample for reference
original_sample_min <- min(data_1b)

# Further analysis: How many T_star values are equal to the original sample minimum
count_equal_to_original_min <- sum(T_star_values == original_sample_min)
percentage_equal_to_original_min <- (count_equal_to_original_min / B) * 100

# Print analysis results
cat("--- Non-Parametric Bootstrap Analysis ---\n")
cat(sprintf("Minimum of the original sample: %.4f\n", original_sample_min))
cat(sprintf("Number of T* values equal to original sample min: %d\n", count_equal_to_original_min))
cat(sprintf("Percentage of T* values equal to original sample min: %.2f%%\n", percentage_equal_to_original_min))


# --- 4. VISUALIZATION ---

# Create a data frame for ggplot
T_star_df <- data.frame(T_star = T_star_values)

# Create the histogram
histogram_T_star <- ggplot(T_star_df, aes(x = T_star)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.04,
                 fill = "#8E3C2E", color = "black", alpha = 0.7) +
  geom_density(alpha = 0.2, fill = "#F084C1", color = "#FF00B9", linewidth = 1) +
  geom_vline(aes(xintercept = original_sample_min, linetype = "Original Sample Min"), 
             color = "#4F8E38", linewidth = 1.5) +
  scale_linetype_manual(name = "", values = c("Original Sample Min" = "dashed")) +
  labs(
    title = "Histogram of Non-Parametric Bootstrap Estimates for T = min(X_i)",
    subtitle = paste("B =", B, "Bootstrap samples; n =", n_obs),
    x = "T* (Minimum of Bootstrap Sample)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top",
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the histogram
ggsave("plots/04_bootstrap_non_parametric_min.png", plot = histogram_T_star, width = 8, height = 6)