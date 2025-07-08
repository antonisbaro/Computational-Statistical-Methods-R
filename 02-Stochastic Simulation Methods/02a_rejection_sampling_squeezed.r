################################################################################
# 02a_rejection_sampling_squeezed.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script implements the "Squeezed" Rejection Sampling method to generate
# samples from a standard normal distribution, N(0,1). It uses a Laplace(0,1)
# distribution as the proposal and a specific squeezing function to reduce
# the number of expensive target density calculations.
#
# The script:
# 1. Defines the target (Normal), proposal (Laplace), and squeezing functions.
# 2. Implements the inversion method to sample from the Laplace proposal.
# 3. Runs the main Squeezed Rejection Sampling loop to generate 10,000 samples.
# 4. Performs a theoretical and empirical analysis of the algorithm's
#    efficiency (acceptance rate, cost savings).
# 5. Visualizes the distribution of the generated samples against the
#    theoretical N(0,1) density.
#
# This script corresponds to Exercise 2(a) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)


# --- 2. DEFINE FUNCTIONS FOR SAMPLING ---

# Target density function f(x): Standard Normal N(0,1)
f_target <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

# Proposal density function g(x): Laplace(0,1)
g_proposal <- function(x) {
  0.5 * exp(-abs(x))
}

# Function to generate a sample from g_proposal(x) using the inversion method
generate_from_g_inversion <- function() {
  u <- runif(1)
  if (u < 0.5) {
    x_sample <- log(2 * u)
  } else {
    x_sample <- -log(2 * (1 - u))
  }
  return(x_sample)
}

# Envelope constant M such that f(x) <= M*g(x)
M_constant <- sqrt(2 * exp(1) / pi)

# Squeezing function s(x) such that 0 <= s(x) <= f(x)
s_squeeze <- function(x) {
  term_in_max <- 1 - (x^2 / 2)
  (1 / sqrt(2 * pi)) * pmax(0, term_in_max)
}


# --- 3. SQUEEZED REJECTION SAMPLING LOOP ---

# Parameters
num_samples_target <- 10000
set.seed(12345) # For reproducibility

# Initialize storage and counters
generated_samples <- numeric(num_samples_target)
count_accepted <- 0
total_attempts <- 0
f_calculations_avoided <- 0

cat("Starting Squeezed Rejection Sampling...\n")
while (count_accepted < num_samples_target) {
  total_attempts <- total_attempts + 1
  
  # Step 1: Generate a candidate from the proposal distribution g(x)
  Y_cand <- generate_from_g_inversion()
  
  # Step 2: Generate a uniform random number for the test
  U_rand <- runif(1)
  
  # Calculate required values for the tests
  G_Y_cand <- M_constant * g_proposal(Y_cand) # M*g(y)
  s_Y_cand <- s_squeeze(Y_cand) # s(y)
  
  # Step 3: "Cheap" Squeeze Test
  if (U_rand <= s_Y_cand / G_Y_cand) {
    # Accept based on the squeeze test
    count_accepted <- count_accepted + 1
    generated_samples[count_accepted] <- Y_cand
    f_calculations_avoided <- f_calculations_avoided + 1
  } else {
    # Step 4: "Expensive" Full Test (only if squeeze test fails)
    f_Y_cand <- f_target(Y_cand) # f(y)
    
    if (U_rand <= f_Y_cand / G_Y_cand) {
      # Accept based on the full test
      count_accepted <- count_accepted + 1
      generated_samples[count_accepted] <- Y_cand
    }
    # Else: Reject the candidate and loop again
  }
}
cat("Sampling finished.\n\n")


# --- 4. PERFORMANCE ANALYSIS ---

cat("--- Squeezed Rejection Sampling Process Summary ---\n")
cat(sprintf("Target number of samples: %d\n", num_samples_target))
cat(sprintf("Constant M used: %.4f\n", M_constant))
cat(sprintf("Total proposals generated (attempts): %d\n", total_attempts))
cat(sprintf("Empirical overall acceptance probability: %.4f\n", count_accepted / total_attempts))
cat(sprintf("Theoretical overall acceptance probability (1/M): %.4f\n", 1 / M_constant))
cat(sprintf("Times f(x) calculation was avoided (successful squeezes): %d\n", f_calculations_avoided))
cat(sprintf("Proportion of attempts that were successful squeezes: %.4f\n", f_calculations_avoided / total_attempts))
cat(sprintf("Proportion of accepted samples from successful squeezes: %.4f\n\n", f_calculations_avoided / count_accepted))


# --- 5. VISUALIZATION ---

# Create a data frame for ggplot
samples_df <- data.frame(x_values = generated_samples)

# Create the histogram with the theoretical N(0,1) density overlaid
histogram_squeezed <- ggplot(samples_df, aes(x = x_values)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.25,
                 fill = "#8E3C2E", color = "black", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "#F084C1", linewidth = 1.1, linetype = "dashed") +
  labs(
    title = "Histogram of Samples from N(0,1) via Squeezed Rejection",
    subtitle = paste(num_samples_target, "samples generated. M =", round(M_constant, 4)),
    x = "Generated Value X ~ N(0,1)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the histogram
ggsave("plots/06_rejection_sampling_squeezed.png", plot = histogram_squeezed, width = 8, height = 6)    