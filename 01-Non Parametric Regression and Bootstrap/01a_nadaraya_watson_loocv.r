################################################################################
# 01a_nadaraya_watson_loocv.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script performs bandwidth selection for the Nadaraya-Watson kernel
# regression estimator using brute-force Leave-One-Out Cross-Validation (LOOCV).
# It visualizes the MSE curve and then uses the optimal bandwidth to plot
# the final regression fit over the original data.
#
# This script corresponds to Exercise 1(a)(i) of the
# "Computational Statistics & Stochastic Optimization" project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

# Load necessary libraries for plotting
library(ggplot2)

# Load the dataset
# Assumes the script is run from the project's root directory.
tryCatch({
  data_pairs <- readRDS("data/data1pairs.rds")
}, error = function(e) {
  stop("Error loading data. Make sure 'data/data1pairs.rds' exists and the working directory is correct.")
})

# For convenience, assign X and Y to separate vectors
X <- data_pairs$X
Y <- data_pairs$Y
n <- length(X)


# --- 2. LOOCV GRID SEARCH FOR OPTIMAL BANDWIDTH ---

# Define a sequence of candidate bandwidths (hx_values)
hx_values <- seq(0.1, 2.5, by = 0.05)

# Initialize a vector to store the Cross-Validation MSE for each hx
cv_mse <- numeric(length(hx_values))

# Loop through each candidate bandwidth
cat("Starting LOOCV grid search...\n")
for (j in 1:length(hx_values)) {
  h_current <- hx_values[j]
  squared_errors <- numeric(n) # To store squared errors for this h_current

  # LOOCV loop: for each data point i, use it as the test point
  for (i in 1:n) {
    # Data for training (all points except i)
    X_train <- X[-i]
    Y_train <- Y[-i]

    # Point to predict
    x_test_point <- X[i]

    # Use ksmooth to get the prediction for the left-out point
    ksmooth_fit_loocv <- ksmooth(
      x = X_train,
      y = Y_train,
      kernel = "normal",
      bandwidth = h_current,
      x.points = x_test_point
    )

    y_pred_loocv <- ksmooth_fit_loocv$y
    squared_errors[i] <- (Y[i] - y_pred_loocv)^2
  }

  # Calculate the mean squared error for the current hx
  cv_mse[j] <- mean(squared_errors, na.rm = TRUE)

  # Optional: Print progress
  # cat(sprintf("Bandwidth: %.2f | LOOCV MSE: %.6f\n", h_current, cv_mse[j]))
}
cat("LOOCV grid search finished.\n")


# --- 3. PROCESS AND DISPLAY RESULTS ---

# Find the optimal hx that minimizes CV(MSE)
optimal_hx_index <- which.min(cv_mse)
optimal_hx <- hx_values[optimal_hx_index]
min_cv_mse <- cv_mse[optimal_hx_index]

# Print the results to the console
cat("\n--- Naive LOOCV Results ---\n")
cat(sprintf("Optimal bandwidth (hx): %.2f\n", optimal_hx))
cat(sprintf("Minimum CV(MSE): %.6f\n", min_cv_mse))


# --- 4. VISUALIZE LOOCV MSE vs. BANDWIDTH ---

# Create a data frame for ggplot
plot_data_gg <- data.frame(Bandwidth = hx_values, MSE = cv_mse)

# Define the label for the optimal hx line in the legend
optimal_hx_legend_label <- paste("Optimal hx =", format(round(optimal_hx, 2), nsmall = 2))

# Create the ggplot object
loocv_gg_plot <- ggplot(plot_data_gg, aes(x = Bandwidth, y = MSE)) +
  # MSE Curve
  geom_line(color = "#8E3C2E", linewidth = 1.2) +
  # Vertical line for optimal hx
  geom_vline(aes(xintercept = optimal_hx, linetype = "Optimal"), color = "#F084C1", linewidth = 0.8) +
  # Point marking the optimal hx and minimum MSE
  geom_point(aes(x = optimal_hx, y = min_cv_mse), color = "#F084C1", size = 2.5, shape = 19) +
  # Define the linetype and label for the legend
  scale_linetype_manual(
    name = "",
    values = c("Optimal" = "dashed"),
    labels = c("Optimal" = optimal_hx_legend_label)
  ) +
  # Labels and Title
  labs(
    title = "LOOCV MSE for Nadaraya-Watson",
    x = "Bandwidth (hx)",
    y = "LOOCV MSE"
  ) +
  # Theme and styling
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.6), colour = "grey70", linewidth = 0.5),
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_blank()
  )

# Save the 1plot
ggsave("plots/01_loocv_mse_vs_bandwidth.png", plot = loocv_gg_plot, width = 8, height = 6)


# --- 5. VISUALIZE FINAL NADARAYA-WATSON FIT ---

cat("Generating final regression fit plot...\n")
# Generate a sequence of x values for plotting the smooth curve
x_grid <- seq(min(X), max(X), length.out = 200)

# Get the Nadaraya-Watson estimate using the optimal hx
nw_estimate_optimal <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = optimal_hx, x.points = x_grid)

# Prepare data for ggplot
plot_data_scatter <- data.frame(X_orig = X, Y_orig = Y)
plot_data_line <- data.frame(X_smooth = nw_estimate_optimal$x, Y_smooth = nw_estimate_optimal$y)

# Create the final regression plot
final_plot <- ggplot() +
  # Scatter plot of the original data
  geom_point(data = plot_data_scatter, aes(x = X_orig, y = Y_orig), color = "#F084C1", alpha = 0.7, size = 2) +
  # Add the Nadaraya-Watson smoothed line
  geom_line(data = plot_data_line, aes(x = X_smooth, y = Y_smooth), color = "#8E3C2E", linewidth = 1.2) +
  # Add titles and labels
  labs(
    title = paste("Nadaraya-Watson Regression (Gaussian Kernel, hx =", round(optimal_hx, 2), ")"),
    x = "X",
    y = "Y",
    caption = "Data: data1pairs.rds"
  ) +
  # Apply a theme
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, face = "italic"),
    panel.background = element_rect(fill = "#FBFAF7")
  )

# Save the second plot
ggsave("plots/02_nadaraya_watson_fit.png", plot = final_plot, width = 8, height = 6)