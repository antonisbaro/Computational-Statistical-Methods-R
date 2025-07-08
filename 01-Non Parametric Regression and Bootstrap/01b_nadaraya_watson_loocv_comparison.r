################################################################################
# 01b_nadaraya_watson_loocv_comparison.R
#
# Author: Antonios Barotsakis
# Date: 2024-05-21
#
# Description:
# This script compares three different methods for calculating the
# Leave-One-Out Cross-Validation (LOOCV) error for the Nadaraya-Watson
# kernel estimator: Naive (ksmooth), Efficient (manual), and Manual Naive.
# It then generates a detailed comparison plot, replicating the final
# visualization from the project report, including annotations for optimal
# bandwidths for each method.
#
# This script corresponds to Exercise 1(a)(ii) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(ggplot2)

# Load the dataset
tryCatch({
  data_pairs <- readRDS("data/data1pairs.rds")
}, error = function(e) {
  stop("Error loading data. Make sure 'data/data1pairs.rds' exists and the working directory is correct.")
})

X <- data_pairs$X
Y <- data_pairs$Y
n <- length(X)

# Define a sequence of candidate bandwidths
hx_values <- seq(0.1, 2.5, by = 0.05)


# --- 2. DEFINE KERNEL FUNCTION ---

gaussian_kernel <- function(u) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)
}
K_0 <- gaussian_kernel(0)


# --- 3. RUN ALL THREE LOOCV METHODS ---

# Initialize vectors to store results
cv_mse_naive_ksmooth <- numeric(length(hx_values))
cv_mse_efficient <- numeric(length(hx_values))
cv_mse_manual_naive <- numeric(length(hx_values))

cat("Starting LOOCV calculations for all three methods...\n")

for (j in 1:length(hx_values)) {
  h_current <- hx_values[j]
  sq_err_naive_ksmooth <- numeric(n)
  sq_err_efficient <- numeric(n)
  sq_err_manual_naive <- numeric(n)
  
  for (i in 1:n) {
    # Method 1: Naive ksmooth
    ksmooth_fit <- ksmooth(x = X[-i], y = Y[-i], kernel = "normal", bandwidth = h_current, x.points = X[i])
    sq_err_naive_ksmooth[i] <- (Y[i] - ksmooth_fit$y)^2
    
    # Method 2: Efficient
    u_values_eff <- (X[i] - X) / h_current; weights_ij_eff <- gaussian_kernel(u_values_eff)
    S_i_h <- sum(weights_ij_eff * Y); D_i_h <- sum(weights_ij_eff)
    y_pred_eff <- (S_i_h - K_0 * Y[i]) / (D_i_h - K_0)
    sq_err_efficient[i] <- (Y[i] - y_pred_eff)^2
    
    # Method 3: Manual Naive
    u_values_man <- (X[i] - X[-i]) / h_current; kernel_weights_man <- gaussian_kernel(u_values_man)
    y_pred_man <- sum(kernel_weights_man * Y[-i]) / sum(kernel_weights_man)
    sq_err_manual_naive[i] <- (Y[i] - y_pred_man)^2
  }
  
  cv_mse_naive_ksmooth[j] <- mean(sq_err_naive_ksmooth, na.rm = TRUE)
  cv_mse_efficient[j] <- mean(sq_err_efficient, na.rm = TRUE)
  cv_mse_manual_naive[j] <- mean(sq_err_manual_naive, na.rm = TRUE)
}
cat("All calculations finished.\n\n")


# --- 4. PREPARE DATA & PLOT SETTINGS ---

# Find optimal hx values for each method
optimal_hx_naive_ksmooth <- hx_values[which.min(cv_mse_naive_ksmooth)]
optimal_hx_manual_naive <- hx_values[which.min(cv_mse_manual_naive)]
optimal_hx_efficient <- hx_values[which.min(cv_mse_efficient)]

# Create the combined data frame for plotting
plot_comparison_df <- data.frame(
  Bandwidth = rep(hx_values, 3),
  MSE = c(cv_mse_naive_ksmooth, cv_mse_manual_naive, cv_mse_efficient),
  Method = factor(
    rep(c("Naive (ksmooth)", "Manual Naive (by hand)", "Efficient (by hand)"), each = length(hx_values)),
    levels = c("Naive (ksmooth)", "Manual Naive (by hand)", "Efficient (by hand)")
  )
)

# Define customizable colors and linetypes exactly as in the report
color_palette <- c("Naive (ksmooth)"= "#8E3C2E", "Manual Naive (by hand)"= "#F084C1", "Efficient (by hand)"= "#4F8E38")
linetype_palette <- c("Naive (ksmooth)"= "solid", "Manual Naive (by hand)"= "dashed", "Efficient (by hand)"= "dotted")

# Create a list of optimal hx values for easier access in ggplot
optimal_values <- list(
  "Naive (ksmooth)" = optimal_hx_naive_ksmooth,
  "Manual Naive (by hand)" = optimal_hx_manual_naive,
  "Efficient (by hand)" = optimal_hx_efficient
)

# --- 5. GENERATE THE FINAL COMPARISON PLOT ---

cat("Generating final comparison plot...\n")
mse_comparison_plot <- ggplot(plot_comparison_df, aes(x = Bandwidth, y = MSE, color = Method, linetype = Method)) +
  geom_line(linewidth = 1) +
  
  # Add vertical lines for optimal hx values, without showing them in the legend
  geom_vline(aes(xintercept = optimal_values[["Naive (ksmooth)"]]), color = color_palette["Naive (ksmooth)"], linetype = linetype_palette["Naive (ksmooth)"], linewidth = 0.8, show.legend = FALSE) +
  geom_vline(aes(xintercept = optimal_values[["Manual Naive (by hand)"]]), color = color_palette["Manual Naive (by hand)"], linetype = linetype_palette["Manual Naive (by hand)"], linewidth = 0.8, show.legend = FALSE) +
  geom_vline(aes(xintercept = optimal_values[["Efficient (by hand)"]]), color = color_palette["Efficient (by hand)"], linetype = linetype_palette["Efficient (by hand)"], linewidth = 0.8, show.legend = FALSE) +

  # Add text annotations for each optimal hx
  annotate("text", x = optimal_values[["Naive (ksmooth)"]], y = max(plot_comparison_df$MSE, na.rm = TRUE) * 0.95, 
           label = paste("ksmooth opt=", format(round(optimal_values[["Naive (ksmooth)"]], 2), nsmall = 2)),
           color = color_palette["Naive (ksmooth)"], hjust = -0.05, vjust = 0, size = 3.0, fontface = "bold") +
  annotate("text", x = optimal_values[["Manual Naive (by hand)"]], y = max(plot_comparison_df$MSE, na.rm = TRUE) * 0.85, 
           label = paste("Man. Naive opt=", format(round(optimal_values[["Manual Naive (by hand)"]], 2), nsmall = 2)),
           color = color_palette["Manual Naive (by hand)"], hjust = -0.05, vjust = 0, size = 3.0, fontface = "bold") +
  annotate("text", x = optimal_values[["Efficient (by hand)"]], y = max(plot_comparison_df$MSE, na.rm = TRUE) * 0.75, 
           label = paste("Man. Eff. opt=", format(round(optimal_values[["Efficient (by hand)"]], 2), nsmall = 2)),
           color = color_palette["Efficient (by hand)"], hjust = -0.05, vjust = 0, size = 3.0, fontface = "bold") +
  
  # Apply manual scales for color and linetype
  scale_color_manual(values = color_palette, name = "Method:") +
  scale_linetype_manual(values = linetype_palette, name = "Method:") +
  
  # Set axis limits
  coord_cartesian(
    ylim = c(min(plot_comparison_df$MSE, na.rm = TRUE) * 0.95, max(plot_comparison_df$MSE, na.rm = TRUE) * 1.1),
    xlim = range(plot_comparison_df$Bandwidth, na.rm = TRUE), expand = TRUE
  ) +
  
  # Set labels and theme
  labs(title = "Comparison of LOOCV MSE: ksmooth vs. Manual Implementations", x = "Bandwidth (hx)", y = "LOOCV MSE") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.background = element_rect(fill = "#FBFAF7", colour = NA),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 1),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.8), colour = "grey70", linewidth = 0.5),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.2, "cm")
  )

# Save the plot
ggsave("plots/03_loocv_comparison_all_methods.png", plot = mse_comparison_plot, width = 9, height = 6, dpi = 300)