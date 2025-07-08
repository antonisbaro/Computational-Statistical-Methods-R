################################################################################
# 04b_model_selection_lasso.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script applies the LASSO regression method for automatic variable
# selection on the Boston housing dataset. It uses k-fold cross-validation
# to find the optimal regularization parameter (lambda) and saves the
# resulting cross-validation plot.
#
# This script corresponds to Exercise 4(b) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

library(glmnet)
library(MASS) # For the Boston dataset
data(Boston)


# --- 2. PREPARE DATA FOR GLMNET ---

x_matrix_full <- model.matrix(medv ~ ., data = Boston)[, -1]
y_vector <- Boston$medv


# --- 3. PERFORM CROSS-VALIDATION FOR LASSO ---

set.seed(123)
cat("Performing 10-fold cross-validation for LASSO...\n")
cv_lasso_fit <- cv.glmnet(x = x_matrix_full, y = y_vector, alpha = 1, nfolds = 10)
cat("Cross-validation finished.\n\n")


# --- 4. VISUALIZE CV RESULTS AND SAVE PLOT ---

# Define the output path for the plot
plot_path <- "plots/08_lasso_cv_plot.png"
cat(sprintf("Saving LASSO cross-validation plot to: %s\n", plot_path))

# Open a PNG device
png(plot_path, width = 8, height = 6, units = "in", res = 300)

# Create the plot
plot(cv_lasso_fit)
title("LASSO Cross-Validation: MSE vs. log(Lambda)", line = 2.5, cex.main = 1.2)

# Close the device, which saves the file
dev.off()

# Also display the plot in the R session (optional)
# plot(cv_lasso_fit)
# title("LASSO Cross-Validation: MSE vs. log(Lambda)", line = 2.5, cex.main = 1.2)


# --- 5. EXTRACT LAMBDAS AND ANALYZE MODEL ---

lambda_min <- cv_lasso_fit$lambda.min
lambda_1se <- cv_lasso_fit$lambda.1se

cat("--- LASSO Cross-Validation Results ---\n")
cat(sprintf("Lambda min (gives minimum CV-MSE): %.8f\n", lambda_min))
cat(sprintf("Lambda 1se (more parsimonious model): %.8f\n\n", lambda_1se))

# Get coefficients for the lambda.1se model
coefficients_lambda_1se <- coef(cv_lasso_fit, s = "lambda.1se")

cat("--- Analysis of LASSO Model with lambda.1se ---\n")
cat("Coefficients:\n")
print(round(coefficients_lambda_1se, 6))

# Identify selected predictors
selected_predictors_lasso_indices <- which(coefficients_lambda_1se != 0)
selected_predictors_lasso <- rownames(coefficients_lambda_1se)[selected_predictors_lasso_indices]
selected_predictors_lasso <- setdiff(selected_predictors_lasso, "(Intercept)")

cat("\nPredictors selected by LASSO (lambda.1se):\n")
cat(paste(selected_predictors_lasso, collapse = ", "), "\n")
cat(sprintf("Number of predictors selected (excluding intercept): %d\n", length(selected_predictors_lasso)))       