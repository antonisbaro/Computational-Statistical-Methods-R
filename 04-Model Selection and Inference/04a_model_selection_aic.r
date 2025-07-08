################################################################################
# 04a_model_selection_aic.R
#
# Author: Antonios Barotsakis
#
# Description:
# This script performs an exhaustive search (best subset selection) to find
# the optimal linear regression model for predicting the 'medv' variable in
# the Boston housing dataset. It evaluates all 2^13 possible models and
# selects the one that minimizes the Akaike Information Criterion (AIC).
#
# The script:
# 1. Loads the Boston dataset from the MASS package.
# 2. Iterates through all possible combinations of the 13 predictors using
#    a binary representation.
# 3. For each combination, it fits a linear model (lm).
# 4. Calculates the AIC for the model.
# 5. Identifies the model with the minimum AIC and prints a detailed
#    summary of its formula, predictors, and performance.
#
# This script corresponds to Exercise 4(a) of the project.
#
################################################################################


# --- 1. SETUP & LIBRARIES ---

# The Boston dataset is in the MASS package
library(MASS)
data(Boston)


# --- 2. INITIALIZE VARIABLES ---

# Separate predictors (X) and response (y)
response_variable_name <- "medv"
predictor_names <- setdiff(colnames(Boston), response_variable_name)

n_total_predictors <- length(predictor_names)
n_models_to_evaluate <- 2^n_total_predictors

# Create a list to store detailed results for each model
model_results_list <- vector("list", n_models_to_evaluate)


# --- 3. EXHAUSTIVE MODEL SEARCH (BEST SUBSET SELECTION) ---

cat(sprintf("Starting exhaustive search for %d models...\n", n_models_to_evaluate))
# Loop through all possible combinations of predictors
# We use the integer i (from 0 to 2^p - 1) as a bitmask to select predictors.
for (i in 0:(n_models_to_evaluate - 1)) {
  # The model counter is i + 1
  model_counter <- i + 1
  
  # Convert i to a binary vector to select predictors
  binary_selection <- as.integer(intToBits(i))[1:n_total_predictors]
  selected_predictors <- predictor_names[binary_selection == 1]
  
  # Construct the model formula string
  if (length(selected_predictors) == 0) {
    # Model with intercept only
    current_formula_str <- paste(response_variable_name, "~ 1")
  } else {
    current_formula_str <- paste(response_variable_name, "~", paste(selected_predictors, collapse = " + "))
  }
  current_formula <- as.formula(current_formula_str)
  
  # Fit the linear model
  current_lm <- lm(current_formula, data = Boston)
  
  # Store all relevant information in the list
  model_results_list[[model_counter]] <- list(
    formula_str = current_formula_str,
    predictors = selected_predictors,
    num_predictors = length(selected_predictors),
    aic = AIC(current_lm),
    model_object = current_lm
  )
  
  # Optional: Print progress
  if (model_counter %% 500 == 0) {
    cat(sprintf("... Evaluated %d / %d models.\n", model_counter, n_models_to_evaluate))
  }
}
cat("Exhaustive search finished.\n")


# --- 4. IDENTIFY AND DISPLAY THE BEST MODEL ---

# Extract AIC values from the list of results
all_aics <- sapply(model_results_list, function(x) x$aic)

# Find the index of the model with the minimum AIC
min_aic_index <- which.min(all_aics)
best_model_info <- model_results_list[[min_aic_index]]

# Output a summary of the best model found
cat("\n--- Full Enumeration Results (AIC Criterion) ---\n")
cat(sprintf("Total models evaluated: %d\n", n_models_to_evaluate))
cat(sprintf("Best model formula based on AIC:\n%s\n", best_model_info$formula_str))
cat(sprintf("Minimum AIC value: %.3f\n", best_model_info$aic))
cat(sprintf("Number of predictors in the best model: %d\n", best_model_info$num_predictors))
cat(sprintf("Predictors in the best model: %s\n\n", paste(best_model_info$predictors, collapse = ", ")))

# Print the detailed summary of the best lm object
cat("--- Summary of the Best Model (M1) ---\n")
print(summary(best_model_info$model_object))   