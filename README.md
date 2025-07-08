# Computational and Stochastic Optimization Methods in R

This repository contains a collection of projects developed for the "Computational Statistics and Stochastic Optimization" course at NTUA. It serves as a practical, in-depth exploration of fundamental statistical methods, implemented from scratch in **R**. The focus is on bridging the gap between statistical theory and practical, efficient code.

---

## 🚀 Core Algorithms & Concepts Explored

This project is structured around four key exercises, each tackling a different area of computational statistics:

1.  **Non-Parametric Regression**:
    -   Implementation of the **Nadaraya-Watson kernel estimator** with a Gaussian kernel.
    -   Selection of the optimal bandwidth (`h_x`) using **Leave-One-Out Cross-Validation (LOOCV)**.
    -   Developed both a naive (loop-based) and a mathematically **efficient (vectorized) LOOCV** implementation and compared their results, highlighting the importance of efficient computation.

2.  **Resampling and Simulation Methods**:
    -   **Bootstrap**: Implemented both non-parametric and parametric Bootstrap to estimate the sampling distribution of the `min(X)` statistic. This exercise critically evaluates the limitations of non-parametric Bootstrap for extreme value statistics.
    -   **Rejection Sampling**: Implemented the **Squeezed Rejection Sampling** algorithm to generate samples from a standard Normal distribution using a Laplace proposal distribution, demonstrating a powerful technique for efficient simulation.
    -   **Importance Sampling**: Used Importance Sampling to reduce the variance of a Monte Carlo estimator for a definite integral, showcasing a ~87% reduction in standard error compared to the classic Monte Carlo method.

3.  **Parameter Estimation with Latent Variables**:
    -   Implementation of the **Expectation-Maximization (EM) algorithm** from scratch to estimate the mixing probability (`p`) in a mixture of two exponential distributions.
    -   Visualized the convergence of the parameter estimate, demonstrating the iterative nature of the algorithm.

4.  **Model and Variable Selection**:
    -   Applied model selection techniques to the "Boston Housing" dataset.
    -   **Exhaustive Search**: Performed an exhaustive search over all 2¹³ possible models to find the optimal model based on the **Akaike Information Criterion (AIC)**.
    -   **LASSO Regression**: Utilized LASSO (via `glmnet`) for automatic variable selection and regularization, comparing the resulting parsimonious model with the one selected by AIC.
    -   **Residual Bootstrap**: Constructed a 95% confidence interval for a specific regression coefficient (`rm`) using the residual bootstrap method, providing a robust inference without assuming normality of errors.

---

## 🛠️ Key Skills Demonstrated

-   **Algorithmic Implementation**: Translating complex statistical formulas (LOOCV, EM, Bootstrap) into functional and efficient R code.
-   **Statistical Theory**: Deep understanding of the assumptions, advantages, and limitations of each method.
-   **Model Validation & Selection**: Mastery of techniques like Cross-Validation, AIC, and LASSO.
-   **Stochastic Simulation**: Practical application of Monte Carlo methods for simulation and integration.
-   **Data Visualization**: Using `ggplot2` to create insightful plots for model diagnostics, convergence analysis, and distribution visualization.

---

## 💻 Technology Stack

-   **Language**: R
-   **Core Packages**: `ggplot2` (for visualization), `glmnet` (for LASSO), `fitdistrplus` (for distribution fitting).

---

## ✍️ Author

-   Antonios Barotsakis
