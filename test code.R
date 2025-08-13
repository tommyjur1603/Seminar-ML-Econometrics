
# test

# Install if not already
# install.packages("grf")
# install.packages("ggplot2")
# install.packages("dplyr")

library(grf)
library(ggplot2)
library(dplyr)

set.seed(123)

# --- Simulation function ---
run_simulation <- function(n = 2000, d = 3, num.trees = 2000, honesty = TRUE) {
  
  # Covariates: first is signal, rest are noise
  X <- matrix(runif(n * d), n, d)
  
  # True CATE depends only on X1
  tau_true <- 2 * (X[,1] - 0.5)
  
  # Baseline outcome f0
  f0 <- sin(2 * pi * X[,1]) + X[,2]  # depends a bit on X2 to make baseline nontrivial
  
  # Propensity score for overlap
  ps <- plogis(0.5 * X[,1] - 0.25 * X[,2] + 0.1 * rowMeans(X))
  
  # Treatment assignment
  W <- rbinom(n, 1, ps)
  
  # Outcomes
  Y <- f0 + tau_true * W + rnorm(n)
  
  # Fit causal forest
  cf <- causal_forest(X, Y, W, honesty = honesty, num.trees = num.trees)
  
  # Predictions with variance estimates
  pred <- predict(cf, estimate.variance = TRUE)
  tau_hat <- pred$predictions
  se_hat <- sqrt(pred$variance.estimates)
  
  # CI bounds
  ci_lower <- tau_hat - 1.96 * se_hat
  ci_upper <- tau_hat + 1.96 * se_hat
  
  # Metrics
  rmse <- sqrt(mean((tau_hat - tau_true)^2))
  bias <- mean(tau_hat - tau_true)
  coverage <- mean(tau_true >= ci_lower & tau_true <= ci_upper)
  
  list(
    rmse = rmse,
    bias = bias,
    coverage = coverage,
    df = data.frame(
      tau_true = tau_true,
      tau_hat = tau_hat,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      X1 = X[,1],
      d = d
    )
  )
}

# --- Dimensions to test ---
dims <- c(3, 6, 12, 24)
results <- list()

# Run simulation for each dimension
for (d in dims) {
  cat("Running d =", d, "\n")
  results[[as.character(d)]] <- run_simulation(d = d)
}

# --- Summary table ---
summary_table <- data.frame(
  d = dims,
  rmse = sapply(results, function(res) res$rmse),
  bias = sapply(results, function(res) res$bias),
  coverage = sapply(results, function(res) res$coverage)
)

print(summary_table)

# --- Plot 1: RMSE and Coverage vs d ---
df_metrics <- summary_table %>%
  tidyr::pivot_longer(cols = c(rmse, coverage), names_to = "metric", values_to = "value")

ggplot(df_metrics, aes(x = d, y = value, color = metric)) +
  geom_line() +
  geom_point(size = 2) +
  labs(title = "Effect of Dimensionality on RMSE and Coverage",
       x = "Dimension (d)", y = "Metric Value",
       color = "Metric") +
  theme_minimal()

# --- Plot 2: Distribution of estimates for each d ---
df_all <- do.call(rbind, lapply(results, function(res) res$df))

ggplot(df_all, aes(x = tau_hat, fill = as.factor(d))) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 20) +
  geom_vline(aes(xintercept = mean(tau_true)), color = "black", linetype = "dashed") +
  labs(title = "Distribution of Estimated CATEs by Dimension",
       x = "Estimated Treatment Effect", fill = "Dimension") +
  theme_minimal()

# --- Plot 3: True vs Estimated for each d ---
ggplot(df_all, aes(x = tau_true, y = tau_hat, color = as.factor(d))) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "True vs Estimated CATE by Dimension",
       x = "True CATE", y = "Estimated CATE", color = "Dimension") +
  theme_minimal()

