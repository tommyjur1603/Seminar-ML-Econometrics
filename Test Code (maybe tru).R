
library(grf)
library(ggplot2)
library(dplyr)

set.seed(123)

# --- Simulation function ---
run_simulation <- function(n = 2000, d = 3, num.trees = 2000, honesty = TRUE) {
  
  X <- matrix(runif(n * d), n, d)  # Covariates: first is signal, rest noise
  
  tau_true <- 2 * (X[,1] - 0.5)  # True CATE depends only on X1
  
  f0 <- sin(2 * pi * X[,1]) + X[,2]  # Baseline
  
  ps <- plogis(0.5 * X[,1] - 0.25 * X[,2] + 0.1 * rowMeans(X))
  
  W <- rbinom(n, 1, ps)
  
  Y <- f0 + tau_true * W + rnorm(n)
  
  cf <- causal_forest(X, Y, W, honesty = honesty, num.trees = num.trees)
  
  pred <- predict(cf, estimate.variance = TRUE)
  tau_hat <- pred$predictions
  se_hat <- sqrt(pred$variance.estimates)
  
  ci_lower <- tau_hat - 1.96 * se_hat
  ci_upper <- tau_hat + 1.96 * se_hat
  
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

# --- Higher dimensions ---
dims <- c(3, 10, 50, 100, 200)
results <- list()

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
  scale_x_log10(breaks = dims) +  # log scale to handle big jumps in d
  labs(title = "Effect of Dimensionality on RMSE and Coverage",
       x = "Dimension (d)", y = "Metric Value",
       color = "Metric") +
  theme_minimal()

# --- Plot 2: Histogram of estimates ---
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


## MAYBE BETTER PLOTS?

# --- Improved Histogram: Faceted view ---
ggplot(df_all, aes(x = tau_hat)) +
  geom_histogram(fill = "#2C7BB6", color = "white", bins = 20, alpha = 0.7) +
  geom_vline(aes(xintercept = mean(tau_true)), color = "red", linetype = "dashed") +
  facet_wrap(~ d, scales = "free_y") +
  labs(title = "Distribution of Estimated CATEs by Dimension",
       subtitle = "Dashed red line = mean of true CATE",
       x = "Estimated Treatment Effect",
       y = "Count") +
  theme_minimal(base_size = 14)

# --- Optional: Histogram with both true and estimated (facets) ---
# Create a long-format dataset for true and estimated
df_hist_long <- data.frame(
  value = c(df_all$tau_true, df_all$tau_hat),
  type = rep(c("True", "Estimated"), each = nrow(df_all)),
  d = rep(df_all$d, 2)
)

ggplot(df_hist_long, aes(x = value, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20, color = "white") +
  scale_fill_manual(values = c("True" = "#D7191C", "Estimated" = "#2C7BB6")) +
  facet_wrap(~ d, scales = "free_y") +
  labs(title = "True vs Estimated CATE Distributions by Dimension",
       x = "Treatment Effect", y = "Count") +
  theme_minimal(base_size = 14)

# --- Improved Scatter: Faceted with smoothing ---
ggplot(df_all, aes(x = tau_true, y = tau_hat)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#2C7BB6") +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ d) +
  labs(title = "True vs Estimated CATE by Dimension",
       subtitle = "Red line = fitted slope; Black dashed = 45° line",
       x = "True CATE", y = "Estimated CATE") +
  theme_minimal(base_size = 14)
