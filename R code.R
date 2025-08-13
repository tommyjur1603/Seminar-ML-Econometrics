
library(grf)
library(ggplot2)

set.seed(123)

# --- Simulation parameters ---
n <- 1000
p <- 3

# Covariates
X <- matrix(runif(n * p), n, p)

# True heterogeneous treatment effect
tau_true <- 2 * (X[,1] - 0.5)  # ranges from -1 to 1

# Baseline function f0
f0 <- X[,2] + sin(2 * pi * X[,3])

# Propensity score model for overlap
ps <- plogis(0.5 * X[,1] - 0.25 * X[,2] + 0.1 * X[,3])

# Treatment assignment
W <- rbinom(n, 1, ps)

# Outcomes
epsilon <- rnorm(n)
Y <- f0 + tau_true * W + epsilon

# --- Fit causal forest ---
cf <- causal_forest(X, Y, W, num.trees = 2000, honesty = TRUE)

# Predictions and CIs
pred <- predict(cf, estimate.variance = TRUE)
tau_hat <- pred$predictions
se_hat <- sqrt(pred$variance.estimates)

# 95% confidence intervals
ci_lower <- tau_hat - 1.96 * se_hat
ci_upper <- tau_hat + 1.96 * se_hat

# --- Evaluation metrics ---
rmse <- sqrt(mean((tau_hat - tau_true)^2))
coverage <- mean(tau_true >= ci_lower & tau_true <= ci_upper)
cat("RMSE:", rmse, "\n")
cat("Coverage:", coverage, "\n")

# --- Visualization 1: True vs Estimated over X1 ---
df_plot <- data.frame(
  X1 = X[,1],
  tau_true = tau_true,
  tau_hat = tau_hat,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

ggplot(df_plot, aes(x = X1)) +
  geom_point(aes(y = tau_true), color = "black", alpha = 0.5, size = 1) +
  geom_line(aes(y = tau_hat), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "blue") +
  labs(y = "Treatment Effect", title = "True vs Estimated CATE",
       subtitle = "Black: True, Blue: Estimated with 95% CI") +
  theme_minimal()

# --- Visualization 2: Histogram of treatment effects ---
df_hist <- data.frame(
  tau_true = tau_true,
  tau_hat = tau_hat
)

ggplot(df_hist) +
  geom_histogram(aes(x = tau_true), fill = "black", alpha = 0.4, bins = 20) +
  geom_histogram(aes(x = tau_hat), fill = "blue", alpha = 0.4, bins = 20) +
  labs(x = "Treatment Effect",
       y = "Count",
       title = "Distribution of True vs Estimated CATE",
       subtitle = "Black = True, Blue = Estimated") +
  theme_minimal()
