rm(list = ls())
while (!is.null(dev.list())) dev.off()

library(multDM)
library(sandwich)
library(R.matlab)

# dpmsvSVo2022_12.mat
# dpmsv16992022_12.mat

data <- readMat("dpmsvSVo2022_12.mat")
data1 <- readMat("dpmsv16992022_12.mat")

# Extract forecasts, true values, and CRPS from the datasets
dpmsv_y <- data$dpmsv.y
dpmsv_yhat <- data$dpmsv.yhat
dpmsv_crps <- data$dpmsv.crps

# SVo_yhat <- data$SVo.yhat
# SVo_crps <- data$SVo.crps

SVo_yhat <- data1$dpmsv.yhat1
SVo_crps <- data1$dpmsv.crps1

t1 <- 423  # Start of evaluation window
t2 <- 434  # End of evaluation window

nn <- 16 # choosing variable

# Install necessary packages

# Initialize vectors to store DMW test statistics and p-values (length 4 each)
dmw_statistics_point <- numeric(4)  # Store test statistics for point forecasts
dmw_p_values_point <- numeric(4)    # Store p-values for point forecasts

dmw_statistics_crps <- numeric(4)   # Store test statistics for CRPS
dmw_p_values_crps <- numeric(4)     # Store p-values for CRPS

# Define the forecast horizons of interest
horizons <- c(1, 3, 12, 24)

# Loop over the forecast horizons
for (i in seq_along(horizons)) {
  h <- horizons[i]  # Get the current forecast horizon
  
  # === Point Forecasts ===
  # Extract forecasts and true values for the current horizon and evaluation window
  y_true <- dpmsv_y[nn, h, t1:t2]
  yhat1 <- dpmsv_yhat[nn, h, t1:t2]
  yhat2 <- SVo_yhat[nn, h, t1:t2]
  
  # Filter valid cases (remove NaNs)
  valid_idx <- complete.cases(y_true, yhat1, yhat2)
  y_true_vec <- y_true[valid_idx]
  yhat1_vec <- yhat1[valid_idx]
  yhat2_vec <- yhat2[valid_idx]
  
  if (length(y_true_vec) > 0) {
    # Compute squared errors for both models
    errors1 <- (yhat1_vec - y_true_vec)^2
    errors2 <- (yhat2_vec - y_true_vec)^2
    
    # Compute the loss differential
    loss_diff <- errors1 - errors2
    
    # Fit a constant model for the loss differential
    model <- lm(loss_diff ~ 1)
    
    # Compute Newey-West standard error with lag = h + 1
    nw_cov <- NeweyWest(model, lag = h + 1, prewhite = FALSE)  # Get the variance-covariance matrix
    nw_se <- sqrt(nw_cov[1, 1])  # Extract the standard error from the [1,1] element
    
    
    # Calculate the DMW test statistic
    dmw_stat <- mean(loss_diff) / nw_se
    
    # Calculate the p-value
    p_value <- 2 * (1 - pnorm(abs(dmw_stat)))
    
    # Store the results
    dmw_statistics_point[i] <- dmw_stat
    dmw_p_values_point[i] <- p_value
  } else {
    # Store NA if no valid data
    dmw_statistics_point[i] <- NA
    dmw_p_values_point[i] <- NA
  }
  
  # === CRPS Forecasts ===
  # Extract CRPS values for both models for the current horizon
  crps1 <- dpmsv_crps[nn, h, t1:t2]
  crps2 <- SVo_crps[nn, h, t1:t2]
  
  # Filter valid cases (remove NaNs)
  valid_idx_crps <- complete.cases(crps1, crps2)
  crps1_vec <- crps1[valid_idx_crps]
  crps2_vec <- crps2[valid_idx_crps]
  
  if (length(crps1_vec) > 0) {
    # Compute the loss differential for CRPS
    loss_diff_crps <- crps1_vec - crps2_vec
    
    # Fit a constant model for the loss differential
    model_crps <- lm(loss_diff_crps ~ 1)
    
    # Compute Newey-West standard error with lag = h + 1
    nw_cov <- NeweyWest(model_crps, lag = h + 1, prewhite = FALSE)  # Get variance-covariance matrix
    nw_se_crps <- sqrt(nw_cov[1, 1])  # Extract the standard error
    
    # Calculate the DMW test statistic
    dmw_stat_crps <- mean(loss_diff_crps) / nw_se_crps
    
    # Calculate the p-value
    p_value_crps <- 2 * (1 - pnorm(abs(dmw_stat_crps)))
    
    # Store the results
    dmw_statistics_crps[i] <- dmw_stat_crps
    dmw_p_values_crps[i] <- p_value_crps
  } else {
    # Store NA if no valid data
    dmw_statistics_crps[i] <- NA
    dmw_p_values_crps[i] <- NA
  }
}

# Print the results for point forecasts
cat("DMW Test Statistics (Point Forecasts):\n", dmw_statistics_point, "\n")
cat("DMW Test P-values (Point Forecasts):\n", dmw_p_values_point, "\n")

# Print the results for CRPS forecasts
cat("DMW Test Statistics (CRPS):\n", dmw_statistics_crps, "\n")
cat("DMW Test P-values (CRPS):\n", dmw_p_values_crps, "\n")

# Function to assign asterisks based on p-value
get_significance <- function(p) {
  if (is.na(p)) {
    return("")  # Handle NA case
  } else if (p < 0.001) {
    return("***")  # Highly significant
  } else if (p < 0.01) {
    return("**")  # Strongly significant
  } else if (p < 0.05) {
    return("*")  # Statistically significant
  } else if (p < 0.1) {
    return(".")  # Marginally significant (optional)
  } else {
    return("")  # Not significant
  }
}

# Print results with asterisks for Point Forecasts
cat("DMW Test Results (Point Forecasts):\n")
for (i in seq_along(dmw_statistics_point)) {
  stat <- dmw_statistics_point[i]
  p <- dmw_p_values_point[i]
  significance <- get_significance(p)
  cat(sprintf("Horizon %d: Statistic = %.3f, p-value = %.3f%s\n", 
              horizons[i], stat, p, significance))
}

# Print results with asterisks for CRPS Forecasts
cat("\nDMW Test Results (CRPS Forecasts):\n")
for (i in seq_along(dmw_statistics_crps)) {
  stat <- dmw_statistics_crps[i]
  p <- dmw_p_values_crps[i]
  significance <- get_significance(p)
  cat(sprintf("Horizon %d: Statistic = %.3f, p-value = %.3f%s\n", 
              horizons[i], stat, p, significance))
}

