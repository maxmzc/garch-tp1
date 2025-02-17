library(PerformanceAnalytics)

load("indices.rda")  # loads "prices"

# Extract prices since January 2005, compute log returns, and drop first line (NA)
SP500 <- prices[,"SP500"]["2005-01-01/"] 
FTSE  <- prices[,"FTSE100"]["2005-01-01/"] 
rets_SP500 <- CalculateReturns(SP500, method = "log")[-1,] 
rets_FTSE  <- CalculateReturns(FTSE, method = "log")[-1,] 

f_forecast_var <- function(y, level) {
  ### Compute the VaR forecast of a GARCH(1,1) model with Normal errors at the desired risk level
  #  INPUTS:
  #   y     : vector (T x 1) of observations (log-returns)
  #   level : scalar risk level (e.g. 0.95 for a VaR at the 95% risk level)
  #  OUTPUTS:
  #   VaR   : scalar VaR forecast 
  #   sig2  : vector (T+1 x 1) conditional variances
  #   theta : vector GARCH parameters
  
  # Starting values and bounds 
  theta0 <- c(0.1 * var(y), 0.1, 0.8)  # (omega, alpha, beta)
  LB     <- c(1e-6, 1e-6, 1e-6)            # lower bounds: omega >= 0, alpha > 0, beta > 0
  
  # Stationarity condition: alpha + beta <= 0.999
  A <- c(0, 1, 1)  # selects alpha and beta
  b <- 0.999
  ui <- matrix(A, nrow = 1)
  ci = -b
  
  # Run the optimization using optim (L-BFGS-B) with a penalty for violating stationarity
  opt_result <- constrOptim(theta = theta0,
                            f = function(theta) f_nll(theta, y),
                            grad = NULL, # Numerical gradient will be computed automatically
                            ui = ui,
                            ci = ci,
                            control = list())
  
  print(opt_result)  # Print the optimization result
  
  theta <- opt_result$par
  
  # Recompute the conditional variance
  sig2 <- f_ht(theta, y)
  
  # Compute the next-day forecast volatility (last element of sig2)
  sigma_next <- sqrt(sig2[length(sig2)])
  
  # Compute the VaR using the specified risk level.
  VaR <- qnorm(level, mean = 0, sd = sigma_next)
  
  out <- list(VaR_Forecast = VaR, 
              ConditionalVariances = sig2, 
              GARCH_param = theta)
  
  out
}

f_nll <- function(theta, y, A = c(0, 1, 1), b = 0.999) {
  T <- length(y)
  sig2 <- f_ht(theta, y)[1:T]
  
  # Avoid non-positive variances
  if (any(sig2 <= 0)) return(1e20)
  
  ll <- sum(dnorm(y, mean = 0, sd = sqrt(sig2), log = TRUE))
  
  nll <- -ll
  nll
}


f_ht <- function(theta, y) {
  ### Computes the vector of conditional variances for a GARCH(1,1) model
  #  INPUTS:
  #   theta: vector (omega, alpha, beta)
  #   y    : vector (T x 1) of log-returns
  #  OUTPUT:
  #   sig2 : vector (T+1 x 1) of conditional variances
  
  a0 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  
  T <- length(y)
  sig2 <- rep(NA, T + 1)
  
  # Initialize with unconditional variance
  sig2[1] <- a0 / (1 - a1 - b1)
  
  # Recursively compute the conditional variances
  for (t in 1:T) {
    sig2[t + 1] <- a0 + a1 * y[t]^2 + b1 * sig2[t]
  }
  
  sig2
}


# ----- Rolling Backtest and Plotting -----

# Settings for the rolling backtest
window_size <- 1000  # Rolling window of 1000 days
n_forecast <- 1000   # Number of forecasts (days) after the window

# Convert the returns to numeric vectors
rets_SP500_num <- as.numeric(rets_SP500)
rets_FTSE_num  <- as.numeric(rets_FTSE)

# Preallocate vectors to store VaR forecasts
VaR_SP500_series <- rep(NA, n_forecast)
VaR_FTSE_series  <- rep(NA, n_forecast)

# Loop over the forecasting period
for (i in 1:n_forecast) {
  # Define the rolling window indices
  window_idx <- i:(window_size + i - 1)
  
  # Forecast VaR for SP500
  forecast_SP500 <- f_forecast_var(rets_SP500_num[window_idx], 0.95)
  VaR_SP500_series[i] <- forecast_SP500$VaR_Forecast
  
  # Forecast VaR for FTSE100
  forecast_FTSE <- f_forecast_var(rets_FTSE_num[window_idx], 0.95)
  VaR_FTSE_series[i] <- forecast_FTSE$VaR_Forecast
}

# Realized returns following the rolling window period
realized_SP500 <- rets_SP500_num[(window_size + 1):(window_size + n_forecast)]
realized_FTSE  <- rets_FTSE_num[(window_size + 1):(window_size + n_forecast)]
time_index <- 1:n_forecast

# Save the plot to a PNG file
png("VaR_backtest.png", width = 800, height = 600)
par(mfrow = c(2,1), mar = c(4,4,2,1))

# Plot for SP500
plot(time_index, realized_SP500, type = "l", col = "black", lwd = 2,
     ylab = "Returns", xlab = "Time", main = "SP500: Realized Returns & VaR Forecast (95%)")
lines(time_index, VaR_SP500_series, col = "red", lwd = 2)
legend("bottomleft", legend = c("Realized Returns", "VaR (95%)"), 
       col = c("black", "red"), lwd = 2, bty = "n")

# Plot for FTSE100
plot(time_index, realized_FTSE, type = "l", col = "black", lwd = 2,
     ylab = "Returns", xlab = "Time", main = "FTSE100: Realized Returns & VaR Forecast (95%)")
lines(time_index, VaR_FTSE_series, col = "red", lwd = 2)
legend("bottomleft", legend = c("Realized Returns", "VaR (95%)"), 
       col = c("black", "red"), lwd = 2, bty = "n")

dev.off()
cat("The backtest plots have been saved to 'VaR_backtest.png'.\n")

backtest_results <- list(
  window_size         = window_size,
  n_forecast          = n_forecast,
  VaR_SP500_series    = VaR_SP500_series,
  VaR_FTSE_series     = VaR_FTSE_series,
  realized_SP500      = realized_SP500,
  realized_FTSE       = realized_FTSE
)

save(backtest_results, file = "backtest_results.rda")
cat("The backtest results have been saved to 'backtest_results.rda'.\n")