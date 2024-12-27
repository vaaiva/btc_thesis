################################################
## --- ROLLING BOOTSTRAPPED HMM FORECASTS --- ##
################################################
# This section compares the optimal MCMC sets 
# and two benchmark sets: 
covariates_benchmark1 <- c("CNYUSD", "DJI", "GC_F", "VOL")     # a mix of macroeconomic and BTC-specific variables
covariates_benchmark2 <- c("MINER","HASH", "VOL", "halving")   # only BTC-specific variables
# with, a common benchmark transition covariate:
trans_benchmark <- c("VIX")


# Optimal observed MCMC sets:
covariates4 <- c("VOL")                       # early sample - emission
trans1 <- c("EURUSD", "VIX")                  # early sample - transition
covariates5 <- c("CNYUSD", "DJI", "EURUSD")   # recent sample - emission
trans2 <- c("CL_F", "EURUSD")                 # recent sample - transition
covariates6 <- c("GBPUSD", "HASH")            # full sample - emission
trans3 <- c("VIX")                            # full sample - transition


###########################################
## -------------- MAIN FORECASTING FUNCTION

rolling_window_nh_hmm <- function(df, response_var, covariates, covariates_trans, states, window_size, forecast_horizon, n_bootstrap = 1000) {
  
  # Constants and storage
  n <- nrow(df)
  mae_results <- c()
  rmse_results <- c()
  mape_results <- c()
  coverage_prob_results <- c()
  all_mean_forecasts <- c()
  all_ci_lower <- c()
  all_ci_upper <- c()
  all_actual_values <- c()
  skipped_iterations <- c()
  
  # ----------- Main function
  for (i in 1:(n - window_size - forecast_horizon)) {
    print(paste("Iteration:", i))
    
    # Defining the rolling window, and the remaining testing window for each iteration
    df_train <- df[i:(i + window_size - 1), ]
    df_test <- df[(i + window_size):(i + window_size + forecast_horizon - 1), ]
    
    # Fitting the NH-HMM model on the training set
    fitted_model <- tryCatch({
      fit_nh_hmm(df_train, response_var, covariates, covariates_trans)
    }, error = function(e) {
      skipped_iterations <<- c(skipped_iterations, i) # if the model fails to converge, store that as a skipped iteration
      return(NULL)
    })
    
    # Skip if the model fails to converge, move on further
    if (is.null(fitted_model)) next
    
    
    # Extracting emission coefficients and standard deviations
    coefficients <- list()
    sigma <- list()
    for (state in 1:states) {
      response <- fitted_model@response[[state]][[1]]
      coefficients[[state]] <- response@parameters$coefficients
      sigma[[state]] <- response@parameters$sd
    }
    
    # Extracting transition coefficients
    betas <- list()
    for (state in 1:states) {
      betas[[state]] <- fitted_model@transition[[state]]@parameters$coefficients[, state]
    }
    
    # Preparing matrices for emissions and transition, adding intercept as =1
    emission_covariates_matrix <- cbind(1, as.matrix(df_test[, covariates]))
    transition_covariates_matrix <- cbind(1, as.matrix(df_test[, covariates_trans]))
    
    # Initializing from the last observed states
    posterior_states <- posterior(fitted_model)$state
    last_state <- tail(posterior_states, 1)
    
    # Bootstrapping forecasts
    bootstrap_forecasts <- matrix(NA, nrow = n_bootstrap, ncol = forecast_horizon)
    for (b in 1:n_bootstrap) {
      current_state <- last_state
      for (h in 1:forecast_horizon) {
        mu_t <- emission_covariates_matrix[h, ] %*% coefficients[[current_state]]
        bootstrap_forecasts[b, h] <- rnorm(1, mean = mu_t, sd = sigma[[current_state]])  # forecast
        
        # Calculating dynamic transition probabilities using the covariates, to initialize the next step ahead STATE
        eta <- transition_covariates_matrix[h, ] %*% betas[[current_state]]
        if (!is.finite(eta)) {
          next_state_prob <- c(0.5, 0.5)  # Assign equal probabilities if eta is invalid
        } else {
          next_state_prob <- pmax(pmin(exp(eta) / (1 + exp(eta)), 1 - 1e-6), 1e-6)  # 1e-6 is used to ensure non-zero values
        }
        
        next_state_prob <- if (current_state == 1) c(next_state_prob, 1 - next_state_prob) else c(1 - next_state_prob, next_state_prob)
        
        
        # Adding a fallback in case probabilities are still invalid, not due to eta not being finite
        if (any(!is.finite(next_state_prob)) || sum(next_state_prob) != 1) {
          next_state_prob <-  c(0.5, 0.5) # assigning equal probabilities in such case
        }
        # Updating the next state
        current_state <- sample(1:states, size = 1, prob = next_state_prob) 
      }
    }
    
    # Calculating mean forecasts and confidence intervals 
    mean_forecasts <- colMeans(bootstrap_forecasts, na.rm = TRUE)
    ci_lower <- apply(bootstrap_forecasts, 2, quantile, probs = 0.025, na.rm = TRUE)
    ci_upper <- apply(bootstrap_forecasts, 2, quantile, probs = 0.975, na.rm = TRUE)
    
    # Calculating accuracy metrics
    actual_values <- df_test$BTC_USD[1:forecast_horizon]
    mae <- mean(abs(actual_values - mean_forecasts))
    rmse <- sqrt(mean((actual_values - mean_forecasts)^2))
    mape <- mean(abs((actual_values - mean_forecasts) / actual_values)) * 100
    
    # Calculating coverage probability (percentage of actual values within the 95% CI)
    coverage_prob <- mean((actual_values >= ci_lower) & (actual_values <= ci_upper)) * 100
    
    # Storing everything
    mae_results <- c(mae_results, mae)
    rmse_results <- c(rmse_results, rmse)
    mape_results <- c(mape_results, mape)
    coverage_prob_results <- c(coverage_prob_results, coverage_prob)
    all_mean_forecasts <- c(all_mean_forecasts, mean_forecasts)
    all_ci_lower <- c(all_ci_lower, ci_lower)
    all_ci_upper <- c(all_ci_upper, ci_upper)
    all_actual_values <- c(all_actual_values, actual_values)
  }
  
  return(list(
    mae = mae_results,
    rmse = rmse_results,
    mape = mape_results,
    coverage_prob = coverage_prob_results,
    mean_forecasts = all_mean_forecasts,
    ci_lower = all_ci_lower,
    ci_upper = all_ci_upper,
    actual_values = all_actual_values,
    skipped_iterations = skipped_iterations
  ))
}


rolling_nh_hmm_results <- rolling_window_nh_hmm(
  df = df_scaled,                      #TODO update the subsample used
  response_var = "BTC_USD",          
  covariates = covariates_benchmark2,  #TODO update the emission covariates used
  covariates_trans = trans_benchmark,  #TODO update the transition covariates used
  states = 2,
  window_size = 100,                 
  forecast_horizon = 1,                #TODO update the forecast horizons (1-step, 5-step, 30-steps ahead) as needed
  n_bootstrap = 1000                   # number of bootstrap simulations
)

# Accuracy metrics
mean(rolling_nh_hmm_results$mae)
mean(rolling_nh_hmm_results$rmse)
mean(rolling_nh_hmm_results$mape)
mean(rolling_nh_hmm_results$coverage_prob)


#TODO 
# repeat the procedure for each subsample:
# 1. with the corresponding optimal MCMC set, and it's optimal transition covariate set
# 2. with bencmark Set 1 and Set 2, using benchmark transition set (VIX)
# 3. for 1, 5, and 30 steps ahead.

