
################################################
## --- ROLLING BOOTSTRAPPED HMM FORECASTS --- ##
################################################
# This section compares the optimal MCMC sets 
# and two benchmark sets: 
covariates_benchmark1 <- c("CNYUSD", "DJI", "GC_F", "VOL")     # a mix of macroeconomic and BTC-specific variables
covariates_benchmark2 <- c("MINER","HASH", "VOL", "halving")   # only BTC-specific variables

# Optimal observed MCMC sets:
covariates1 <- c("GC_F", "INFL", "halving")    # early sample
covariates2 <- c("CNYUSD", "GBPUSD", "HASH")   # recent sample
covariates3 <- c("CNYUSD", "INFL", "VOL")      # full sample



###########################################
## -------------- MAIN FORECASTING FUNCTION

rolling_window_hmm <- function(df, response_var, covariates, states, window_size, forecast_horizon, n_bootstrap) {
  
  # Constants and storage:
  n <- nrow(df)  
  mae_results <- c()
  rmse_results <- c()
  mape_results <- c()
  coverage_prob_results <- c()  
  all_mean_forecasts <- c()
  all_ci_lower <- c()
  all_ci_upper <- c()
  all_actual_values <- c()
  coefficient_series <- vector("list", states)     # storage for the rolled coefficients
  for (state in 1:states) {
    coefficient_series[[state]] <- data.frame(matrix(ncol = length(covariates) + 1, nrow = 0))
  }

  # ----------- Main function
  for (i in 1:(n - window_size - forecast_horizon)) {
    print(paste("Iteration:", i)) 
    
    # Defining the rolling window, and the remaining testing window for each iteration
    df_train <- df[i:(i + window_size - 1), ]
    df_test <- df[(i + window_size):(i + window_size + forecast_horizon - 1), ]
    
    # Fitting the HMM model on the training set
    fitted_model <- fit_hmm_c(df_train, response_var, covariates, states)
    
    # Initializing storage for the bootstrap forecasts
    bootstrap_forecasts <- matrix(NA, nrow = n_bootstrap, ncol = forecast_horizon)
    
    # Extracting emission coefficients and sd
    coefficients <- list()
    sigma <- list()
    for (state in 1:states) {
      response <- fitted_model@response[[state]][[1]]
      coefficients[[state]] <- response@parameters$coefficients
      sigma[[state]] <- response@parameters$sd
      coefficient_series[[state]] <- rbind(coefficient_series[[state]], coefficients[[state]])
    }
    
    # Initializing from the last observed state
    posterior_states <- posterior(fitted_model)$state
    last_state <- tail(posterior_states, 1)
    covariates_matrix <- as.matrix(df_train[, covariates])
    covariates_matrix <- cbind(1, covariates_matrix)  # adding intercept as =1
    
    # Bootstrapping the forecasts
    for (b in 1:n_bootstrap) {
      current_state <- last_state
      for (h in 1:forecast_horizon) {
        mu_t <- covariates_matrix[nrow(df_train), ] %*% coefficients[[current_state]]
        bootstrap_forecasts[b, h] <- rnorm(1, mean = mu_t, sd = sigma[[current_state]])  # forecast
        
        # Simulating the next state using the transition probabilities
        next_state_prob <- fitted_model@transition[[current_state]]@parameters$coefficients
        next_state <- sample(1:states, size = 1, prob = next_state_prob)
        current_state <- next_state
      }
    }
    
    # Calculating mean forecasts and confidence intervals 
    mean_forecasts <- colMeans(bootstrap_forecasts)
    ci_lower <- apply(bootstrap_forecasts, 2, quantile, probs = 0.025)
    ci_upper <- apply(bootstrap_forecasts, 2, quantile, probs = 0.975)
    
    
    # Calculating accuracy metrics
    actual_values <- df_test$BTC_USD[1:forecast_horizon] # actual values for comparison
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
  
  # Naming the rolled coefficients
  for (state in 1:states) {
    colnames(coefficient_series[[state]]) <- c("(Intercept)", covariates)
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
    coefficients_over_time = coefficient_series
  ))
}


rolling_hmm_c_results <- rolling_window_hmm(
  df = df_scaled,              #TODO update the subsample used
  response_var = "BTC_USD",          
  covariates = covariates_benchmark2,    #TODO update the list of covariates as needed          
  states = 2,                        
  window_size = 100,                 
  forecast_horizon = 1,        #TODO update the forecast horizons (1-step, 5-step, 30-steps ahead) as needed
  n_bootstrap = 1000           # number of bootstrap simulations
)

# Accuraccy metrics
mean(rolling_hmm_c_results$mae)
mean(rolling_hmm_c_results$rmse)
mean(rolling_hmm_c_results$mape)
mean(rolling_hmm_c_results$coverage_prob)

#TODO 
# repeat the procedure for each subsample:
# 1. with the corresponding optimal MCMC set
# 2. with bencmark Set 1 and Set 2
# 3. for 1, 5, and 30 steps ahead.


##### ---------------- PLOTS:

#----  1. Actual vs forecasted values
#TODO update the "hmm_final1" with the final model (with the corresponding list of covariates)
# this is for shading the State 1 periods in grey
hmm_final1 <- fit_hmm_c(d = df_scaled, response_var = "BTC_USD",  covariates = covariates_benchmark2, states = states)
state_sequence <- posterior(hmm_final1)$state  # the obtained state sequence is for grey shaded bars for State 1


plot_actual_vs_forecast <- function(df, actual_values, forecasted_values, ci_lower, ci_upper, state_sequence, start_point, forecast_horizon) {
  
  n_forecast <- length(forecasted_values)
  
  # Adjusting the length of forecast and actual data to fit the available date range
  df <- df[1:(start_point + n_forecast - 1), ]
  
  # Creating a data frame for the actual and forecasted part
  df_plot <- data.frame(
    date = df$date,
    actual = c(df$BTC_USD[1:(start_point - 1)], actual_values[1:n_forecast]),  # actual values
    forecast = c(rep(NA, start_point - 1), forecasted_values),                 # starting forecast after the window
    ci_lower = c(rep(NA, start_point - 1), ci_lower),                          # confidence intervals start after the window
    ci_upper = c(rep(NA, start_point - 1), ci_upper),
    state = c(state_sequence[1:(start_point + n_forecast - 1)]) 
  )
  
  # Detecting where the State 1 periods are
  state1_ranges <- df_plot %>%
    mutate(is_state1 = state == 1,
           block = cumsum(c(1, diff(is_state1)) != 0)) %>%  # creating blocks for when state changes
    filter(is_state1) %>%
    group_by(block) %>%
    summarise(start = min(date), end = max(date), .groups = "drop")
  
  # Plotting 
  plot <- ggplot(df_plot, aes(x = date)) +
    geom_rect(data = state1_ranges, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey", alpha = 0.2) +
    geom_line(aes(y = actual, color = "Actual"), size = 1) +  # Actual values
    geom_line(aes(y = forecast, color = "Forecast"), size = 0.8) +  # Forecasted values
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "blue", alpha = 0.2) +  # Confidence intervals
    labs(title = paste0("Actual vs Forecasted BTC Price with ", forecast_horizon, "-step Ahead Forecast"),
         y = "BTC Price", x = "Date") +
    scale_color_manual(values = c("Actual" = "black", "Forecast" = "lightblue")) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  print(plot)
}

plot_actual_vs_forecast(
  df = df_scaled,                                         #TODO adjust the subsample as needed
  actual_values = rolling_hmm_c_results$actual_values, 
  forecasted_values = rolling_hmm_c_results$mean_forecasts, 
  ci_lower = rolling_hmm_c_results$ci_lower, 
  ci_upper = rolling_hmm_c_results$ci_upper, 
  state_sequence = state_sequence,
  start_point = 100 + 1,                                 #TODO adjust the start of the forecast: 100 + forecast horizon
  forecast_horizon = 1                                   #TODO adjust based forecasting horizon
)


#---- 2. Dynamic rolled coefficients by state
plot_coefficients_by_state <- function(rolling_hmm_results, dates, start_index) {
  
  # Extracting coefficients for each state
  coefficients_series <- rolling_hmm_results$coefficients_over_time
  states <- length(coefficients_series)
  
  # Adjusting dates to match the rolling window
  plot_dates <- dates[start_index:(start_index + nrow(coefficients_series[[1]]) - 1)]
  
  # Combining coefficients into one data frame, reshaping the data
  plot_data <- bind_rows(
    lapply(1:states, function(state) {
      data.frame(Date = plot_dates, coefficients_series[[state]], State = paste("State", state))
    })
  )
  plot_data_long <- plot_data %>%
    pivot_longer(cols = -c(Date, State), names_to = "Coefficient", values_to = "Value")
  
  # Plot 
  ggplot(plot_data_long, aes(x = Date, y = Value, color = Coefficient)) +
    geom_line(size = 1) +
    facet_wrap(~ State, scales = "free_y") +
    labs(title = "Coefficients Over Time by State",
         x = "Date", y = "Coefficient Value") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Coefficient"))
}

plot_coefficients_by_state(rolling_hmm_c_results, df_scaled$date, 102) 
#TODO update the "dates" as needed


# -------- 3. Dynamic MAPE over time
#TODO update the dates as needed 
rolling_hmm_c_results$coefficients_over_time

plot_data <- data.frame(
  Date = df_scaled$date[102:nrow(df_scaled)], #TODO update as needed
  MAPE = rolling_hmm_c_results$mape
)
# Plot:
ggplot(plot_data, aes(x = Date, y = MAPE)) +
  geom_line(color = "blue", size = 1) +  
  labs(
    x = "",  
    y = "MAPE, %",  
    title = "MAPE over time, 2016-06:2024-12"  # TODO update the title as needed
  ) +
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5),  
    axis.title.y = element_text(size = 12), 
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )

