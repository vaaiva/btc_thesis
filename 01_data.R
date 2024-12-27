
# ---------------------------------
# -- DATA PREPARATION
# ---------------------------------
#This part downloads, merges, and transorms the data employed for the analysis. 


################################## -----  Loading the data

tidyquant_df <- tq_get(c("BTC-USD", "ETH-USD", 
                          "EURUSD=X", "GBPUSD=X", "JPYUSD=X","CNYUSD=X",
                          "^SPX", "^DJI", "^IXIC",
                          "CL=F","GC=F",
                          "^VIX"),
                          from = "2016-06-01", to = "2024-12-02") %>% dplyr::select(symbol, date, close) %>% arrange(date) %>% spread(symbol, close) #selecting 2024-12-02 instead of 12-01 as the final day, since the 1st was a Sunday
tidyquant_df$date <- as.Date(tidyquant_df$date)

# Other variables:
#---- inflation
df_yields<- tq_get(c("DGS5","DFII5"), get = "economic.data", from = "2016-05-31") %>%
  arrange(date) %>% spread(symbol, price) %>% mutate(INFL = DGS5 - DFII5) %>% dplyr::select(date, INFL)
df_yields$date <- as.Date(df_yields$date)

#---- btc volume only
btc_volume_df <- tq_get("BTC-USD", from = "2016-06-01") %>% 
  dplyr::select(date, volume)
btc_volume_df$date <- as.Date(btc_volume_df$date)

#---- days untill halving
halving_dates <- as.Date(c("2016-07-09", "2020-05-11", "2024-04-19", "2028-03-26"))

days_until_halving <- function(date, halving_dates) {
  future_halving <- halving_dates[halving_dates >= date]
  if (length(future_halving) == 0) {
    return(NA)  # No future halving
  }
  return(as.numeric(future_halving[1] - date))
}
tidyquant_df$halving <- sapply(tidyquant_df$date, days_until_halving, halving_dates)

#-- Downloading bitcoin hash rate and block size data
# API
urls <- list(
  hash_rate = "https://api.blockchain.info/charts/hash-rate?timespan=all&format=json",
  block_size = "https://api.blockchain.info/charts/avg-block-size?timespan=all&format=json", 
  miners_rev = "https://api.blockchain.info/charts/miners-revenue?timespan=all&format=json"
)

data_list <- lapply(urls, function(url) {
  response <- GET(url)
  data <- fromJSON(content(response, "text"))
  

  data_values <- data$values
  data_frame <- data.frame(
    date = as.Date(as.POSIXct(data_values$x, origin = "1970-01-01", tz = "UTC")),  # converting to Date format is necessary for merging later
    value = data_values$y
  )
  return(data_frame)
})

# Merging
tidyquant_df <- merge(tidyquant_df, data_list[[1]], by = "date", all.x = TRUE)  # with hash rate [[1]]
tidyquant_df <- merge(tidyquant_df, data_list[[2]], by = "date", all.x = TRUE)  # with block size [[2]]
tidyquant_df <- merge(tidyquant_df, data_list[[3]], by = "date", all.x = TRUE)  # with miner renevue [[3]]
tidyquant_df <- merge(tidyquant_df, df_yields, by = "date", all.x = TRUE) 
tidyquant_df <- merge(tidyquant_df, btc_volume_df, by = "date", all.x = TRUE)



######################### ------ Reshaping, transformations
# Renaming the variables for better handling; lagging exchange rates by one due to the timing of observation in Yahoo Finance
df1 <- tidyquant_df %>%
  rename(
    BTC_USD = `BTC-USD`,
    ETH_USD = `ETH-USD`,
    EURUSD = `EURUSD=X`,
    GBPUSD = `GBPUSD=X`,
    JPYUSD = `JPYUSD=X`,
    CNYUSD = `CNYUSD=X`,
    SPX = `^SPX`,
    DJI = `^DJI`,
    IXIC = `^IXIC`,
    CL_F = `CL=F`,
    GC_F = `GC=F`,
    VIX = `^VIX`, 
    HASH = 'value.x', 
    BLOCK = 'value.y',
    MINER = 'value', 
    VOL = 'volume'
  ) %>%
  mutate(
    EURUSD = lag(EURUSD, 1),
    GBPUSD = lag(GBPUSD, 1),
    JPYUSD = lag(JPYUSD, 1),
    CNYUSD = lag(CNYUSD, 1)
  )


# Additionaly interpolating exchange rate data cause it is sometimes not observed, also filling the gaps for any remaining variables just in case there are observations missing
interpl_variables <- setdiff(names(df1), c("date", "BTC_USD", "ETH_USD", "halving"))

for (variable_name in interpl_variables) {
  
  ts_variable <- xts(df1[[variable_name]], order.by = as.Date(df1$date))
  interpolated_ts <- na.approx(ts_variable, rule = 2)
  df1[[variable_name]] <- coredata(interpolated_ts)
}

# Manual imputation: the price of crude oil futures was negative (-37) on 2020-04-20. I'm setting it to 1 as otherwise we can't take a log
df1$CL_F[1419:1421] <- 1

# Log transforming the entire dataset
df_log <- df1 %>%
  mutate(across(-c(date, halving), log))     %>%      #FIXME
  filter(row_number() <= n()-1) #THE LAST BIT CAN BE COMMENTED, it's a quick fix for when BTC is not reported as promptly

# Creating a week's lag
lags <- setdiff(names(df_log), c("date", "BTC_USD", "ETH_USD"))
df_log_lag <- df_log
for (variable_name in lags){
    df_log_lag[[variable_name]] <- lag(df_log_lag[[variable_name]], 7)
  }
# filtering out the first 7 days that are now missing due to lagging
df_log_lag <- df_log_lag %>% 
  filter(!is.na(IXIC)) 

# A few subsets of the data, according to date
df_log_16_19 <- df_log %>% filter (date < "2019-06-01")
df_log_19_24 <- df_log %>% filter (date >= "2019-06-01")  

df_log_lag_16_19 <- df_log_lag %>% filter (date < "2019-06-01") 
df_log_lag_19_24 <- df_log_lag %>% filter (date >= "2019-06-01")


# Scaling 
df_scaled <- df_log_lag %>% 
  mutate(across(-date, scale))
df_scaled_16_19 <- df_log_lag_16_19 %>%
  mutate(across(-date, scale))
df_scaled_19_24 <- df_log_lag_19_24 %>%
  mutate(across(-date, scale))

describe(df_scaled)  # a quick look into the summary of the data



##################### ------ PLOTS:

macro_variables <- c("EURUSD", "GBPUSD", "JPYUSD", "CNYUSD", "SPX", "DJI", "IXIC", "CL_F", "GC_F", "VIX", "INFL")
btc_variables <- c("BLOCK", "HASH", "MINER", "VOL", "halving")


# Reshaping data for ggplot
df_long_macro <- melt(df_scaled[, c("date", "BTC_USD", macro_variables)], id.vars = "date", 
                      variable.name = "Variable", value.name = "Value")

df_long_btc <- melt(df_scaled[, c("date", "BTC_USD", btc_variables)], id.vars = "date", 
                    variable.name = "Variable", value.name = "Value")

# Plot 1: BTC price vs macro variables 
ggplot(df_long_macro, aes(x = date)) +
  geom_line(aes(y = Value, color = Variable), alpha = 0.7) +
  geom_line(data = df_scaled, aes(y = BTC_USD, x = date), color = "black", size = 1) +
  labs(title = "Bitcoin Price vs. Macroeconomic Variables",
       x = "Date", y = "Value") +
  scale_color_manual(values = c("black", rainbow(length(macro_variables)))) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")

# Plot 2: BTC price vs Bitcoin-specific Variables
ggplot(df_long_btc, aes(x = date)) +
  geom_line(aes(y = Value, color = Variable), alpha = 0.7) +
  geom_line(data = df_scaled, aes(y = BTC_USD, x = date), color = "black", size = 1) +
  labs(title = "Bitcoin Price vs. Bitcoin-Specific Variables",
       x = "Date", y = "Value") +
  scale_color_manual(values = c("black", rainbow(length(btc_variables)))) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")

# Plot 3 (pdf): each variable vs BTC price (pdf!)
create_time_series_comparison_plots <- function(df, response_var, output_file = "btc_time_series_comparison_plots.pdf") {

  pdf(output_file, width = 10, height = 6)

  df_melted <- melt(df, id.vars = "date", variable.name = "variable", value.name = "value")
  
  variables <- setdiff(names(df), c("date", response_var))
  
  for (var in variables) {
    df_filtered <- df_melted[df_melted$variable %in% c(var, response_var), ]
    
    # the plot:
    p <- ggplot(df_filtered, aes(x = date, y = value, color = variable, group = variable)) +
      geom_line() +
      labs(title = paste("Time Series Comparison:", response_var, "vs", var),
           x = "Date", y = "Value") +
      theme_minimal()

    print(p)
  }
  dev.off()
}

create_time_series_comparison_plots(df_scaled, "BTC_USD", "btc_time_series_comparison_plots.pdf")  # using the full data

# Plot 4: correlation matrix of the original dataset
cor_matrix <- cor(df1[, c("EURUSD", "GBPUSD", "JPYUSD", "CNYUSD", "SPX", "DJI", "IXIC", "CL_F", "GC_F", "VIX", "BLOCK", "HASH", "INFL", "MINER", "VOL", "halving")], use = "complete.obs")
ggcorrplot(cor_matrix, method = "circle", hc.order = TRUE, type = "lower", 
           lab = TRUE, lab_size = 3, colors = c("red", "white", "blue"))
# high correlations: EUR, GBP, CNY and JPY ex rates
# SPX, IXIC, DJI (+ HASH rate and GOLD?)


#################### -------------- Alternative procedures for dimension reduction [reported, but not useful]:
# these are just tryouts

## ----------------- Lasso

# subsetting the data:
y_lasso <- df_scaled$BTC_USD               

covariates0 <- c("EURUSD", "GBPUSD", "JPYUSD", "CNYUSD", "SPX", "DJI", "IXIC", "CL_F", "GC_F", "VIX", "BLOCK", "HASH", "INFL", "MINER", "VOL", "halving")
X_lasso <- as.matrix(df_scaled[, covariates0])  


#  Lasso model:
lasso_model <- glmnet(X_lasso, y_lasso, alpha = 1)  # alpha = 1 specifies Lasso (L1 regularization)

# Cross-validation to find the best lambda (regularization strength)
cv_lasso <- cv.glmnet(X_lasso, y_lasso, alpha = 1)

# Choosing the best lambda, which selects the largest value such that the CV error
# is within 1 SE of the minimum. Because the one that minimizes CV, doesn't
# shrink any of the coefficients
best_lambda <- cv_lasso$lambda.1se                

# Fiting the Lasso model again using the best lambda
lasso_best <- glmnet(X_lasso, y_lasso, alpha = 1, lambda = best_lambda)
lasso_coefficients <- coef(lasso_best)

# Which are the non-zero coefficients? they are potential covariates
selected_covariates <- which(lasso_coefficients != 0)
colnames(X_lasso) #numbers of covariates: 1  2  3  4  5  6  7  8  9 11 12 13 14 15 16 17
# still not good enough, because of multicollinearity


# ------- Elastic Net (alpha = 0.5 gives equal weight to Lasso and Ridge penalties)
elastic_net_model <- glmnet(X_lasso, y_lasso, alpha = 0.5)

# Cross-validation to find the best lambda for Elastic Net
cv_elastic_net <- cv.glmnet(X_lasso, y_lasso, alpha = 0.5)
best_lambda_enet <- cv_elastic_net$lambda.min

elastic_net_best <- glmnet(X_lasso, y_lasso, alpha = 0.5, lambda = best_lambda_enet) # Elastic Net 
elastic_net_coefficients <- coef(elastic_net_best)

# The non-zero coefficients:
selected_covariates_enet <- which(elastic_net_coefficients != 0)
colnames(X_lasso) #numbers of covariates: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
# still not good enough, doesn't shrink any coefficients


#--------- Random Forest model
rf_model <- randomForest(X_lasso, y_lasso)

importance_rf <- importance(rf_model) #variable importance
sorted_importance <- sort(importance_rf[,1], decreasing = TRUE) # ranking covariates by importance

# Ranked importance of each covariate:
print(names(sorted_importance))


