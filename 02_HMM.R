
# ---------------------------------
# -- HOMOGENEOUS HMM
# ---------------------------------
# this section fits the homogeneous version of the HMM with covariates
# and extracts the corresponding output parameters from the depmixs4 summary



############################################## ---------- MAIN FUNCTION (with covariates)

fit_hmm_c <- function(df, response_var, covariates, states) {
  
  covariates<- paste(covariates, collapse = " + ")   
  
  model <- depmix(
    response = as.formula(paste(response_var, "~", covariates)),  
    data = df,
    nstates = states,  # Number of hidden states
    family = gaussian()  # Using Gaussian distribution for emissions
  )
  
  fitted_model <- fit(model)
  
  print(summary(fitted_model))
  return(fitted_model)
}
##############################################

# Setting constants:
#TODO - update the number of states if needed
states <- 2 
# Full list of covariates:
covariates0 <- c("EURUSD", "GBPUSD", "JPYUSD", "CNYUSD", "SPX", "DJI", "IXIC", "CL_F", "GC_F", "VIX", "BLOCK", "HASH", "INFL", "MINER", "VOL", "halving")

##### Fitting the initial models:
hmm_c_0_1 <- fit_hmm_c(d = df_scaled, response_var = "BTC_USD",  covariates = covariates0, states = states)
hmm_c_0_2 <- fit_hmm_c(d = df_scaled_16_19, response_var = "BTC_USD",  covariates = covariates0, states = states)
hmm_c_0_3 <- fit_hmm_c(d = df_scaled_19_24, response_var = "BTC_USD",  covariates = covariates0, states = states)

# Convergence, log likelihood, AIC, BIC
hmm_c_0_1
hmm_c_0_2
hmm_c_0_3


#adjustable for any quick tryouts such as:
#covariates7 <- c("CNYUSD", "VIX", "HASH", "VOL", "halving")
#covariates6 <- c("MINER")
#covariates8 <- c("HASH", "VOL", "halving")
#covariates9 <- c("MINER","HASH", "VOL", "halving" )
covariates10 <- c("MINER","HASH",  "halving")
#covariates11 <- c("CNYUSD", "VIX","MINER")
hmm_c_x_3 <- fit_hmm_c(d = df_scaled_19_24, response_var = "BTC_USD",  covariates = covariates10, states = states)

# posterior state sequence:
state_sequence <- posterior(hmm_c_0_1)$state
state_runs <- rle(state_sequence)                                        # lengths of consecutive runs for each state
average_duration <- tapply(state_runs$lengths, state_runs$values, mean)  #  average duration for each state



#########################################################################
#
#-------------------------------------- EXTRACTING PRIORS
#
# For the next procedure, which will be MCMC, we need to extract
# resulting parameters, since they will be used as priors:


extract_hmm_c_priors <- function(fitted_model, df, response_var, covariates, num_states) {

  y <- as.numeric(df[[response_var]])
  
  # Initial state probabilities and state sequence
  init_probs <- fitted_model@prior@parameters$coefficients
  posterior_states <- posterior(fitted_model)$state
  
  # Storage for emission parameters
  trans_probs <- matrix(NA, nrow = num_states, ncol = num_states)   
  coefficients <- list()    # emission coefficients 
  sigma <- list()           # emission standard deviations 
  emission_probs <- matrix(NA, nrow = nrow(df), ncol = num_states)  # emission probabilities for each state
  var_prior <- list()      # coefficient variance
  
  # Preparing covariates
  covariates <- as.matrix(df[, covariates])               # matrix of covariates (X_t)
  covariates <- cbind(1, covariates)                      # adding intercept column 
  colnames(covariates)[1] <- "(Intercept)"
  

  for (state in 1:num_states) {
    # Transition probabilities
    trans_probs[state, ] <- fitted_model@transition[[state]]@parameters$coefficients
  
    # ---- From emission model:
    response <- fitted_model@response[[state]][[1]]   # specifying the model by the S4 object. Unable to extract the coefficients otherwise
    
    # Emission coefficients 
    coefficients[[state]] <- response@parameters$coefficients
    
    # Emission standard deviation (sigma)
    sigma[[state]] <- response@parameters$sd
    
    # Dynamic mean based on covariates
    mu_t <- covariates %*% coefficients[[state]]    
   
    # Emission probabilities 
    emission_probs[, state] <- dnorm(y, mean = mu_t, sd = sigma[[state]])
    
    var_prior[[state]] <- var(coefficients[[state]])   # coefficient variance
  }
  
  return(list(
    init_probs = init_probs,
    posterior_states = posterior_states,
    trans_probs = trans_probs,
    coefficients = coefficients,
    var_prior = var_prior,
    sigma = sigma,
    emission_probs = emission_probs,
    y = y,
    covariates = covariates,
    covariate_names = colnames(covariates)
  ))
}

hmm_c_priors <- extract_hmm_c_priors(fitted_model = hmm_c_0_1, df = df_scaled, 
                                     response_var = "BTC_USD", covariates = covariates0, 
                                     num_states = 2)
#TODO
# Adjust the "fitted_model" and "df" as needed, dependent on which prior parameters are needed
# respectively:
# hmm_c_0_1 - df_scaled
# hmm_c_0_2 - df_scaled_16_19
# hmm_c_0_3 - df_scaled_19_24