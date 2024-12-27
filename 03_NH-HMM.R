
# ---------------------------------
# -- NON - HOMOGENEOUS HMM
# ---------------------------------
# this section fits the non-homogeneous version of the HMM with and without covariates
# and extracts the corresponding output parameters from the depmixs4 summary


############################################## ---------- MAIN FUNCTION

fit_nh_hmm <- function(df, response_var, covariates, covariates_trans) {
  
  covariates<- paste(covariates, collapse = " + ")
  covariates_trans <- paste(covariates_trans, collapse = "+")
  
  model <- depmix(
    response = as.formula(paste(response_var, "~", covariates)),  
    data = df,
    nstates = 2,  
    transition = as.formula(paste("~", covariates_trans)),  
    family = gaussian()  
  )
  
  fitted_model <- fit(model)
  
  print(summary(fitted_model))
  return(fitted_model)
}
################################################## 



# Full list of covariates:
covariates0 <- c("EURUSD", "GBPUSD", "JPYUSD", "CNYUSD", "SPX", "DJI", "IXIC", "CL_F", "GC_F", "VIX", "BLOCK", "HASH", "INFL", "MINER", "VOL", "halving")

##### Fitting the initial models:
nh_hmm0_1 <- fit_nh_hmm(d = df_scaled, response_var = "BTC_USD",  covariates = covariates0, covariates_trans = covariates0)
nh_hmm0_2 <- fit_nh_hmm(d = df_scaled_16_19, response_var = "BTC_USD",  covariates = covariates0, covariates_trans = covariates0)
nh_hmm0_3 <- fit_nh_hmm(d = df_scaled_19_24, response_var = "BTC_USD",  covariates = covariates0, covariates_trans = covariates0)

# Convergence, log likelihood, AIC, BIC
nh_hmm0_1
nh_hmm0_2
nh_hmm0_3



############ A useful plot whenever if it's needed to inspect the states:
#prstates <- apply(posterior(nh_hmm0_3)[,c("S1","S2")], 
#                  1, which.max)
#plot( df_log_19_24$date, prstates, type="b", xlab="Time", ylab="State")

#adjustable for any quick tryouts such as:
covariates12 <- c("DJI", "VOL")
trans1 <- c("MINER")
nh_hmmx_3 <- fit_nh_hmm(d = df_scaled_19_24, response_var = "BTC_USD",  covariates = covariates12, covariates_trans = trans1)
nh_hmmx_3 



#########################################################################
#
#-------------------------------------- EXTRACTING PRIORS
#
# For the next procedure, which will be MCMC, we need to extract
# resulting transition and emission parameters, since they will be used as priors:

extract_nh_hmm_priors <- function(fitted_model, df, covariates) {
  
  y <- df$BTC_USD    
  
  # Initial state probabilities, state sequence
  init_probs <- fitted_model@prior@parameters$coefficients
  posterior_states <- posterior(fitted_model)$state    
  
  # storage for the other priors
  betas <- list()           # logistic regression coefficients for transitions into state 1 and state 2
  trans_probs <- list()     # transition probabilities from the obtained betas
  var_prior_beta <- list()  # variance of transition coefficients
  
  coefficients <- list()    # emision coefficients
  sigma <- list()           # emission sd's
  var_prior <- list()       # emision coefficient variance
  emission_probs <- matrix(NA, nrow = nrow(df), ncol = 2)  # emission probabilities for each state
  
  # Preparing covariates
  covariates <- as.matrix(df[, covariates])                # matrix of X_t
  covariates <- cbind(1, covariates) 
  colnames(covariates)[1] <- "(Intercept)" # adding a first column to the matrix equal to 1. 
  # ^ this represents the intercept column, and is important to match the dimensions
  # in the upcomming logistic regression for the transition probability extraction
  

  for (state in 1:2) {
    
    # ---- From transition model: 
    
    # Logistic regression coefficients for the current state
    betas[[state]] <- fitted_model@transition[[state]]@parameters$coefficients[, state]
    
    eta <- covariates %*% betas[[state]]   # calculating the linear combination of covariates and new betas (etas)
    
    # Transition probabilities
    trans_probs[[state]] <- exp(eta) / (1 + exp(eta))
    
    var_prior_beta[[state]] <- var(betas[[state]])         # coefficient variance
    
    #
    # ---- From emission model:
    response <- fitted_model@response[[state]][[1]]         # specifying the model by the S4 object. Unable to extract the coefficients otherwise
    
    # Emission coefficients
    coefficients[[state]] <- response@parameters$coefficients
    
    # Emission standard deviation (sigma)
    sigma[[state]] <- response@parameters$sd
    
    # Dynamic mean
    mu_t <- covariates %*% coefficients[[state]]      
    
    # Emission probabilities
    emission_probs[, state] <- dnorm(y, mean = mu_t, sd = sigma[[state]])

    var_prior[[state]] <- var(coefficients[[state]])        # coefficient variance
  }
      return(list(init_probs = init_probs,
              posterior_states = posterior_states,
              betas = betas,
              beta_state1 = betas[[1]], 
              beta_state2 = betas[[2]],
              trans_state1 = trans_probs[[1]],
              trans_state2 = trans_probs[[2]], 
              var_prior_beta = var_prior_beta,
              mu = coefficients,
              sigma = sigma,
              var_prior = var_prior,
              #mu_state1 = coefficients[[1]],
              #mu_state2 = coefficients[[2]],
              #sigma_state1 = sigma[[1]],
              #sigma_state2 = sigma[[2]],
              emission_probs = emission_probs, 
              y = y, 
              covariates = covariates, #the updated covariates with an additional column of ones
              covariate_names = colnames(covariates)
              )
             )
}

nh_hmm_priors <- extract_nh_hmm_priors(nh_hmm0_1, df_scaled, covariates0)
#TODO
# Adjust the "fitted_model" and "df" as needed, dependent on which prior parameters are needed
# respectively:
# nh_hmm0_1 - df_scaled
# nh_hmm0_2 - df_scaled_16_19
# nh_hmm0_3 - df_scaled_19_24