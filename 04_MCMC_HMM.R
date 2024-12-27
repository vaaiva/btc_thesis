
##############################
# ------ MCMC- HMM procedure
##############################

# The procedure uses the previously obtained initial state probabilities,
# transition probabilities and  coefficients from the regular HMM. 
# They're used as priors for the MCMC process.
hmm_c_priors <- extract_hmm_c_priors(fitted_model = hmm_c_0_1, df = df_scaled, 
                                     response_var = "BTC_USD", covariates = covariates0, 
                                     num_states = 2) 
#TODO
# Adjust the "fitted_model" and "df" as needed, dependent on which prior parameters are needed
# respectively:
# hmm_c_0_1 - df_scaled
# hmm_c_0_2 - df_scaled_16_19
# hmm_c_0_3 - df_scaled_19_24

# It's useful to adjust initial probabilities
# by a small epsilon to avoid division by 0:
hmm_c_priors$init_probs <- hmm_c_priors$init_probs+0.0001                          # choosing epsilon = 1e-5
hmm_c_priors$init_probs <- (hmm_c_priors$init_probs)/sum(hmm_c_priors$init_probs)  # renorming such that they add up to 1




## ----------------------------------
# Step 1: FORWARD-BACKWARD ALGORITHM
## ----------------------------------
# 

forward_backward_homogeneous<- function(init_probs, trans_probs, emission_probs) {
  
  # Constants and storage:
  T <- nrow(emission_probs)           # number of observations (T)
  K <- 2                              # number of states  
  alpha <- matrix(0, T, K)            # forward probabilities (alpha)
  beta <- matrix(0, T, K)             # backward probabilities (beta)
  gamma <- matrix(0, T, K)            # posterior state probabilities (gamma)
  latent_states <- integer(T)         # latent state sequence (z_t)
  
  
  # Forward probabilities (alpha) for the first observation
  alpha[1, 1] <- init_probs[1] * emission_probs[1, 1]
  alpha[1, 2] <- init_probs[2] * emission_probs[1, 2]
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])    # normalizing
  
  # Forward recursion (t=2 to T)
  for (t in 2:T) {
    alpha[t, 1] <- (alpha[t-1, 1] * trans_probs[1, 1] + alpha[t-1, 2] * trans_probs[2, 1]) * emission_probs[t, 1]
    alpha[t, 2] <- (alpha[t-1, 1] * trans_probs[1, 2] + alpha[t-1, 2] * trans_probs[2, 2]) * emission_probs[t, 2]
    
    alpha[t, ] <- alpha[t, ] / sum(alpha[t, ]) # normalizing
  }
  
  # Initializing backward probabilities (beta)
  beta[T, 1] <- 1
  beta[T, 2] <- 1
  
  # Backward recursion (t=T-1 to 1)
  for (t in (T-1):1) {
    beta[t, 1] <- (beta[t+1, 1] * trans_probs[1, 1] * emission_probs[t+1, 1] +
                     beta[t+1, 2] * trans_probs[1, 2] * emission_probs[t+1, 2])
    beta[t, 2] <- (beta[t+1, 1] * trans_probs[2, 1] * emission_probs[t+1, 1] +
                     beta[t+1, 2] * trans_probs[2, 2] * emission_probs[t+1, 2])
    
    beta[t, ] <- beta[t, ] / sum(beta[t, ]) # normalizing
  }
  
  # Posterior state probabilities (gamma)
  for (t in 1:T) {
    gamma[t, 1] <- alpha[t, 1] * beta[t, 1]
    gamma[t, 2] <- alpha[t, 2] * beta[t, 2]
    
    gamma[t, ] <- gamma[t, ] / sum(gamma[t, ])   # normalizing
  }
  
  # Sampling the latent states (z_t) based on posterior state probabilities
  for (t in 1:T) {
    latent_states[t] <- sample(1:K, size = 1, prob = gamma[t, ])
  }
  
  return(list(alpha = alpha, beta = beta, gamma = gamma, latent_states = latent_states))
}

forward_backward_c_results <- forward_backward_homogeneous(
  init_probs = hmm_c_priors$init_probs,
  trans_probs = hmm_c_priors$trans_probs,
  emission_probs = hmm_c_priors$emission_probs
)


#
## ----------------------------------
# Step 2: GIBBS SAMPLER
## ----------------------------------

gibbs_sampler <- function(latent_states, mu, sigma, var_prior, y, covariates) {
  
  # Constants and storage:
  K <- 2                                                  # number of states
  n_covariates <- ncol(covariates)                        # number of covariates (including intercept)
  mu_updated <- list()                                    # updated emission coefficients as a list for each state
  sigma_updated <- numeric(K)                             # updated sd
  

  for (state in 1:K) {
    
    # Filtering the data for the current state based on latent states
    y_state <- y[latent_states == state]
    covariates_state <- covariates[latent_states == state, ]
    n_state <- length(y_state)                             # number of observations in the current state
    
    # Initializing a vector to store the updated coefficients for each state
    mu_updated[[state]] <- setNames(numeric(n_covariates), colnames(covariates))  #alsi setting the same colnames
    
    # Gibbs update for EACH covariate coefficient in the emission model
    for (j in 1:n_covariates) {
      covariate_j <- covariates_state[, j]
      
      # Updating the mean for EACH covariate j using normal distribution (conjugate prior)
      mean_post_mu <- (sum(covariate_j * y_state) + mu[[state]][j] * var_prior[[state]]) / (sum(covariate_j^2) + var_prior[[state]])
      var_post_mu <- sigma[[state]] / (sum(covariate_j^2) + var_prior[[state]])
      mu_updated[[state]][j] <- rnorm(1, mean = mean_post_mu, sd = sqrt(var_post_mu))
    }
    
    # Updating the variance (sigma) for the current state using inverse-gamma distribution (not needed but could be useful?)
    shape_post_sigma <- (n_state / 2) + 0.001
    scale_post_sigma <- sum((y_state - covariates_state %*% mu_updated[[state]])^2) / 2 + 0.001
    sigma_updated[state] <- sqrt(1 / rgamma(1, shape = shape_post_sigma, rate = scale_post_sigma))
  }
  return(mu_updated) 
}

gibbs_mu_c <- gibbs_sampler(
  latent_states = forward_backward_c_results$latent_states, 
  mu = hmm_c_priors$coefficients, 
  sigma = hmm_c_priors$sigma, 
  var_prior = hmm_c_priors$var_prior, 
  y = hmm_c_priors$y, 
  covariates = hmm_c_priors$covariates
)


## ----------------------
# Step 3: DOUBLE REVERSIBLE JUMP
## ----------------------
# 1. Deciding which covariates to keep in the model (propose_covariate_change)
# 2. Calculating the likelihood (calculate_emission_likelihood)
# 3. Combining into the actual function of a single step of the reversible jump (reversible_jump_step_c)

propose_covariate_change <- function(current_covariates, all_covariates) {
  current_covariates <- c("(Intercept)", setdiff(current_covariates, "(Intercept)"))  # ensuring that the intercept is always included
  
  
  # If there are only two covariates (Intercept + one other), force to add a covariate
  if (length(current_covariates) == 2) {
    available_covariates <- setdiff(all_covariates, current_covariates)
    if (length(available_covariates) > 0) {
      selected_covariate <- sample(available_covariates, 1)
      return(c(current_covariates, selected_covariate)) 
    }
  }
  
  # Then, randomly add or remove one covariate, based on a random draw. Add if draw is <0.5, remove otherwise
  if (runif(1) < 0.5 && length(setdiff(all_covariates, current_covariates)) > 0) {
    # Add a covariate 
    available_covariates <- setdiff(all_covariates, current_covariates)
    if (length(available_covariates) > 0) {  # additional check if there are available covariates to add
      selected_covariate <- sample(available_covariates, 1)
      return(c(current_covariates, selected_covariate))
    }
  } else if (length(setdiff(current_covariates, "(Intercept)")) > 0) {
    # Remove a covariate 
    removable_covariates <- setdiff(current_covariates, "(Intercept)")
    if (length(removable_covariates) > 0) {   # additional check if there are available covariates to remove
      selected_covariate <- sample(removable_covariates, 1)
      return(setdiff(current_covariates, selected_covariate))
    }
  }
  
  return(current_covariates)
}

# Helper function to calculate the emission likelihood
calculate_emission_likelihood <- function(y, latent_states, covariates, coefficients) {
  likelihood <- 0   # initial value
  
  for (state in 1:2) {
    coefficients_state <- coefficients[[state]]
    mu <- covariates %*% coefficients_state
    
    # Log likelihood for observations belonging to the current state
    likelihood <- likelihood + sum(dnorm(y[latent_states == state], mean = mu[latent_states == state], log = TRUE))
  }
  return(likelihood)
}


# Main reversible jump step in the MCMC flow 
reversible_jump_step_c <- function(y, latent_states, current_emission_covariates,
                                 current_mu, var_prior, sigma, 
                                 covariate_names, covariates, prior_mu) {
  
  # 1: Proposing new covariate set
  proposed_emission_names <- propose_covariate_change(colnames(current_emission_covariates), covariate_names)
  proposed_emission_covariates <- covariates[, proposed_emission_names, drop=FALSE]    

  # updating the current mu's (storage of coefficients, for comparing both likelihoods)
  prior_current_mu <- list(
    prior_mu[[1]][proposed_emission_names],
    prior_mu[[2]][proposed_emission_names]
  )

  # 2: Simulating new coefficients based on proposed covariates
  proposed_mu <- gibbs_sampler(latent_states, prior_current_mu, sigma, var_prior, y, proposed_emission_covariates)
 
  
  # Initializing a matrix that would have small values for the unsampled covariates 
  full_proposed_mu <- list(
    rep(0.00001, length(covariate_names)),
    rep(0.00001, length(covariate_names))
  )

  # Filling in the coefficients only for selected covariates
  full_proposed_mu[[1]][match(proposed_emission_names, covariate_names)] <- proposed_mu[[1]]
  full_proposed_mu[[2]][match(proposed_emission_names, covariate_names)] <- proposed_mu[[2]]

  
  # 3: Calculating likelihood for both models 
  emission_likelihood_current <- calculate_emission_likelihood(y, latent_states, current_emission_covariates, current_mu)
  emission_likelihood_proposed <- calculate_emission_likelihood(y, latent_states, proposed_emission_covariates, proposed_mu)
  
  # 4: Calculate acceptance ratios
  penalty_add <- 0.1        # favouring simpler models
  penalty_remove <- 0.05
  penalty <- if (length(proposed_emission_names) > length(current_emission_covariates)) penalty_add else -penalty_remove
  emission_acceptance_ratio <- exp(emission_likelihood_proposed - emission_likelihood_current - penalty)
  #emission_acceptance_ratio <- exp(emission_likelihood_proposed - emission_likelihood_current)  # without penalty
 
  # 5: Metropolis-Hastings acceptance step. The random draw ensures that some degree of lower acceptance ratios are still, by chance, accepted
  if (runif(1) < emission_acceptance_ratio) {
    final_emission_covariates <- proposed_emission_covariates
    final_mu <- proposed_mu
  } else {
    final_emission_covariates <- current_emission_covariates
    final_mu <- current_mu
  }

  return(list(emission_covariates = final_emission_covariates,
              mu = final_mu,
              full_proposed_mu = full_proposed_mu))
}



## -------------------------- ##
# ---- THE MAIN MCMC LOOP ---- #
## -------------------------- ##
initial_covariates <- c("(Intercept)", "VIX", "HASH") # starting with a shorter list of covariates

mcmc_hmm_c <- function(y, init_probs, trans_probs, emission_probs, hmm_c_priors, num_iterations, burn_in) {
  
  # Initialization for the first iteration:
  current_latent_states <- forward_backward_c_results$latent_states
  current_emission_covariates <- hmm_c_priors$covariates[, initial_covariates, drop=FALSE]
  current_mu <- list(
    hmm_c_priors$coefficients[[1]][initial_covariates],
    hmm_c_priors$coefficients[[2]][initial_covariates]
  )
  var_prior <- hmm_c_priors$var_prior
  sigma <- hmm_c_priors$sigma
  covariate_names <- hmm_c_priors$covariate_names
  covariates <- hmm_c_priors$covariates
  prior_mu <- hmm_c_priors$coefficients

  
  # Initializing matrices to store the coefficients for each state for coda (mcmc objects)
  post_burn_in_iterations <- num_iterations - burn_in
  all_mu_state1 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  all_mu_state2 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  colnames(all_mu_state1) <- colnames(all_mu_state2) <- covariate_names
  
  # Storage for tracking covariates selection (posterior inclusion prob)
  emission_covariate_count <- matrix(0, ncol = ncol(covariates), nrow = post_burn_in_iterations)
  colnames(emission_covariate_count) <- covariate_names
  
  ##
  # -----------------  Main MCMC loop
  
  for (iter in 1:num_iterations) {
    cat(paste("Iteration:", iter, "\n"))
    
    # 1: Resampling latent states using Forward-Backward algorithm
    forward_backward_c_results <- forward_backward_homogeneous(init_probs, trans_probs, emission_probs)
    current_latent_states <- forward_backward_c_results$latent_states
    
    # 2: Reversible jump step (covariate change, gibbs sampling, likelihood calculation, acceptance)
    model_update <- reversible_jump_step_c(y = y,
                                         latent_states = current_latent_states,
                                         current_emission_covariates = current_emission_covariates,  
                                         current_mu = current_mu,                                    
                                         var_prior = var_prior,
                                         sigma = sigma,
                                         covariate_names = covariate_names,
                                         covariates = covariates, 
                                         prior_mu = prior_mu)  # after the step, the model is updated (or not) on the "decided" new parameters
                                                               # so every iteration is conditional on the last one
    
    
    # 3: Storing the sampled coefficients only post-burn-in period
    if (iter > burn_in) {
      index <- iter - burn_in
      all_mu_state1[index, ] <- model_update$full_proposed_mu[[1]]
      all_mu_state2[index, ] <- model_update$full_proposed_mu[[2]]
      
      # 4: Tracking which covariates were selected in this iteration
      selected_emission_covariates <- colnames(current_emission_covariates)
      emission_covariate_count[index, selected_emission_covariates] <- 1
    }

    # 5: Updating the "current" variables for the next iteration
    current_emission_covariates <- model_update$emission_covariates
    current_mu <- model_update$mu
  }
  
  # Converting to mcmc objects for diagnostics
  mcmc_mu_state1 <- as.mcmc(all_mu_state1)
  mcmc_mu_state2 <- as.mcmc(all_mu_state2)
  mcmc_coefs <- mcmc.list(mcmc_mu_state1, mcmc_mu_state2)
  
  # Calculating posterior inclusion probabilities for each covariate
  emission_inclusion_prob <- colMeans(emission_covariate_count)
  
  return(list(mcmc_coefs = mcmc_coefs, 
              final_mu = current_mu,
              emission_inclusion_prob = emission_inclusion_prob))
}


mcmc_hmm_c_results <- mcmc_hmm_c(y = hmm_c_priors$y,
                            init_probs = hmm_c_priors$init_probs,
                            trans_probs = hmm_c_priors$trans_probs,
                            emission_probs = hmm_c_priors$emission_probs,
                            hmm_c_priors = hmm_c_priors,
                            num_iterations = 100000, 
                            burn_in = 30000)

# Posterior inclusion
mcmc_hmm_c_results$emission_inclusion_prob

# Some diagnostics:
ess_state1 <- effectiveSize(mcmc_hmm_c_results$mcmc_coefs[[1]])
ess_state2 <- effectiveSize(mcmc_hmm_c_results$mcmc_coefs[[2]])
print(ess_state1)
print(ess_state2)
gelman.diag(mcmc_hmm_c_results$mcmc_coefs)
gelman.plot(mcmc_hmm_c_results$mcmc_coefs)

#TODO
############# IMPORTANT - to initialize the procedure for the other subsamples as well, 
# update the "hmm_c_priors" through "extract_hmm_priors_c" function, using, respectively:
# hmm_c_0_1 - df_scaled
# hmm_c_0_2 - df_scaled_16_19
# hmm_c_0_3 - df_scaled_19_24





#summary(mcmc_hmm_c_results$mcmc_coefs)
#plot(mcmc_hmm_c_results$mcmc_coefs) # too large
#plot(mcmc_hmm_c_results$mcmc_coefs[[1]])
#plot(mcmc_hmm_c_results$mcmc_coefs[[2]])


