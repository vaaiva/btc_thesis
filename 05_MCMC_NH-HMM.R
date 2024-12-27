
##############################
# ------ MCMC-NH-HMM procedure
##############################

# The procedure uses the previously obtained initial state probabilities,
# transition probabilities and logistic regression coefficients from the regular NH-HMM. 
# They're used as priors for the MCMC process.
nh_hmm_priors <- extract_nh_hmm_priors(nh_hmm0_2, df_scaled_16_19, covariates0) 
#TODO
# Adjust the "fitted_model" and "df" as needed, dependent on which prior parameters are needed
# respectively:
# nh_hmm0_1 - df_scaled
# nh_hmm0_2 - df_scaled_16_19
# nh_hmm0_3 - df_scaled_19_24

# It's useful to adjust initial probabilities
# by a small epsilon to avoid division by 0:
nh_hmm_priors$init_probs <- nh_hmm_priors$init_probs+0.0001                          # choosing epsilon = 1e-5
nh_hmm_priors$init_probs <- (nh_hmm_priors$init_probs)/sum(nh_hmm_priors$init_probs) # renorming such that they add up to 1




## ----------------------------------
# Step 1: FORWARD-BACKWARD ALGORITHM
## ----------------------------------
# 

forward_backward<- function(init_probs, trans_state1, trans_state2, emission_probs) {
  
  # Constants and storage
  T <- nrow(emission_probs)           # number of observations (T)
  K <- 2                              # number of states  
  alpha <- matrix(0, T, K)            # forward probabilities (alpha)
  beta <- matrix(0, T, K)             # backward probabilities (beta)
  gamma <- matrix(0, T, K)            # posterior state probabilities (gamma)
  latent_states <- integer(T)         # latent state sequence (z_t)
  
  
  # Forward probabilities (alpha) for the first observation
  alpha[1, 1] <- init_probs[1] * emission_probs[1, 1]
  alpha[1, 2] <- init_probs[2] * emission_probs[1, 2]
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ]) # normalizing
  
  # Forward recursion (t=2 to T)
  for (t in 2:T) {
    alpha[t, 1] <- (alpha[t-1, 1] * trans_state1[t-1] + alpha[t-1, 2] * (1 - trans_state2[t-1])) * emission_probs[t, 1]
    alpha[t, 2] <- (alpha[t-1, 1] * (1 - trans_state1[t-1]) + alpha[t-1, 2] * trans_state2[t-1]) * emission_probs[t, 2]
    alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])# normalizing
  }
  
  # Initializing backward probabilities (beta)
  beta[T, 1] <- 1
  beta[T, 2] <- 1
  
  # Backward Recursion (t=T-1 to 1)
  for (t in (T-1):1) {
    beta[t, 1] <- (beta[t+1, 1] * trans_state1[t] * emission_probs[t+1, 1] +
                     beta[t+1, 2] * (1 - trans_state2[t]) * emission_probs[t+1, 2])
    beta[t, 2] <- (beta[t+1, 1] * (1 - trans_state1[t]) * emission_probs[t+1, 1] +
                     beta[t+1, 2] * trans_state2[t] * emission_probs[t+1, 2])
    beta[t, ] <- beta[t, ] / sum(beta[t, ])  # normalizing
  }
  
  # Posterior state probabilities (gamma)
  for (t in 1:T) {
    gamma[t, 1] <- alpha[t, 1] * beta[t, 1]
    gamma[t, 2] <- alpha[t, 2] * beta[t, 2]
    gamma[t, ] <- gamma[t, ] / sum(gamma[t, ]) # normalizing
  }
  
  # Sampling the latent states (z_t) based on posterior state probabilities
  for (t in 1:T) {
    latent_states[t] <- sample(1:K, size = 1, prob = gamma[t, ])
  }
  return(list(alpha = alpha, beta = beta, gamma = gamma, latent_states = latent_states))
}

forward_backward_results <- forward_backward(
  init_probs = nh_hmm_priors$init_probs, 
  trans_state1 = nh_hmm_priors$trans_state1, 
  trans_state2 = nh_hmm_priors$trans_state2, 
  emission_probs = nh_hmm_priors$emission_probs
)

## ----------------------
# Step 2: GIBBS SAMPLING 
## ----------------------
# sampling the state-dependent regression coefficients

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
    mu_updated[[state]] <- setNames(numeric(n_covariates), colnames(covariates))  #also setting the same colnames
    
    # Gibbs update for EACH covariate coefficient in the emission model
    for (j in 1:n_covariates) {
      covariate_j <- covariates_state[, j]
      
      # Updating the mean for EACH covariate j using normal distribution (conjugate prior)
      mean_post_mu <- (sum(covariate_j * y_state) + mu[[state]][j] * var_prior[[state]]) / (sum(covariate_j^2) + var_prior[[state]])
      var_post_mu <- sigma[[state]] / (sum(covariate_j^2) + var_prior[[state]])
      mu_updated[[state]][j] <- rnorm(1, mean = mean_post_mu, sd = sqrt(var_post_mu))
    }
    
    # Updating the variance (sigma) for the current state using inverse-gamma distribution
    shape_post_sigma <- (n_state / 2) + 0.001
    scale_post_sigma <- sum((y_state - covariates_state %*% mu_updated[[state]])^2) / 2 + 0.001
    sigma_updated[state] <- sqrt(1 / rgamma(1, shape = shape_post_sigma, rate = scale_post_sigma))
  }
  return(mu_updated) # also sigma_updated?
}

gibbs_mu <- gibbs_sampler(
  latent_states = forward_backward_results$latent_states, 
  mu = nh_hmm_priors$mu, 
  sigma = nh_hmm_priors$sigma, 
  var_prior = nh_hmm_priors$var_prior, 
  y = nh_hmm_priors$y, 
  covariates = nh_hmm_priors$covariates
)

## --------------------------------------------
# Step 3: POLYA GAMMA AUGMENTATION AND SAMPLING
## --------------------------------------------
# sampling the logistic regression coefficients

polya_gamma_sampler <- function(latent_states, covariates, betas, var_prior_beta) {
  
  # Constants and storage:
  T <- nrow(covariates)   # number of observations
  K <- 2                  # number of states
  p <- ncol(covariates)   # number of covariates (including intercept)
  beta_updated <- list()  # updated beta coefficients for each state
  var_prior_beta[[1]] <- 0.0001  # manually setting, to prevent a 0 value for state 1


  for (state in 1:K) {

    omega <- numeric(T)    # storage for polya gamma variables
    
    # The linear predictor: covariate * current beta for the current state
    eta <- covariates %*% betas[[state]]  
    
    # Sampling Polya-Gamma latent variables for each observation
    for (t in 1:T) {
      omega[t] <- rpg(1, 1, eta[t])  # ~ PG(1, eta_t) 
      
    ##
    #--- Updating betas iteratively for each covariate ---
  
    beta_updated[[state]] <- numeric(p)                    # storage
    names(beta_updated[[state]]) <- names(betas[[state]])  # preserving coefficient names
    
    prior_var_j <- var_prior_beta[[state]]
    
    # for each covariate j for the current state
    for (j in 1:p) {
      covariate_j <- covariates[, j]      # the j-th covariate
      prior_mean_j <- betas[[state]][j]   # betas as the prior mean for the current coefficient beta_j
      
      # Posterior variance and mean for beta_j using PG augmentation
      V_beta_inv_j <- sum(omega * covariate_j^2) + 1 / prior_var_j    # posterior precision for beta_j
      V_beta_j <- 1 / V_beta_inv_j                                    # posterior variance for beta_j
      
      m_beta_j <- V_beta_j * (sum(covariate_j * as.numeric(latent_states == state)) + prior_mean_j / prior_var_j)  # posterior mean for beta_j
      
      # Drawing the new beta coefficient for the j-th covariate from a normal distribution:
      beta_updated[[state]][j] <- rnorm(1, mean = m_beta_j, sd = sqrt(V_beta_j))
    }
    }
  }
  return(beta_updated)
}

polya_gamma_beta <- polya_gamma_sampler(
  latent_states = forward_backward_results$latent_states, 
  covariates = nh_hmm_priors$covariates, 
  betas = nh_hmm_priors$betas, 
  var_prior_beta = nh_hmm_priors$var_prior_beta
)


## ----------------------
# Step 4: DOUBLE REVERSIBLE JUMP
## ----------------------
# 1. Deciding which covariates to keep in the model (propose_covariate_change)
# 2. Calculating the likelihood (calculate_emission_likelihood, calculate_transition_likelihood)
# 3. Combining into the actual function of a single step of the reversible jump (reversible_jump_step)

propose_covariate_change <- function(current_covariates, all_covariates) {
  current_covariates <- c("(Intercept)", setdiff(current_covariates, "(Intercept)"))  # ensuring that the intercept is always included

  # If there are only two covariates (Intercept + one other), force to add a covariate
  if (length(current_covariates) == 2) {
    available_covariates <- setdiff(all_covariates, current_covariates)
    if (length(available_covariates) > 0) {
      selected_covariate <- sample(available_covariates, 1)
      return(c(current_covariates, selected_covariate)) 
  }}
  
  # Then, randomly add or remove one covariate, based on a random draw. Add if draw is <0.5, remove otherwise
  if (runif(1) < 0.5 && length(setdiff(all_covariates, current_covariates)) > 0) {
    # Add a covariate 
    available_covariates <- setdiff(all_covariates, current_covariates)
    if (length(available_covariates) > 0) {        # additional check if there are available covariates to add
      selected_covariate <- sample(available_covariates, 1)
      return(c(current_covariates, selected_covariate))
    }
  } else if (length(setdiff(current_covariates, "(Intercept)")) > 0) {
    # Remove a covariate 
    removable_covariates <- setdiff(current_covariates, "(Intercept)")
    if (length(removable_covariates) > 0) {       # additional check if there are available covariates to remove
      selected_covariate <- sample(removable_covariates, 1)
      return(setdiff(current_covariates, selected_covariate))
    }
  }
  return(current_covariates)
  }
  
  
# Helper function to calculate the emission likelihood
calculate_emission_likelihood <- function(y, latent_states, covariates, coefficients) {
  likelihood <- 0 #initial value
  
  for (state in 1:2) {
    coefficients_state <- coefficients[[state]]
    mu <- covariates %*% coefficients_state

    # Log likelihood for observations belonging to the current state
    likelihood <- likelihood + sum(dnorm(y[latent_states == state], mean = mu[latent_states == state], log = TRUE))
  }
  return(likelihood)
}

# Helper function to calculate the transition likelihood
calculate_transition_likelihood <- function(latent_states, covariates, coefficients) {
  likelihood <- 0  #initial value
  for (t in 2:length(latent_states)) {                              #starting at t=2, because at t=1 we don't have a state to transition from:)
    beta_state <- coefficients[[latent_states[t-1]]]                #using coefficients of that latent state
    beta_state[1] <- 0

    eta <- covariates[t-1, ] %*% beta_state
    prob <- exp(eta) / (1 + exp(eta))
    
    # Bounding the probabilities to avoid extreme values
    prob <- max(min(prob, 0.999), 0.001)
    
    # Checking if latent_states[t] == 1 and prob is valid
    if (latent_states[t] == 1 && !is.na(prob) && prob > 0 && prob < 1) {
      likelihood <- likelihood + log(prob)
    } else if (latent_states[t] == 2 && !is.na(prob) && prob > 0 && prob < 1) {
      likelihood <- likelihood + log(1 - prob)
    } else {
      print(paste("!!!!!!!!! Invalid at t =", t))
    }
  }
  return(likelihood)
}


# Main reversible jump step in the MCMC flow
reversible_jump_step <- function(y, latent_states, current_emission_covariates, current_transition_covariates, 
                                 current_mu, current_betas, var_prior, sigma, 
                                 covariate_names, covariates, var_prior_beta, prior_mu, prior_betas) {
  
  # 1: Proposing new covariate sets for both emission and transition models
  proposed_emission_names <- propose_covariate_change(colnames(current_emission_covariates), covariate_names)
  proposed_transition_names <- propose_covariate_change(colnames(current_transition_covariates), covariate_names)
  
  proposed_emission_covariates <- covariates[, proposed_emission_names, drop=FALSE]    
  proposed_transition_covariates <- covariates[, proposed_transition_names, drop=FALSE]    
  
  #updating the mu's and the betas
  prior_current_mu <- list(
    prior_mu[[1]][proposed_emission_names],
    prior_mu[[2]][proposed_emission_names]
  )
  prior_current_betas <- list(
    prior_betas[[1]][proposed_transition_names],
    prior_betas[[2]][proposed_transition_names]
  )
  
  # 2: Simulating new coefficients based on proposed covariates
  proposed_mu <- gibbs_sampler(latent_states, prior_current_mu, sigma, var_prior, y, proposed_emission_covariates)
  proposed_betas <- polya_gamma_sampler(latent_states, proposed_transition_covariates, prior_current_betas, var_prior_beta)
  
  # Initializing a matrix that would have small values for the unsampled covariates
   full_proposed_mu <- list(
    rep(0.00001, length(covariate_names)),
    rep(0.00001, length(covariate_names))
  )
  full_proposed_betas <- list(
    rep(0.00001, length(covariate_names)),
    rep(0.00001, length(covariate_names))
  )
  
  # Filling in the coefficients only for selected covariates
  full_proposed_mu[[1]][match(proposed_emission_names, covariate_names)] <- proposed_mu[[1]]
  full_proposed_mu[[2]][match(proposed_emission_names, covariate_names)] <- proposed_mu[[2]]
  full_proposed_betas[[1]][match(proposed_transition_names, covariate_names)] <- proposed_betas[[1]]
  full_proposed_betas[[2]][match(proposed_transition_names, covariate_names)] <- proposed_betas[[2]]
  
  
  
  # 3: Calculating likelihoods for both models, current and proposed
  emission_likelihood_current <- calculate_emission_likelihood(y, latent_states, current_emission_covariates, current_mu)
  emission_likelihood_proposed <- calculate_emission_likelihood(y, latent_states, proposed_emission_covariates, proposed_mu)
  
  transition_likelihood_current <- calculate_transition_likelihood(latent_states, current_transition_covariates, current_betas)
  transition_likelihood_proposed <- calculate_transition_likelihood(latent_states, proposed_transition_covariates, proposed_betas)
  
  # 4: Calculating acceptance ratios
  penalty_add <- 0.1           # favouring simpler models
  penalty_remove <- 0.05
  penalty <- if (length(proposed_emission_names) > length(current_emission_covariates)) 
    penalty_add else -penalty_remove
  emission_acceptance_ratio <- exp(emission_likelihood_proposed - emission_likelihood_current - penalty)
  penalty <- if (length(proposed_transition_names) > length(current_transition_covariates)) 
    penalty_add else -penalty_remove
  transition_acceptance_ratio <- exp(transition_likelihood_proposed - transition_likelihood_current - penalty)
  
  # 5: Metropolis-Hastings acceptance. The random draw ensures that some degree of lower acceptance ratios are still, by chance, accepted
  if (runif(1) < emission_acceptance_ratio) {
    final_emission_covariates <- proposed_emission_covariates
    final_mu <- proposed_mu
  } else {
    final_emission_covariates <- current_emission_covariates
    final_mu <- current_mu
  }
  
  if (runif(1) < transition_acceptance_ratio) {
    final_transition_covariates <- proposed_transition_covariates
    final_betas <- proposed_betas
  } else {
    final_transition_covariates <- current_transition_covariates
    final_betas <- current_betas
  }

  return(list(emission_covariates = final_emission_covariates,
              transition_covariates = final_transition_covariates,
              mu = final_mu,
              betas = final_betas, 
              full_proposed_mu = full_proposed_mu,
              full_proposed_betas = full_proposed_betas))
}

## -------------------------- ##
# ---- THE MAIN MCMC LOOP ---- #
## -------------------------- ##
initial_covariates <- c("(Intercept)", "VIX", "HASH") # starting with a shorter list of covariates


mcmc_nh_hmm <- function(y, init_probs, trans_state1, trans_state2, emission_probs, nh_hmm_priors, num_iterations, burn_in) {
  
  # Initialization:
  current_latent_states <- nh_hmm_priors$latent_states
  current_emission_covariates <-  nh_hmm_priors$covariates[, initial_covariates, drop=FALSE]
  current_transition_covariates <- nh_hmm_priors$covariates[, initial_covariates, drop=FALSE]
  current_mu <- list(
    nh_hmm_priors$mu[[1]][initial_covariates],
    nh_hmm_priors$mu[[2]][initial_covariates]
  )
  current_betas <- list(
    nh_hmm_priors$betas[[1]][initial_covariates],
    nh_hmm_priors$betas[[2]][initial_covariates]
  )
  var_prior <- nh_hmm_priors$var_prior
  sigma <- nh_hmm_priors$sigma
  covariate_names <- nh_hmm_priors$covariate_names
  covariates <- nh_hmm_priors$covariates
  var_prior_beta <- nh_hmm_priors$var_prior_beta
  prior_mu <- nh_hmm_priors$mu
  prior_betas <- nh_hmm_priors$betas
  
  
  # Initializing matrices to store the coefficients for each state for coda
  post_burn_in_iterations <- num_iterations - burn_in
  all_mu_state1 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  all_mu_state2 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  all_betas_state1 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  all_betas_state2 <- matrix(NA, nrow = post_burn_in_iterations, ncol = length(covariate_names))
  colnames(all_mu_state1) <- colnames(all_mu_state2) <- covariate_names
  colnames(all_betas_state1) <- colnames(all_betas_state2) <- covariate_names
  
  # Storage for tracking covariates selection
  emission_covariate_count <- matrix(0, ncol = ncol(covariates), nrow = post_burn_in_iterations)
  transition_covariate_count <- matrix(0, ncol = ncol(covariates), nrow = post_burn_in_iterations)
  colnames(emission_covariate_count) <- covariate_names
  colnames(transition_covariate_count) <- covariate_names
  
  
  ##
  # -----------------  Main MCMC loop
  
  for (iter in 1:num_iterations) {
    cat(paste("Iteration:", iter, "\n"))
    
    # 1: Resampling latent states using forward-backward algorithm
    forward_backward_results <- forward_backward(init_probs, trans_state1, trans_state2, emission_probs)
    current_latent_states <- forward_backward_results$latent_states
    
    # 2: Performing a reversible jump step
    model_update <- reversible_jump_step(y = y,
                                         latent_states = current_latent_states,
                                         current_emission_covariates = current_emission_covariates,
                                         current_transition_covariates = current_transition_covariates,
                                         current_mu = current_mu,
                                         current_betas = current_betas,
                                         var_prior_beta = var_prior_beta,
                                         var_prior = var_prior,
                                         sigma = sigma,
                                         covariate_names = covariate_names,
                                         covariates = covariates, 
                                         prior_mu = prior_mu,
                                         prior_betas = prior_betas)
    
    # 3: Storing the sampled coefficients after the burn-in period

    if (iter > burn_in) {
      index <- iter - burn_in
    all_mu_state1[index, ] <- model_update$full_proposed_mu[[1]]
    all_mu_state2[index, ] <- model_update$full_proposed_mu[[2]]
    all_betas_state1[index, ] <- model_update$full_proposed_betas[[1]]
    all_betas_state2[index, ] <- model_update$full_proposed_betas[[2]]
    
    # 4: Tracking which covariates were selected in this iteration
    selected_emission_covariates <- colnames(current_emission_covariates)
    selected_transition_covariates <- colnames(current_transition_covariates)
    
    emission_covariate_count[index, selected_emission_covariates] <- 1
    transition_covariate_count[index, selected_transition_covariates] <- 1
    }
    
    # 5: Updating the "current" variables for the next iteration
    current_emission_covariates <- model_update$emission_covariates
    current_transition_covariates <- model_update$transition_covariates
    current_mu <- model_update$mu
    current_betas <- model_update$betas
  }
  
# Converting to mcmc objects for diagnostics
mcmc_mu_state1 <- as.mcmc(all_mu_state1)
mcmc_mu_state2 <- as.mcmc(all_mu_state2)
mcmc_betas_state1 <- as.mcmc(all_betas_state1)
mcmc_betas_state2 <- as.mcmc(all_betas_state2)
mcmc_coefs <- mcmc.list(mcmc_mu_state1, mcmc_mu_state2, mcmc_betas_state1, mcmc_betas_state2)

# Calculating inclusion probabilities for each covariate
emission_inclusion_prob <- colMeans(emission_covariate_count)
transition_inclusion_prob <- colMeans(transition_covariate_count)
  
  return(list(mcmc_coefs = mcmc_coefs, 
              final_mu = current_mu,
              final_betas = current_betas,
              emission_inclusion_prob = emission_inclusion_prob,
              transition_inclusion_prob = transition_inclusion_prob))
}


nh_mcmc_results <- mcmc_nh_hmm(y = nh_hmm_priors$y,
                            init_probs = nh_hmm_priors$init_probs,
                            trans_state1 = nh_hmm_priors$trans_state1,
                            trans_state2 = nh_hmm_priors$trans_state2,
                            emission_probs = nh_hmm_priors$emission_probs,
                            nh_hmm_priors = nh_hmm_priors,
                            num_iterations = 100000, 
                            burn_in = 30000)

# Posterior inclusion
nh_mcmc_results$emission_inclusion_prob
nh_mcmc_results$transition_inclusion_prob

# Some diagnostics:
ess_state1mu <- effectiveSize(nh_mcmc_results$mcmc_coefs[[1]])
ess_state2mu <- effectiveSize(nh_mcmc_results$mcmc_coefs[[2]])
ess_state1beta <- effectiveSize(nh_mcmc_results$mcmc_coefs[[3]])
ess_state2beta <- effectiveSize(nh_mcmc_results$mcmc_coefs[[4]])
print(ess_state1mu)
print(ess_state2mu)
print(ess_state1beta)
print(ess_state2beta)

gelman.diag(nh_mcmc_results$mcmc_coefs)
gelman.plot(nh_mcmc_results$mcmc_coefs)

#TODO
############# IMPORTANT - to initialize the procedure for the other subsamples as well, 
# update the "nh_hmm_priors" through "extract_nh_hmm_priors" function, using, respectively:
# nh_hmm0_1 - df_scaled
# nh_hmm0_2 - df_scaled_16_19
# nh_hmm0_3 - df_scaled_19_24


#summary(mcmc_results$mcmc_coefs)
#plot(mcmc_results$mcmc_coefs[[1]])
#plot(mcmc_results$mcmc_coefs[[2]])
#plot(mcmc_results$mcmc_coefs[[3]])
#plot(mcmc_results$mcmc_coefs[[4]])



