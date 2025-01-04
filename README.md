## --------------
## ----- SCRIPTS:
# 00_main: 
loading the required libraries, sourcing the subsequent scripts if needed.
# 01_data:
importing the data directly from Yahoo Finance and
Blockchain.com. Initial preparation - merging, mutating, renaming,
interpolating missing values. Data transformation - log transforming,
scaling. Initial data plots and summaries. Alternative variable
selection procedure tryouts (Lasso, EN, RF).
# 02_HMM: 
defining the framework for homogeneous HMM’s. Function for
extracting model parameters, to be used as priors in ’04_MCMC_HMM’.
# 03_NH-HMM: 
defining the framework for non - homogeneous HMM’s. Function
for extracting model parameters, to be used as priors in ’05_MCMC_NH-HMM’.
# 04_MCMC_HMM: 
defining interim functions and the main loop for the
Bayesian MCMC covariate selection with homogeneous HMM’s. Obtaining
posterior inclusion results, diagnostics measures.
# 05_MCMC_NH-HMM: 
defining interim functions and the main loop for the
Bayesian MCMC covariate selection with non-homogeneous HMM’s. Obtaining
posterior inclusion results, diagnostics measures.
# 06_HMM_rolling_forecast: 
specifying the benchmarks sets. Defining the
bootstrapped rolling windows forecasting procedure with a homogeneous HMM for a given forecast
horizon. Obtaining accuraccy metrics. Plots - Actual vs. Forecasted
values, dynamic rolled coefficients over time, dynamic MAPE over time.
# 07_NH_HMM_rolling_forecast: 
specifying the benchmarks sets. Defining
the bootstrapped rolling windows forecasting procedure with a non-homogeneous HMM for a given
forecast horizon. Obtaining accuraccy metrics.
## -----------
# 
The current version in github.com specifies the full sample
calculations in all the functions. For calculations with the
subsamples, ensure the function parameters are adjusted accordingly.
The parameters to be adjusted are marked in the code and highlighted by
’#TODO’, with descriptive comments on what exactly has to be changed in
each case.
