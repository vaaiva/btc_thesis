
rm(list=ls())

set.seed(1207)

# Loading libraries
library(tidyverse)    # important to load BEFORE 'tidyquant' due to confilts for select()
library(tidyquant)
library(xts)
library(depmixS4)     # for HMM
library(httr)         # for loading hash rate and block size data 
library(jsonlite)     # for loading hash rate and block size data 
library(ggcorrplot)
library(BayesLogit)   # for polya gamma augmentation
library(coda)         # for storing mcmc coefficients
library(forecast)     
library(glmnet)       # for variable selection with lasso and elastic net
library(randomForest) # for variable selection with RF
library(reshape2)
library(zoo)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(psych)

# FLOW:
source("01_data.R")
source("02_HMM.R")
source("03_NH-HMM.R")
source("04_MCMC_HMM.R")
source("05_MCMC_NH-HMM.R")
source("06_HMM_rolling_forecast.R")
source("07_NH_HMM_rolling_forecast.R")



