# Set bounds on inputs and step size for input slider and NLL plots
min_lambda = 0.98
max_lambda = 1.08
step_lambda = 0.01

# Set simulation parameters for basic scenario
exp.Ns <- 2169 # Expected superpopulation size
phi <- 0.95 # Individual survival rate
alpha <- 8 # Age of sexual maturity
hist.len <- 80 # Duration of simulation in years
srvy.yrs <- c(1995:1998, 2006:2009, 2020) # Survey years
p <- 0.1 # Capture probability
stch.bths <- T # Set whether numbers of births stochastic
clvng.p <- 0 # Additional capture probability when calving
tmp.emgn <- 0 # Probability of males being away from survey area

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 

# Find useful global variables
k <- length(srvy.yrs) # Number of surveys per study (> 2 for POPAN)
f.year <- srvy.yrs[k] # Final survey year
srvy.gaps <- as.integer(diff(srvy.yrs)) # Survey gaps
stdy.len <- sum(srvy.gaps) # Length of study
n.srvy.prs <- choose(k, 2) # Number of pairs of surveys

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to work sometimes!
library(TMB)
compile("TMB_objective_functions/CloseKinNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CloseKinNLL"))