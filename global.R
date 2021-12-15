# Set bounds on inputs and step size for input slider and NLL plots
min_rho = 0.03
max_rho = 0.13
step_rho = 0.005
min_phi = 0.91
max_phi = 0.99
step_phi = 0.005

# Set simulation parameters for basic scenario
exp.Ns <- 2169 # Expected superpopulation size
alpha <- 8 # Age of sexual maturity
hist.len <- 80 # Duration of simulation in years
p <- 0.1 # Capture probability
stch.bths <- T # Set whether numbers of births stochastic
clvng.p <- 0 # Additional capture probability when calving
tmp.emgn <- 0 # Probability of males being away from survey area

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to work sometimes!
library(TMB)
compile("TMB_objective_functions/POPANNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/POPANNLL"))
compile("TMB_objective_functions/CloseKinNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CloseKinNLL"))