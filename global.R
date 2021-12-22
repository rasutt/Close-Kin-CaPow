# Set bounds on inputs and step size for input slider and NLL plots
min_rho = 0
max_rho = 0.16
step_rho = 0.005
min_phi = 0.9
max_phi = 1
step_phi = 0.005

# Model choices
mod_choices = c("Popan", "Close-kin")

# Set simulation parameters for basic scenario
alpha <- 8 # Age of sexual maturity
p <- 0.1 # Capture probability
stch.bths <- T # Set whether numbers of births stochastic
clvng.p <- 0 # Additional capture probability when calving
tmp.emgn <- 0 # Probability of males being away from survey area

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 