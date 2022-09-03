# Set bounds on inputs and step size for input slider and NLL plots
min_rho = 0
max_rho = 0.16
step_rho = 0.005
min_phi = 0.9
max_phi = 1
step_phi = 0.005

# Model choices
mod.choices = c("Popan", "Close-kin")

# Set simulation parameters for basic scenario
stch.bths <- T # Set whether numbers of births stochastic

# Probability of males permanently emigrating each year. Note, code is not yet
# consistent with non-zero value, see readme.
pmt.emgn <- 0 