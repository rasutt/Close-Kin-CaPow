# Set parameters and simulate studies.  ~10 seconds for 1000 studies

# Set simulation name
# sim.name <- "test"
# sim.name <- "det_N"
# sim.name <- "0.4_clvng.p"
# sim.name <- "0.5_tmp.emgn"
sim.name <- "0.1_pmt.emgn"

# Set number of studies to simulate
n.stds.sim <- 1000

# Load simulation and POPAN functions
source("Functions/SimPopStud.R")
source("Functions/PopanFuncs.R")

# Set simulation parameters for basic scenario
exp.Ns <- 2169 # Expected superpopulation size
lambda <- 1.03 # Population growth rate
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

# Change parameters for alternative scenarios
# stch.bths <- F # Set whether numbers of births stochastic
# clvng.p <- 0.4 # Additional capture probability when calving
# tmp.emgn <- 0.5 # Probability of males being away from survey area

# Find useful global variables
rho <- lambda - phi # "Population birthrate"
k <- length(srvy.yrs) # Number of surveys per study (> 2 for POPAN)
f.year <- srvy.yrs[k] # Final survey year
srvy.gaps <- as.integer(diff(srvy.yrs)) # Survey gaps
stdy.len <- sum(srvy.gaps) # Length of study
n.srvy.prs <- choose(k, 2) # Number of pairs of surveys

# Expected final population size.  gaps variables only used here
lambda.gaps <- lambda^srvy.gaps # Population growth rate between surveys
phi.gaps <- phi^srvy.gaps # Individual survival rate between surveys
exp.N.fin <- sum(pent_func(lambda.gaps, phi.gaps) * exp.Ns * 
                   prod(phi.gaps) / cumprod(c(1, phi.gaps)))

# Expected population size over time
exp.N.t <- exp.N.fin / lambda^((hist.len - 1):0)
N.init <- round(exp.N.t[1]) # Initial population size

# Create list for population and capture histories
hists.lst <- vector("list", n.stds.sim)

# Display progress
cat("Simulating study: ")

# Loop over histories
for (hist.ind in 1:n.stds.sim) {
  # Display progress
  if (hist.ind %% 100 == 1) cat(hist.ind, "")
  
  # Simulate family and capture histories of population of animals over time
  hists.lst[[hist.ind]] <- SimPopStud()
}

# Display progress
cat(hist.ind, "")

# Save parameters and histories
save(list = objects(),
     file = paste0("Simulated data/", sim.name, "_sim.Rdata"))
