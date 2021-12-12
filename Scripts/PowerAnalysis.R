# Load populations and studies, and find and analyze maximum likelihood
# estimates for POPAN, CKMR, and combined models

# Set sim name if loading previously simulated data
sim.name <- "test"
# sim.name <- "det_N"
# sim.name <- "0.4_clvng.p"
# sim.name <- "0.5_tmp.emgn"
# sim.name <- "0.1_pmt.emgn"

# Load parameters and simulated populations and studies
load(file = paste0("Simulated data/", sim.name, "_sim.Rdata"))
sim.name <- "test2"

# Set number of studies to fit models to (reduce for faster testing)
n.stds.fit <- 1000

# Load functions
funcs <- list.files("Functions")
for (i in 1:length(funcs)) source(paste0("Functions/",funcs[i]))

# Load TMB library and likelihood functions.  Have to restart R for compile and
# dyn.load to work sometimes!
library(TMB)
compile("TMB_objective_functions/POPANNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/POPANNLL"))
compile("TMB_objective_functions/CloseKinNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CloseKinNLL"))
compile("TMB_objective_functions/CombNLL.cpp")
dyn.load(dynlib("TMB_objective_functions/CombNLL"))

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start <- c(rho, phi, NA)
ck.lwr <- c(0, 0.75, NA)
ck.upr <- c(0.35, 1, Inf)
ppn.start <- cbd.start <- c(ck.start, rep(p, k))
ppn.lwr <- cbd.lwr <- c(ck.lwr, rep(0, k))
ppn.upr <- cbd.upr <- c(ck.upr, rep(1, k))

# Create vectors for superpopulation and final population sizes
Ns.vec <- N.fin.vec <- numeric(n.stds.fit)

# Create matrices for estimates
ppn.ests <- ck.ests <- cbd.ests <- ppn.tmb.ests <- ck.tmb.ests <- 
  cbd.tmb.ests <- matrix(nrow = n.stds.fit, ncol = 5 + k, dimnames = list(
    NULL, c("lambda", "phi", "N_final", "Ns", paste0("p", 1:k), "cnvg")))

# Loop over histories
for (hist.ind in 1:n.stds.fit) {
  # Display progress
  cat("History:", hist.ind, "\n")
  
  # Get simulated family and capture histories of population of animals over
  # time
  pop.cap.hist <- hists.lst[[hist.ind]]
  
  # Store superpopulation and final population size
  Ns.vec[hist.ind] <- attributes(pop.cap.hist)$Ns
  N.fin.vec[hist.ind] <- attributes(pop.cap.hist)$N.t.vec[hist.len]
  
  # Get numbers of animals captured in study and each survey
  n.cap.hists <- nrow(pop.cap.hist)
  ns.caps <- attributes(pop.cap.hist)$ns.caps
  
  # Summarise data for POPAN model
  pop.sum <- FindPopSum()
  
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs()

  # Update optimiser starting-values and bounds
  ppn.start[3] <- attributes(pop.cap.hist)$Ns
  ppn.lwr[3] <- n.cap.hists
  ck.start[3] <- cbd.start[3] <- N.fin.vec[hist.ind]
  ck.lwr[3] <- cbd.lwr[3] <- ns.caps[k]
  
  # Try to fit models
  # ppn.ests[hist.ind, -3] <- TryPOPAN()
  ppn.tmb.ests[hist.ind, ] <- TryPOPANTMB()
  # ck.ests[hist.ind, -(4:(4 + k))] <- TryCloseKin()
  ck.tmb.ests[hist.ind, -(5:(4 + k))] <- TryCloseKinTMB()
  # cbd.ests[hist.ind, -4] <- TryComb()
  cbd.tmb.ests[hist.ind, ] <- TryCombTMB()
}

# # Find POPAN estimates of N_2020
# ppn.ests[, 3] <- N_fin_vec_func(ppn.ests[, 4], ppn.ests[, 1], ppn.ests[, 2])
# 
# # Find close kin and combined model estimates of Ns
# ck.ests[, 4] <- Ns_vec_func(ck.ests[, 3], ck.ests[, 1], ck.ests[, 2])
# cbd.ests[, 4] <- Ns_vec_func(cbd.ests[, 3], cbd.ests[, 1], cbd.ests[, 2])

# Combine model estimates as list
mod.ests.lst <- list(
  # ppn = ppn.ests,
  ppn_tmb = ppn.tmb.ests,
  # ck = ck.ests,
  ck_tmb = ck.tmb.ests, 
  # cbd = cbd.ests,
  cbd_tmb = cbd.tmb.ests
)  

# Compare model results
CompMods()