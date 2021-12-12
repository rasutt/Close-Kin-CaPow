# Function to try to fit POPAN model
TryPOPAN <- function() {
  # Run optimiser starting from true values
  pop.opt <- nlminb(
    start = ppn.start,
    obj = popan_nll,
    # scale = c(0.1, 1, 1000, rep(0.1, k)),
    control = list(iter.max = 200),
    lower = ppn.lwr,
    upper = ppn.upr,
  )
  
  # lambda = rho + phi
  pop.opt$par[1] <- sum(pop.opt$par[1:2])
  
  # Show results
  if (pop.opt$convergence == 0) {
    cat("Optimiser reports success for POPAN model \n")
    cat("Estimates:", round(pop.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for POPAN model \n")
  }
  
  # Return results
  c(pop.opt$par, pop.opt$convergence)
}