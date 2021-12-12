# Function to try to fit combined model
TryComb <- function() {
  # Run optimiser starting from true values
  comb.opt <- nlminb(
    start = cbd.start,
    obj = CombNLL,
    # scale = c(0.1, 1, 1000, rep(0.1, k)),
    control = list(iter.max = 200),
    lower = cbd.lwr,
    upper = cbd.upr,
  )
  
  # lambda = rho + phi
  comb.opt$par[1] <- sum(comb.opt$par[1:2])
  
  # Show results
  if (comb.opt$convergence == 0) {
    cat("Optimiser reports success for combined model \n")
    cat("Estimates:", round(comb.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for combined model \n")
  }
  
  # Return results
  c(comb.opt$par, comb.opt$convergence)
}