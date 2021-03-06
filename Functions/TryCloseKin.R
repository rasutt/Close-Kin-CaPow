# Try to fit close kin model
TryCloseKin <- function() {
  # Run optimiser starting from true values
  ck.opt <- nlminb(
    start = ck.start,
    obj = CloseKinNLL,
    # scale = c(0.1, 1, 1000, rep(0.1, k)),
    control = list(eval.max = 400, iter.max = 1000),
    lower = ck.lwr, 
    upper = ck.upr,
  )
  
  # lambda = rho + phi
  ck.opt$par[1] <- sum(ck.opt$par[1:2])
  
  # Show results
  if (ck.opt$convergence == 0) {
    cat("Optimiser reports success for close kin model \n")
    cat("Estimates:", round(ck.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for close kin model \n")
    cat("Message:", ck.opt$message, "\n")
  }
  
  # Return results
  c(ck.opt$par, ck.opt$convergence)
}