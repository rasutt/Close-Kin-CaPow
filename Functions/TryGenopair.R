# Try to fit genopair model in R
TryGenopair <- function(
    gpps, smp.yr.ind.prs, k, f.year, srvy.yrs, alpha, ck.start, ck.lwr, ck.upr
) {
  # Run optimiser starting from true values
  ck.opt <- nlminb(
    start = ck.start,
    objective = GenopairNLLR, 
    gradient = NULL, hessian = NULL,
    k, f.year, srvy.yrs, alpha, gpps, smp.yr.ind.prs,
    # scale = c(0.1, 1, 1000, rep(0.1, k)),
    control = list(eval.max = 400, iter.max = 1000),
    lower = ck.lwr, 
    upper = ck.upr
  )
  
  # lambda = rho + phi
  ck.opt$par[1] <- sum(ck.opt$par[1:2])
  
  # Show results
  if (ck.opt$convergence == 0) {
    cat("Optimiser reports success for genopair model in R\n")
    cat("Estimates:", round(ck.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for genopair model in R\n")
    cat("Message:", ck.opt$message, "\n")
  }
  
  # Return results
  c(ck.opt$par, ck.opt$convergence)
}