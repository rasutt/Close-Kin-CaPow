# Function to try to fit POPAN model
TryPOPANTMB <- function(
  k, srvy.gaps, n.cap.hists, pop.sum, ppn.start, ppn.lwr, ppn.upr
) {
  # Create TMB function
  data <- list(k = k, srvygaps = srvy.gaps, ncaphists = n.cap.hists, 
               firsttab = pop.sum$first.tab, 
               lasttab = pop.sum$last.tab, caps = pop.sum$caps, 
               noncaps = pop.sum$non.caps, survives = pop.sum$survives[-k])
  parameters <- list(pars = ppn.start)
  obj <- MakeADFun(data, parameters, DLL = "POPANNLL", silent = T)

  # Run optimiser starting from true values
  ppn.opt <- try(
    nlminb(
      start = obj$par,
      obj = obj$fn,
      gradient = obj$gr,
      hessian = obj$he,
      # scale = c(0.1, 1, 1000, rep(0.1, k)),
      control = list(iter.max = 200),
      lower = ppn.lwr,
      upper = ppn.upr,
    )
  )
  
  # If optimiser hit error
  if (inherits(ppn.opt, "try-error")) {
    cat("Optimiser reports error for POPAN model using TMB \n")
    return(NA)
  }
  
  # Get estimates and standard errors from TMB, replacing rho with lambda and
  # inserting estimate for final population size
  est.se.df = summary(sdreport(obj))[c(k + 4, 2, k + 5, 3:(k + 3)), ]

  # Show results
  if (ppn.opt$convergence == 0) {
    cat("Optimiser reports success for POPAN model using TMB \n")
    cat("Estimates:", round(est.se.df[, 1], 3), "\n")
  } else {
    cat("Optimiser reports failure for POPAN model using TMB \n")
  }
  
  # Return results, including derived estimate of final population size
  list(est.se.df = est.se.df, cnvg = ppn.opt$convergence)
}