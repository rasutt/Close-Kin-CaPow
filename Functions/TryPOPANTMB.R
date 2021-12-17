# Function to try to fit POPAN model
TryPOPANTMB <- function() {
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
  
  # Get estimates and standard errors for all parameters requested
  sd.rep.summ = summary(sdreport(obj))
  
  # Replace estimate for rho with that for lambda, as code assumes first
  # estimated parameter is lambda for now
  ppn.opt$par[1] <- sd.rep.summ[k + 4, 1]
  
  # Show results
  if (ppn.opt$convergence == 0) {
    cat("Optimiser reports success for POPAN model using TMB \n")
    cat("Estimates:", round(ppn.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for POPAN model using TMB \n")
  }
  
  print("PPN sdrep")
  print(sd.rep.summ)
  
  # Return results, including derived estimate of final population size
  list(
    c(ppn.opt$par[1:2], tail(sd.rep.summ, 1)[1], 
      ppn.opt$par[3:(k + 3)], ppn.opt$convergence),
    c(sd.rep.summ[k + 4, 2], sd.rep.summ[2, 2], tail(sd.rep.summ, 1)[2],
      sd.rep.summ[3:(k + 3), 2])
  )
}