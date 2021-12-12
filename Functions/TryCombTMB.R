# Function to try to fit combined model
TryCombTMB <- function() {
  # Create TMB function
  data <- list(k = k, srvygaps = srvy.gaps, ncaphists = n.cap.hists, 
               firsttab = pop.sum$first.tab, 
               lasttab = pop.sum$last.tab, caps = pop.sum$caps, 
               noncaps = pop.sum$non.caps, survives = pop.sum$survives[-k],
               nsPOPswtn = ns.kps.lst$ns.POPs.wtn,
               nsHSPswtn = ns.kps.lst$ns.HSPs.wtn,
               nsPOPsbtn = ns.kps.lst$ns.POPs.btn,
               nsSPsbtn = ns.kps.lst$ns.SPs.btn,
               fyear = f.year, srvyyrs = srvy.yrs, 
               nscaps = ns.caps, alpha = alpha)
  parameters <- list(pars = cbd.start)
  obj <- MakeADFun(data, parameters, DLL = "CombNLL", silent = T)
  
  # Run optimiser starting from true values
  comb.opt <- try(
    nlminb(
      start = obj$par,
      obj = obj$fn,
      gradient = obj$gr,
      hessian = obj$he,
      # scale = c(0.1, 1, 1000, rep(0.1, k)),
      control = list(iter.max = 200),
      lower = cbd.lwr,
      upper = cbd.upr,
    )
  )
  
  # If optimiser hit error
  if (inherits(comb.opt, "try-error")) {
    cat("Optimiser reports error for combined model using TMB \n")
    return(NA)
  }
  
  # lambda = rho + phi
  comb.opt$par[1] <- sum(comb.opt$par[1:2])
  
  # Show results
  if (comb.opt$convergence == 0) {
    cat("Optimiser reports success for combined model using TMB \n")
    cat("Estimates:", round(comb.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for combined model using TMB \n")
  }
  
  # Return results
  c(comb.opt$par[1:3], summary(sdreport(obj))[length(cbd.start) + 1, 1],
    comb.opt$par[c(4:(k + 3))], comb.opt$convergence)
}