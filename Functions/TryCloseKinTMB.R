# Try to fit close kin model with TMB
TryCloseKinTMB <- function() {
  # Create TMB function
  data <- list(
    nsPOPswtn = ns.kps.lst$ns.POPs.wtn, nsHSPswtn = ns.kps.lst$ns.HSPs.wtn, 
    nsPOPsbtn = ns.kps.lst$ns.POPs.btn,
    nsSPsbtn = ns.kps.lst$ns.SPs.btn, k = k, srvygaps = srvy.gaps, 
    fyear = f.year, srvyyrs = srvy.yrs, nscaps = ns.caps, alpha = alpha)
  parameters <- list(pars = ck.start)
  obj <- MakeADFun(data, parameters, DLL = "CloseKinNLL", silent = T)
  
  # Run optimiser starting from true values
  ck.opt <- try(
    nlminb(
      start = obj$par,
      obj = obj$fn,
      grad = obj$gr,
      hess = obj$he,
      # scale = c(0.1, 1, 1000),
      control = list(iter.max = 400),
      lower = ck.lwr, 
      upper = ck.upr,
    )
  )
  
  # If optimiser hit error
  if (inherits(ck.opt, "try-error")) {
    cat("Optimiser reports error for close kin model using TMB \n")
    return(NA)
  }
  
  # lambda = rho + phi
  ck.opt$par[1] <- sum(ck.opt$par[1:2])
  
  # Show results
  if (ck.opt$convergence == 0) {
    cat("Optimiser reports success for close kin model using TMB \n")
    cat("Estimates:", round(ck.opt$par, 3), "\n")
  } else {
    cat("Optimiser reports failure for close kin model using TMB \n")
    cat("Message:", ck.opt$message, "\n")
  }
  
  sd.rep.summ = summary(sdreport(obj))
  print("CK sdrep")
  print(sd.rep.summ)
  
  # Return results
  list(
    c(ck.opt$par[1:3], sd.rep.summ[4, 1], ck.opt$convergence),
    c(tail(sd.rep.summ, 1)[2], sd.rep.summ[2:4, 2])
  )
}