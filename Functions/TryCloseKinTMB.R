# Try to fit close kin model with TMB
TryCloseKinTMB <- function(
  k, srvy.gaps, f.year, srvy.yrs, ns.caps, ns.kps.lst, ck.start, ck.lwr, ck.upr,
  alpha
) {
  # Create TMB function
  data <- list(
    nsPOPswtn = ns.kps.lst$ns.POPs.wtn, nsHSPswtn = ns.kps.lst$ns.HSPs.wtn, 
    nsPOPsbtn = ns.kps.lst$ns.POPs.btn, nsSPsbtn = ns.kps.lst$ns.SPs.btn, k = k, 
    srvygaps = srvy.gaps, fyear = f.year, srvyyrs = srvy.yrs, nscaps = ns.caps, 
    alpha = alpha
  )
  obj <- MakeADFun(data, list(pars = ck.start), DLL = "CloseKinNLL", silent = T)
  
  # Run optimiser starting from true values
  ck.opt <- try(
    nlminb(
      start = obj$par, obj = obj$fn, grad = obj$gr, hess = obj$he,
      # scale = c(0.1, 1, 1000),
      control = list(iter.max = 400), lower = ck.lwr, upper = ck.upr,
    )
  )
  
  # If optimiser hit error
  if (inherits(ck.opt, "try-error")) {
    cat("Optimiser reports error for close kin model using TMB \n")
    return(NA)
  }
  
  # Get estimates and standard errors from TMB, replacing rho with lambda
  est.se.df = summary(sdreport(obj))[c(4, 2:3, 5), ]

  # Show results
  if (ck.opt$convergence == 0) {
    cat("Optimiser reports success for close kin model using TMB \n")
    cat("Estimates:", round(est.se.df[, 1], 3), "\n")
  } else {
    cat("Optimiser reports failure for close kin model using TMB \n")
    cat("Message:", ck.opt$message, "\n")
  }
  
  # Return results
  list(est.se.df = est.se.df, cnvg = ck.opt$convergence)
}