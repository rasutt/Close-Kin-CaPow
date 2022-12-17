# Make TMB objective function by providing starting parameter values, model
# type, and data required.
MakeTMBObj <- function(
    ck.start, mdltp = c("true kinship", "genopair"),
    k = NA, srvygaps = NA, fyear = NA, srvyyrs = NA, alpha = NA, 
    nsSPsbtn = NA, nsPOPsbtn = NA, nsPOPswtn = NA, # nsHSPswtn = NA, 
    nscaps = NA,
    gpprobs = matrix(NA, 1, 1), sampyrinds = matrix(NA, 1, 1)
) {
  # Create TMB function
  data <- list(
    mdltp = mdltp,
    k = k, srvygaps = srvygaps, fyear = fyear, srvyyrs = srvyyrs, 
    alpha = alpha, 
    nsSPsbtn = nsSPsbtn, nsPOPsbtn = nsPOPsbtn, nsPOPswtn = nsPOPswtn, 
    # nsHSPswtn = nsHSPswtn, 
    nscaps = nscaps,
    gpprobs = gpprobs, sampyrinds = sampyrinds, npairs = nrow(gpprobs)
  )
  MakeADFun(data, list(pars = ck.start), DLL = "UnifiedNLL", silent = T)
}

# Try to fit close kin model with TMB
TryModelTMB <- function(obj, lwr, upr, mdltp = c("true kinship", "genopair")) {
  # Run optimiser starting from true values
  opt <- try(
    nlminb(
      start = obj$par, obj = obj$fn, grad = obj$gr, hess = obj$he,
      scale = 1 / obj$par,
      control = list(iter.max = 400), lower = lwr, upper = upr
    )
  )
  
  # If optimiser hit error
  if (inherits(opt, "try-error")) {
    cat("Optimiser reports error for", mdltp, "model using TMB \n")
    return(NA)
  }
  
  # Get estimates and standard errors from TMB, replacing rho with lambda
  est.se.df = summary(sdreport(obj))[c(4, 2:3, 5), ]
  
  # Show results
  if (opt$convergence == 0) {
    cat("Optimiser reports success for", mdltp, "model using TMB \n")
    cat("Estimates:", round(est.se.df[, 1], 3), "\n")
  } else {
    cat("Optimiser reports failure for", mdltp, "model using TMB \n")
    cat("Message:", opt$message, "\n")
  }
  
  # Return results
  list(est.se.df = est.se.df, cnvg = opt$convergence)
}
