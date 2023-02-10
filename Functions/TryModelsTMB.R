# Make TMB objective function by providing starting parameter values, model
# type, and data required.
MakeTMBObj <- function(
    start, mdl_tp = c(
      "popan", "full true kinship", "offset true kinship", "genopair"
    ),
    k = NA, srvy_gaps = NA, fnl_year = NA, srvy_yrs = NA, 
    n_cap_hists = NA, first_tab = NA, last_tab = NA, caps = NA, non_caps = NA, 
    survives = NA,
    alpha = NA, knshp_st_bool = NA,
    ns_SPs_btn = NA, ns_POPs_wtn = NA, ns_POPs_btn = NA, ns_HSPs_wtn = NA, 
    ns_HSPs_btn = NA, 
    ns_caps = NA,
    ns_pairs_wtn = NA, ns_pairs_btn = NA,
    gp_probs = matrix(NA, 1, 1), smp_yr_ind_prs = matrix(NA, 1, 1)
) {
  # Create TMB function
  data <- list(
    mdl_tp = mdl_tp,
    k = k, srvy_gaps = srvy_gaps, fnl_year = fnl_year, srvy_yrs = srvy_yrs, 
    n_cap_hists = n_cap_hists, first_tab = first_tab, last_tab = last_tab, 
    caps = caps, non_caps = non_caps, survives = survives,
    alpha = alpha, knshp_st_bool = knshp_st_bool,
    ns_SPs_btn = ns_SPs_btn, ns_POPs_wtn = ns_POPs_wtn, 
    ns_POPs_btn = ns_POPs_btn, ns_HSPs_wtn = ns_HSPs_wtn, 
    ns_HSPs_btn = ns_HSPs_btn, ns_caps = ns_caps,
    ns_pairs_wtn = ns_pairs_wtn, ns_pairs_btn = ns_pairs_btn,
    gp_probs = gp_probs, smp_yr_ind_prs = smp_yr_ind_prs, 
    n_pairs = nrow(gp_probs)
  )
  MakeADFun(data, list(pars = start), DLL = "UnifiedNLL", silent = T)
}

# Try to fit close kin model with TMB
TryModelTMB <- function(
    obj, lwr, upr, 
    mdl.tp = c(
      "popan", "full true kinship", "offset true kinship", "genopair"
    )
) {
  # cat(mdl.tp, "\n")
  # cat("Optimiser trace: \n")
  
  # Run optimiser starting from true values
  opt <- try(
    suppressWarnings(
      nlminb(
        start = obj$par, obj = obj$fn, grad = obj$gr, hess = obj$he,
        scale = 1 / obj$par, lower = lwr, upper = upr,
        control = list(
          trace = 1,
          iter.max = 2000,
          rel.tol = 1e-4,
          x.tol = 1.2e-2
        )
      )
    )
  )
  
  cat("Reported values from TMB: \n")
  print(summary(sdreport(obj)))
  
  # If optimiser hit error
  if (inherits(opt, "try-error")) {
    cat("Optimiser reports error for", mdl.tp, "model using TMB \n")
    return(NA)
  }
  
  # If Popan model
  if (mdl.tp == "popan") {
    # Get estimates and standard errors from TMB, replacing rho with lambda and
    # inserting estimate for final population size
    k = length(lwr) - 3
    est.se.df = summary(sdreport(obj))[c(k + 4, 2, k + 5, 3:(k + 3)), ]
  } else {
    # Get estimates and standard errors from TMB, replacing rho with lambda
    est.se.df = summary(sdreport(obj))[c(4, 2:3, 5), ]
  }
  
  # Show results
  if (opt$convergence == 0) {
    cat("Optimiser reports success for", mdl.tp, "model using TMB \n")
    cat("Estimates:", round(est.se.df[, 1], 3), "\n")
  } else {
    cat("Optimiser reports failure for", mdl.tp, "model using TMB \n")
    cat("Message:", opt$message, "\n")
  }
  
  # Return results
  list(est.se.df = est.se.df, cnvg = opt$convergence)
}
