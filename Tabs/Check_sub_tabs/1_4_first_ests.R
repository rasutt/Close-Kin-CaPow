# Find estimates for first study

# Popan model
first.ppn.ests = reactive({
  # Try to fit popan likelihood model
  rslt = TryModelTMB(ppn.obj(), ppn.lwr(), ppn.upr(), "popan")
  
  # If no error return results, otherwise return NA and non-convergence
  if (!all(is.na(rslt))) {
    c(rslt$est.se.df[, 1], rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# True kinship model
first.ck.ests = reactive({
  # Try to fit true kinship likelihood model
  rslt = TryModelTMB(tk.obj(), ck.lwr(), ck.upr(), "true kinship")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Full genopair model 
first.FGP.ests = reactive({
  # Try to fit genopair likelihood model
  rslt = TryModelTMB(GPP.obj.fll(), ck.lwr(), ck.upr(), "genopair")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Offset model
first.OGP.ests = reactive({
  # Try to fit genopair likelihood model
  rslt = TryModelTMB(GPP.obj.offst(), ck.lwr(), ck.upr(), "genopair")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Print results for first study
output$firstResults = renderTable({
  # True parameter values
  res.mat = matrix(c(par.vals(), NA), nrow = 1)
  res.mat[1, 3:4] = c(sim.lst()$N.fin.vec[1], sim.lst()$Ns.vec[1])
  colnames(res.mat) = c(est.par.names(), "Convergence")
  
  # Add results
  res.mat = rbind(
    res.mat, first.ppn.ests(), first.ck.ests(),
    first.FGP.ests(), first.OGP.ests()
  )
  
  # Add model names
  res.df = data.frame(
    model = c(
      "True values", "Popan", "True kinship", "Full genopair", "Offset genopair"
    ),
    res.mat
  )
  res.df[, ncol(res.df)] = res.df[, ncol(res.df)] == 0
  
  res.df
}, digits = 3)

# Kinpair probabilities for first study after optimisation of true kinships
# likelihood
frst.KP.prbs.mat = reactive({
  matrix(
    summary(sdreport(tk.obj()))[-(1:5), "Estimate"], 
    nrow = k()^2, dimnames = list(srvy.pr = NULL, kp.tp = c("SP", "POP", "HSP"))
  )
})
output$firstHSPPrbsTK = renderTable(
  frmt.kp.prbs(frst.KP.prbs.mat()[, "HSP"]), rownames = T
)
output$firstPOPPrbsTK = renderTable(
  frmt.kp.prbs(frst.KP.prbs.mat()[, "POP"]), rownames = T
)
output$firstSPPrbsTK = renderTable(
  frmt.kp.prbs(frst.KP.prbs.mat()[, "SP"]), rownames = T
)
