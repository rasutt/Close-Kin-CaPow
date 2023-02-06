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
  rslt = TryModelTMB(ftk.obj(), ck.lwr(), ck.upr(), "full true kinship")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Offset true kinship model
first.otk.ests = reactive({
  # Try to fit true kinship likelihood model
  rslt = TryModelTMB(otk.obj(), ck.lwr(), ck.upr(), "offset true kinship")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Full genopair model 
first.fg.ests = reactive({
  # Try to fit genopair likelihood model
  rslt = TryModelTMB(fg.obj(), ck.lwr(), ck.upr(), "genopair")
  
  # If no error
  if (!all(is.na(rslt))) {
    # Combine with missing values for capture probabilities
    c(rslt$est.se.df[, 1], rep(NA, k()), rslt$cnvg)
  } else {
    c(rep(NA, k() + 4), 1)
  }
})

# Offset model
first.og.ests = reactive({
  # Try to fit genopair likelihood model
  rslt = TryModelTMB(og.obj(), ck.lwr(), ck.upr(), "genopair")
  
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
    res.mat, first.ppn.ests(), first.ck.ests(), first.otk.ests(),
    first.fg.ests(), first.og.ests()
  )
  
  # Add model names
  res.df = data.frame(
    model = c(
      "True values", "Popan", "True kinship", "Offset true kinship", 
      "Full genopair", "Offset genopair"
    ),
    res.mat
  )
  res.df[, ncol(res.df)] = res.df[, ncol(res.df)] == 0
  
  res.df
}, digits = 3)

# Kinpair probabilities for first study after optimisation of true kinships
# likelihood
frst.KP.prbs.TKs.lst = reactive({
  ftk.obj()$report()
})

output$firstKPPrbsTK = renderTable({
  frmt.kp.prbs(frst.KP.prbs.TKs.lst())
}, rownames = T)

