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

# Kinpair probabilities for true parameter values from TMB objective function
frst.KP.prbs.TVs.lst = reactive({
  ftk.obj()$report(c(rho(), phi(), exp.N.fin()))
})

# Function to format kinpair probabilities for display
frmt.kp.prbs = function(kp.prbs.lst) {
  mat = t(
    sapply(
      kp.prbs.lst[c("prbs_SPs", "prbs_POPs", "prbs_HSPs")], 
      function(kp.prbs.mat) {
        # Pull out and combine values within and between surveys
        c(diag(kp.prbs.mat), kp.prbs.mat[upper.tri(kp.prbs.mat)])
      }
    )
  )
  mat = rbind(c(kp.prbs.lst[["exp_N_s_yrs"]], rep(NA, n.srvy.prs())), mat)
  df = data.frame(format(mat, digits = 3, scientific = T))
  colnames(df) = c(srvy.yrs(), srvy.prs())
  rownames(df) = c("Population size", knshp.chcs)
  df
}

# Kinpair probabilities for first study given true parameter values, from TMB
output$firstKPPrbsTMB = renderTable({
  frmt.kp.prbs(frst.KP.prbs.TVs.lst())
}, rownames = T)

# Kinpair probabilities for first study given true parameter values, from R
output$firstKPPrbsR = renderTable({
  wtn = pred.ns.kps.pop()$wtn
  btn = pred.ns.kps.pop()$btn
  
  prbs.wtn = cbind(wtn[, "N.s.yrs"], 0, wtn[, c("POPs", "HSPs")] / wtn[, "APs"])
  prbs.btn = cbind(NA, btn[, c("SPs", "POPs", "HSPs")] / btn[, "APs"])
  
  df = data.frame(format(t(rbind(prbs.wtn, prbs.btn)), dig = 3, sci = T))
  
  colnames(df) = c(srvy.yrs(), srvy.prs())
  rownames(df) = c("Population size", knshp.chcs)
  
  df
}, rownames = T)

# Kinpair probabilities for first study after optimisation of true kinships
# likelihood
frst.KP.prbs.TKs.lst = reactive({
  ftk.obj()$report()
})

output$firstKPPrbsTK = renderTable({
  frmt.kp.prbs(frst.KP.prbs.TKs.lst())
}, rownames = T)

