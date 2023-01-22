# Code and outputs for first study estimates sub-tab of check-tab

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start = reactive({
  vec = c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  names(vec) = c("rho", "phi", "N.final")
  vec
})
ppn.start = reactive(c(ck.start()[-3], FS.atts()$Ns, rep(p(), k())))

# Get numbers of animals captured in study
n.cap.hists = reactive(nrow(fst.std()))

# Popan likelihood TMB objective function
ppn.obj = reactive({
  # Create TMB function
  MakeTMBObj(
    ppn.start(), "popan",
    k(), srvy.gaps(), 
    n_cap_hists = n.cap.hists(), first_tab = frst.pop.sum()$first.tab, 
    last_tab = frst.pop.sum()$last.tab, caps = frst.pop.sum()$caps, 
    non_caps = frst.pop.sum()$non.caps, survives = frst.pop.sum()$survives[-k()]
  )
})

# True kinship likelihood TMB objective function - True kinships 
tk.obj = reactive({
  # Get numbers of animals captured in each survey
  ns.caps = FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst = FindNsKinPairs(k(), n.srvy.prs(), fst.std())

  # Create TMB function
  MakeTMBObj(
    ck.start(), "true kinship",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ], 
    ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
    ns_HSPs_btn = ns.kps.lst$btn[3, ], ns_caps = ns.caps
  )
})

# Genopair likelihood TMB objective function
GPP.obj.fll = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    gp_probs = GPPs.fll(), smp_yr_ind_prs = frst.SYIPs.fll()
  )
})
GPP.obj.offst = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln, 
    gp_probs = GPPs.offst(), smp_yr_ind_prs = frst.SYIPs.offst()
  )
})

n.pts = 200

# Create general optimizer bounds
ck.lwr = reactive(c(0, 0.75, FS.atts()$ns.caps[k()]))
ck.upr = reactive(c(0.35, 1, Inf))
ppn.lwr = reactive(c(ck.lwr()[1:2], n.cap.hists(), rep(0, k())))
ppn.upr = reactive(c(ck.upr(), rep(1, k())))

# Matrix of values to plot NLL at for each parameter, n_points x 3
ck.par.mat = reactive({
  ck.upr = ck.upr()
  ck.upr[3] = 1e4
  sapply(1:3, function(p) {
    seq(ck.lwr()[p], ck.upr[p], len = n.pts)
  })
})
ppn.par.mat = reactive({
  cbind(
    ck.par.mat()[, 1:2],
    seq(ppn.lwr()[3], 1e4, len = n.pts)
  )
})

# Negative log-likelihood surfaces for first study, n_points x n_pars x n_models
frst.nll.srfcs = reactive({
  obj.par = array(
    0, c(n.pts, 3, 4), dimnames = list(
      NULL, NULL, 
      c("Popan", "True kinship", "Full genopair", "Offset genopair")
    )
  )
  
  # Loop over parameters
  for (p in 1:3) {
    # Reset all parameter values
    ck.par.vals = ck.start()
    ppn.par.vals = ppn.start()
    
    # Loop over points for plotting
    for (i in 1:n.pts) {
      # Set value of current parameter
      ck.par.vals[p] = ck.par.mat()[i, p]
      ppn.par.vals[p] = ppn.par.mat()[i, p]
      
      # Find NLL values at current parameter values
      obj.par[i, p, "Popan"] = ppn.obj()$fn(ppn.par.vals)
      obj.par[i, p, "True kinship"] = tk.obj()$fn(ck.par.vals)
      obj.par[i, p, "Full genopair"] = GPP.obj.fll()$fn(ck.par.vals)
      obj.par[i, p, "Offset genopair"] = GPP.obj.offst()$fn(ck.par.vals)
    }
  }
  
  obj.par
})

# Function to plot NLL over each parameter while holding others at true values,
# for one model
plot.nll = function(mdl.nm) {
  p.nms = c("Rho", "Phi", if (mdl.nm == "Popan") "Ns" else "N.final")
  
  # Loop over parameters
  for (p in 1:3) {
    plot(
      if (mdl.nm == "Popan") ppn.par.mat()[, p] else ck.par.mat()[, p], 
      frst.nll.srfcs()[, p, mdl.nm], 
      main = paste(p.nms[p], "-", mdl.nm), 
      xlab = p.nms[p], 
      ylab = mdl.nm, 
      type = 'l', col = 1
    )
    abline(
      v = if (mdl.nm == "Popan") ppn.start()[p] else ck.start()[p], col = 2
    )
    # legend(
    #   "topright", legend = c("Likelihood", "True value"), col = 1:2, lty = 1
    # )
  }
}

# Plot NLL surfaces for each parameter with others held at true values, for each
# model
output$firstNLLSurfs = renderPlot({
  par(mfrow = c(4, 3), mar = rep(2, 4))
  plot.nll("Popan")
  plot.nll("True kinship")
  plot.nll("Full genopair")
  plot.nll("Offset genopair")
})

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

# Print kinpair probabilities for first study under true parameters
output$firstKPPrbs = renderTable({
  sumsd = summary(sdreport(tk.obj()))[, "Estimate"]
  srvy.pr.inds = 1:n.srvy.prs()^2
  prbs_SPs = sumsd[5 + srvy.pr.inds]

  # True parameter values
  res.mat = data.frame(matrix(
    format(prbs_SPs, digits = 3, scientific = T), nrow = k()
  ))
  colnames(res.mat) = rownames(res.mat) = srvy.yrs()
  
  res.mat
}, rownames = T)

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


