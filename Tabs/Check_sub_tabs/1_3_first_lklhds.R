# Code and outputs for first study estimates sub-tab of check-tab

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start = reactive({
  vec = c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  names(vec) = c("rho", "phi", "N.final")
  vec
})

# Create general optimizer starting-values and bounds, NAs filled in below
ppn.start = reactive(c(ck.start()[-3], FS.atts()$Ns, rep(p(), k())))

# Get numbers of animals captured in study
n.cap.hists = reactive(nrow(frst.std()))

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
ftk.obj = reactive({
  # Get numbers of animals captured in each survey
  ns.caps = FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst = FindNsKinPairs(k(), n.srvy.prs(), frst.std())
  
  # Create TMB function
  MakeTMBObj(
    ck.start(), "full true kinship",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ], 
    ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
    ns_HSPs_btn = ns.kps.lst$btn[3, ], ns_caps = ns.caps
  )
})

# Offset true kinship objective function from TMB
otk.obj = reactive({
  # Get numbers of pairs for each pair of surveys
  ns.pairs = table(data.frame(frst.osyips()))

  # Find numbers of kin pairs
  ns.kps.lst = FindNsOKPs(
    k(), n.srvy.prs(), frst.std(), frst.osiips(), frst.osyips()
  )

  # Create TMB function
  MakeTMBObj(
    ck.start(), "offset true kinship",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(),
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    ns_SPs_btn = ns.kps.lst$btn[1, ], ns_POPs_wtn = ns.kps.lst$wtn[1, ],
    ns_POPs_btn = ns.kps.lst$btn[2, ], ns_HSPs_wtn = ns.kps.lst$wtn[2, ],
    ns_HSPs_btn = ns.kps.lst$btn[3, ], 
    ns_pairs_wtn = diag(ns.pairs), 
    ns_pairs_btn = t(ns.pairs)[lower.tri(ns.pairs)]
  )
})

# Full and offset genopair objective functions from TMB
fg.obj = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln,
    gp_probs = frst.fgps(), smp_yr_ind_prs = frst.fsyips()
  )
})
og.obj = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), knshp_st_bool = all.knshps.bln, 
    gp_probs = frst.ogps(), smp_yr_ind_prs = frst.osyips()
  )
})

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
    0, c(n.pts, 3, 5), dimnames = list(
      NULL, NULL, 
      c("Popan", "True kinship", "Offset true kinship", "Full genopair", 
        "Offset genopair")
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
      obj.par[i, p, "True kinship"] = ftk.obj()$fn(ck.par.vals)
      obj.par[i, p, "Offset true kinship"] = otk.obj()$fn(ck.par.vals)
      obj.par[i, p, "Full genopair"] = fg.obj()$fn(ck.par.vals)
      obj.par[i, p, "Offset genopair"] = og.obj()$fn(ck.par.vals)
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
  par(mfrow = c(5, 3), mar = rep(2, 4))
  plot.nll("Popan")
  plot.nll("True kinship")
  plot.nll("Offset true kinship")
  plot.nll("Full genopair")
  plot.nll("Offset genopair")
})


