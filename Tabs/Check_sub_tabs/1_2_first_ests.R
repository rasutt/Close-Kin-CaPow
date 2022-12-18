# Code and outputs for first study estimates sub-tab of check-tab

# Genopair probabilities for optimization
GPPs.fll = reactive(FindGPPs(frst.LGPPs.KP.fll()[, -2]))
GPPs.offst = reactive(FindGPPs(frst.LGPPs.KP.offst()[, -2]))

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start = reactive(c(rho(), phi(), FS.atts()$N.t.vec[hist.len()]))

# Genopair likelihood TMB objective function
GPP.obj.fll = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), 
    gpprobs = GPPs.fll(), sampyrinds = frst.SYIPs.fll()
  )
})
GPP.obj.offst = reactive({
  MakeTMBObj(
    ck.start(), "genopair",
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), 
    alpha = alpha(), 
    gpprobs = GPPs.offst(), sampyrinds = frst.SYIPs.offst()
  )
})

# Close-kin likelihood TMB objective function - True kinships 
ck.obj = reactive({
  # Get numbers of animals captured in each survey
  ns.caps <- FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), fst.std())
  
  # Create TMB function
  MakeTMBObj(
    ck.start = ck.start(), mdltp = "true kinship",
    k = k(), srvygaps = srvy.gaps(), fyear = fnl.year(), 
    srvyyrs = srvy.yrs(), 
    alpha = alpha(), 
    nsSPsbtn = ns.kps.lst$btn[1, ], nsPOPsbtn = ns.kps.lst$btn[2, ],
    nsPOPswtn = ns.kps.lst$wtn[1, ], # nsHSPswtn = ns.kps.lst$wtn[3, ],
    nscaps = ns.caps
  )
})

# # Popan likelihood TMB objective function 
# ppn.obj = reactive({
#   # Get numbers of animals captured in each survey
#   ns.caps <- FS.atts()$ns.caps
#   
#   # Find numbers of kin pairs
#   ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), fst.std())
#   
#   # Create TMB function
#   MakeTMBObj(
#     ck.start = ck.start(), mdltp = "true kinship",
#     k = k(), srvygaps = srvy.gaps(), fyear = fnl.year(), 
#     srvyyrs = srvy.yrs(), 
#     alpha = alpha(), 
#     nsSPsbtn = ns.kps.lst$btn[1, ], nsPOPsbtn = ns.kps.lst$btn[2, ],
#     nsPOPswtn = ns.kps.lst$wtn[1, ], # nsHSPswtn = ns.kps.lst$wtn[3, ],
#     nscaps = ns.caps
#   )
# })

n.pts = 200

# Create general optimizer bounds
ck.lwr <- reactive(c(0, 0.75, FS.atts()$ns.caps[k()]))
ck.upr <- reactive(c(0.35, 1, Inf))

# Matrix of values to plot NLL at for each parameter, n_points x 3
par.mat = reactive({
  ck.upr = ck.upr()
  ck.upr[3] = 1e4
  sapply(1:3, function(p) {
    seq(ck.lwr()[p], ck.upr[p], len = n.pts)
  })
})

# Negative log-likelihood surfaces for first study, n_points x n_pars x n_models
frst.nll.srfcs = reactive({
  obj.par = array(
    0, c(n.pts, 3, 3), dimnames = list(
      NULL, NULL, c("Full genopair", "Offset genopair", "True kinship")
    )
  )
  
  # Loop over parameters
  for (p in 1:3) {
    # Reset all parameter values
    par.vals = ck.start()
    
    # Loop over points for plotting
    for (i in 1:n.pts) {
      # Set value of current parameter
      par.vals[p] = par.mat()[i, p]
      
      # Find NLL values at current parameter values
      obj.par[i, p, "Full genopair"] = GPP.obj.fll()$fn(par.vals)
      obj.par[i, p, "Offset genopair"] = GPP.obj.offst()$fn(par.vals)
      obj.par[i, p, "True kinship"] = ck.obj()$fn(par.vals)
    }
  }
  
  obj.par
})

# Function to plot NLL over each parameter while holding others at true values,
# for one model
plot.nll = function(mdl.nm) {
  par(mfrow = c(1, 3))
  p.nms = c("Rho", "Phi", "N.final")
  
  # Loop over parameters
  for (p in 1:3) {
    plot(
      par.mat()[, p], frst.nll.srfcs()[, p, mdl.nm], main = p.nms[p], 
      xlab = p.nms[p], ylab = "NLL", type = 'l', col = 1
    )
    abline(v = ck.start()[p], col = 2)
    legend(
      "topright", legend = c("Likelihood", "True value"), col = 1:2, lty = 1
    )
  }
}

# Plot NLL surfaces for each parameter with others held at true values, for each
# model
output$firstFGPNLLSurfs = renderPlot(plot.nll("Full genopair"))
output$firstOGPNLLSurfs = renderPlot(plot.nll("Offset genopair"))
output$firstCKNLLSurfs = renderPlot(plot.nll("True kinship"))

# Find estimates for full genopair model for first study
first.FGP.ests = reactive({
  # Try to fit genopair likelihood model
  print(table(frst.SYIPs.fll()[, 1], frst.SYIPs.fll()[, 2]))
  print(str(GPPs.fll()))
  gp.tmb = TryModelTMB(GPP.obj.fll(), ck.lwr(), ck.upr(), "genopair")
  
  # Combine with missing values for capture probabilities
  c(gp.tmb$est.se.df[, 1], rep(NA, k()), gp.tmb$cnvg)
})
first.OGP.ests = reactive({
  # Try to fit genopair likelihood model
  print(table(frst.SYIPs.offst()[, 1], frst.SYIPs.offst()[, 2]))
  print(str(GPPs.offst()))
  gp.tmb = TryModelTMB(GPP.obj.offst(), ck.lwr(), ck.upr(), "genopair")

  # Combine with missing values for capture probabilities
  c(gp.tmb$est.se.df[, 1], rep(NA, k()), gp.tmb$cnvg)
})

# Find estimates for full genopair model for first study
first.ck.ests = reactive({
  # Get numbers of animals captured in each survey
  ns.caps <- FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), fst.std())
  
  # Try to fit close-kin likelihood model
  ck.tmb = TryModelTMB(ck.obj(), ck.lwr(), ck.upr(), "true kinship")
  
  c(ck.tmb$est.se.df[, 1], rep(NA, k()), ck.tmb$cnvg)
})

# Print results for first study
output$firstResults <- renderTable({
  # True parameter values
  res.mat = matrix(c(par.vals(), NA), nrow = 1)
  res.mat[1, 3:4] = c(sim.lst()$N.fin.vec[1], sim.lst()$Ns.vec[1])
  colnames(res.mat) = c(est.par.names(), "Convergence")

  # Add results
  res.mat = rbind(res.mat, first.FGP.ests(), first.OGP.ests(), first.ck.ests())

  # Add model names
  res.df = data.frame(
    model = 
      c("True values", "Full genopair", "Offset genopair", "True kinship"),
    res.mat
  )
  res.df[, ncol(res.df)] = res.df[, ncol(res.df)] == 0

  res.df
}, digits = 3)


