# Code and outputs for first study estimates sub-tab of check-tab

# Genopair probabilities for optimization
gpp.opt = reactive({
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  lg.gpp.slct = frst.LGPPs.KP.fll()[, -2]
  gpp.slct = exp(lg.gpp.slct)
  all_undrflw = rowSums(gpp.slct) == 0
  
  # If there is underflow adjust log-probabilities by factor giving equal
  # weight to smallest and largest values to avoid both under and overflow
  if (any(all_undrflw)) {
    cat("Proportion of pairs for which probabilities given all kinships 
        underflow to zero:", mean(all_undrflw), "\n")
    
    # Want smallest maximum kinship probability and largest probability to be
    # equally far from one
    adj = mean(c(min(apply(lg.gpp.slct, 1, max)), max(lg.gpp.slct)))
    lg.gpp.adj = lg.gpp.slct - adj
    gpp.adj = exp(lg.gpp.adj)
    
    # Show adjustment and results
    cat("Probabilities adjusted by factor of exp(", adj, ")\n", sep = "")
    print("Adjusted log-probabilities and probabilities of genopairs:")
    print(summary(lg.gpp.adj))
    print(summary(gpp.adj))
    cat("Proportion of pairs for which adjusted probabilities given all 
        kinships underflow to zero:", mean(rowSums(gpp.adj) == 0), "\n")
    
    return(gpp.adj)
  } 
  else {
    return(gpp.slct)
  }
})

# Genopair likelihood TMB objective function
gpp.obj = reactive({
  ck.start <- c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  
  MakeGPObj(
    gpp.opt(), frst.SYIPs.fll(),
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), alpha(), ck.start
  )
})

# Close-kin likelihood TMB objective function - True kinships 
ck.obj = reactive({
  ck.start <- c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  
  # Get numbers of animals captured in each survey
  ns.caps <- FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), fst.std())
  
  # Create TMB function
  data <- list(
    nsPOPswtn = ns.kps.lst$wtn[1, ], 
    nsSPsbtn = ns.kps.lst$btn[1, ], nsPOPsbtn = ns.kps.lst$btn[2, ], 
    k = k(), srvygaps = srvy.gaps(), fyear = fnl.year(), srvyyrs = srvy.yrs(), 
    nscaps = ns.caps, alpha = alpha()
  )
  obj <- MakeADFun(data, list(pars = ck.start), DLL = "CloseKinNLL", silent = T)
})

# Create general optimizer starting-values and bounds, NAs filled in below
ck.start = reactive({
  c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
})

n.pts = 200

# Matrix of values to plot NLL at for each parameter, n_points x 3
par.mat = reactive({
  ck.lwr <- c(0, 0.75, FS.atts()$ns.caps[k()])
  ck.upr <- c(0.35, 1, 1e4)
  sapply(1:3, function(p) {
    seq(ck.lwr[p], ck.upr[p], len = n.pts)
  })
})

# Negative log-likelihood surfaces for first study, n_points x n_pars x n_models
frst.nll.srfcs = reactive({
  obj.par = array(
    0, c(n.pts, 3, 2), 
    dimnames = list(NULL, NULL, c("Full genopair", "True kinship"))
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
      obj.par[i, p, 1] = gpp.obj()$fn(par.vals)
      obj.par[i, p, 2] = ck.obj()$fn(par.vals)
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
output$firstGPNLLSurfs = renderPlot(plot.nll("Full genopair"))
output$firstCKNLLSurfs = renderPlot(plot.nll("True kinship"))

# Find estimates for full genopair model for first study - The estimates for the
# final population size don't seem to follow the likelihood surface like they do
# for the true kinships model... Specifying the hessian and parameter scale in
# nlminb seems to help somewhat
first.gp.ests = reactive({
  # Create general optimizer starting-values and bounds
  ck.start <- c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  ck.lwr <- c(0, 0.75, FS.atts()$ns.caps[k()])
  ck.upr <- c(0.35, 1, Inf)
  
  # Try to fit genopair likelihood model
  print(table(frst.SYIPs.fll()[, 1], frst.SYIPs.fll()[, 2]))
  print(str(gpp.opt()))
  gp.tmb = TryGenopairTMB(
    gpp.opt(), frst.SYIPs.fll(),
    k(), srvy.gaps(), fnl.year(), srvy.yrs(), ck.start, ck.lwr, ck.upr, alpha()
  )
  
  # # Checked that it give the same results in R
  # gp.r = TryGenopair(
  #   if (any(all_undrflw)) gpp.adj else gpp.slct, 
  #   smp.yr.ind.prs() + 1,
  #   k(), fnl.year(), srvy.yrs(), alpha(), 
  #   ck.start, ck.lwr, ck.upr
  # )
  # print(gp.r)
  
  # Combine with missing values for capture probabilities
  rbind(gp.tmb$est.se.df, matrix(NA, k(), 2))
})

# Find estimates for full genopair model for first study
first.ck.ests = reactive({
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  ck.lwr <- c(0, 0.75, FS.atts()$ns.caps[k()])
  ck.upr <- c(0.35, 1, Inf)
  
  # Get numbers of animals captured in each survey
  ns.caps <- FS.atts()$ns.caps
  
  # Find numbers of kin pairs
  ns.kps.lst <- FindNsKinPairs(k(), n.srvy.prs(), fst.std())
  
  # Try to fit close-kin likelihood model
  ck.tmb = TryCloseKinTMB(
    ns.kps.lst, k(), srvy.gaps(), fnl.year(), srvy.yrs(), ns.caps, 
    ck.start, ck.lwr, ck.upr, alpha()
  )
  
  rbind(ck.tmb$est.se.df, matrix(NA, k(), 2))
})

# Print results for first study
output$firstResults <- renderTable({
  # True parameter values
  res.mat = matrix(par.vals(), nrow = 1)
  res.mat[1, 3:4] = c(sim.lst()$N.fin.vec[1], sim.lst()$Ns.vec[1])
  colnames(res.mat) = est.par.names()

  # Add results
  res.mat = rbind(
    res.mat, 
    t(first.gp.ests()[, 1]),
    t(first.ck.ests()[, 1])
  )

  # Add model names
  res.df = data.frame(
    model = c("True values", "Genopairs", "True kinships"), res.mat
  )

  res.df
}, digits = 3)


