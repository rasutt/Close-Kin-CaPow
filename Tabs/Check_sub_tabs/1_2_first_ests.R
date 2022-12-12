# Code and outputs for first study estimates sub-tab of check-tab

# Genopair probabilities for optimization
gpp.opt = reactive({
  # Get genopair probabilities (by excluding probabilities giveb half-sibs for
  # now) and check for pairs where all probabilities underflow to zero
  lg.gpp.slct = frst.lg.gp.prbs.KP()[, -2]
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
    gpp.opt(), frst.smp.yr.ind.prs(),
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

ck.start = reactive({
  N.fnl = FS.atts()$N.t.vec[hist.len()]
  
  # Create general optimizer starting-values and bounds, NAs filled in below
  c(rho(), phi(), N.fnl)
})

n.pts = 200

par.vec = reactive(function(p) {
  ck.lwr <- c(0, 0.75, FS.atts()$ns.caps[k()])
  ck.upr <- c(0.35, 1, 1e4)
  par.vec = seq(ck.lwr[p], ck.upr[p], len = n.pts)
})

frst.nll.srfcs = reactive({
  obj.par = array(0, c(n.pts, 2, 3))
  for (p in 1:3) {
    p.strt = ck.start()
    
    for (i in 1:n.pts) {
      p.strt[p] = par.vec()(p)[i]
      obj.par[i, 1, p] = gpp.obj()$fn(p.strt)
      obj.par[i, 2, p] = ck.obj()$fn(p.strt)
    }
  }
  obj.par
})

plot.nll = function(mdl.ind) {
  par(mfrow = c(1, 3))
  p.nms = c("Rho", "Phi", "N.final")
  for (p in 1:3) {
    plot(
      par.vec()(p), frst.nll.srfcs()[, mdl.ind, p], main = p.nms[p], 
      xlab = p.nms[p], 
      ylab = "NLL", type = 'l', col = 1:2
    )
    abline(v = ck.start()[p], col = 3)
    if (p == 3) {
      legend(
        "topright", legend = c("Likelihood", "True value"), col = 1:2
      )
    }
  }
}

# Negative log-likelihood surfaces for each parameter with others held at true
# values
output$firstGPNLLSurfs = renderPlot({
  plot.nll(1)
})
output$firstCKNLLSurfs = renderPlot({
  plot.nll(2)
})

first.gp.ests = reactive({
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), FS.atts()$N.t.vec[hist.len()])
  ck.lwr <- c(0, 0.75, FS.atts()$ns.caps[k()])
  ck.upr <- c(0.35, 1, Inf)
  
  # Try to fit genopair likelihood model
  print(table(frst.smp.yr.ind.prs()[, 1], frst.smp.yr.ind.prs()[, 2]))
  print(str(gpp.opt()))
  gp.tmb = TryGenopairTMB(
    gpp.opt(), frst.smp.yr.ind.prs(),
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
  
  rbind(gp.tmb$est.se.df, matrix(NA, k(), 2))
})

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

  res.mat = rbind(
    res.mat, 
    t(first.gp.ests()[, 1]),
    t(first.ck.ests()[, 1])
  )

  res.df = data.frame(
    model = c("True values", "Genopairs", "True kinships"), 
    res.mat  
  )

  # Formatting
  res.df[, 5] = as.integer(res.df[, 5])
  res.df[, 6] = as.integer(res.df[, 6])
  res.df
}, digits = 3)


