# Code and outputs for first study estimates sub-tab of check-tab

# Table of estimates from first study optimising genopair likelihood
output$firstGPEsts = renderTable({
  # Create general optimizer starting-values and bounds, NAs filled in below
  ck.start <- c(rho(), phi(), attributes(fst.std())$N.t.vec[hist.len()])
  ck.lwr <- c(0, 0.75, attributes(fst.std())$ns.caps[k()])
  ck.upr <- c(0.35, 1, Inf)
  
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
  }
  
  # Try to fit genopair likelihood model
  print(table(frst.smp.yr.ind.prs()[, 1], frst.smp.yr.ind.prs()[, 2]))
  print(str(gpp.slct))
  gp.tmb = TryGenopairTMB(
    if (any(all_undrflw)) gpp.adj else gpp.slct, 
    frst.smp.yr.ind.prs(),
    # smp.yr.ind.prs()[wtn.prs.inds, ], 
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
  
  gp.tmb$est.se.df
})
