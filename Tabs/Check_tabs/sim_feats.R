# Outputs for simulation features sub-tab of check-tab

# Display current simulation parameter values
output$currParsSltd <- renderTable({
  # Make data frame for display
  par.vals.df(
    c(phi(), rho(), exp.N.base(), base.yr(), paste(srvy.yrs(), collapse = ", "), 
      p(), L()),
    c("Survival rate", "Per capita birth rate", "Expected population size in 
      base year", "Base year for expected population size", "Survey years",
      "Base level capture probability", "Number of loci"),
    NULL
  )
}, digits = 3)

# Display current implied parameter values
output$currParsImpld <- renderTable({
  frmt.pars.impld(lambda(), beta(), N.init(), exp.N.t()[hist.len()], exp.Ns())
}, digits = 3)

# Display current simulation values
output$currSimOpts = renderTable(sim.opts.bio.scen()[1:2])

# Display current simulation values
output$currBioScen = renderTable(sim.opts.bio.scen()[-(1:2)])

# Display numbers of datasets removed and reasons
output$nsDtstsRmvd = renderTable({
  par.vals.df(
    c(sim.lst()$n.extinct, sim.lst()$n.no.pairs, n.sims()), 
    c("Number of extinct populations", "Number of studies with < 2 samples", 
      "Number of datasets retained"),
    NULL
  )
})