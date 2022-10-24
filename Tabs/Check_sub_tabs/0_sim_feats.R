# Outputs for simulation features sub-tab of check-tab

# Display simulation parameter values
output$lastParsSltd <- renderTable({
  # Make data frame for display
  par.vals.df(
    c(phi(), rho(), exp.N.base(), base.yr(), paste(srvy.yrs(), collapse = ", "), 
      p()),
    c("Survival rate", "Per capita birth rate", "Expected population size in 
      base year", "Base year for expected population size", "Survey years",
      "Base level capture probability"),
    NULL
  )
}, digits = 3)

# Display implied parameter values
output$lastParsImpld <- renderTable(
  frmt.pars.impld(lambda(), exp.N.t()[hist.len()], exp.Ns()), digits = 3
)

# Display last simulation values
output$lastSimOpts = renderTable(sim.opts.bio.scen()[1:2])

# Display last simulation values
output$lastBioScen = renderTable(sim.opts.bio.scen()[-(1:2)])

