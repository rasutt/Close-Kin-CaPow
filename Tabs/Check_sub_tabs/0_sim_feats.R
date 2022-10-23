# Outputs for simulation features sub-tab of check-tab

# Display simulation parameter values
output$lastParVals <- renderTable({
  # Make data frame for display
  par.vals.df(par.vals(), par.names())
}, digits = 3)

# Display last simulation values
output$lastSimOpts = renderTable(sim.opts.bio.scen()[1:2])

# Display last simulation values
output$lastBioScen = renderTable(sim.opts.bio.scen()[-(1:2)])

