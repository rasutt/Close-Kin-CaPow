# Simulation and outputs bound to simulate button

load("sims_and_checks.Rdata")
sim.lst = reactiveVal(sims)

# Simulate population and capture histories
observeEvent(input$simulate, {
  # Initial population size
  N.init = round(exp.N.t()[1])
  
  # Create list for population and capture histories
  hists.lst <- vector("list", n_sims())
  # Create vectors for final and super-population sizes
  N.fin.vec <- Ns.vec <- numeric(n_sims())
  
  # Loop over histories
  withProgress({
    for (hist.ind in 1:n_sims()) {
      # Simulate family and capture histories of population of animals over
      # time
      hists.lst[[hist.ind]] <- SimPopStud(
        phi(), lambda(), N.init, hist.len(), srvy.yrs(), k(), f.year(), p(),
        clvng.p(), tmp.emgn(), alpha(), clvng.ints()
      )
      # Collect final and super-population sizes
      N.fin.vec[hist.ind] <- tail(attributes(hists.lst[[hist.ind]])$N.t.vec, 1)
      Ns.vec[hist.ind] <- attributes(hists.lst[[hist.ind]])$Ns
      
      # Update progress. Unexplained "Error in as.vector: object 'x' not
      # found" seen 19/12/2021 coming from incProgress...
      incProgress(1/n_sims())
    }
  }, value = 0, message = "Simulating populations")
  
  sim.lst(list(hists.lst = hists.lst, N.fin.vec = N.fin.vec, Ns.vec = Ns.vec))
})

## Last simulation

# Display simulation parameter values
output$simParVals <- renderTable({
  # Force simulation when starting app in sim tab
  # req(sim.lst())
  
  # Make data frame for display
  par_vals_df(sim.par.vals(), sim.par.names())
}, digits = 3)

# Display last simulation values
output$lastSimVals = renderTable(sim.vals())

# Plot expected population size over time
output$simExpPop <- renderPlot({
  plot(
    (f.year() - input$hist.len + 1):f.year(), exp.N.t(), 
    col = 'red', lwd = 2, t = 'l', ylim = c(0, max(exp.N.t())),
    xlab = 'Year', ylab = 'Exp_N_t', 
    main = "Expected population size over time"
  )
  # Base year
  abline(v = base.yr(), h = exp.N.base(), col = 2)
  # Surveys
  abline(v = srvy.yrs(), lty = 2)
  # Add legend
  legend(
    "topleft", legend = c("Over time", "In base year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
})

## Next simulation

# Display selected parameter values
output$selParVals <- renderTable({
  # Make data frame for display
  par_vals_df(sim.par.vals.rct(), sim.par.names.rct())
}, digits = 3)

# Display next simulation values
output$nextSimVals = renderTable(sim.vals.rct())

# Plot expected population size over time
output$selExpPop <- renderPlot({
  plot(
    (f.year.rct() - input$hist.len + 1):f.year.rct(), exp.N.t.rct(), 
    col = 'red', lwd = 2, t = 'l', ylim = c(0, max(exp.N.t.rct())),
    xlab = 'Year', ylab = 'Exp_N_t', 
    main = "Expected population size over time"
  )
  # Base year
  abline(v = input$base.yr, h = input$exp.N.base, col = 2)
  # Surveys
  abline(v = srvy.yrs.rct(), lty = 2)
  # Add legend
  legend(
    "topleft", legend = c("Over time", "In base year", "Survey years"),
    col = c(2, 2, 1), lwd = c(2, 1, 1), lty = c(1, 1, 2)
  )
})