# Output for populations sub-tab in checks tab

# Plot population sizes over time
output$checkExpPop = renderPlot({
  # Plot population trajectories
  matplot(
    sim.yrs(), t(checks.lst()$N.t.mat), type = 'l',
    col = rgb(0, 0, 0, alpha = 0.1), lty = 1, 
    ylim = c(0, max(checks.lst()$N.t.mat)),
    xlab = 'Year', ylab = 'Nt', main = "Population sizes over time"
  )
  # Expected value
  lines(sim.yrs(), exp.N.t(), col = 2, lwd = 2)
  # Average value
  lines(sim.yrs(), colMeans(checks.lst()$N.t.mat), col = 4, lwd = 2)
  # Surveys
  abline(v = srvy.yrs(), lty = 2)
  # Add legend
  legend(
    "topleft", legend = c("Expected", "Simulated", "Mean", "Survey years"),
    col = c(2, rgb(0, 0, 0, alpha = 0.1), 4, 1), lwd = c(2, 1, 2, 1), 
    lty = c(1, 1, 1, 2)
  )
})
