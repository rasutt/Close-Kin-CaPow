# Reactive implied parameter outputs

# Population growth rate
output$lambda = renderText({
  paste("Population growth rate (lambda):", input$rho + input$phi)
})
# Expected super-population size
output$exp.Ns = renderText({
  paste("Expected super-population size (Ns):", round(exp.Ns.rct()))
})
# Plot expected population size over time
output$expPopPlot <- renderPlot({
  plot(
    (f.year.rct() - input$hist.len + 1):f.year.rct(), exp.N.t.rct(), 
    col = 'red', lwd = 2, t = 'l', ylim = c(0, max(exp.N.t.rct())),
    xlab = 'Year', ylab = 'E(N[t])', main = "Expected population size over time"
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