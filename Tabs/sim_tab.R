# Reactive implied parameter outputs

# Population growth rate
output$lambda = renderText({
  paste("Population growth rate (lambda):", input$rho + input$phi)
})

# Plot expected population size over time
output$expPopPlot <- renderPlot({
  plot(
    (f.year.rct() - input$hist.len + 1):f.year.rct(), exp.N.sim(), 
    col = 'red', lwd = 2, t = 'l', ylim = c(0, max(exp.N.sim())),
    xlab = 'Year', ylab = 'Nt', main = "Expected population size over time"
  )
  abline(v = input$base.yr, h = input$exp.N.base, col = 2)
  # Surveys
  abline(v = srvy.yrs.rct(), lty = 2)
  # Add legend
  legend(
    "topleft", 
    legend = c(
      "Expected population size", 
      "Population size in base year", 
      "Survey years"
    ),
    col = c(2, 2, 1),
    lwd = c(2, 1, 1),
    lty = c(1, 1, 2)
  )
})