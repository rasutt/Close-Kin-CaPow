# Reactive implied parameter outputs

# Population growth rate
output$lambda = renderText({
  paste("Population growth rate (lambda):", input$rho + input$phi)
})
# Number of surveys
output$k = renderText(paste("Number of surveys (k):", k.rct()))
# Final survey year
output$f.year = renderText(paste("Final survey year:", f.year.rct()))

# Plot expected population size over time
output$expPopPlot <- renderPlot({
  plot(
    (f.year.rct() - input$hist.len + 1):f.year.rct(), exp.N.t.rct(), 
    col = 'red', lwd = 2, t = 'l',
    xlab = 'Year', ylab = 'Nt', main = "Expected population size over time"
  )
  # Surveys
  abline(v = srvy.yrs.rct(), lty = 2)
})