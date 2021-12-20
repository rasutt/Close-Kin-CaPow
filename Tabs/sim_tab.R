# Reactive implied parameter outputs

# Population growth rate
output$lambda = renderText({
  paste("Population growth rate (lambda):", input$rho + input$phi)
})
# Number of surveys
output$k = renderText(paste("Number of surveys (k):", k.rct()))
# Final survey year
output$f.year = renderText(paste("Final survey year:", f.year.rct()))