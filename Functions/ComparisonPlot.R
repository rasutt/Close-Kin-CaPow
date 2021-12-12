# Function to compare estimates from POPAN and close kin models
ComparisonPlot <- function(ests.lst, par, true.val) {
  # Find empirical CVs and biases
  ests.cv <- sapply(ests.lst, sd) * 100
  ests.mean <- sapply(ests.lst, mean)
  ests.bias <- (ests.mean - true.val) * 100

  # If not proportional differences
  if (!true.val == 0) {
    ests.cv <- ests.cv / ests.mean
    ests.bias <- ests.bias / true.val
  }
  
  # Plot estimates
  boxplot(
    ests.lst,
    main = par,
    sub = paste0(
      "Bias (1DP): ", paste0(round(ests.bias, 1), "%", collapse = ", "), "\n",
      "CV (1DP): ", paste0(round(ests.cv, 1), "%", collapse = ", ")
    )
  )
  abline(h = true.val, col = 'red')
}