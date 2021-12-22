# Function to compare estimates from POPAN and close kin models
ComparisonPlot <- function(ests.lst, par, true.val) {
  # Find empirical CVs and biases
  ests.mean = sapply(ests.lst, mean)
  ests.sd = sapply(ests.lst, sd)
  if (true.val == 0) {
    # For proportional differences
    ests.bias <- ests.mean * 100
    ests.cv <- ests.sd * 100
  } else {
    # For standard estimates
    ests.bias <- (ests.mean - true.val) / true.val * 100
    ests.cv <- ests.sd / ests.mean * 100
  }
  
  # Plot estimates
  boxplot(
    ests.lst, main = par, show.names = T,
    sub = paste0(
      "Bias (1DP): ", paste0(round(ests.bias, 1), "%", collapse = ", "), "\n",
      "CV (1DP): ", paste0(round(ests.cv, 1), "%", collapse = ", ")
    )
  )
  abline(h = true.val, col = 'red')
}