# Compare model results
CompMods <- function() {
  # View estimates for first few populations
  print(lapply(mod.ests.lst, function(mod.ests) head(round(mod.ests, 3))))
  
  # View convergence codes.  1: iteration limit, 52: error, see message
  # component
  print(lapply(mod.ests.lst, function(mod.ests) table(mod.ests[, "cnvg"])))
  
  # Find optimisation attempts that converged
  cvgd.ests.lst <- lapply(mod.ests.lst, function(ests.mat)
    ests.mat[!ests.mat[, "cnvg"], ])
  
  # Plot estimates from all models side-by-side
  
  # Open file to save plots in
  pdf(file = paste0("Plots/", sim.name, "_comparison_plot.pdf"))
  
  # Set four plots per page
  par(mfrow = c(2, 2))
  
  # Plot estimates for lambda
  ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 1]),
                 "Lambda", lambda)
  
  # Plot estimates for Phi
  ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 2]),
                 "Phi", phi)
  
  # Plot estimates of N_2020
  ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 3]),
                 "N_2020", exp.N.t[hist.len])
  
  # Plot estimates of superpopulation size
  ComparisonPlot(lapply(cvgd.ests.lst, function(ests.mat) ests.mat[, 4]),
                 "Ns", exp.Ns)
  
  # Close file with plots saved
  dev.off()
  
  # Check correlation between estimates from POPAN and close kin models when
  # both converged
  bth.cvgd <- !ppn.tmb.ests[, "cnvg"] & !ck.tmb.ests[, "cnvg"]
  ppn.bth.cvgd <- ppn.tmb.ests[bth.cvgd, ]
  ck.bth.cvgd <- ck.tmb.ests[bth.cvgd, ]
  
  plot(ppn.bth.cvgd[, 1], ck.bth.cvgd[, 1], main = "Lambda", xlab = "POPAN", 
       ylab = "close kin", asp = 1,
       sub = paste("R2 = ", round(cor(ppn.bth.cvgd[, 1], ck.bth.cvgd[, 1]), 2)))
  plot(ppn.bth.cvgd[, 2], ck.bth.cvgd[, 2], main = "Phi", xlab = "POPAN", 
       ylab = "close kin", asp = 1,
       sub = paste("R2 = ", round(cor(ppn.bth.cvgd[, 2], ck.bth.cvgd[, 2]), 2)))
  
  # Check POPAN model estimates of capture probabilities
  boxplot(cvgd.ests.lst$ppn_tmb[, 5:(k + 4)], 
          main = "Capture probability estimates")
}
